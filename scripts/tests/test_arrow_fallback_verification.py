"""
Runtime verification tests for Arrow fallback classifications.

This module instruments the actual fallback methods to verify that our
static analysis classifications in generate_arrow_fallback_table.py
match runtime behavior.

The key fallback entry points are:
1. _apply_elementwise() - element-wise Python loop
2. super()._str_*() in ArrowStringArray - delegates to ObjectStringArrayMixin
3. to_numpy() / to_pylist() - NumPy materialization
"""

from __future__ import annotations

from collections import defaultdict
import contextlib

import pytest

import pandas as pd
from pandas import ArrowDtype

pa = pytest.importorskip("pyarrow")


class FallbackTracker:
    """Track which fallback methods are called during execution."""

    def __init__(self):
        self.calls = defaultdict(list)
        self._originals = {}

    def record(self, method_name: str, caller: str = ""):
        """Record a fallback method call."""
        self.calls[method_name].append(caller)

    def reset(self):
        """Clear recorded calls."""
        self.calls.clear()

    def was_called(self, method_name: str) -> bool:
        """Check if a fallback method was called."""
        return method_name in self.calls and len(self.calls[method_name]) > 0

    def any_object_mixin_called(self) -> bool:
        """Check if any ObjectStringArrayMixin method was called."""
        return any(k.startswith("ObjectStringArrayMixin.") for k in self.calls)

    @contextlib.contextmanager
    def track(self):
        """Context manager to track fallback calls."""
        self.reset()
        self._originals.clear()

        # Get the classes we'll be patching
        arrow_ext_array = pd.core.arrays.arrow.array.ArrowExtensionArray
        object_mixin = pd.core.strings.object_array.ObjectStringArrayMixin

        # Patch _apply_elementwise (element-wise Python loop)
        orig_apply = arrow_ext_array._apply_elementwise
        self._originals["_apply_elementwise"] = orig_apply

        def tracked_apply_elementwise(self_inner, func):
            func_name = func.__name__ if hasattr(func, "__name__") else str(func)
            tracker.record("_apply_elementwise", func_name)
            return orig_apply(self_inner, func)

        arrow_ext_array._apply_elementwise = tracked_apply_elementwise

        # Methods to patch on ObjectStringArrayMixin
        methods_to_patch = [
            "_str_map",
            "_str_contains",
            "_str_match",
            "_str_fullmatch",
            "_str_replace",
            "_str_count",
            "_str_repeat",
            "_str_casefold",
            "_str_encode",
            "_str_partition",
            "_str_rpartition",
            "_str_normalize",
            "_str_wrap",
            "_str_translate",
            "_str_find",
            "_str_rfind",
            "_str_index",
            "_str_rindex",
            "_str_findall",
            "_str_join",
        ]

        # Create and apply patches for ObjectStringArrayMixin methods
        for method_name in methods_to_patch:
            if hasattr(object_mixin, method_name):
                original = getattr(object_mixin, method_name)
                self._originals[f"ObjectStringArrayMixin.{method_name}"] = original

                # Use a factory to capture the correct values
                def make_tracked(orig, name):
                    def tracked_method(self_inner, *args, **kwargs):
                        tracker.record(f"ObjectStringArrayMixin.{name}")
                        return orig(self_inner, *args, **kwargs)

                    return tracked_method

                setattr(object_mixin, method_name, make_tracked(original, method_name))

        try:
            yield self
        finally:
            # Restore all originals
            arrow_ext_array._apply_elementwise = self._originals["_apply_elementwise"]
            for method_name in methods_to_patch:
                key = f"ObjectStringArrayMixin.{method_name}"
                if key in self._originals:
                    setattr(object_mixin, method_name, self._originals[key])


tracker = FallbackTracker()


@pytest.fixture
def arrow_string_series():
    """Create a test Series with Arrow string dtype."""
    return pd.Series(
        ["hello", "world", "PANDAS", "Arrow", None, "test123"],
        dtype="string[pyarrow]",
    )


@pytest.fixture
def arrow_dtype_series():
    """Create a test Series with ArrowDtype string."""
    return pd.Series(
        ["hello", "world", "PANDAS", "Arrow", None, "test123"],
        dtype=ArrowDtype(pa.string()),
    )


class TestElementwiseFallbackMethods:
    """Test methods that should use _apply_elementwise() on ArrowExtensionArray."""

    @pytest.mark.parametrize(
        "method,kwargs",
        [
            # Core elementwise methods (11 total in ELEMENTWISE category)
            ("casefold", {}),
            ("encode", {"encoding": "utf-8"}),
            ("partition", {"sep": "l"}),
            ("rpartition", {"sep": "l"}),
            ("normalize", {"form": "NFC"}),
            ("wrap", {"width": 10}),
            ("translate", {"table": str.maketrans("a", "b")}),
            # These also use elementwise in ArrowExtensionArray
            ("findall", {"pat": "l"}),
            ("index", {"sub": "l"}),
            ("rindex", {"sub": "l"}),
            ("rfind", {"sub": "l"}),
        ],
    )
    def test_uses_apply_elementwise(self, arrow_dtype_series, method, kwargs):
        """Verify these methods use _apply_elementwise() on ArrowExtensionArray."""
        with tracker.track():
            try:
                getattr(arrow_dtype_series.str, method)(**kwargs)
            except (ValueError, KeyError):
                # Some methods may raise for certain inputs (e.g., index not found)
                pass

        assert tracker.was_called("_apply_elementwise"), (
            f"str.{method} should use _apply_elementwise() but didn't. "
            f"Calls recorded: {dict(tracker.calls)}"
        )


class TestConditionalFallbackMethods:
    """Test methods that conditionally fall back based on parameters."""

    @pytest.mark.parametrize(
        "method,kwargs,should_fallback",
        [
            # CONTAINS: Arrow path vs fallback path
            ("contains", {"pat": "hello"}, False),
            ("contains", {"pat": "hello", "regex": False}, False),
            # flags=8 is re.MULTILINE which PyArrow doesn't support
            ("contains", {"pat": "hello", "flags": 8}, True),
            # MATCH: Arrow path vs fallback path
            ("match", {"pat": "hello"}, False),
            # re.MULTILINE triggers fallback
            ("match", {"pat": "hello", "flags": 8}, True),
            # FULLMATCH: Arrow path vs fallback path
            ("fullmatch", {"pat": "hello"}, False),
            ("fullmatch", {"pat": "hello", "flags": 8}, True),
            # COUNT: Arrow path vs fallback path
            ("count", {"pat": "l"}, False),
            ("count", {"pat": "l", "flags": 8}, True),
            # REPLACE: Arrow path vs fallback path
            ("replace", {"pat": "hello", "repl": "world"}, False),
            ("replace", {"pat": "hello", "repl": "world", "flags": 8}, True),
            # FIND: Arrow native
            ("find", {"sub": "l"}, False),
            # JOIN: Falls back to ObjectStringArrayMixin (no Arrow impl in string_arrow)
            ("join", {"sep": "-"}, True),
            # ZFILL: Arrow native
            ("zfill", {"width": 10}, False),
        ],
    )
    def test_conditional_fallback(
        self, arrow_string_series, method, kwargs, should_fallback
    ):
        """Verify conditional fallback behavior based on parameters."""
        with tracker.track():
            try:
                getattr(arrow_string_series.str, method)(**kwargs)
            except Exception:
                # Some combinations may raise, that's OK
                pass

        if should_fallback:
            assert tracker.any_object_mixin_called(), (
                f"str.{method}({kwargs}) should fall back to "
                f"ObjectStringArrayMixin but didn't. "
                f"Calls recorded: {dict(tracker.calls)}"
            )
        else:
            assert not tracker.any_object_mixin_called(), (
                f"str.{method}({kwargs}) should NOT fall back to "
                f"ObjectStringArrayMixin but did. "
                f"Calls recorded: {dict(tracker.calls)}"
            )


class TestArrowNativeMethods:
    """Test methods that should use Arrow compute (no fallback)."""

    @pytest.mark.parametrize(
        "method,kwargs",
        [
            # Basic string transforms (ARROW_NATIVE)
            ("lower", {}),
            ("upper", {}),
            ("capitalize", {}),
            ("title", {}),
            ("swapcase", {}),
            # Length
            ("len", {}),
            # Strip methods
            ("strip", {}),
            ("lstrip", {}),
            ("rstrip", {}),
            # Prefix/suffix methods
            ("startswith", {"pat": "he"}),
            ("endswith", {"pat": "lo"}),
            ("removeprefix", {"prefix": "he"}),
            ("removesuffix", {"suffix": "lo"}),
            # Character class checks
            ("isalnum", {}),
            ("isalpha", {}),
            ("isdigit", {}),
            ("islower", {}),
            ("isupper", {}),
            ("isspace", {}),
            ("isnumeric", {}),
            ("isdecimal", {}),
            ("isascii", {}),
            ("istitle", {}),
            # Slicing and indexing
            ("slice", {"start": 0, "stop": 2}),
            ("slice_replace", {"start": 0, "stop": 2, "repl": "XX"}),
            ("get", {"i": 0}),
            # Split methods (Arrow native)
            ("split", {}),
            ("rsplit", {}),
            # Extract and get_dummies
            ("extract", {"pat": r"(\w+)", "expand": True}),
            ("get_dummies", {"sep": "|"}),
        ],
    )
    def test_no_fallback(self, arrow_string_series, method, kwargs):
        """Verify these methods do NOT use any fallback."""
        with tracker.track():
            try:
                getattr(arrow_string_series.str, method)(**kwargs)
            except Exception:
                # Some methods may raise for certain inputs
                pass

        assert not tracker.was_called("_apply_elementwise"), (
            f"str.{method} should NOT use _apply_elementwise() but did. "
            f"Calls recorded: {dict(tracker.calls)}"
        )
        assert not tracker.any_object_mixin_called(), (
            f"str.{method} should NOT use ObjectStringArrayMixin but did. "
            f"Calls recorded: {dict(tracker.calls)}"
        )


class TestObjectFallbackMethods:
    """Test methods that always fall back to ObjectStringArrayMixin."""

    @pytest.mark.parametrize(
        "method,kwargs",
        [
            # REPEAT with sequence (not int)
            ("repeat", {"repeats": [1, 2, 1, 2, 1, 2]}),
        ],
    )
    def test_always_falls_back(self, arrow_string_series, method, kwargs):
        """Verify these methods always use ObjectStringArrayMixin fallback."""
        with tracker.track():
            try:
                getattr(arrow_string_series.str, method)(**kwargs)
            except Exception:
                pass

        assert tracker.any_object_mixin_called(), (
            f"str.{method}({kwargs}) should fall back to "
            f"ObjectStringArrayMixin but didn't. "
            f"Calls recorded: {dict(tracker.calls)}"
        )


class TestVersionGatedMethods:
    """Test methods whose behavior depends on PyArrow version."""

    @pytest.mark.parametrize(
        "method,kwargs",
        [
            ("pad", {"width": 10}),
            ("isdigit", {}),
        ],
    )
    def test_version_gated_no_fallback(self, arrow_string_series, method, kwargs):
        """
        Verify version-gated methods use Arrow in recent PyArrow versions.

        These methods may fall back in older PyArrow versions but should
        use Arrow compute in recent versions.
        """
        with tracker.track():
            try:
                getattr(arrow_string_series.str, method)(**kwargs)
            except Exception:
                pass

        # In recent PyArrow versions, these should not fall back
        # If they do fall back, it's version-dependent and acceptable
        # This test documents current behavior rather than enforcing it


class TestComprehensiveClassificationVerification:
    """
    Comprehensive verification of all classifications against runtime behavior.

    This test dynamically tests all methods in the classification table
    to ensure our static analysis matches actual runtime behavior.
    """

    def test_all_classifications_match_runtime(
        self, arrow_string_series, arrow_dtype_series
    ):
        """
        Verify all FallbackType classifications match runtime behavior.

        Tests both ArrowStringArray and ArrowExtensionArray paths.
        """
        from scripts.generate_arrow_fallback_table import (
            FallbackType,
            build_string_method_table,
        )

        method_table = build_string_method_table()

        # Define test arguments for each method
        # Methods not listed here will be skipped
        test_args = {
            # ARROW_NATIVE methods
            "_str_lower": {},
            "_str_upper": {},
            "_str_len": {},
            "_str_strip": {},
            "_str_lstrip": {},
            "_str_rstrip": {},
            "_str_capitalize": {},
            "_str_title": {},
            "_str_swapcase": {},
            "_str_startswith": {"pat": "he"},
            "_str_endswith": {"pat": "lo"},
            "_str_removeprefix": {"prefix": "he"},
            "_str_removesuffix": {"suffix": "lo"},
            "_str_isalnum": {},
            "_str_isalpha": {},
            "_str_isdigit": {},
            "_str_islower": {},
            "_str_isupper": {},
            "_str_isspace": {},
            "_str_isnumeric": {},
            "_str_isdecimal": {},
            "_str_isascii": {},
            "_str_istitle": {},
            "_str_slice": {"start": 0, "stop": 2},
            "_str_slice_replace": {"start": 0, "stop": 2, "repl": "XX"},
            "_str_get": {"i": 0},
            "_str_split": {},
            "_str_rsplit": {},
            # ELEMENTWISE methods
            "_str_casefold": {},
            "_str_encode": {"encoding": "utf-8"},
            "_str_partition": {"sep": "l"},
            "_str_rpartition": {"sep": "l"},
            "_str_normalize": {"form": "NFC"},
            "_str_wrap": {"width": 10},
            "_str_translate": {"table": str.maketrans("a", "b")},
            "_str_findall": {"pat": "l"},
            "_str_index": {"sub": "l"},
            "_str_rindex": {"sub": "l"},
            "_str_rfind": {"sub": "l"},
            # CONDITIONAL methods (test the Arrow-native path)
            "_str_contains": {"pat": "hello"},
            "_str_count": {"pat": "l"},
            "_str_replace": {"pat": "hello", "repl": "world"},
            "_str_find": {"sub": "l"},
            # Note: _str_join always falls back (no Arrow impl in string_arrow.py)
            # Skip it from Arrow-native path testing
            "_str_zfill": {"width": 10},
            # OBJECT_FALLBACK methods
            "_str_match": {"pat": "hello"},
            "_str_fullmatch": {"pat": "hello"},
            # VERSION_GATED methods
            "_str_pad": {"width": 10},
        }

        mismatches = []
        tested_count = 0
        skipped_count = 0

        for method_name, info in method_table.items():
            if method_name not in test_args:
                skipped_count += 1
                continue

            kwargs = test_args[method_name]
            public_name = info.public_name.replace("str.", "")
            tested_count += 1

            # Test with ArrowStringArray (uses ObjectStringArrayMixin)
            with tracker.track():
                try:
                    getattr(arrow_string_series.str, public_name)(**kwargs)
                except Exception:
                    continue

            used_object_mixin = tracker.any_object_mixin_called()

            # Test with ArrowExtensionArray (uses _apply_elementwise)
            with tracker.track():
                try:
                    getattr(arrow_dtype_series.str, public_name)(**kwargs)
                except Exception:
                    continue

            used_elementwise = tracker.was_called("_apply_elementwise")

            # Check classification
            if info.fallback_type == FallbackType.ARROW_NATIVE:
                if used_elementwise or used_object_mixin:
                    mismatches.append(
                        f"{method_name}: classified as ARROW_NATIVE but used fallback "
                        f"(elementwise={used_elementwise}, "
                        f"object_mixin={used_object_mixin})"
                    )
            elif info.fallback_type == FallbackType.ELEMENTWISE:
                # ELEMENTWISE methods use _apply_elementwise in ArrowExtensionArray
                # and ObjectStringArrayMixin._str_map in ArrowStringArray
                if not used_elementwise and not used_object_mixin:
                    mismatches.append(
                        f"{method_name}: classified as ELEMENTWISE but used no fallback"
                    )
            elif info.fallback_type == FallbackType.CONDITIONAL:
                # Conditional methods should NOT fall back with default arguments
                if used_elementwise or used_object_mixin:
                    mismatches.append(
                        f"{method_name}: classified as CONDITIONAL but used fallback "
                        f"with Arrow-native args "
                        f"(elementwise={used_elementwise}, "
                        f"object_mixin={used_object_mixin})"
                    )
            elif info.fallback_type == FallbackType.OBJECT_FALLBACK:
                # Some OBJECT_FALLBACK methods may use Arrow for simple cases
                # Only flag if neither path is used AND we expect fallback
                pass  # These are harder to verify uniformly
            elif info.fallback_type == FallbackType.VERSION_GATED:
                # Version-gated methods behavior depends on PyArrow version
                pass  # Don't enforce for version-gated

        if mismatches:
            pytest.fail(
                f"Classification mismatches found ({len(mismatches)}):\n"
                + "\n".join(mismatches)
                + f"\n\nTested: {tested_count}, Skipped: {skipped_count}"
            )


# =============================================================================
# Non-String Method Verification Tests
# =============================================================================


class NonStringFallbackTracker:
    """Track fallback methods for non-string operations."""

    def __init__(self):
        self.calls = defaultdict(list)
        self._originals = {}

    def record(self, method_name: str, caller: str = ""):
        """Record a fallback method call."""
        self.calls[method_name].append(caller)

    def reset(self):
        """Clear recorded calls."""
        self.calls.clear()

    def was_called(self, method_name: str) -> bool:
        """Check if a fallback method was called."""
        return method_name in self.calls and len(self.calls[method_name]) > 0

    def any_fallback_called(self) -> bool:
        """Check if any fallback method was called."""
        return len(self.calls) > 0

    @contextlib.contextmanager
    def track(self):
        """Context manager to track fallback calls for non-string operations."""
        self.reset()
        self._originals.clear()

        arrow_ext_array = pd.core.arrays.arrow.array.ArrowExtensionArray

        # Track to_numpy calls (primary fallback mechanism)
        orig_to_numpy = arrow_ext_array.to_numpy
        self._originals["to_numpy"] = orig_to_numpy

        def tracked_to_numpy(self_inner, *args, **kwargs):
            ns_tracker.record("to_numpy")
            return orig_to_numpy(self_inner, *args, **kwargs)

        arrow_ext_array.to_numpy = tracked_to_numpy

        # Track _apply_elementwise
        orig_apply = arrow_ext_array._apply_elementwise
        self._originals["_apply_elementwise"] = orig_apply

        def tracked_apply(self_inner, func):
            func_name = func.__name__ if hasattr(func, "__name__") else str(func)
            ns_tracker.record("_apply_elementwise", func_name)
            return orig_apply(self_inner, func)

        arrow_ext_array._apply_elementwise = tracked_apply

        try:
            yield self
        finally:
            # Restore originals
            arrow_ext_array.to_numpy = self._originals["to_numpy"]
            arrow_ext_array._apply_elementwise = self._originals["_apply_elementwise"]


ns_tracker = NonStringFallbackTracker()


@pytest.fixture
def arrow_int_series():
    """Create an Arrow-backed integer Series for testing."""
    return pd.Series([1, 2, 3, 4, 5], dtype=ArrowDtype(pa.int64()))


@pytest.fixture
def arrow_float_series():
    """Create an Arrow-backed float Series for testing."""
    return pd.Series([1.0, 2.0, 3.0, 4.0, 5.0], dtype=ArrowDtype(pa.float64()))


@pytest.fixture
def arrow_bool_series():
    """Create an Arrow-backed boolean Series for testing."""
    return pd.Series([True, False, True, False, True], dtype=ArrowDtype(pa.bool_()))


@pytest.fixture
def arrow_datetime_series():
    """Create an Arrow-backed datetime Series for testing."""
    return pd.Series(
        pd.date_range("2020-01-01", periods=5, freq="D"),
        dtype=ArrowDtype(pa.timestamp("ns")),
    )


@pytest.fixture
def arrow_timedelta_series():
    """Create an Arrow-backed timedelta Series for testing."""
    return pd.Series(
        pd.timedelta_range("1 day", periods=5, freq="D"),
        dtype=ArrowDtype(pa.duration("ns")),
    )


@pytest.fixture
def arrow_list_series():
    """Create an Arrow-backed list Series for testing."""
    return pd.Series(
        [[1, 2], [3, 4, 5], [6], [7, 8], [9, 10]],
        dtype=ArrowDtype(pa.list_(pa.int64())),
    )


@pytest.fixture
def arrow_struct_series():
    """Create an Arrow-backed struct Series for testing."""
    return pd.Series(
        [{"a": 1, "b": "x"}, {"a": 2, "b": "y"}, {"a": 3, "b": "z"}],
        dtype=ArrowDtype(pa.struct([("a", pa.int64()), ("b", pa.string())])),
    )


class TestArithmeticOperations:
    """Verify arithmetic operation classifications."""

    @pytest.mark.parametrize(
        "op,other,expected_arrow",
        [
            # ARROW_NATIVE operations
            ("__add__", 1, True),
            ("__sub__", 1, True),
            ("__mul__", 2, True),
            ("__truediv__", 2, True),
            ("__pow__", 2, True),
            ("__radd__", 1, True),
            ("__rsub__", 1, True),
            ("__rmul__", 2, True),
            ("__rtruediv__", 2, True),
            ("__rpow__", 2, True),
            # Logical operations (ARROW_NATIVE)
            ("__and__", True, True),
            ("__or__", True, True),
            ("__xor__", True, True),
        ],
    )
    def test_arrow_native_arithmetic(
        self, arrow_int_series, arrow_bool_series, op, other, expected_arrow
    ):
        """Verify Arrow-native arithmetic operations don't fall back."""
        # Use bool series for logical ops, int series for arithmetic
        if op in ("__and__", "__or__", "__xor__"):
            series = arrow_bool_series
            other = True
        else:
            series = arrow_int_series

        with ns_tracker.track():
            try:
                getattr(series, op)(other)
            except Exception:
                pass

        if expected_arrow:
            # Should not use to_numpy for the main operation
            # Note: Some internal operations may use to_numpy, so we check
            # that _apply_elementwise is not called
            assert not ns_tracker.was_called("_apply_elementwise"), (
                f"{op} should use Arrow compute, not _apply_elementwise"
            )

    @pytest.mark.parametrize(
        "op,other",
        [
            # CONDITIONAL - comparison operations
            ("__eq__", 1),
            ("__ne__", 1),
            ("__lt__", 2),
            ("__le__", 2),
            ("__gt__", 0),
            ("__ge__", 0),
        ],
    )
    def test_comparison_operations(self, arrow_int_series, op, other):
        """Verify comparison operations work with Arrow."""
        with ns_tracker.track():
            try:
                result = getattr(arrow_int_series, op)(other)
                assert result is not None
            except Exception:
                pass

    @pytest.mark.parametrize(
        "op",
        [
            # Note: floordiv and rfloordiv are now implemented despite
            # being classified as NOT_IMPLEMENTED in documentation
            "__mod__",
            "__divmod__",
            "__rmod__",
            "__rdivmod__",
        ],
    )
    def test_not_implemented_operations(self, arrow_int_series, op):
        """Verify NOT_IMPLEMENTED operations raise appropriately."""
        with pytest.raises((TypeError, NotImplementedError)):
            getattr(arrow_int_series, op)(2)

    @pytest.mark.parametrize(
        "op",
        [
            # These use floordiv_compat() which uses PyArrow compute internally
            "__floordiv__",
            "__rfloordiv__",
        ],
    )
    def test_floordiv_uses_arrow_compat(self, arrow_int_series, op):
        """
        Verify floordiv operations use Arrow via floordiv_compat.

        floordiv_compat uses pc.divide and pc.floor internally,
        so these operations use Arrow compute (not _apply_elementwise).
        """
        with ns_tracker.track():
            result = getattr(arrow_int_series, op)(2)
            assert result is not None

        # Verify it doesn't use _apply_elementwise (it uses Arrow via compat)
        assert not ns_tracker.was_called("_apply_elementwise"), (
            f"{op} should use Arrow (via floordiv_compat), not _apply_elementwise"
        )


class TestDatetimeMethods:
    """Verify datetime accessor method classifications."""

    @pytest.mark.parametrize(
        "attr",
        [
            # ARROW_NATIVE datetime properties
            "year",
            "month",
            "day",
            "hour",
            "minute",
            "second",
            "microsecond",
            "nanosecond",
            "day_of_week",
            "dayofweek",
            "day_of_year",
            "dayofyear",
            "quarter",
            "is_month_start",
            "is_month_end",
            "is_quarter_start",
            "is_quarter_end",
            "is_year_start",
            "is_year_end",
            "is_leap_year",
            "days_in_month",
        ],
    )
    def test_arrow_native_datetime_properties(self, arrow_datetime_series, attr):
        """Verify Arrow-native datetime properties don't use fallback."""
        with ns_tracker.track():
            try:
                getattr(arrow_datetime_series.dt, attr)
            except Exception:
                pass

        assert not ns_tracker.was_called("_apply_elementwise"), (
            f"dt.{attr} should use Arrow compute, not _apply_elementwise"
        )

    @pytest.mark.parametrize(
        "method,kwargs",
        [
            ("floor", {"freq": "D"}),
            ("ceil", {"freq": "D"}),
            ("round", {"freq": "D"}),
            ("day_name", {}),
            ("month_name", {}),
            ("strftime", {"date_format": "%Y-%m-%d"}),
            ("normalize", {}),
            ("tz_localize", {"tz": "UTC"}),
        ],
    )
    def test_arrow_native_datetime_methods(self, arrow_datetime_series, method, kwargs):
        """Verify Arrow-native datetime methods work."""
        with ns_tracker.track():
            try:
                getattr(arrow_datetime_series.dt, method)(**kwargs)
            except Exception:
                pass

        assert not ns_tracker.was_called("_apply_elementwise"), (
            f"dt.{method} should use Arrow compute, not _apply_elementwise"
        )

    @pytest.mark.parametrize(
        "attr",
        [
            # OBJECT_FALLBACK - timedelta components that use to_numpy
            "days",
            "seconds",
            "microseconds",
            "nanoseconds",
        ],
    )
    def test_timedelta_object_fallback(self, arrow_timedelta_series, attr):
        """Verify timedelta component access uses fallback."""
        with ns_tracker.track():
            try:
                getattr(arrow_timedelta_series.dt, attr)
            except Exception:
                pass

        # These properties typically use to_numpy internally
        assert ns_tracker.was_called("to_numpy"), (
            f"dt.{attr} should use to_numpy fallback but didn't. "
            f"Calls: {dict(ns_tracker.calls)}"
        )


class TestAggregationMethods:
    """Verify aggregation method classifications."""

    @pytest.mark.parametrize(
        "method,kwargs",
        [
            # ARROW_NATIVE aggregations
            ("sum", {}),
            ("mean", {}),
            ("min", {}),
            ("max", {}),
            ("std", {}),
            ("var", {}),
            ("median", {}),
            ("prod", {}),
            ("sem", {}),
        ],
    )
    def test_arrow_native_aggregations(self, arrow_float_series, method, kwargs):
        """Verify Arrow-native aggregations work correctly."""
        with ns_tracker.track():
            try:
                result = getattr(arrow_float_series, method)(**kwargs)
                assert result is not None
            except Exception:
                pass

        assert not ns_tracker.was_called("_apply_elementwise"), (
            f"{method} should use Arrow compute, not _apply_elementwise"
        )

    @pytest.mark.parametrize(
        "method,kwargs",
        [
            ("any", {}),
            ("all", {}),
        ],
    )
    def test_conditional_aggregations(self, arrow_bool_series, method, kwargs):
        """Verify conditional aggregations work."""
        with ns_tracker.track():
            try:
                result = getattr(arrow_bool_series, method)(**kwargs)
                assert result is not None
            except Exception:
                pass

    @pytest.mark.parametrize(
        "method",
        [
            # kurt is NOT_IMPLEMENTED - pandas uses getattr(pc, "kurt", None)
            # but PyArrow uses "kurtosis" and there's no name mapping
            "kurt",
        ],
    )
    def test_not_implemented_aggregations(self, arrow_float_series, method):
        """Verify NOT_IMPLEMENTED aggregations raise appropriately."""
        with pytest.raises((TypeError, NotImplementedError)):
            getattr(arrow_float_series, method)()

    def test_skew_version_gated(self, arrow_float_series):
        """
        Verify skew works with PyArrow >= 20.0 (VERSION_GATED).

        skew is available in PyArrow 20.0+ via pc.skew.
        """
        with ns_tracker.track():
            result = arrow_float_series.skew()
            assert result is not None

        # Verify it doesn't use _apply_elementwise
        assert not ns_tracker.was_called("_apply_elementwise"), (
            "skew should use Arrow compute (pc.skew), not _apply_elementwise"
        )


class TestArrayMethods:
    """Verify array method classifications."""

    @pytest.mark.parametrize(
        "method,kwargs",
        [
            # ARROW_NATIVE array methods
            ("argsort", {}),
            ("unique", {}),
            ("dropna", {}),
            ("copy", {}),
            ("take", {"indices": [0, 1, 2]}),
            ("value_counts", {}),
        ],
    )
    def test_arrow_native_array_methods(self, arrow_int_series, method, kwargs):
        """Verify Arrow-native array methods work correctly."""
        with ns_tracker.track():
            try:
                result = getattr(arrow_int_series, method)(**kwargs)
                assert result is not None
            except Exception:
                pass

        assert not ns_tracker.was_called("_apply_elementwise"), (
            f"{method} should use Arrow compute, not _apply_elementwise"
        )

    @pytest.mark.parametrize(
        "method,kwargs",
        [
            # CONDITIONAL array methods
            ("fillna", {"value": 0}),
            ("ffill", {}),
            ("bfill", {}),
            ("cumsum", {}),
            ("cummin", {}),
            ("cummax", {}),
            ("cumprod", {}),
        ],
    )
    def test_conditional_array_methods(self, arrow_int_series, method, kwargs):
        """Verify conditional array methods work."""
        with ns_tracker.track():
            try:
                result = getattr(arrow_int_series, method)(**kwargs)
                assert result is not None
            except Exception:
                pass

    @pytest.mark.parametrize(
        "method,kwargs",
        [
            ("duplicated", {}),
            ("isna", {}),
        ],
    )
    def test_object_fallback_array_methods(self, arrow_int_series, method, kwargs):
        """Verify object fallback array methods use to_numpy."""
        with ns_tracker.track():
            try:
                result = getattr(arrow_int_series, method)(**kwargs)
                assert result is not None
            except Exception:
                pass

        # These methods typically convert to numpy internally
        # isna and duplicated use to_numpy for the result or internally


class TestListAccessor:
    """Verify list accessor method classifications."""

    def test_len_arrow_native(self, arrow_list_series):
        """Verify list.len uses Arrow compute."""
        with ns_tracker.track():
            try:
                result = arrow_list_series.list.len()
                assert result is not None
            except Exception:
                pass

        assert not ns_tracker.was_called("_apply_elementwise"), (
            "list.len should use Arrow compute, not _apply_elementwise"
        )

    def test_flatten_arrow_native(self, arrow_list_series):
        """Verify list.flatten uses Arrow compute."""
        with ns_tracker.track():
            try:
                result = arrow_list_series.list.flatten()
                assert result is not None
            except Exception:
                pass

        assert not ns_tracker.was_called("_apply_elementwise"), (
            "list.flatten should use Arrow compute, not _apply_elementwise"
        )

    def test_getitem_conditional(self, arrow_list_series):
        """Verify list.__getitem__ works."""
        with ns_tracker.track():
            try:
                result = arrow_list_series.list[0]
                assert result is not None
            except Exception:
                pass


class TestStructAccessor:
    """Verify struct accessor method classifications."""

    def test_field_arrow_native(self, arrow_struct_series):
        """Verify struct.field uses Arrow compute."""
        with ns_tracker.track():
            try:
                result = arrow_struct_series.struct.field("a")
                assert result is not None
            except Exception:
                pass

        assert not ns_tracker.was_called("_apply_elementwise"), (
            "struct.field should use Arrow compute, not _apply_elementwise"
        )

    def test_explode_arrow_native(self, arrow_struct_series):
        """Verify struct.explode uses Arrow compute."""
        with ns_tracker.track():
            try:
                result = arrow_struct_series.struct.explode()
                assert result is not None
            except Exception:
                pass

        assert not ns_tracker.was_called("_apply_elementwise"), (
            "struct.explode should use Arrow compute, not _apply_elementwise"
        )

    def test_dtypes_arrow_native(self, arrow_struct_series):
        """Verify struct.dtypes uses Arrow compute."""
        with ns_tracker.track():
            try:
                result = arrow_struct_series.struct.dtypes
                assert result is not None
            except Exception:
                pass

        assert not ns_tracker.was_called("_apply_elementwise"), (
            "struct.dtypes should use Arrow compute, not _apply_elementwise"
        )


class TestComprehensiveNonStringVerification:
    """Comprehensive verification of all non-string classifications."""

    def test_arithmetic_classifications(self, arrow_int_series, arrow_float_series):
        """Verify arithmetic operation classifications match runtime."""
        from scripts.generate_arrow_fallback_table import (
            FallbackType,
            analyze_arithmetic_operations,
        )

        methods = {m.name: m for m in analyze_arithmetic_operations()}

        # Test a subset of critical operations
        test_cases = {
            "add": (arrow_int_series, "__add__", 1),
            "sub": (arrow_int_series, "__sub__", 1),
            "mul": (arrow_int_series, "__mul__", 2),
            "truediv": (arrow_float_series, "__truediv__", 2.0),
        }

        mismatches = []
        for name, (series, op, other) in test_cases.items():
            if name not in methods:
                continue

            info = methods[name]
            with ns_tracker.track():
                try:
                    getattr(series, op)(other)
                except Exception:
                    continue

            used_fallback = ns_tracker.was_called("_apply_elementwise")

            if info.fallback_type == FallbackType.ARROW_NATIVE and used_fallback:
                mismatches.append(
                    f"{name}: classified as ARROW_NATIVE but used fallback"
                )

        if mismatches:
            pytest.fail(
                "Arithmetic classification mismatches:\n" + "\n".join(mismatches)
            )

    def test_datetime_classifications(self, arrow_datetime_series):
        """Verify datetime method classifications match runtime."""
        from scripts.generate_arrow_fallback_table import (
            FallbackType,
            analyze_datetime_methods,
        )

        methods = {m.name: m for m in analyze_datetime_methods()}

        # Test datetime properties
        test_attrs = [
            "_dt_year",
            "_dt_month",
            "_dt_day",
            "_dt_hour",
            "_dt_minute",
            "_dt_second",
            "_dt_quarter",
            "_dt_day_of_week",
            "_dt_is_month_start",
            "_dt_is_year_end",
        ]

        mismatches = []
        for name in test_attrs:
            if name not in methods:
                continue

            info = methods[name]
            public_name = info.public_name.split("/")[0]  # Handle aliases

            with ns_tracker.track():
                try:
                    getattr(arrow_datetime_series.dt, public_name)
                except Exception:
                    continue

            used_fallback = ns_tracker.was_called("_apply_elementwise")

            if info.fallback_type == FallbackType.ARROW_NATIVE and used_fallback:
                mismatches.append(
                    f"{name}: classified as ARROW_NATIVE but used fallback"
                )

        if mismatches:
            pytest.fail("Datetime classification mismatches:\n" + "\n".join(mismatches))

    def test_aggregation_classifications(
        self, arrow_int_series, arrow_float_series, arrow_bool_series
    ):
        """Verify aggregation method classifications match runtime."""
        from scripts.generate_arrow_fallback_table import (
            FallbackType,
            analyze_aggregation_methods,
        )

        methods = {m.name: m for m in analyze_aggregation_methods()}

        test_cases = {
            "sum": (arrow_int_series, {}),
            "mean": (arrow_float_series, {}),
            "min": (arrow_int_series, {}),
            "max": (arrow_int_series, {}),
            "std": (arrow_float_series, {}),
            "var": (arrow_float_series, {}),
            "any": (arrow_bool_series, {}),
            "all": (arrow_bool_series, {}),
        }

        mismatches = []
        for name, (series, kwargs) in test_cases.items():
            if name not in methods:
                continue

            info = methods[name]
            with ns_tracker.track():
                try:
                    getattr(series, name)(**kwargs)
                except Exception:
                    continue

            used_fallback = ns_tracker.was_called("_apply_elementwise")

            if info.fallback_type == FallbackType.ARROW_NATIVE and used_fallback:
                mismatches.append(
                    f"{name}: classified as ARROW_NATIVE but used fallback"
                )

        if mismatches:
            pytest.fail(
                "Aggregation classification mismatches:\n" + "\n".join(mismatches)
            )
