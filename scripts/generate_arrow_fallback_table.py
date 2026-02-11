"""
Generate Arrow fallback documentation via runtime observation.

This script runs all operations on all Arrow-backed dtypes and observes:
1. Whether the operation succeeds or raises an error
2. Whether the result dtype is Arrow-backed or falls back to NumPy
3. Whether to_numpy() or _apply_elementwise() was called (instrumented)

This is more accurate than AST analysis because it observes actual behavior.

Usage:
    # Generate RST documentation
    python scripts/generate_arrow_fallback_table.py

    # Generate detailed JSON (grouped by category)
    python scripts/generate_arrow_fallback_table.py --format json

    # Generate JSON lookup table for downstream libs: method -> dtype -> status
    python scripts/generate_arrow_fallback_table.py --format json-lookup

    # Generate JSON by dtype for downstream libs: dtype -> method -> status
    python scripts/generate_arrow_fallback_table.py --format json-by-dtype

    # Check if docs are up to date
    python scripts/generate_arrow_fallback_table.py --check
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
import sys
from typing import Any

import numpy as np

import pandas as pd

# =============================================================================
# Configuration: All Arrow dtypes to test
# =============================================================================

ARROW_DTYPES: dict[str, str] = {
    # String types
    "string": "string[pyarrow]",
    "large_string": "large_string[pyarrow]",
    # Integer types
    "int8": "int8[pyarrow]",
    "int16": "int16[pyarrow]",
    "int32": "int32[pyarrow]",
    "int64": "int64[pyarrow]",
    "uint8": "uint8[pyarrow]",
    "uint16": "uint16[pyarrow]",
    "uint32": "uint32[pyarrow]",
    "uint64": "uint64[pyarrow]",
    # Float types
    "float32": "float[pyarrow]",
    "float64": "double[pyarrow]",
    # Boolean
    "bool": "bool[pyarrow]",
    # Temporal types
    "timestamp_us": "timestamp[us][pyarrow]",
    "timestamp_ns": "timestamp[ns][pyarrow]",
    "timestamp_us_tz": "timestamp[us, UTC][pyarrow]",
    "date32": "date32[pyarrow]",
    "date64": "date64[pyarrow]",
    "duration_us": "duration[us][pyarrow]",
    "duration_ns": "duration[ns][pyarrow]",
    "time64_us": "time64[us][pyarrow]",
    # Binary
    "binary": "binary[pyarrow]",
    "large_binary": "large_binary[pyarrow]",
}

# Simplified dtype groups for documentation
DTYPE_GROUPS: dict[str, list[str]] = {
    "string": ["string", "large_string"],
    "integer": [
        "int8",
        "int16",
        "int32",
        "int64",
        "uint8",
        "uint16",
        "uint32",
        "uint64",
    ],
    "float": ["float32", "float64"],
    "bool": ["bool"],
    "timestamp": ["timestamp_us", "timestamp_ns", "timestamp_us_tz"],
    "date": ["date32", "date64"],
    "duration": ["duration_us", "duration_ns"],
    "time": ["time64_us"],
    "binary": ["binary", "large_binary"],
}


def validate_dtypes() -> set[str]:
    """
    Validate that all configured Arrow dtypes are supported.

    Returns a set of dtype names that failed validation (should be skipped).
    """
    import datetime
    import warnings

    # Test values appropriate for each dtype category
    test_values: dict[str, list] = {
        "string": ["a"],
        "large_string": ["a"],
        "binary": [b"a"],
        "large_binary": [b"a"],
        "bool": [True],
    }

    invalid_dtypes: set[str] = set()
    for dtype_name, dtype_str in ARROW_DTYPES.items():
        # Get appropriate test value
        if dtype_name in test_values:
            data = test_values[dtype_name]
        elif "timestamp" in dtype_name or "duration" in dtype_name:
            data = [1]  # integers work for these
        elif "date" in dtype_name:
            data = [datetime.date(2024, 1, 1)]
        elif dtype_name.startswith("time"):  # time64, not timestamp
            data = [datetime.time(10, 0, 0)]
        else:
            data = [1]  # numeric default

        try:
            pd.Series(data, dtype=dtype_str)
        except Exception as e:
            warnings.warn(
                f"Skipping dtype '{dtype_name}' ({dtype_str}): {e}",
                stacklevel=2,
            )
            invalid_dtypes.add(dtype_name)
    return invalid_dtypes


# =============================================================================
# Configuration: Methods to test with their required arguments
# =============================================================================

# String methods (Series.str.*)
STRING_METHODS: dict[str, dict[str, Any]] = {
    # No arguments needed
    "len": {},
    "lower": {},
    "upper": {},
    "capitalize": {},
    "title": {},
    "swapcase": {},
    "casefold": {},
    "isalnum": {},
    "isalpha": {},
    "isascii": {},
    "isdecimal": {},
    "isdigit": {},
    "islower": {},
    "isnumeric": {},
    "isspace": {},
    "istitle": {},
    "isupper": {},
    # With arguments
    "strip": {},
    "lstrip": {},
    "rstrip": {},
    "startswith": {"pat": "a"},
    "endswith": {"pat": "a"},
    "contains": {"pat": "a"},
    "match": {"pat": "a"},
    "fullmatch": {"pat": "a"},
    "find": {"sub": "a"},
    "rfind": {"sub": "a"},
    "index": {"sub": "a"},
    "rindex": {"sub": "a"},
    "count": {"pat": "a"},
    "replace": {"pat": "a", "repl": "b"},
    "repeat": {"repeats": 2},
    "pad": {"width": 10},
    "center": {"width": 10},
    "ljust": {"width": 10},
    "rjust": {"width": 10},
    "zfill": {"width": 10},
    "slice": {"start": 0, "stop": 2},
    "slice_replace": {"start": 0, "stop": 1, "repl": "X"},
    "get": {"i": 0},
    "join": {"sep": "-"},
    "split": {"pat": " "},
    "rsplit": {"pat": " "},
    "partition": {"sep": " "},
    "rpartition": {"sep": " "},
    "removeprefix": {"prefix": "a"},
    "removesuffix": {"suffix": "a"},
    "wrap": {"width": 10},
    "translate": {"table": str.maketrans("a", "b")},
    "encode": {"encoding": "utf-8"},
    "normalize": {"form": "NFC"},
    "findall": {"pat": r"\w+"},
    "extract": {"pat": r"(\w)(\w)"},
    "extractall": {"pat": r"(\w)"},
    "get_dummies": {"sep": "|"},
}

# Methods with parameter-dependent fallback behavior
# These are tested with multiple parameter combinations
STRING_PARAMETER_VARIANTS: dict[str, list[tuple[str, dict[str, Any]]]] = {
    # split: expand=True returns DataFrame (Arrow), expand=False returns list Series
    "split": [
        ("split(expand=False)", {"pat": " ", "expand": False}),
        ("split(expand=True)", {"pat": " ", "expand": True}),
    ],
    "rsplit": [
        ("rsplit(expand=False)", {"pat": " ", "expand": False}),
        ("rsplit(expand=True)", {"pat": " ", "expand": True}),
    ],
    # extract: expand=True returns DataFrame, expand=False returns Series
    # ArrowExtensionArray requires named groups, so use them for consistent testing
    "extract": [
        ("extract(expand=False)", {"pat": r"(?P<first>\w)", "expand": False}),
        (
            "extract(expand=True)",
            {"pat": r"(?P<first>\w)(?P<second>\w)", "expand": True},
        ),
    ],
    # replace: case=False falls back to object dtype, callable repl also falls back
    "replace": [
        ("replace(case=True)", {"pat": "a", "repl": "b", "case": True}),
        ("replace(case=False)", {"pat": "a", "repl": "b", "case": False}),
        # callable repl requires regex=True
        ("replace(repl=callable)", {"pat": "a", "repl": lambda _: "b", "regex": True}),
    ],
    # contains: flags!=0 falls back to object dtype
    "contains": [
        ("contains(flags=0)", {"pat": "a", "flags": 0}),
        ("contains(flags=re.I)", {"pat": "a", "flags": 2}),  # re.IGNORECASE = 2
    ],
    # match: flags!=0 falls back to object dtype
    "match": [
        ("match(flags=0)", {"pat": "a", "flags": 0}),
        ("match(flags=re.I)", {"pat": "a", "flags": 2}),
    ],
    # fullmatch: flags!=0 falls back to object dtype
    "fullmatch": [
        ("fullmatch(flags=0)", {"pat": ".*", "flags": 0}),
        ("fullmatch(flags=re.I)", {"pat": ".*", "flags": 2}),
    ],
    # count: flags!=0 falls back to object dtype
    "count": [
        ("count(flags=0)", {"pat": "a", "flags": 0}),
        ("count(flags=re.I)", {"pat": "a", "flags": 2}),
    ],
    # findall: flags parameter changes behavior
    "findall": [
        ("findall(flags=0)", {"pat": r"\w+", "flags": 0}),
        ("findall(flags=re.I)", {"pat": r"\w+", "flags": 2}),
    ],
    # pad: side='both' needs PyArrow 17+ for utf8_center
    "pad": [
        ("pad(side=left)", {"width": 10, "side": "left"}),
        ("pad(side=right)", {"width": 10, "side": "right"}),
        ("pad(side=both)", {"width": 10, "side": "both"}),
    ],
}

# Datetime methods (Series.dt.*)
DATETIME_METHODS: dict[str, dict[str, Any]] = {
    # Properties (no args)
    "year": {},
    "month": {},
    "day": {},
    "hour": {},
    "minute": {},
    "second": {},
    "microsecond": {},
    "nanosecond": {},
    "dayofweek": {},
    "day_of_week": {},
    "weekday": {},
    "dayofyear": {},
    "day_of_year": {},
    "quarter": {},
    "is_month_start": {},
    "is_month_end": {},
    "is_quarter_start": {},
    "is_quarter_end": {},
    "is_year_start": {},
    "is_year_end": {},
    "is_leap_year": {},
    "days_in_month": {},
    "daysinmonth": {},
    "date": {},
    "time": {},
    "tz": {},
    "unit": {},
    # Methods
    "normalize": {},
    "strftime": {"date_format": "%Y-%m-%d"},
    "round": {"freq": "h"},
    "floor": {"freq": "h"},
    "ceil": {"freq": "h"},
    "day_name": {},
    "month_name": {},
    "tz_localize": {"tz": "UTC"},
    "tz_convert": {"tz": "US/Eastern"},  # needs tz-aware input
    "isocalendar": {},
    "to_pydatetime": {},
    "as_unit": {"unit": "ms"},
}

# Timedelta-specific methods (Series.dt.* for duration)
TIMEDELTA_METHODS: dict[str, dict[str, Any]] = {
    "days": {},
    "seconds": {},
    "microseconds": {},
    "nanoseconds": {},
    "total_seconds": {},
    "to_pytimedelta": {},
    "as_unit": {"unit": "ms"},
}

# Aggregation methods
AGGREGATION_METHODS: dict[str, dict[str, Any]] = {
    "sum": {},
    "mean": {},
    "median": {},
    "min": {},
    "max": {},
    "std": {},
    "var": {},
    "sem": {},
    "prod": {},
    "any": {},
    "all": {},
    "count": {},
    "skew": {},
    "kurt": {},
}

# Array methods (direct Series methods)
ARRAY_METHODS: dict[str, dict[str, Any]] = {
    "unique": {},
    "dropna": {},
    "fillna": {"value": 0},  # value depends on dtype
    "isna": {},
    "notna": {},
    "argsort": {},
    "sort_values": {},
    "duplicated": {},
    "drop_duplicates": {},
    "searchsorted": {"value": 0},  # value depends on dtype
    "round": {"decimals": 2},
    "diff": {},
    "shift": {"periods": 1},
    "cumsum": {},
    "cumprod": {},
    "cummin": {},
    "cummax": {},
    "ffill": {},
    "bfill": {},
    "abs": {},
    "clip": {"lower": 0, "upper": 10},
}

# Algorithmic methods - separated due to behavior that may vary across versions
# These have explicit kwargs for deterministic results
ALGORITHMIC_METHODS: dict[str, dict[str, Any]] = {
    "value_counts": {"dropna": True, "sort": True},
    "factorize": {"sort": False, "use_na_sentinel": True},
    "rank": {"method": "average", "na_option": "keep"},
}

# Array methods with parameter-dependent fallback behavior
ARRAY_PARAMETER_VARIANTS: dict[str, list[tuple[str, dict[str, Any]]]] = {
    # fillna: limit=None uses Arrow, limit=N falls back
    "fillna": [
        ("fillna(limit=None)", {"value": 0}),  # value adjusted per dtype later
        ("fillna(limit=1)", {"value": 0, "limit": 1}),
    ],
    # interpolate: method="linear" with defaults uses Arrow, others fall back
    "interpolate": [
        ("interpolate(method=linear)", {"method": "linear"}),
        ("interpolate(method=pad)", {"method": "pad"}),
    ],
}

# Datetime methods with parameter-dependent fallback behavior
DATETIME_PARAMETER_VARIANTS: dict[str, list[tuple[str, dict[str, Any]]]] = {
    # tz_localize: ambiguous="raise" is required, others raise NotImplementedError
    "tz_localize": [
        ("tz_localize(ambiguous=raise)", {"tz": "UTC", "ambiguous": "raise"}),
        ("tz_localize(ambiguous=NaT)", {"tz": "UTC", "ambiguous": "NaT"}),
    ],
}

# =============================================================================
# Version-gated methods: behavior changes based on PyArrow version
# =============================================================================

# Maps method -> {accessor, min_version, description, below_behavior,
#                 at_or_above_behavior}
# These methods have different fallback behavior depending on PyArrow version
VERSION_GATED_METHODS: dict[str, dict[str, Any]] = {
    # _str_pad with side="both" needs PyArrow 17.0 for utf8_center
    "center": {
        "accessor": "str",
        "min_version": "17.0",
        "description": "Uses pc.utf8_center for side='both' padding",
        "below_behavior": "object",  # Falls back to object dtype
        "at_or_above_behavior": "arrow",  # Native Arrow
    },
    "pad(side=both)": {
        "accessor": "str",
        "min_version": "17.0",
        "description": "Uses pc.utf8_center for side='both' padding",
        "below_behavior": "object",
        "at_or_above_behavior": "arrow",
    },
    # _str_isdigit needs PyArrow 21.0 for correct utf8_is_digit behavior
    "isdigit": {
        "accessor": "str",
        "min_version": "21.0",
        "description": "PyArrow < 21.0 has incorrect utf8_is_digit for some digits",
        "below_behavior": "elementwise",  # Uses _apply_elementwise
        "at_or_above_behavior": "arrow",
    },
    # _str_zfill needs PyArrow 21.0 for utf8_zfill
    "zfill": {
        "accessor": "str",
        "min_version": "21.0",
        "description": "pc.utf8_zfill was added in PyArrow 21.0",
        "below_behavior": "elementwise",
        "at_or_above_behavior": "arrow",
    },
}

# Arithmetic operations
ARITHMETIC_OPS: dict[str, tuple[str, Any]] = {
    "add": ("__add__", 1),
    "sub": ("__sub__", 1),
    "mul": ("__mul__", 2),
    "truediv": ("__truediv__", 2),
    "floordiv": ("__floordiv__", 2),
    "mod": ("__mod__", 2),
    "pow": ("__pow__", 2),
}

# Comparison operations
COMPARISON_OPS: dict[str, tuple[str, Any]] = {
    "eq": ("__eq__", 1),
    "ne": ("__ne__", 1),
    "lt": ("__lt__", 1),
    "le": ("__le__", 1),
    "gt": ("__gt__", 1),
    "ge": ("__ge__", 1),
}


# =============================================================================
# Result classification
# =============================================================================


class ResultType(Enum):
    """Classification of operation result."""

    ARROW_NATIVE = "arrow"  # Result is Arrow-backed, no fallback detected
    NUMPY_FALLBACK = "numpy"  # Result fell back to NumPy dtype
    ELEMENTWISE = "elementwise"  # Used _apply_elementwise
    OBJECT_FALLBACK = "object"  # Converted to object dtype
    NOT_IMPLEMENTED = "not_implemented"  # Raises NotImplementedError
    TYPE_ERROR = "type_error"  # Raises TypeError (not supported for dtype)
    OTHER_ERROR = "error"  # Other exception


@dataclass
class OperationResult:
    """Result of running an operation."""

    method: str
    dtype: str
    result_type: ResultType
    result_dtype: str | None = None
    error_message: str | None = None
    used_to_numpy: bool = False
    used_elementwise: bool = False


# =============================================================================
# Instrumentation for detecting fallback paths
# =============================================================================


class FallbackTracker:
    """Track fallback calls during operation execution."""

    def __init__(self):
        self.to_numpy_called = False
        self.elementwise_called = False

    def reset(self):
        self.to_numpy_called = False
        self.elementwise_called = False


@contextmanager
def track_fallbacks():
    """Context manager to track fallback calls."""
    tracker = FallbackTracker()

    # Try to import ArrowExtensionArray for instrumentation
    # If pandas layout changes or Arrow is unavailable, skip instrumentation
    try:
        from pandas.core.arrays.arrow import ArrowExtensionArray
    except (ImportError, AttributeError):
        # Can't instrument, just yield tracker without patching
        yield tracker
        return

    original_to_numpy = ArrowExtensionArray.to_numpy
    original_elementwise = getattr(ArrowExtensionArray, "_apply_elementwise", None)

    # If _apply_elementwise doesn't exist, skip that instrumentation
    if original_elementwise is None:
        yield tracker
        return

    def patched_to_numpy(self, *args, **kwargs):
        tracker.to_numpy_called = True
        return original_to_numpy(self, *args, **kwargs)

    def patched_elementwise(self, *args, **kwargs):
        tracker.elementwise_called = True
        return original_elementwise(self, *args, **kwargs)

    # Manual monkey-patching (avoiding unittest.mock which is banned)
    ArrowExtensionArray.to_numpy = patched_to_numpy
    ArrowExtensionArray._apply_elementwise = patched_elementwise
    try:
        yield tracker
    finally:
        ArrowExtensionArray.to_numpy = original_to_numpy
        ArrowExtensionArray._apply_elementwise = original_elementwise


# =============================================================================
# Test data creation
# =============================================================================


def create_test_series(dtype_name: str) -> pd.Series | None:
    """Create a test Series for the given dtype."""
    dtype_str = ARROW_DTYPES.get(dtype_name)
    if dtype_str is None:
        return None

    try:
        if dtype_name in ["string", "large_string"]:
            data = ["hello world", "test string", "abc def", None, "xyz"]
        elif dtype_name in ["binary", "large_binary"]:
            data = [b"hello", b"world", b"test", None, b"xyz"]
        elif dtype_name == "bool":
            data = [True, False, True, None, False]
        elif "int" in dtype_name or "uint" in dtype_name:
            data = [1, 2, 3, None, 5]
        elif "float" in dtype_name:
            data = [1.5, 2.5, 3.5, None, 5.5]
        elif "timestamp" in dtype_name:
            data = pd.to_datetime(
                ["2024-01-01", "2024-01-02", "2024-01-03", None, "2024-01-05"]
            )
            if "tz" in dtype_name or "UTC" in dtype_str:
                data = data.tz_localize("UTC")
        elif "date" in dtype_name:
            import datetime

            data = [
                datetime.date(2024, 1, 1),
                datetime.date(2024, 1, 2),
                datetime.date(2024, 1, 3),
                None,
                datetime.date(2024, 1, 5),
            ]
        elif "duration" in dtype_name:
            data = pd.to_timedelta(["1 day", "2 days", "3 days", None, "5 days"])
        elif "time" in dtype_name:
            import datetime

            data = [
                datetime.time(10, 0, 0),
                datetime.time(11, 30, 0),
                datetime.time(12, 45, 30),
                None,
                datetime.time(14, 0, 0),
            ]
        else:
            data = [1, 2, 3, None, 5]

        return pd.Series(data, dtype=dtype_str)
    except Exception:
        return None


def create_test_series_no_na(dtype_name: str) -> pd.Series | None:
    """Create a test Series without NA values for methods that require sorted data."""
    dtype_str = ARROW_DTYPES.get(dtype_name)
    if dtype_str is None:
        return None

    try:
        if dtype_name in ["string", "large_string"]:
            data = ["apple", "banana", "cherry", "date", "elderberry"]
        elif dtype_name in ["binary", "large_binary"]:
            data = [b"apple", b"banana", b"cherry", b"date", b"elderberry"]
        elif dtype_name == "bool":
            data = [False, False, True, True, True]
        elif "int" in dtype_name or "uint" in dtype_name:
            data = [1, 2, 3, 4, 5]
        elif "float" in dtype_name:
            data = [1.5, 2.5, 3.5, 4.5, 5.5]
        elif "timestamp" in dtype_name:
            data = pd.to_datetime(
                ["2024-01-01", "2024-01-02", "2024-01-03", "2024-01-04", "2024-01-05"]
            )
            if "tz" in dtype_name or "UTC" in dtype_str:
                data = data.tz_localize("UTC")
        elif "date" in dtype_name:
            import datetime

            data = [
                datetime.date(2024, 1, 1),
                datetime.date(2024, 1, 2),
                datetime.date(2024, 1, 3),
                datetime.date(2024, 1, 4),
                datetime.date(2024, 1, 5),
            ]
        elif "duration" in dtype_name:
            data = pd.to_timedelta(["1 day", "2 days", "3 days", "4 days", "5 days"])
        elif "time" in dtype_name:
            import datetime

            data = [
                datetime.time(10, 0, 0),
                datetime.time(11, 0, 0),
                datetime.time(12, 0, 0),
                datetime.time(13, 0, 0),
                datetime.time(14, 0, 0),
            ]
        else:
            data = [1, 2, 3, 4, 5]

        return pd.Series(data, dtype=dtype_str)
    except Exception:
        return None


def get_fillna_value(dtype_name: str) -> Any:
    """Get appropriate fillna value for dtype."""
    if dtype_name in ["string", "large_string"]:
        return "filled"
    elif dtype_name in ["binary", "large_binary"]:
        return b"filled"
    elif dtype_name == "bool":
        return False
    elif "timestamp" in dtype_name or "date" in dtype_name:
        return pd.Timestamp("2024-01-01")
    elif "duration" in dtype_name:
        return pd.Timedelta("1 day")
    else:
        return 0


def get_searchsorted_value(dtype_name: str) -> Any:
    """Get appropriate searchsorted value for dtype."""
    if dtype_name in ["string", "large_string"]:
        return "m"
    elif dtype_name in ["binary", "large_binary"]:
        return b"m"
    elif dtype_name == "bool":
        return True
    elif "timestamp" in dtype_name or "date" in dtype_name:
        return pd.Timestamp("2024-01-02")
    elif "duration" in dtype_name:
        return pd.Timedelta("2 days")
    else:
        return 2


# =============================================================================
# Run operations and classify results
# =============================================================================


def _is_arrow_backed(dtype) -> bool:
    """Check if a dtype is Arrow-backed."""
    # ArrowDtype has pyarrow_dtype attribute
    if hasattr(dtype, "pyarrow_dtype"):
        return True
    # StringDtype with pyarrow storage
    if hasattr(dtype, "storage") and dtype.storage == "pyarrow":
        return True
    # Check dtype name for [pyarrow] suffix
    dtype_str = str(dtype)
    if "[pyarrow]" in dtype_str:
        return True
    return False


def classify_result(
    series: pd.Series,
    result: Any,
    tracker: FallbackTracker,
) -> ResultType:
    """Classify the result of an operation."""
    # Check if elementwise was used (definitive fallback indicator)
    if tracker.elementwise_called:
        return ResultType.ELEMENTWISE

    # For Series/DataFrame results, prioritize checking the result dtype
    # This avoids false positives from incidental to_numpy calls in pandas internals
    # (e.g., DataFrame construction may call to_numpy on Index objects)
    if isinstance(result, (pd.Series, pd.DataFrame)):
        if isinstance(result, pd.DataFrame):
            # For DataFrames, check if all columns are Arrow-backed
            has_arrow = any(
                _is_arrow_backed(result[col].dtype) for col in result.columns
            )
        else:
            has_arrow = _is_arrow_backed(result.dtype)

        if has_arrow:
            return ResultType.ARROW_NATIVE
        elif isinstance(result, pd.Series) and result.dtype == object:
            return ResultType.OBJECT_FALLBACK
        else:
            # Result is numpy-backed; to_numpy may or may not have been called
            return ResultType.NUMPY_FALLBACK

    # Handle tuple results (e.g., factorize returns (codes, uniques))
    # If any element is a numpy array, classify as numpy fallback
    # This is because the numpy array indicates data was converted from Arrow
    if isinstance(result, tuple):
        for item in result:
            if isinstance(item, np.ndarray):
                return ResultType.NUMPY_FALLBACK
        return ResultType.ARROW_NATIVE

    # For scalar/array results, use instrumentation as indicator
    if tracker.to_numpy_called:
        return ResultType.NUMPY_FALLBACK

    return ResultType.ARROW_NATIVE


def run_operation(
    series: pd.Series,
    method_name: str,
    kwargs: dict[str, Any],
    accessor: str | None = None,
    display_name: str | None = None,
) -> OperationResult:
    """Run an operation and return the result classification."""
    dtype_name = str(series.dtype)
    # Use display_name for results if provided (for parameter variants)
    result_name = display_name if display_name else method_name

    with track_fallbacks() as tracker:
        try:
            # Get the method
            if accessor:
                obj = getattr(series, accessor)
            else:
                obj = series

            method = getattr(obj, method_name)

            # Handle special cases for kwargs
            actual_kwargs = kwargs.copy()

            # Call the method
            if callable(method):
                result = method(**actual_kwargs)
            else:
                # It's a property
                result = method

            # Classify the result
            result_type = classify_result(series, result, tracker)

            # Get result dtype
            if isinstance(result, pd.Series):
                result_dtype = str(result.dtype)
            elif isinstance(result, pd.DataFrame):
                result_dtype = "DataFrame"
            elif isinstance(result, np.ndarray):
                result_dtype = str(result.dtype)
            else:
                result_dtype = type(result).__name__

            return OperationResult(
                method=result_name,
                dtype=dtype_name,
                result_type=result_type,
                result_dtype=result_dtype,
                used_to_numpy=tracker.to_numpy_called,
                used_elementwise=tracker.elementwise_called,
            )

        except NotImplementedError as e:
            return OperationResult(
                method=result_name,
                dtype=dtype_name,
                result_type=ResultType.NOT_IMPLEMENTED,
                error_message=str(e),
            )
        except TypeError as e:
            return OperationResult(
                method=result_name,
                dtype=dtype_name,
                result_type=ResultType.TYPE_ERROR,
                error_message=str(e),
            )
        except Exception as e:
            return OperationResult(
                method=result_name,
                dtype=dtype_name,
                result_type=ResultType.OTHER_ERROR,
                error_message=f"{type(e).__name__}: {e}",
            )


def run_arithmetic_op(
    series: pd.Series,
    op_name: str,
    dunder: str,
    operand: Any,
) -> OperationResult:
    """Run an arithmetic operation."""
    dtype_name = str(series.dtype)

    with track_fallbacks() as tracker:
        try:
            method = getattr(series, dunder)
            result = method(operand)

            result_type = classify_result(series, result, tracker)

            return OperationResult(
                method=op_name,
                dtype=dtype_name,
                result_type=result_type,
                result_dtype=str(result.dtype)
                if isinstance(result, pd.Series)
                else type(result).__name__,
                used_to_numpy=tracker.to_numpy_called,
                used_elementwise=tracker.elementwise_called,
            )
        except NotImplementedError as e:
            return OperationResult(
                method=op_name,
                dtype=dtype_name,
                result_type=ResultType.NOT_IMPLEMENTED,
                error_message=str(e),
            )
        except TypeError as e:
            return OperationResult(
                method=op_name,
                dtype=dtype_name,
                result_type=ResultType.TYPE_ERROR,
                error_message=str(e),
            )
        except Exception as e:
            return OperationResult(
                method=op_name,
                dtype=dtype_name,
                result_type=ResultType.OTHER_ERROR,
                error_message=f"{type(e).__name__}: {e}",
            )


# =============================================================================
# Run all tests
# =============================================================================


def run_all_tests() -> dict[str, list[OperationResult]]:
    """Run all operations on all dtypes and collect results."""
    # Validate dtypes upfront and warn about any that won't work
    invalid_dtypes = validate_dtypes()
    if invalid_dtypes:
        import warnings

        warnings.warn(
            f"Skipping {len(invalid_dtypes)} unsupported dtypes: {invalid_dtypes}",
            stacklevel=2,
        )

    results: dict[str, list[OperationResult]] = {
        "string_methods": [],
        "datetime_methods": [],
        "timedelta_methods": [],
        "aggregations": [],
        "array_methods": [],
        "algorithmic_methods": [],
        "arithmetic": [],
        "comparison": [],
    }

    # String methods - only on string dtypes
    for dtype_name in DTYPE_GROUPS["string"]:
        series = create_test_series(dtype_name)
        if series is None:
            continue

        for method, kwargs in STRING_METHODS.items():
            # Skip methods that have parameter variants (tested separately below)
            if method in STRING_PARAMETER_VARIANTS:
                continue
            result = run_operation(series, method, kwargs, accessor="str")
            results["string_methods"].append(result)

        # Test parameter variants
        for method, variants in STRING_PARAMETER_VARIANTS.items():
            for display_name, kwargs in variants:
                result = run_operation(
                    series, method, kwargs, accessor="str", display_name=display_name
                )
                results["string_methods"].append(result)

    # Datetime methods - on timestamp dtypes
    for dtype_name in DTYPE_GROUPS["timestamp"]:
        series = create_test_series(dtype_name)
        if series is None:
            continue

        for method, kwargs in DATETIME_METHODS.items():
            # Skip methods that have parameter variants (tested separately below)
            if method in DATETIME_PARAMETER_VARIANTS:
                continue
            # Skip tz_convert for non-tz-aware
            if (
                method == "tz_convert"
                and "tz" not in dtype_name
                and "UTC" not in ARROW_DTYPES[dtype_name]
            ):
                continue
            # Skip tz_localize for tz-aware
            if method == "tz_localize" and (
                "tz" in dtype_name or "UTC" in ARROW_DTYPES[dtype_name]
            ):
                continue

            result = run_operation(series, method, kwargs, accessor="dt")
            results["datetime_methods"].append(result)

        # Test datetime parameter variants (only on non-tz-aware for tz_localize)
        for method, variants in DATETIME_PARAMETER_VARIANTS.items():
            # Skip tz_localize for tz-aware
            if method == "tz_localize" and (
                "tz" in dtype_name or "UTC" in ARROW_DTYPES[dtype_name]
            ):
                continue
            for display_name, kwargs in variants:
                result = run_operation(
                    series, method, kwargs, accessor="dt", display_name=display_name
                )
                results["datetime_methods"].append(result)

    # Timedelta methods - on duration dtypes
    for dtype_name in DTYPE_GROUPS["duration"]:
        series = create_test_series(dtype_name)
        if series is None:
            continue

        for method, kwargs in TIMEDELTA_METHODS.items():
            result = run_operation(series, method, kwargs, accessor="dt")
            results["timedelta_methods"].append(result)

    # Aggregations - on numeric dtypes
    for dtype_name in DTYPE_GROUPS["integer"] + DTYPE_GROUPS["float"]:
        series = create_test_series(dtype_name)
        if series is None:
            continue

        for method, kwargs in AGGREGATION_METHODS.items():
            result = run_operation(series, method, kwargs)
            results["aggregations"].append(result)

    # Array methods - on all dtypes
    for dtype_name, dtype_str in ARROW_DTYPES.items():
        series = create_test_series(dtype_name)
        if series is None:
            continue

        # Define dtype categories for skip rules
        non_numeric = [
            "string",
            "large_string",
            "binary",
            "large_binary",
            "bool",
        ]
        unsigned = ["uint8", "uint16", "uint32", "uint64"]

        for method, kwargs in ARRAY_METHODS.items():
            # Skip methods that have parameter variants (tested separately below)
            if method in ARRAY_PARAMETER_VARIANTS:
                continue

            # Adjust kwargs for dtype
            actual_kwargs = kwargs.copy()

            # Skip rules for incompatible dtype/method combinations
            if method == "fillna":
                actual_kwargs["value"] = get_fillna_value(dtype_name)
            elif method == "searchsorted":
                # searchsorted requires sorted data without NAs
                series_no_na = create_test_series_no_na(dtype_name)
                if series_no_na is None:
                    continue
                actual_kwargs["value"] = get_searchsorted_value(dtype_name)
                result = run_operation(series_no_na, method, actual_kwargs)
                results["array_methods"].append(result)
                continue  # Skip the normal run_operation below
            elif method == "clip":
                if dtype_name in non_numeric:
                    continue  # clip doesn't make sense for these
                if "timestamp" in dtype_name or "date" in dtype_name:
                    actual_kwargs = {
                        "lower": pd.Timestamp("2024-01-02"),
                        "upper": pd.Timestamp("2024-01-04"),
                    }
                elif "duration" in dtype_name:
                    actual_kwargs = {
                        "lower": pd.Timedelta("2 days"),
                        "upper": pd.Timedelta("4 days"),
                    }
            elif method == "round":
                if dtype_name not in DTYPE_GROUPS["float"]:
                    continue  # round only for float
            elif method == "diff":
                if dtype_name in ["string", "large_string", "binary", "large_binary"]:
                    continue  # diff doesn't make sense for string/binary
            elif method in ("cumsum", "cumprod"):
                if dtype_name in ["string", "large_string", "binary", "large_binary"]:
                    continue  # cumulative ops don't make sense for string/binary
            elif method == "abs":
                # abs doesn't make sense for non-numeric or unsigned types
                if dtype_name in non_numeric or dtype_name in unsigned:
                    continue

            result = run_operation(series, method, actual_kwargs)
            results["array_methods"].append(result)

        # Test array parameter variants
        for method, variants in ARRAY_PARAMETER_VARIANTS.items():
            # Skip interpolate for non-numeric types (always errors)
            if method == "interpolate" and dtype_name in non_numeric:
                continue

            for display_name, kwargs in variants:
                actual_kwargs = kwargs.copy()
                # Adjust fillna value for dtype
                if method == "fillna" and "value" in actual_kwargs:
                    actual_kwargs["value"] = get_fillna_value(dtype_name)
                result = run_operation(
                    series, method, actual_kwargs, display_name=display_name
                )
                results["array_methods"].append(result)

    # Algorithmic methods - on all dtypes (may error on some)
    for dtype_name in ARROW_DTYPES:
        series = create_test_series(dtype_name)
        if series is None:
            continue

        for method, kwargs in ALGORITHMIC_METHODS.items():
            result = run_operation(series, method, kwargs)
            results["algorithmic_methods"].append(result)

    # Arithmetic - on numeric dtypes
    for dtype_name in DTYPE_GROUPS["integer"] + DTYPE_GROUPS["float"]:
        series = create_test_series(dtype_name)
        if series is None:
            continue

        for op_name, (dunder, operand) in ARITHMETIC_OPS.items():
            result = run_arithmetic_op(series, op_name, dunder, operand)
            results["arithmetic"].append(result)

    # Comparison - on numeric dtypes
    for dtype_name in DTYPE_GROUPS["integer"] + DTYPE_GROUPS["float"]:
        series = create_test_series(dtype_name)
        if series is None:
            continue

        for op_name, (dunder, operand) in COMPARISON_OPS.items():
            result = run_arithmetic_op(series, op_name, dunder, operand)
            results["comparison"].append(result)

    return results


# =============================================================================
# Format output
# =============================================================================


def format_rst_table(all_results: dict[str, list[OperationResult]]) -> str:
    """Format results as RST tables."""
    title = "Arrow Fallbacks"
    lines = [
        ".. _arrow-fallbacks:",
        "",
        "{{ header }}",
        "",
        "*" * len(title),
        title,
        "*" * len(title),
        "",
        "This document shows the runtime behavior of pandas methods on",
        "Arrow-backed arrays. Results are determined by actually running",
        "each operation and observing the outcome.",
        "",
        ".. note::",
        "",
        "   Results may vary based on method parameters. For example, some string",
        "   methods like ``split(expand=True)`` use Arrow, while",
        "   ``split(expand=False)`` may fall back to NumPy. This table shows",
        "   results for default or common parameter values only.",
        "",
        ".. note::",
        "",
        "   Some methods have version-gated behavior that depends on your PyArrow",
        "   version. See :ref:`arrow-version-gated` below for details.",
        "",
        "Legend",
        "======",
        "",
        ".. list-table::",
        "   :widths: 20 80",
        "   :header-rows: 1",
        "",
        "   * - Status",
        "     - Description",
        "   * - |arrow|",
        "     - Operation returns Arrow-backed result",
        "   * - |numpy|",
        "     - Operation falls back to NumPy",
        "   * - |elementwise|",
        "     - Uses element-wise Python processing",
        "   * - |object|",
        "     - Returns object dtype",
        "   * - |notimpl|",
        "     - Raises NotImplementedError",
        "   * - |typeerror|",
        "     - Not supported for this dtype",
        "   * - |error|",
        "     - Other error",
        "   * - |mixed|",
        "     - Mixed results within dtype group",
        "",
        ".. |arrow| replace:: ✓ Arrow",
        ".. |numpy| replace:: → NumPy",
        ".. |elementwise| replace:: ⟳ Elem",
        ".. |object| replace:: → Object",
        ".. |notimpl| replace:: ✗ N/I",
        ".. |typeerror| replace:: ✗ Type",
        ".. |error| replace:: ✗ Err",
        ".. |mixed| replace:: ~ Mixed",
        "",
    ]

    status_map = {
        ResultType.ARROW_NATIVE: "|arrow|",
        ResultType.NUMPY_FALLBACK: "|numpy|",
        ResultType.ELEMENTWISE: "|elementwise|",
        ResultType.OBJECT_FALLBACK: "|object|",
        ResultType.NOT_IMPLEMENTED: "|notimpl|",
        ResultType.TYPE_ERROR: "|typeerror|",
        ResultType.OTHER_ERROR: "|error|",
    }

    # String methods table
    if all_results["string_methods"]:
        lines.extend(
            [
                "",
                "String Methods (Series.str.*)",
                "=============================",
                "",
            ]
        )
        lines.extend(_format_method_table(all_results["string_methods"], status_map))

    # Datetime methods table
    if all_results["datetime_methods"]:
        lines.extend(
            [
                "",
                "Datetime Methods (Series.dt.*)",
                "==============================",
                "",
            ]
        )
        lines.extend(_format_method_table(all_results["datetime_methods"], status_map))

    # Timedelta methods table
    if all_results["timedelta_methods"]:
        lines.extend(
            [
                "",
                "Timedelta Methods (Series.dt.*)",
                "===============================",
                "",
            ]
        )
        lines.extend(_format_method_table(all_results["timedelta_methods"], status_map))

    # Aggregation methods table
    if all_results["aggregations"]:
        lines.extend(
            [
                "",
                "Aggregation Methods",
                "===================",
                "",
            ]
        )
        lines.extend(_format_method_table(all_results["aggregations"], status_map))

    # Array methods table
    if all_results["array_methods"]:
        lines.extend(
            [
                "",
                "Array Methods",
                "=============",
                "",
            ]
        )
        lines.extend(_format_method_table(all_results["array_methods"], status_map))

    # Algorithmic methods table
    if all_results["algorithmic_methods"]:
        lines.extend(
            [
                "",
                "Algorithmic Methods",
                "===================",
                "",
                "These methods have behavior that may vary across versions.",
                "Explicit kwargs are used for deterministic results.",
                "",
            ]
        )
        lines.extend(
            _format_method_table(all_results["algorithmic_methods"], status_map)
        )

    # Arithmetic operations table
    if all_results["arithmetic"]:
        lines.extend(
            [
                "",
                "Arithmetic Operations",
                "=====================",
                "",
            ]
        )
        lines.extend(_format_method_table(all_results["arithmetic"], status_map))

    # Comparison operations table
    if all_results["comparison"]:
        lines.extend(
            [
                "",
                "Comparison Operations",
                "=====================",
                "",
            ]
        )
        lines.extend(_format_method_table(all_results["comparison"], status_map))

    # Version-gated methods section
    lines.extend(_format_version_gated_section())

    return "\n".join(lines) + "\n"


def _format_version_gated_section() -> list[str]:
    """Format the version-gated methods section for RST."""
    lines = [
        "",
        ".. _arrow-version-gated:",
        "",
        "Version-Gated Methods",
        "=====================",
        "",
        "The following methods have behavior that depends on your installed PyArrow",
        "version. The table above shows results for the current environment; if you",
        "have an older PyArrow version, some methods may fall back to slower paths.",
        "",
        ".. list-table::",
        "   :widths: 25 15 60",
        "   :header-rows: 1",
        "",
        "   * - Method",
        "     - Min PyArrow",
        "     - Notes",
    ]

    for method_name, info in sorted(VERSION_GATED_METHODS.items()):
        accessor = info["accessor"]
        if accessor:
            method_key = f"{accessor}.{method_name}"
        else:
            method_key = method_name

        min_ver = info["min_version"]
        below = info["below_behavior"]
        above = info["at_or_above_behavior"]
        desc = info["description"]

        lines.append(f"   * - ``{method_key}``")
        lines.append(f"     - {min_ver}")
        lines.append(f"     - {desc}. Below {min_ver}: {below}; {min_ver}+: {above}.")

    return lines


def _format_method_table(
    results: list[OperationResult],
    status_map: dict[ResultType, str],
) -> list[str]:
    """Format a single method table with dtype groups for readability."""
    # Build reverse mapping: dtype_str -> group_name
    # Include both the ARROW_DTYPES values and the actual str(dtype) representations
    dtype_to_group: dict[str, str] = {}
    for group_name, dtype_keys in DTYPE_GROUPS.items():
        for dtype_key in dtype_keys:
            dtype_str = ARROW_DTYPES.get(dtype_key)
            if dtype_str:
                dtype_to_group[dtype_str] = group_name
                # Also add the actual str(dtype) representation which may differ
                # e.g., "string[pyarrow]" -> "string", "bool[pyarrow]" -> "bool"
                try:
                    actual_dtype_str = str(pd.Series([None], dtype=dtype_str).dtype)
                    if actual_dtype_str != dtype_str:
                        dtype_to_group[actual_dtype_str] = group_name
                except Exception:
                    pass

    # Group by method
    by_method: dict[str, dict[str, ResultType]] = {}
    groups_seen: set[str] = set()

    for r in results:
        if r.method not in by_method:
            by_method[r.method] = {}
        by_method[r.method][r.dtype] = r.result_type
        group = dtype_to_group.get(r.dtype, "other")
        groups_seen.add(group)

    # Define group order for consistent columns
    group_order = [
        "string",
        "integer",
        "float",
        "bool",
        "timestamp",
        "date",
        "duration",
        "time",
        "binary",
    ]
    # Only include groups that have data
    groups = [g for g in group_order if g in groups_seen]

    def get_group_status(method_results: dict[str, ResultType], group: str) -> str:
        """Get status for a dtype group. Returns unified status or 'mixed'."""
        dtype_keys = DTYPE_GROUPS.get(group, [])
        statuses: set[ResultType] = set()

        for dtype_key in dtype_keys:
            dtype_str = ARROW_DTYPES.get(dtype_key)
            if not dtype_str:
                continue
            # Check both the ARROW_DTYPES value and actual str(dtype) representation
            if dtype_str in method_results:
                statuses.add(method_results[dtype_str])
            else:
                # Try the actual dtype string representation
                try:
                    actual_str = str(pd.Series([None], dtype=dtype_str).dtype)
                    if actual_str in method_results:
                        statuses.add(method_results[actual_str])
                except Exception:
                    pass

        if not statuses:
            return ""  # No data for this group
        if len(statuses) == 1:
            # All dtypes in group have same status
            return status_map.get(statuses.pop(), "|error|")
        # Mixed results within group
        return "|mixed|"

    lines = [
        ".. list-table::",
        f"   :widths: 25 {' '.join(['10'] * len(groups))}",
        "   :header-rows: 1",
        "",
        "   * - Method",
    ]

    lines.extend(f"     - {g}" for g in groups)

    for method in sorted(by_method.keys()):
        lines.append(f"   * - ``{method}``")
        for g in groups:
            status = get_group_status(by_method[method], g)
            lines.append(f"     - {status}")

    return lines


def format_json(all_results: dict[str, list[OperationResult]]) -> str:
    """Format results as JSON (detailed, grouped by category)."""
    import json

    output = {}
    for category, results in all_results.items():
        output[category] = [
            {
                "method": r.method,
                "dtype": r.dtype,
                "result_type": r.result_type.value,
                "result_dtype": r.result_dtype,
                "error_message": r.error_message,
                "used_to_numpy": r.used_to_numpy,
                "used_elementwise": r.used_elementwise,
            }
            for r in results
        ]

    return json.dumps(output, indent=2)


def _get_current_pyarrow_version() -> str:
    """Get current PyArrow version as string."""
    try:
        import pyarrow as pa

        return pa.__version__
    except ImportError:
        return "unknown"


def _build_version_notes() -> dict[str, dict[str, Any]]:
    """Build version notes for version-gated methods."""
    version_notes = {}
    for method_name, info in VERSION_GATED_METHODS.items():
        accessor = info["accessor"]
        if accessor:
            method_key = f"{accessor}.{method_name}"
        else:
            method_key = method_name

        version_notes[method_key] = {
            "min_pyarrow_version": info["min_version"],
            "description": info["description"],
            "below_version_behavior": info["below_behavior"],
            "at_or_above_version_behavior": info["at_or_above_behavior"],
        }
    return version_notes


def format_json_lookup(all_results: dict[str, list[OperationResult]]) -> str:
    """
    Format results as JSON lookup table for downstream consumption.

    Structure: {method: {dtype: status}}
    Allows O(1) lookup: data["str.lower"]["string[pyarrow]"] -> "arrow"
    """
    import json

    # Map category to accessor prefix
    accessor_map = {
        "string_methods": "str",
        "datetime_methods": "dt",
        "timedelta_methods": "dt",
        "aggregations": None,
        "array_methods": None,
        "algorithmic_methods": None,
        "arithmetic": None,
        "comparison": None,
    }

    output: dict[str, Any] = {
        "schema_version": 2,  # Bumped for version_gated addition
        "pandas_version": pd.__version__,
        "pyarrow_version": _get_current_pyarrow_version(),
        "methods": {},
        "version_gated": _build_version_notes(),
    }

    for category, results in all_results.items():
        accessor = accessor_map.get(category)
        for r in results:
            # Build method key with accessor prefix
            if accessor:
                method_key = f"{accessor}.{r.method}"
            else:
                method_key = r.method

            if method_key not in output["methods"]:
                output["methods"][method_key] = {}

            output["methods"][method_key][r.dtype] = r.result_type.value

    return json.dumps(output, indent=2)


def format_json_by_dtype(all_results: dict[str, list[OperationResult]]) -> str:
    """
    Format results as JSON grouped by dtype for downstream consumption.

    Structure: {dtype: {method: status}}
    Useful for checking all methods for a specific dtype.
    """
    import json

    # Map category to accessor prefix
    accessor_map = {
        "string_methods": "str",
        "datetime_methods": "dt",
        "timedelta_methods": "dt",
        "aggregations": None,
        "array_methods": None,
        "algorithmic_methods": None,
        "arithmetic": None,
        "comparison": None,
    }

    output: dict[str, Any] = {
        "schema_version": 2,  # Bumped for version_gated addition
        "pandas_version": pd.__version__,
        "pyarrow_version": _get_current_pyarrow_version(),
        "dtypes": {},
        "version_gated": _build_version_notes(),
    }

    for category, results in all_results.items():
        accessor = accessor_map.get(category)
        for r in results:
            # Build method key with accessor prefix
            if accessor:
                method_key = f"{accessor}.{r.method}"
            else:
                method_key = r.method

            if r.dtype not in output["dtypes"]:
                output["dtypes"][r.dtype] = {}

            output["dtypes"][r.dtype][method_key] = r.result_type.value

    return json.dumps(output, indent=2)


# =============================================================================
# Main
# =============================================================================


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--format",
        choices=["rst", "json", "json-lookup", "json-by-dtype"],
        default="rst",
        help=(
            "Output format: rst (documentation), json (detailed), "
            "json-lookup (method->dtype->status), json-by-dtype (dtype->method->status)"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output file path (default: stdout)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Check if generated docs match existing file",
    )
    args = parser.parse_args()

    # Run all tests
    print("Running all operations...", file=sys.stderr)
    all_results = run_all_tests()

    # Count results
    total = sum(len(r) for r in all_results.values())
    print(f"Collected {total} results", file=sys.stderr)

    # Format output
    if args.format == "rst":
        output = format_rst_table(all_results)
    elif args.format == "json":
        output = format_json(all_results)
    elif args.format == "json-lookup":
        output = format_json_lookup(all_results)
    elif args.format == "json-by-dtype":
        output = format_json_by_dtype(all_results)
    else:
        output = format_rst_table(all_results)

    # Check mode
    if args.check:
        target = Path("doc/source/user_guide/arrow_fallbacks.rst")
        if not target.exists():
            print(f"ERROR: {target} does not exist.", file=sys.stderr)
            sys.exit(1)
        current = target.read_text()
        if current != output:
            print(f"ERROR: {target} is out of date.", file=sys.stderr)
            sys.exit(1)
        print("Arrow fallback documentation is up to date.")
        sys.exit(0)

    # Write output
    if args.output:
        args.output.write_text(output)
        print(f"Written to {args.output}", file=sys.stderr)
    else:
        # Output already has trailing newline, avoid double newline
        sys.stdout.write(output)


if __name__ == "__main__":
    main()
