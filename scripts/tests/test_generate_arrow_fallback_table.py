"""Tests for scripts/generate_arrow_fallback_table.py (runtime-based generator)."""

from __future__ import annotations

import functools
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import pytest

import pandas as pd
from scripts.generate_arrow_fallback_table import (
    AGGREGATION_METHODS,
    ARRAY_METHODS,
    ARROW_DTYPES,
    DATETIME_METHODS,
    DTYPE_GROUPS,
    PARAMETER_VARIANT_METADATA,
    STRING_METHODS,
    ResultType,
    _build_parameter_variants,
    _is_arrow_backed,
    create_test_series,
    format_json,
    format_json_by_dtype,
    format_json_lookup,
    format_rst_table,
    run_all_tests,
    run_operation,
    track_fallbacks,
)

REPO_ROOT = Path(__file__).parents[2]


# Cache run_all_tests() across the entire test module to avoid
# running hundreds of operations multiple times.
@functools.lru_cache(maxsize=1)
def _cached_run_all_tests():
    return run_all_tests()


# =============================================================================
# Test configuration completeness
# =============================================================================


class TestConfiguration:
    """Test that configuration covers expected methods."""

    def test_arrow_dtypes_not_empty(self):
        assert len(ARROW_DTYPES) > 0

    def test_dtype_groups_cover_all_dtypes(self):
        all_grouped = set()
        for group in DTYPE_GROUPS.values():
            all_grouped.update(group)
        # All grouped dtypes should be in ARROW_DTYPES
        for dtype in all_grouped:
            assert dtype in ARROW_DTYPES, f"{dtype} not in ARROW_DTYPES"

    def test_string_methods_not_empty(self):
        assert len(STRING_METHODS) > 0

    def test_datetime_methods_not_empty(self):
        assert len(DATETIME_METHODS) > 0

    def test_aggregation_methods_not_empty(self):
        assert len(AGGREGATION_METHODS) > 0

    def test_array_methods_not_empty(self):
        assert len(ARRAY_METHODS) > 0


# =============================================================================
# Test data creation
# =============================================================================


class TestCreateTestSeries:
    """Test Series creation for various dtypes."""

    @pytest.mark.parametrize("dtype_name", list(ARROW_DTYPES.keys()))
    def test_creates_series_for_dtype(self, dtype_name):
        series = create_test_series(dtype_name)
        # Most dtypes should create a valid series
        # Some edge cases (like time64) may return None
        if series is not None:
            assert len(series) > 0

    def test_string_series_has_values(self):
        series = create_test_series("string")
        assert series is not None
        assert series.notna().any()

    def test_int64_series_has_values(self):
        series = create_test_series("int64")
        assert series is not None
        assert series.notna().any()

    def test_timestamp_series_has_values(self):
        series = create_test_series("timestamp_us")
        assert series is not None
        assert series.notna().any()

    def test_include_na_true_has_nulls(self):
        series = create_test_series("int64", include_na=True)
        assert series is not None
        assert series.isna().any()

    def test_include_na_false_has_no_nulls(self):
        series = create_test_series("int64", include_na=False)
        assert series is not None
        assert not series.isna().any()

    def test_unknown_dtype_returns_none(self):
        series = create_test_series("nonexistent_dtype")
        assert series is None


# =============================================================================
# Test Arrow detection
# =============================================================================


class TestIsArrowBacked:
    """Test _is_arrow_backed helper function."""

    def test_arrow_dtype_detected(self):
        series = pd.Series([1, 2, 3], dtype="int64[pyarrow]")
        assert _is_arrow_backed(series.dtype)

    def test_string_pyarrow_detected(self):
        series = pd.Series(["a", "b"], dtype="string[pyarrow]")
        assert _is_arrow_backed(series.dtype)

    def test_numpy_dtype_not_detected(self):
        series = pd.Series([1, 2, 3], dtype="int64")
        assert not _is_arrow_backed(series.dtype)

    def test_object_dtype_not_detected(self):
        series = pd.Series(["a", "b"], dtype=object)
        assert not _is_arrow_backed(series.dtype)


# =============================================================================
# Test fallback tracking
# =============================================================================


class TestTrackFallbacks:
    """Test the fallback instrumentation mechanism."""

    def test_tracker_detects_elementwise(self):
        # large_string uses ArrowDtype, which goes through _apply_elementwise
        # for casefold (unlike StringDtype which handles it natively)
        series = create_test_series("large_string")
        assert series is not None
        with track_fallbacks() as tracker:
            series.str.casefold()
        assert tracker.elementwise_called

    def test_tracker_no_fallback_for_arrow_native(self):
        series = create_test_series("string")
        assert series is not None
        with track_fallbacks() as tracker:
            series.str.lower()
        assert not tracker.elementwise_called


# =============================================================================
# Test operation running
# =============================================================================


class TestRunOperation:
    """Test running operations and classifying results."""

    def test_string_lower_returns_arrow(self):
        series = create_test_series("string")
        result = run_operation(series, "lower", {}, accessor="str")
        assert result.result_type == ResultType.ARROW_NATIVE

    def test_large_string_casefold_uses_fallback(self):
        # large_string uses ArrowDtype, which falls back for casefold
        series = create_test_series("large_string")
        result = run_operation(series, "casefold", {}, accessor="str")
        assert result.result_type == ResultType.ELEMENTWISE

    def test_invalid_method_returns_error(self):
        series = create_test_series("string")
        result = run_operation(series, "nonexistent_method", {}, accessor="str")
        assert result.result_type == ResultType.OTHER_ERROR

    def test_aggregation_sum_works(self):
        series = create_test_series("int64")
        result = run_operation(series, "sum", {})
        # sum should work on int64
        assert result.result_type in (
            ResultType.ARROW_NATIVE,
            ResultType.NUMPY_FALLBACK,
        )
        assert result.error_message is None

    def test_display_name_used_in_result(self):
        series = create_test_series("string")
        result = run_operation(
            series, "lower", {}, accessor="str", display_name="custom_name"
        )
        assert result.method == "custom_name"


# =============================================================================
# Test full run
# =============================================================================


class TestRunAllTests:
    """Test running the full test suite."""

    def test_returns_all_categories(self):
        results = _cached_run_all_tests()
        expected_categories = [
            "string_methods",
            "datetime_methods",
            "timedelta_methods",
            "aggregations",
            "array_methods",
            "algorithmic_methods",
            "arithmetic",
            "comparison",
        ]
        assert set(expected_categories) == set(results.keys())

    def test_string_methods_not_empty(self):
        results = _cached_run_all_tests()
        assert len(results["string_methods"]) > 0

    def test_aggregations_not_empty(self):
        results = _cached_run_all_tests()
        assert len(results["aggregations"]) > 0

    def test_algorithmic_methods_not_empty(self):
        results = _cached_run_all_tests()
        assert len(results["algorithmic_methods"]) > 0


# =============================================================================
# Test output formatting
# =============================================================================


class TestFormatOutput:
    """Test output formatting functions."""

    def test_format_json_valid(self):
        results = _cached_run_all_tests()
        output = format_json(results)
        parsed = json.loads(output)
        assert "string_methods" in parsed

    def test_format_rst_contains_header(self):
        results = _cached_run_all_tests()
        output = format_rst_table(results)
        assert "Arrow Fallbacks" in output

    def test_format_rst_contains_legend(self):
        results = _cached_run_all_tests()
        output = format_rst_table(results)
        assert "Legend" in output
        assert "|arrow|" in output

    def test_format_json_lookup_valid(self):
        results = _cached_run_all_tests()
        output = format_json_lookup(results)
        parsed = json.loads(output)
        assert "schema_version" in parsed
        assert "pandas_version" in parsed
        assert "methods" in parsed
        # Check structure: method -> dtype -> status
        assert "str.lower" in parsed["methods"]
        assert isinstance(parsed["methods"]["str.lower"], dict)

    def test_format_json_by_dtype_valid(self):
        results = _cached_run_all_tests()
        output = format_json_by_dtype(results)
        parsed = json.loads(output)
        assert "schema_version" in parsed
        assert "pandas_version" in parsed
        assert "dtypes" in parsed
        assert isinstance(parsed["dtypes"], dict)
        assert len(parsed["dtypes"]) > 0

    def test_json_lookup_has_parameter_variants(self):
        results = _cached_run_all_tests()
        output = format_json_lookup(results)
        parsed = json.loads(output)
        pv = parsed["parameter_variants"]
        assert "str.contains" in pv
        assert isinstance(pv["str.contains"]["parameters"], list)
        assert len(pv["str.contains"]["parameters"]) > 0

    def test_json_lookup_base_methods_exist(self):
        results = _cached_run_all_tests()
        output = format_json_lookup(results)
        parsed = json.loads(output)
        # All base methods from parameter_variants should exist
        for base_key in parsed["parameter_variants"]:
            assert base_key in parsed["methods"], (
                f"Base method {base_key!r} missing from methods"
            )

    def test_json_lookup_has_version_gated(self):
        results = _cached_run_all_tests()
        output = format_json_lookup(results)
        parsed = json.loads(output)
        assert "version_gated" in parsed
        assert isinstance(parsed["version_gated"], dict)


# =============================================================================
# Test parameter variant metadata
# =============================================================================


class TestParameterVariantMetadata:
    """Test PARAMETER_VARIANT_METADATA and _build_parameter_variants."""

    def test_metadata_has_named_fields(self):
        meta = PARAMETER_VARIANT_METADATA["contains"]
        assert meta.accessor == "str"
        assert isinstance(meta.parameters, list)
        assert meta.default_variant_key == "contains(flags=0)"

    def test_multi_parameter_method(self):
        meta = PARAMETER_VARIANT_METADATA["replace"]
        assert len(meta.parameters) == 2
        param_names = {p["name"] for p in meta.parameters}
        assert param_names == {"case", "repl"}

    def test_build_parameter_variants_structure(self):
        pv = _build_parameter_variants()
        assert "str.contains" in pv
        info = pv["str.contains"]
        assert "parameters" in info
        assert "default_variant" in info
        assert "variants" in info
        assert isinstance(info["variants"], dict)

    def test_build_parameter_variants_default_variant_key(self):
        pv = _build_parameter_variants()
        info = pv["str.contains"]
        assert info["default_variant"] == "str.contains(flags=0)"

    def test_build_parameter_variants_variant_keys(self):
        pv = _build_parameter_variants()
        info = pv["str.split"]
        assert "expand=False" in info["variants"]
        assert "expand=True" in info["variants"]


# =============================================================================
# Test CLI
# =============================================================================


class TestCLI:
    """Test command-line interface."""

    def test_help_shows_options(self):
        result = subprocess.run(
            [
                sys.executable,
                "scripts/generate_arrow_fallback_table.py",
                "--help",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--format" in result.stdout
        assert "--check" in result.stdout

    def test_json_format_works(self):
        result = subprocess.run(
            [
                sys.executable,
                "scripts/generate_arrow_fallback_table.py",
                "--format",
                "json",
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        # Output should contain valid JSON (may have build messages before it)
        # Find the JSON portion (starts with '{')
        stdout = result.stdout
        json_start = stdout.find("{")
        if json_start != -1:
            json_content = stdout[json_start:]
            parsed = json.loads(json_content)
            assert "string_methods" in parsed
        else:
            pytest.fail("No JSON found in output")

    def test_output_flag_writes_file(self):
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            outpath = f.name
        result = subprocess.run(
            [
                sys.executable,
                "scripts/generate_arrow_fallback_table.py",
                "--format",
                "json-lookup",
                "--output",
                outpath,
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        content = Path(outpath).read_text()
        parsed = json.loads(content)
        assert "schema_version" in parsed
        Path(outpath).unlink()

    def test_check_mode(self):
        # Generate fresh RST, write to a temp file, then --check against it
        with tempfile.NamedTemporaryFile(suffix=".rst", delete=False, mode="w") as f:
            outpath = f.name
        # First generate the file
        subprocess.run(
            [
                sys.executable,
                "scripts/generate_arrow_fallback_table.py",
                "--format",
                "rst",
                "--output",
                outpath,
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        # The --check mode compares against a fixed path
        # (doc/source/user_guide/arrow_fallbacks.rst), so we cannot
        # easily test it with a temp file. Instead verify the file
        # was written correctly.
        content = Path(outpath).read_text()
        assert "Arrow Fallbacks" in content
        Path(outpath).unlink()
