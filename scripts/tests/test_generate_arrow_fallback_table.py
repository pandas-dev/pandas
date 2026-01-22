"""Tests for scripts/generate_arrow_fallback_table.py (runtime-based generator)."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

from scripts.generate_arrow_fallback_table import (
    ARROW_DTYPES,
    DTYPE_GROUPS,
    STRING_METHODS,
    DATETIME_METHODS,
    AGGREGATION_METHODS,
    ARRAY_METHODS,
    ResultType,
    OperationResult,
    create_test_series,
    run_operation,
    run_all_tests,
    format_json,
    format_rst_table,
    _is_arrow_backed,
)

REPO_ROOT = Path(__file__).parents[2]


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


# =============================================================================
# Test Arrow detection
# =============================================================================


class TestIsArrowBacked:
    """Test _is_arrow_backed helper function."""

    def test_arrow_dtype_detected(self):
        import pandas as pd

        series = pd.Series([1, 2, 3], dtype="int64[pyarrow]")
        assert _is_arrow_backed(series.dtype)

    def test_string_pyarrow_detected(self):
        import pandas as pd

        series = pd.Series(["a", "b"], dtype="string[pyarrow]")
        assert _is_arrow_backed(series.dtype)

    def test_numpy_dtype_not_detected(self):
        import pandas as pd

        series = pd.Series([1, 2, 3], dtype="int64")
        assert not _is_arrow_backed(series.dtype)

    def test_object_dtype_not_detected(self):
        import pandas as pd

        series = pd.Series(["a", "b"], dtype=object)
        assert not _is_arrow_backed(series.dtype)


# =============================================================================
# Test operation running
# =============================================================================


class TestRunOperation:
    """Test running operations and classifying results."""

    def test_string_lower_returns_arrow(self):
        series = create_test_series("string")
        result = run_operation(series, "lower", {}, accessor="str")
        assert result.result_type == ResultType.ARROW_NATIVE

    def test_string_casefold_uses_fallback(self):
        series = create_test_series("string")
        result = run_operation(series, "casefold", {}, accessor="str")
        # casefold uses to_numpy, so it's a fallback
        assert result.result_type in (ResultType.NUMPY_FALLBACK, ResultType.ELEMENTWISE)

    def test_invalid_method_returns_error(self):
        series = create_test_series("string")
        result = run_operation(series, "nonexistent_method", {}, accessor="str")
        assert result.result_type == ResultType.OTHER_ERROR

    def test_aggregation_sum_works(self):
        series = create_test_series("int64")
        result = run_operation(series, "sum", {})
        # sum should work on int64
        assert result.result_type in (ResultType.ARROW_NATIVE, ResultType.NUMPY_FALLBACK)
        assert result.error_message is None


# =============================================================================
# Test full run
# =============================================================================


class TestRunAllTests:
    """Test running the full test suite."""

    def test_returns_all_categories(self):
        results = run_all_tests()
        expected_categories = [
            "string_methods",
            "datetime_methods",
            "timedelta_methods",
            "aggregations",
            "array_methods",
            "arithmetic",
            "comparison",
        ]
        for cat in expected_categories:
            assert cat in results

    def test_string_methods_not_empty(self):
        results = run_all_tests()
        assert len(results["string_methods"]) > 0

    def test_aggregations_not_empty(self):
        results = run_all_tests()
        assert len(results["aggregations"]) > 0


# =============================================================================
# Test output formatting
# =============================================================================


class TestFormatOutput:
    """Test output formatting functions."""

    def test_format_json_valid(self):
        results = run_all_tests()
        output = format_json(results)
        # Should be valid JSON
        parsed = json.loads(output)
        assert "string_methods" in parsed

    def test_format_rst_contains_header(self):
        results = run_all_tests()
        output = format_rst_table(results)
        assert "Arrow Method Support Reference" in output

    def test_format_rst_contains_legend(self):
        results = run_all_tests()
        output = format_rst_table(results)
        assert "Legend" in output
        assert "|arrow|" in output


# =============================================================================
# Test CLI
# =============================================================================


class TestCLI:
    """Test command-line interface."""

    def test_help_shows_options(self):
        result = subprocess.run(
            [sys.executable, "scripts/generate_arrow_fallback_table.py", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "--format" in result.stdout
        assert "--check" in result.stdout

    def test_json_format_works(self):
        result = subprocess.run(
            [sys.executable, "scripts/generate_arrow_fallback_table.py", "--format", "json"],
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
