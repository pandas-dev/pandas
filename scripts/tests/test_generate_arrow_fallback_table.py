"""Tests for scripts/generate_arrow_fallback_table.py."""

import json
from pathlib import Path
import subprocess
import sys

import pytest

from scripts.generate_arrow_fallback_table import (
    FallbackType,
    MethodCategory,
    MethodInfo,
    _extract_arithmetic_funcs,
    _extract_cmp_funcs,
    _extract_logical_funcs,
    _generate_summary_stats,
    analyze_arithmetic_operations,
    build_string_method_table,
    format_json,
    get_pyarrow_versions,
)

REPO_ROOT = Path(__file__).parents[2]
ARROW_ARRAY_PATH = REPO_ROOT / "pandas" / "core" / "arrays" / "arrow" / "array.py"


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def arrow_array_source():
    """Read the actual Arrow array source file."""
    return ARROW_ARRAY_PATH.read_text()


@pytest.fixture
def sample_method_info():
    """Factory fixture for creating MethodInfo objects."""

    def _create(
        name="lower",
        public_name="str.lower",
        category=MethodCategory.STRING,
        fallback_type=FallbackType.ARROW_NATIVE,
        arrow_function="pc.utf8_lower",
    ):
        return MethodInfo(
            name=name,
            public_name=public_name,
            category=category,
            fallback_type=fallback_type,
            arrow_function=arrow_function,
        )

    return _create


@pytest.fixture
def sample_all_methods(sample_method_info):
    """Create sample all_methods dict for testing output formatters."""
    return {
        "string": {
            "lower": sample_method_info(
                name="lower",
                public_name="str.lower",
                fallback_type=FallbackType.ARROW_NATIVE,
                arrow_function="pc.utf8_lower",
            ),
            "encode": sample_method_info(
                name="encode",
                public_name="str.encode",
                fallback_type=FallbackType.OBJECT_FALLBACK,
                arrow_function=None,
            ),
        },
        "arithmetic": {
            "add": sample_method_info(
                name="add",
                public_name="+",
                category=MethodCategory.ARITHMETIC,
                fallback_type=FallbackType.ARROW_NATIVE,
                arrow_function="pc.add_checked",
            ),
        },
    }


@pytest.fixture
def pyarrow_versions():
    """Sample PyArrow versions dict for testing."""
    return {"PYARROW_MIN_VERSION": "13.0.0"}


# =============================================================================
# Test extraction from real source files
# =============================================================================


class TestExtractFromRealSource:
    """Tests that verify extraction from actual pandas source files."""

    def test_extract_arithmetic_funcs_from_real_source(self, arrow_array_source):
        """Verify we can extract arithmetic funcs from the real source."""
        result = _extract_arithmetic_funcs(arrow_array_source)
        # Should find the actual operations defined in pandas
        assert len(result) > 0
        assert "add" in result
        assert "sub" in result
        assert "mul" in result

    def test_extract_cmp_funcs_from_real_source(self, arrow_array_source):
        """Verify we can extract comparison funcs from the real source."""
        result = _extract_cmp_funcs(arrow_array_source)
        assert len(result) > 0
        assert "eq" in result
        assert "lt" in result
        assert "gt" in result

    def test_extract_logical_funcs_from_real_source(self, arrow_array_source):
        """Verify we can extract logical funcs from the real source."""
        result = _extract_logical_funcs(arrow_array_source)
        assert len(result) > 0
        # Should have the base operations
        assert "and_" in result
        assert "or_" in result
        # Should NOT have reverse operations (they're filtered out)
        assert "rand_" not in result
        assert "ror_" not in result

    @pytest.mark.parametrize(
        "op,expected_func",
        [
            ("add", "pc.add_checked"),
            ("sub", "pc.subtract_checked"),
            ("mul", "pc.multiply_checked"),
        ],
    )
    def test_arithmetic_funcs_have_correct_values(
        self, arrow_array_source, op, expected_func
    ):
        """Verify arithmetic operations map to correct PyArrow functions."""
        result = _extract_arithmetic_funcs(arrow_array_source)
        assert result[op] == expected_func

    @pytest.mark.parametrize(
        "op,expected_func",
        [
            ("eq", "pc.equal"),
            ("ne", "pc.not_equal"),
            ("lt", "pc.less"),
            ("gt", "pc.greater"),
        ],
    )
    def test_cmp_funcs_have_correct_values(self, arrow_array_source, op, expected_func):
        """Verify comparison operations map to correct PyArrow functions."""
        result = _extract_cmp_funcs(arrow_array_source)
        assert result[op] == expected_func


# =============================================================================
# Test build_string_method_table
# =============================================================================


class TestBuildStringMethodTable:
    """Tests for build_string_method_table."""

    @pytest.fixture
    def string_method_table(self):
        """Cache the string method table for multiple tests."""
        return build_string_method_table()

    def test_returns_non_empty_dict(self, string_method_table):
        assert isinstance(string_method_table, dict)
        assert len(string_method_table) > 0

    def test_all_values_are_method_info(self, string_method_table):
        for info in string_method_table.values():
            assert isinstance(info, MethodInfo)
            assert info.category == MethodCategory.STRING

    @pytest.mark.parametrize(
        "expected_method",
        ["lower", "upper", "len", "strip", "startswith", "endswith"],
    )
    def test_contains_common_string_methods(self, string_method_table, expected_method):
        public_names = {info.public_name for info in string_method_table.values()}
        assert any(expected_method in name for name in public_names)


# =============================================================================
# Test analyze_arithmetic_operations
# =============================================================================


class TestAnalyzeArithmeticOperations:
    """Tests for analyze_arithmetic_operations."""

    @pytest.fixture
    def arithmetic_operations(self):
        """Cache arithmetic operations for multiple tests."""
        return analyze_arithmetic_operations()

    def test_returns_non_empty_list(self, arithmetic_operations):
        assert isinstance(arithmetic_operations, list)
        assert len(arithmetic_operations) > 0

    def test_all_items_are_method_info(self, arithmetic_operations):
        for info in arithmetic_operations:
            assert isinstance(info, MethodInfo)
            assert info.category == MethodCategory.ARITHMETIC

    @pytest.mark.parametrize(
        "expected_op",
        ["+", "-", "*", "/", "==", "!=", "<", ">"],
    )
    def test_includes_common_operators(self, arithmetic_operations, expected_op):
        public_names = {info.public_name for info in arithmetic_operations}
        assert expected_op in public_names


# =============================================================================
# Test get_pyarrow_versions
# =============================================================================


class TestGetPyarrowVersions:
    """Tests for get_pyarrow_versions."""

    @pytest.fixture
    def versions(self):
        return get_pyarrow_versions()

    def test_returns_dict_with_min_version(self, versions):
        assert isinstance(versions, dict)
        assert "PYARROW_MIN_VERSION" in versions

    def test_min_version_is_valid_semver(self, versions):
        version = versions["PYARROW_MIN_VERSION"]
        parts = version.split(".")
        assert len(parts) >= 2
        assert all(part.isdigit() for part in parts)


# =============================================================================
# Test _generate_summary_stats
# =============================================================================


class TestGenerateSummaryStats:
    """Tests for _generate_summary_stats."""

    @pytest.mark.parametrize(
        "expected_text",
        ["Summary", "String methods", "Arithmetic"],
    )
    def test_contains_expected_sections(self, sample_all_methods, expected_text):
        result = _generate_summary_stats(sample_all_methods)
        content = "\n".join(result)
        assert expected_text in content


# =============================================================================
# Test format_json
# =============================================================================


class TestFormatJson:
    """Tests for format_json."""

    def test_returns_valid_json(self, sample_all_methods, pyarrow_versions):
        result = format_json(sample_all_methods, pyarrow_versions)
        data = json.loads(result)
        assert isinstance(data, dict)

    @pytest.mark.parametrize(
        "key",
        ["summary", "categories", "pyarrow_min_version"],
    )
    def test_has_required_top_level_keys(
        self, sample_all_methods, pyarrow_versions, key
    ):
        result = format_json(sample_all_methods, pyarrow_versions)
        data = json.loads(result)
        assert key in data

    @pytest.mark.parametrize(
        "summary_key",
        ["by_category", "totals", "percentages"],
    )
    def test_summary_has_required_keys(
        self, sample_all_methods, pyarrow_versions, summary_key
    ):
        result = format_json(sample_all_methods, pyarrow_versions)
        data = json.loads(result)
        assert summary_key in data["summary"]

    def test_calculates_correct_totals(self, sample_all_methods, pyarrow_versions):
        result = format_json(sample_all_methods, pyarrow_versions)
        data = json.loads(result)
        totals = data["summary"]["totals"]
        # sample_all_methods has 3 methods: 2 arrow, 1 object
        assert totals["total"] == 3
        assert totals["arrow"] == 2
        assert totals["object"] == 1

    def test_calculates_correct_percentages(self, sample_all_methods, pyarrow_versions):
        result = format_json(sample_all_methods, pyarrow_versions)
        data = json.loads(result)
        percentages = data["summary"]["percentages"]
        # 2 arrow out of 3 total = 66.7%
        assert percentages["arrow_native"] == pytest.approx(66.7, abs=0.1)


# =============================================================================
# Test CLI
# =============================================================================


class TestCLI:
    """Tests for CLI functionality."""

    @pytest.fixture
    def script_path(self):
        return str(REPO_ROOT / "scripts" / "generate_arrow_fallback_table.py")

    def test_check_passes_when_docs_in_sync(self, script_path):
        result = subprocess.run(
            [sys.executable, script_path, "--check"],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0
        assert "up to date" in result.stdout

    @pytest.mark.parametrize(
        "expected_option",
        ["--check", "--format", "--output", "--category"],
    )
    def test_help_shows_all_options(self, script_path, expected_option):
        result = subprocess.run(
            [sys.executable, script_path, "--help"],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0
        assert expected_option in result.stdout

    @pytest.mark.parametrize("format_type", ["rst", "md", "json"])
    def test_all_format_options_produce_output(self, script_path, format_type):
        result = subprocess.run(
            [sys.executable, script_path, "--format", format_type],
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0
        assert len(result.stdout) > 0
