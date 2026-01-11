"""
Generate comprehensive documentation for Arrow fallback behavior in pandas.

This script introspects pandas' Arrow-backed arrays to determine which methods
use native PyArrow compute functions vs. falling back to NumPy/Python implementations.

Coverage includes:
- String methods (Series.str.*)
- Arithmetic operations (+, -, *, /, //, %, **)
- Datetime methods (Series.dt.*)
- Aggregation methods (sum, mean, min, max, etc.)
- Array methods (unique, dropna, fillna, etc.)
- List accessor (Series.list.*)
- Struct accessor (Series.struct.*)

Usage:
    # Generate RST documentation to stdout
    python scripts/generate_arrow_fallback_table.py

    # Generate to a specific file
    python scripts/generate_arrow_fallback_table.py --output path/to/file.rst

    # Generate in different formats
    python scripts/generate_arrow_fallback_table.py --format md
    python scripts/generate_arrow_fallback_table.py --format json

    # Generate only a specific category
    python scripts/generate_arrow_fallback_table.py --category string
    python scripts/generate_arrow_fallback_table.py --category datetime

    # Check if generated docs match existing file (for CI/pre-commit)
    python scripts/generate_arrow_fallback_table.py --check

The output can be:
- RST format (default) for inclusion in pandas docs
- Markdown format for README or other documentation
- JSON format for programmatic use

The --check flag compares generated output with the existing documentation file
and exits with code 1 if they differ. This is used by pre-commit hooks and CI
to ensure documentation stays in sync with the source code.
"""

from __future__ import annotations

import argparse
import ast
from dataclasses import (
    dataclass,
    field,
)
from enum import Enum
import json
from pathlib import Path
import re
import sys

# Repository root for file-based analysis (no pandas import needed)
REPO_ROOT = Path(__file__).parent.parent


class FallbackType(Enum):
    """Classification of how a method handles Arrow vs NumPy."""

    ARROW_NATIVE = "arrow"  # Uses PyArrow compute exclusively
    CONDITIONAL = "conditional"  # Falls back under certain conditions
    ELEMENTWISE = "elementwise"  # Uses _apply_elementwise (Python loop)
    OBJECT_FALLBACK = "object"  # Always falls back to ObjectStringArrayMixin
    VERSION_GATED = "version_gated"  # Behavior depends on PyArrow version
    NOT_IMPLEMENTED = "not_implemented"  # Raises NotImplementedError or NotImplemented


class MethodCategory(Enum):
    """Category of methods being analyzed."""

    STRING = "string"
    ARITHMETIC = "arithmetic"
    DATETIME = "datetime"
    AGGREGATION = "aggregation"
    LIST = "list"
    STRUCT = "struct"
    ARRAY_METHODS = "array_methods"


@dataclass
class MethodInfo:
    """Information about a method's Arrow support."""

    name: str
    public_name: str  # The public API name (e.g., Series.str.X, Series.dt.X)
    category: MethodCategory
    fallback_type: FallbackType
    arrow_function: str | None = None  # The pc.X function used
    min_pyarrow_version: str | None = None  # Minimum version for Arrow path
    fallback_conditions: list[str] = field(default_factory=list)
    performance_notes: str | None = None
    github_issues: list[str] = field(default_factory=list)
    notes: str | None = None  # Additional notes


# Map internal method names to public API names
STRING_METHOD_NAME_MAP = {
    "_str_len": "len",
    "_str_lower": "lower",
    "_str_upper": "upper",
    "_str_strip": "strip",
    "_str_lstrip": "lstrip",
    "_str_rstrip": "rstrip",
    "_str_pad": "pad/center/ljust/rjust",
    "_str_get": "get",
    "_str_slice": "slice",
    "_str_slice_replace": "slice_replace",
    "_str_replace": "replace",
    "_str_capitalize": "capitalize",
    "_str_title": "title",
    "_str_swapcase": "swapcase",
    "_str_removeprefix": "removeprefix",
    "_str_removesuffix": "removesuffix",
    "_str_startswith": "startswith",
    "_str_endswith": "endswith",
    "_str_isalnum": "isalnum",
    "_str_isalpha": "isalpha",
    "_str_isascii": "isascii",
    "_str_isdecimal": "isdecimal",
    "_str_isdigit": "isdigit",
    "_str_islower": "islower",
    "_str_isnumeric": "isnumeric",
    "_str_isspace": "isspace",
    "_str_istitle": "istitle",
    "_str_isupper": "isupper",
    "_str_contains": "contains",
    "_str_match": "match",
    "_str_fullmatch": "fullmatch",
    "_str_find": "find",
    "_str_rfind": "rfind",
    "_str_count": "count",
    "_str_repeat": "repeat",
    "_str_partition": "partition",
    "_str_rpartition": "rpartition",
    "_str_casefold": "casefold",
    "_str_encode": "encode",
    "_str_findall": "findall",
    "_str_split": "split",
    "_str_rsplit": "rsplit",
    "_str_join": "join",
    "_str_wrap": "wrap",
    "_str_translate": "translate",
    "_str_normalize": "normalize",
    "_str_zfill": "zfill",
    "_str_index": "index",
    "_str_rindex": "rindex",
    "_str_get_dummies": "get_dummies",
    "_str_extract": "extract",
}

# Datetime method name mappings
DATETIME_METHOD_NAME_MAP = {
    "_dt_year": "year",
    "_dt_month": "month",
    "_dt_day": "day",
    "_dt_hour": "hour",
    "_dt_minute": "minute",
    "_dt_second": "second",
    "_dt_microsecond": "microsecond",
    "_dt_nanosecond": "nanosecond",
    "_dt_day_of_week": "day_of_week/dayofweek/weekday",
    "_dt_day_of_year": "day_of_year/dayofyear",
    "_dt_quarter": "quarter",
    "_dt_is_month_start": "is_month_start",
    "_dt_is_month_end": "is_month_end",
    "_dt_is_quarter_start": "is_quarter_start",
    "_dt_is_quarter_end": "is_quarter_end",
    "_dt_is_year_start": "is_year_start",
    "_dt_is_year_end": "is_year_end",
    "_dt_is_leap_year": "is_leap_year",
    "_dt_days_in_month": "days_in_month/daysinmonth",
    "_dt_date": "date",
    "_dt_time": "time",
    "_dt_tz": "tz",
    "_dt_unit": "unit",
    "_dt_normalize": "normalize",
    "_dt_strftime": "strftime",
    "_dt_round": "round",
    "_dt_floor": "floor",
    "_dt_ceil": "ceil",
    "_dt_day_name": "day_name",
    "_dt_month_name": "month_name",
    "_dt_tz_localize": "tz_localize",
    "_dt_tz_convert": "tz_convert",
    "_dt_isocalendar": "isocalendar",
    "_dt_to_pydatetime": "to_pydatetime",
    "_dt_to_pytimedelta": "to_pytimedelta",
    "_dt_total_seconds": "total_seconds",
    "_dt_as_unit": "as_unit",
    # Timedelta components
    "_dt_days": "days",
    "_dt_hours": "hours (td)",
    "_dt_minutes": "minutes (td)",
    "_dt_seconds": "seconds (td)",
    "_dt_milliseconds": "milliseconds",
    "_dt_microseconds": "microseconds (td)",
    "_dt_nanoseconds": "nanoseconds (td)",
}


def _format_pyarrow_version(major: int | str, minor: int | str) -> str:
    """Format a PyArrow version as semver (major.minor.0)."""
    return f"{major}.{minor}.0"


def _has_arrow_native_pattern(source: str) -> bool:
    """
    Check if source contains Arrow-native patterns beyond just pc.* calls.

    This detects various ways Arrow compute can be used:
    - pc.* (pyarrow.compute as pc)
    - pa.compute.* (direct pyarrow.compute access)
    - compute.* (imported as compute)
    - Arrow array methods like .dictionary_encode(), .unique(), .value_counts()
    - ChunkedArray methods
    - Direct pa.array / pa.chunked_array operations
    """
    # PyArrow compute patterns
    compute_patterns = [
        "pc.",  # import pyarrow.compute as pc
        "pa.compute.",  # pyarrow.compute direct
        "compute.",  # import pyarrow.compute as compute
    ]

    # Arrow array/ChunkedArray native methods
    arrow_array_methods = [
        ".dictionary_encode(",
        ".unique(",
        ".value_counts(",
        ".cast(",
        ".fill_null(",
        ".drop_null(",
        ".take(",
        ".filter(",
        ".sort_indices(",
        ".slice(",
        ".combine_chunks(",
        ".chunk(",
        ".chunks",
    ]

    # Check compute patterns
    for pattern in compute_patterns:
        if pattern in source:
            return True

    # Check arrow array methods (but be careful about false positives)
    # These should be on Arrow objects, not pandas objects
    arrow_indicators = ["._pa_array", "._ndarray", "pa.array", "pa.chunked_array"]
    has_arrow_context = any(ind in source for ind in arrow_indicators)

    if has_arrow_context:
        for method in arrow_array_methods:
            if method in source:
                return True

    return False


def _dedupe_preserve_order(items: list[str]) -> list[str]:
    """Deduplicate a list while preserving insertion order."""
    seen = set()
    result = []
    for item in items:
        if item not in seen:
            seen.add(item)
            result.append(item)
    return result


def get_pyarrow_versions() -> dict[str, str]:
    """Get PyArrow version constants from pandas.compat.pyarrow source file."""
    pyarrow_compat_path = REPO_ROOT / "pandas" / "compat" / "pyarrow.py"

    if not pyarrow_compat_path.exists():
        return {"PYARROW_MIN_VERSION": "13.0.0"}

    source = pyarrow_compat_path.read_text()
    versions = {}

    # Extract PYARROW_MIN_VERSION
    min_version_match = re.search(
        r'PYARROW_MIN_VERSION\s*=\s*["\']([^"\']+)["\']', source
    )
    if min_version_match:
        versions["PYARROW_MIN_VERSION"] = min_version_match.group(1)
    else:
        versions["PYARROW_MIN_VERSION"] = "13.0.0"

    # Extract version checks like pa_version_under17p0
    for match in re.finditer(r"pa_version_under(\d+)p(\d+)", source):
        version_var = f"pa_version_under{match.group(1)}p{match.group(2)}"
        version_str = _format_pyarrow_version(match.group(1), match.group(2))
        versions[version_var] = version_str

    return versions


# =============================================================================
# String Methods Analysis
# =============================================================================


def analyze_arrow_string_mixin() -> list[MethodInfo]:
    """Analyze ArrowStringArrayMixin to identify Arrow-native methods."""
    methods = []

    mixin_path = REPO_ROOT / "pandas" / "core" / "arrays" / "_arrow_string_mixins.py"
    if not mixin_path.exists():
        return methods

    source = mixin_path.read_text()
    tree = ast.parse(source)

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name.startswith("_str_"):
            method_name = node.name
            public_name = STRING_METHOD_NAME_MAP.get(
                method_name, method_name.replace("_str_", "")
            )

            method_source = ast.get_source_segment(source, node)
            info = _analyze_string_method_source(
                method_name, public_name, method_source or ""
            )
            if info:
                methods.append(info)

    return methods


def _analyze_string_method_source(
    method_name: str, public_name: str, source: str
) -> MethodInfo | None:
    """Analyze a string method's source to determine fallback behavior."""
    arrow_func = _extract_string_arrow_function(source)
    fallback_type = _determine_string_fallback_type(source, arrow_func)
    conditions = _extract_string_fallback_conditions(source)
    notes = _extract_string_method_notes(source, fallback_type)
    min_version = _extract_min_pyarrow_version(source)

    # Skip methods that have no Arrow path at all (return None if no pc. calls)
    if fallback_type is None:
        return None

    return MethodInfo(
        name=method_name,
        public_name=public_name,
        category=MethodCategory.STRING,
        fallback_type=fallback_type,
        arrow_function=arrow_func,
        min_pyarrow_version=min_version,
        fallback_conditions=conditions,
        notes=notes,
        performance_notes=_get_performance_note(fallback_type),
        github_issues=_extract_github_issues(source),
    )


def _determine_string_fallback_type(
    source: str, arrow_func: str | None
) -> FallbackType | None:
    """Determine the fallback type for a string method based on source patterns."""
    has_pc_call = _has_arrow_native_pattern(source)
    has_elementwise = "_apply_elementwise" in source
    has_version_gate = "pa_version_under" in source
    has_not_implemented = "NotImplementedError" in source
    has_super_call = "super()." in source or "super()" in source

    # If method has version gating, it's version-gated
    if has_version_gate:
        return FallbackType.VERSION_GATED

    # If it only uses _apply_elementwise without pc. calls, it's elementwise
    if has_elementwise and not has_pc_call:
        return FallbackType.ELEMENTWISE

    # If it raises NotImplementedError for some cases, it's conditional
    if has_not_implemented:
        return FallbackType.CONDITIONAL

    # If it has both elementwise and pc calls, it's conditional
    if has_elementwise and has_pc_call:
        return FallbackType.CONDITIONAL

    # Handle super() calls based on whether there's also an Arrow path
    if has_super_call:
        if has_pc_call:
            # Both Arrow and fallback paths exist
            return FallbackType.CONDITIONAL
        else:
            # Only super() call, no Arrow path - pure object fallback
            return FallbackType.OBJECT_FALLBACK

    # If it has pc. calls, it's Arrow native
    if has_pc_call or arrow_func:
        return FallbackType.ARROW_NATIVE

    # No Arrow path detected
    return None


def _extract_string_arrow_function(source: str) -> str | None:
    """Extract the PyArrow function used for a string method programmatically."""
    # Extract all pc.X function calls from source
    pc_calls = re.findall(r"pc\.(\w+)\s*\(", source)
    if not pc_calls:
        return None

    # Get unique calls, preserving order
    unique_calls = list(dict.fromkeys(pc_calls))

    # Filter out utility functions that aren't the main operation
    utility_funcs = {"if_else", "is_null", "invert", "or_", "and_"}
    main_funcs = [c for c in unique_calls if c not in utility_funcs]

    if main_funcs:
        # Return the first main function (usually the primary operation)
        return f"pc.{main_funcs[0]}"
    elif unique_calls:
        # Fall back to first call if all are utilities
        return f"pc.{unique_calls[0]}"

    return None


def _extract_string_fallback_conditions(source: str) -> list[str]:
    """Extract conditions that trigger fallback from string method source."""
    conditions = []

    # Check for regex flags condition
    if "flags" in source:
        if "if flags" in source or "flags !=" in source or "flags:" in source:
            conditions.append("regex flags used")

    # Check for compiled regex pattern
    if "isinstance(pat, re.Pattern)" in source:
        conditions.append("compiled regex pattern")

    # Check for callable replacement
    if "callable(repl)" in source:
        conditions.append("callable replacement")

    # Check for case-insensitive matching
    if "not case" in source or "case is False" in source or "case=False" in source:
        conditions.append("case-insensitive matching")

    # Check for named group references
    if r"\g<" in source or "\\g<" in source:
        conditions.append("named group references in replacement")

    # Check for array-like repeats
    if "not isinstance(repeats, int)" in source:
        conditions.append("array-like repeats")

    # Check for start/end parameters triggering fallback
    if "start" in source and "end" in source and "_apply_elementwise" in source:
        conditions.append("start/end parameters with edge cases")

    # Check for version gating
    version_patterns = re.findall(r"pa_version_under(\d+)p(\d+)", source)
    if version_patterns:
        major, minor = version_patterns[0]
        conditions.append(f"PyArrow < {_format_pyarrow_version(major, minor)}")

    return conditions


def _extract_string_method_notes(
    source: str, fallback_type: FallbackType | None = None
) -> str | None:
    """Extract notes for string methods programmatically from source patterns."""
    notes = []

    # Add fallback flavor description for object fallback methods
    if fallback_type == FallbackType.OBJECT_FALLBACK:
        flavor = _determine_fallback_flavor(source)
        notes.append(_get_fallback_flavor_description(flavor))
    elif fallback_type == FallbackType.ELEMENTWISE:
        # Check if it uses _apply_elementwise specifically
        if "_apply_elementwise" in source:
            notes.append("via _apply_elementwise()")
        else:
            notes.append("element-wise loop")

    # Check for element-wise Python operations (conditional case)
    if "_apply_elementwise" in source and fallback_type == FallbackType.CONDITIONAL:
        notes.append("_apply_elementwise() fallback")

    # Check for object dtype conversion
    if "astype(object" in source:
        notes.append("Converts to object dtype")

    # Check for NotImplementedError with specific message
    not_impl_match = re.search(r'NotImplementedError\(\s*["\']([^"\']+)["\']', source)
    if not_impl_match:
        # Extract a short summary from the message
        msg = not_impl_match.group(1)
        if len(msg) > 50:
            msg = msg[:47] + "..."
        notes.append(f"Raises: {msg}")

    # Check for GH issue references
    gh_matches = re.findall(r"GH#?(\d+)", source)
    if gh_matches:
        notes.append(f"Related: GH#{gh_matches[0]}")

    return "; ".join(notes) if notes else None


def _extract_min_pyarrow_version(source: str) -> str | None:
    """Extract the minimum PyArrow version required for Arrow path."""
    # Look for pa_version_under patterns
    version_patterns = re.findall(r"pa_version_under(\d+)p(\d+)", source)
    if version_patterns:
        # Return the highest version found (the minimum required)
        versions = [(int(m), int(n)) for m, n in version_patterns]
        max_version = max(versions)
        return _format_pyarrow_version(max_version[0], max_version[1])
    return None


def analyze_arrow_string_array() -> list[MethodInfo]:
    """Analyze ArrowStringArray for conditional fallbacks to ObjectStringArrayMixin."""
    methods = []

    array_path = REPO_ROOT / "pandas" / "core" / "arrays" / "string_arrow.py"
    if not array_path.exists():
        return methods

    source = array_path.read_text()
    tree = ast.parse(source)

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name.startswith("_str_"):
            method_name = node.name
            public_name = STRING_METHOD_NAME_MAP.get(
                method_name, method_name.replace("_str_", "")
            )

            method_source = ast.get_source_segment(source, node)
            if method_source and "super()" in method_source:
                # Determine fallback type first (needed for notes)
                has_pc_call = _has_arrow_native_pattern(method_source)
                has_super_call = "super()" in method_source

                if has_pc_call and has_super_call:
                    # Both Arrow and fallback paths exist
                    fallback_type = FallbackType.CONDITIONAL
                elif has_super_call and not has_pc_call:
                    # Only super() call, no Arrow path - pure object fallback
                    fallback_type = FallbackType.OBJECT_FALLBACK
                else:
                    # Only Arrow path
                    fallback_type = FallbackType.ARROW_NATIVE

                # Use programmatic extraction
                arrow_func = _extract_string_arrow_function(method_source)
                conditions = _extract_string_fallback_conditions(method_source)
                notes = _extract_string_method_notes(method_source, fallback_type)

                methods.append(
                    MethodInfo(
                        name=method_name,
                        public_name=public_name,
                        category=MethodCategory.STRING,
                        fallback_type=fallback_type,
                        arrow_function=arrow_func,
                        fallback_conditions=conditions,
                        notes=notes,
                        performance_notes=_get_performance_note(fallback_type),
                        github_issues=_extract_github_issues(method_source),
                    )
                )

    return methods


def analyze_object_string_mixin() -> list[MethodInfo]:
    """Analyze ObjectStringArrayMixin for methods that always use object dtype."""
    methods = []

    object_array_path = REPO_ROOT / "pandas" / "core" / "strings" / "object_array.py"
    if not object_array_path.exists():
        return methods

    source = object_array_path.read_text()
    tree = ast.parse(source)

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name.startswith("_str_"):
            name = node.name
            public_name = STRING_METHOD_NAME_MAP.get(name, name.replace("_str_", ""))
            method_source = ast.get_source_segment(source, node)

            # Determine specific fallback flavor from source
            if method_source:
                flavor = _determine_fallback_flavor(method_source)
                notes = _get_fallback_flavor_description(flavor)
            else:
                # Default for ObjectStringArrayMixin methods
                notes = "Object string implementation"

            methods.append(
                MethodInfo(
                    name=name,
                    public_name=public_name,
                    category=MethodCategory.STRING,
                    fallback_type=FallbackType.OBJECT_FALLBACK,
                    notes=notes,
                    performance_notes=_get_performance_note(
                        FallbackType.OBJECT_FALLBACK
                    ),
                )
            )

    return methods


def analyze_arrow_extension_array_string_methods() -> list[MethodInfo]:
    """
    Analyze ArrowExtensionArray for string methods that use _apply_elementwise.

    These are string methods implemented directly in array.py that typically
    use _apply_elementwise() because PyArrow doesn't have a native compute
    function for them.
    """
    methods = []

    array_path = REPO_ROOT / "pandas" / "core" / "arrays" / "arrow" / "array.py"
    if not array_path.exists():
        return methods

    source = array_path.read_text()
    tree = ast.parse(source)

    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name.startswith("_str_"):
            method_name = node.name
            public_name = STRING_METHOD_NAME_MAP.get(
                method_name, method_name.replace("_str_", "")
            )

            method_source = ast.get_source_segment(source, node)
            if method_source:
                # Determine fallback type and flavor
                has_elementwise = "_apply_elementwise" in method_source
                has_pc_call = _has_arrow_native_pattern(method_source)

                if has_elementwise and not has_pc_call:
                    # Pure element-wise fallback
                    fallback_type = FallbackType.ELEMENTWISE
                elif has_elementwise and has_pc_call:
                    # Conditional - uses Arrow when possible
                    fallback_type = FallbackType.CONDITIONAL
                elif has_pc_call:
                    fallback_type = FallbackType.ARROW_NATIVE
                else:
                    # Not relevant to our analysis
                    continue

                # Get notes with fallback flavor
                notes = _extract_string_method_notes(method_source, fallback_type)

                methods.append(
                    MethodInfo(
                        name=method_name,
                        public_name=public_name,
                        category=MethodCategory.STRING,
                        fallback_type=fallback_type,
                        arrow_function=_extract_string_arrow_function(method_source),
                        fallback_conditions=_extract_string_fallback_conditions(
                            method_source
                        ),
                        notes=notes,
                        performance_notes=_get_performance_note(fallback_type),
                        github_issues=_extract_github_issues(method_source),
                    )
                )

    return methods


# =============================================================================
# Arithmetic Operations Analysis
# =============================================================================


def analyze_arithmetic_operations() -> list[MethodInfo]:
    """Analyze ArrowExtensionArray for arithmetic operation support programmatically."""
    methods = []

    array_path = REPO_ROOT / "pandas" / "core" / "arrays" / "arrow" / "array.py"
    if not array_path.exists():
        return methods

    source = array_path.read_text()

    # Extract ARROW_ARITHMETIC_FUNCS dictionary from source
    arithmetic_ops = _extract_arithmetic_funcs(source)
    for op_name, arrow_func in arithmetic_ops.items():
        public_name = _get_arithmetic_public_name(op_name)
        if arrow_func is None:
            methods.append(
                MethodInfo(
                    name=op_name,
                    public_name=public_name,
                    category=MethodCategory.ARITHMETIC,
                    fallback_type=FallbackType.NOT_IMPLEMENTED,
                    notes="Not supported for Arrow arrays, raises TypeError",
                    performance_notes=_get_performance_note(
                        FallbackType.NOT_IMPLEMENTED
                    ),
                )
            )
        else:
            methods.append(
                MethodInfo(
                    name=op_name,
                    public_name=public_name,
                    category=MethodCategory.ARITHMETIC,
                    fallback_type=FallbackType.ARROW_NATIVE,
                    arrow_function=arrow_func,
                    performance_notes=_get_performance_note(FallbackType.ARROW_NATIVE),
                )
            )

    # Extract ARROW_CMP_FUNCS dictionary from source
    comparison_ops = _extract_cmp_funcs(source)
    for op_name, arrow_func in comparison_ops.items():
        public_name = _get_comparison_public_name(op_name)
        methods.append(
            MethodInfo(
                name=op_name,
                public_name=public_name,
                category=MethodCategory.ARITHMETIC,
                fallback_type=FallbackType.CONDITIONAL,
                arrow_function=arrow_func,
                fallback_conditions=[
                    "ArrowNotImplementedError for unsupported type combinations"
                ],
                performance_notes=_get_performance_note(FallbackType.CONDITIONAL),
            )
        )

    # Extract ARROW_LOGICAL_FUNCS dictionary from source
    logical_ops = _extract_logical_funcs(source)
    for op_name, arrow_func in logical_ops.items():
        public_name = _get_logical_public_name(op_name)
        methods.append(
            MethodInfo(
                name=op_name,
                public_name=public_name,
                category=MethodCategory.ARITHMETIC,
                fallback_type=FallbackType.ARROW_NATIVE,
                arrow_function=arrow_func,
                notes="Uses Kleene logic for boolean, bit_wise for integers",
                performance_notes=_get_performance_note(FallbackType.ARROW_NATIVE),
            )
        )

    return methods


def _extract_pc_function_from_node(node: ast.expr, source: str) -> str | None:
    """Extract pc.X function name from an AST node."""
    if isinstance(node, ast.Attribute):
        # Direct pc.X reference like pc.add_checked
        if isinstance(node.value, ast.Name) and node.value.id == "pc":
            return f"pc.{node.attr}"
    elif isinstance(node, ast.Name):
        # NotImplemented or other names
        if node.id == "NotImplemented":
            return None
        # Could be a reference to another function like floordiv_compat
        return node.id
    elif isinstance(node, ast.Lambda):
        # Lambda - extract pc.X from the body or helper function name
        node_source = ast.get_source_segment(source, node)
        if node_source:
            # First try to find direct pc.X calls
            pc_match = re.search(r"pc\.(\w+)", node_source)
            if pc_match:
                return f"pc.{pc_match.group(1)}"
            # Check for helper function calls like floordiv_compat(x, y)
            # These are valid implementations, not NotImplemented
            helper_match = re.search(r"(\w+_compat)\s*\(", node_source)
            if helper_match:
                return helper_match.group(1)
    elif isinstance(node, ast.Call):
        # Function call - check if it's pc.X(...)
        if isinstance(node.func, ast.Attribute):
            if isinstance(node.func.value, ast.Name) and node.func.value.id == "pc":
                return f"pc.{node.func.attr}"
    return None


def _extract_dict_from_assignment(source: str, var_name: str) -> dict[str, str | None]:
    """Extract a dictionary assignment from source using AST."""
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return {}

    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == var_name:
                    if isinstance(node.value, ast.Dict):
                        result = {}
                        for key, value in zip(
                            node.value.keys, node.value.values, strict=True
                        ):
                            if isinstance(key, ast.Constant) and isinstance(
                                key.value, str
                            ):
                                pc_func = _extract_pc_function_from_node(value, source)
                                result[key.value] = pc_func
                        return result
    return {}


def _extract_arithmetic_funcs(source: str) -> dict[str, str | None]:
    """Extract ARROW_ARITHMETIC_FUNCS dictionary from source using AST."""
    return _extract_dict_from_assignment(source, "ARROW_ARITHMETIC_FUNCS")


def _extract_cmp_funcs(source: str) -> dict[str, str]:
    """Extract ARROW_CMP_FUNCS dictionary from source using AST."""
    result = _extract_dict_from_assignment(source, "ARROW_CMP_FUNCS")
    # Filter out None values for comparison funcs (they should all have pc.X)
    return {k: v for k, v in result.items() if v is not None}


def _extract_logical_funcs(source: str) -> dict[str, str]:
    """Extract ARROW_LOGICAL_FUNCS dictionary from source using AST."""
    result = _extract_dict_from_assignment(source, "ARROW_LOGICAL_FUNCS")
    # Filter out None values and reverse operations (rand_, ror_, rxor)
    return {k: v for k, v in result.items() if v is not None and not k.startswith("r")}


def _get_arithmetic_public_name(op_name: str) -> str:
    """Get the public name for an arithmetic operation."""
    mapping = {
        "add": "+",
        "radd": "+ (right)",
        "sub": "-",
        "rsub": "- (right)",
        "mul": "*",
        "rmul": "* (right)",
        "truediv": "/",
        "rtruediv": "/ (right)",
        "floordiv": "//",
        "rfloordiv": "// (right)",
        "mod": "%",
        "rmod": "% (right)",
        "divmod": "divmod",
        "rdivmod": "divmod (right)",
        "pow": "**",
        "rpow": "** (right)",
    }
    return mapping.get(op_name, op_name)


def _get_comparison_public_name(op_name: str) -> str:
    """Get the public name for a comparison operation."""
    mapping = {
        "eq": "==",
        "ne": "!=",
        "lt": "<",
        "gt": ">",
        "le": "<=",
        "ge": ">=",
    }
    return mapping.get(op_name, op_name)


def _get_logical_public_name(op_name: str) -> str:
    """Get the public name for a logical operation."""
    mapping = {
        "and_": "&",
        "or_": "|",
        "xor": "^",
    }
    return mapping.get(op_name, op_name)


# =============================================================================
# Datetime Methods Analysis
# =============================================================================


def analyze_datetime_methods() -> list[MethodInfo]:
    """Analyze ArrowExtensionArray for datetime accessor methods."""
    methods = []

    array_path = REPO_ROOT / "pandas" / "core" / "arrays" / "arrow" / "array.py"
    if not array_path.exists():
        return methods

    source = array_path.read_text()
    tree = ast.parse(source)

    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            if node.name.startswith("_dt_"):
                method_name = node.name
                public_name = DATETIME_METHOD_NAME_MAP.get(
                    method_name, method_name.replace("_dt_", "")
                )

                method_source = ast.get_source_segment(source, node)
                if method_source:
                    arrow_func = _extract_datetime_arrow_function(method_source)
                    fallback_type = _determine_datetime_fallback_type(
                        method_source, arrow_func
                    )
                    notes = _extract_datetime_method_notes(method_source, fallback_type)

                    methods.append(
                        MethodInfo(
                            name=method_name,
                            public_name=public_name,
                            category=MethodCategory.DATETIME,
                            fallback_type=fallback_type,
                            arrow_function=arrow_func,
                            notes=notes,
                            performance_notes=_get_performance_note(fallback_type),
                        )
                    )

    return methods


def _determine_datetime_fallback_type(
    source: str, arrow_func: str | None
) -> FallbackType:
    """Determine the fallback type for a datetime method based on source patterns."""
    # Check for _apply_elementwise which indicates Python element-wise processing
    if "_apply_elementwise" in source:
        return FallbackType.ELEMENTWISE

    # Check for version gating
    if "pa_version_under" in source:
        return FallbackType.VERSION_GATED

    # Check for patterns that indicate conversion to Python/NumPy
    conversion_patterns = [
        "_to_timedeltaarray()",
        "_to_datetimearray()",
        "to_numpy(",
        ".to_numpy()",
        "to_pydatetime(",
        ".to_pydatetime()",
        "to_pylist()",
    ]

    has_conversion = any(pattern in source for pattern in conversion_patterns)
    has_pc_call = _has_arrow_native_pattern(source) or arrow_func is not None

    if has_conversion and has_pc_call:
        # Determine if Arrow compute happens before or after conversion
        # by comparing first occurrence positions
        first_pc_pos = source.find("pc.")
        if first_pc_pos == -1 and arrow_func:
            # arrow_func was extracted but pc. not found directly
            first_pc_pos = 0  # Assume it's at the start

        first_conv_pos = len(source)  # Start with max position
        for pattern in conversion_patterns:
            pos = source.find(pattern)
            if pos != -1 and pos < first_conv_pos:
                first_conv_pos = pos

        if first_pc_pos < first_conv_pos:
            # Arrow compute happens first, then conversion for output
            # This is still a meaningful Arrow path
            return FallbackType.CONDITIONAL
        else:
            # Conversion happens before compute - pure fallback
            return FallbackType.OBJECT_FALLBACK

    if has_conversion and not has_pc_call:
        # Only conversion, no Arrow compute - pure object fallback
        return FallbackType.OBJECT_FALLBACK

    # Check for PyArrow compute function usage
    if has_pc_call:
        return FallbackType.ARROW_NATIVE

    # Default to Arrow native if we can't determine otherwise
    return FallbackType.ARROW_NATIVE


def _extract_datetime_arrow_function(source: str) -> str | None:
    """Extract the PyArrow function used for a datetime method programmatically."""
    # Extract pc.X function calls from source regardless of conversion patterns
    # The fallback type determination will handle whether this is meaningful
    pc_calls = re.findall(r"pc\.(\w+)\s*\(", source)
    if pc_calls:
        # Get unique calls, preserving order
        unique_calls = list(dict.fromkeys(pc_calls))
        # Filter out type-checking functions
        significant = [c for c in unique_calls if not c.startswith("is_")]
        if significant:
            if len(significant) == 1:
                return f"pc.{significant[0]}"
            # Return most significant (first non-cast function, or first)
            non_cast = [c for c in significant if c != "cast"]
            if non_cast:
                return f"pc.{non_cast[0]}"
            return f"pc.{significant[0]}"
        elif unique_calls:
            return f"pc.{unique_calls[0]}"

    # Check for ChunkedArray method calls
    if ".cast(" in source:
        return "ChunkedArray.cast"

    return None


def _extract_datetime_method_notes(
    source: str, fallback_type: FallbackType | None = None
) -> str | None:
    """Extract notes for datetime methods programmatically from source patterns."""
    notes = []

    # Add fallback flavor description for object fallback methods
    if fallback_type == FallbackType.OBJECT_FALLBACK:
        flavor = _determine_fallback_flavor(source)
        notes.append(_get_fallback_flavor_description(flavor))
    elif fallback_type == FallbackType.ELEMENTWISE:
        notes.append("element-wise loop")

    # Fallback conversion patterns (additional detail beyond flavor)
    fallback_patterns = {
        "_to_timedeltaarray()": "via TimedeltaArray",
        "_to_datetimearray()": "via DatetimeArray",
        "to_pylist()": "via pylist",
    }
    for pattern, note in fallback_patterns.items():
        if pattern in source:
            notes.append(note)

    # Check for to_numpy (only if not already noted as materialization)
    if "to_numpy(" in source or ".to_numpy()" in source:
        if fallback_type != FallbackType.OBJECT_FALLBACK:
            notes.append("via numpy")

    # Check for element-wise Python operations (conditional case)
    if "_apply_elementwise" in source and fallback_type == FallbackType.CONDITIONAL:
        notes.append("_apply_elementwise() fallback")

    # Check for version gating
    version_patterns = re.findall(r"pa_version_under(\d+)p(\d+)", source)
    if version_patterns:
        major, minor = version_patterns[0]
        notes.append(f"PyArrow >= {_format_pyarrow_version(major, minor)}")

    # Check for specific computation patterns
    if "components." in source:
        # Extract which component is accessed
        comp_match = re.search(r"\.components\.(\w+)", source)
        if comp_match:
            notes.append(f".components.{comp_match.group(1)}")

    return "; ".join(notes) if notes else None


# =============================================================================
# Aggregation Methods Analysis (Series-level only)
# =============================================================================


def analyze_aggregation_methods() -> list[MethodInfo]:
    """
    Analyze ArrowExtensionArray for aggregation method support programmatically.

    NOTE: These are Series-level aggregations (e.g., Series.sum()).
    As of this pandas version, GroupBy aggregations (e.g., df.groupby(...).sum())
    use the standard pandas implementation rather than Arrow compute kernels.
    This may change in future versions.
    """
    methods = []

    array_path = REPO_ROOT / "pandas" / "core" / "arrays" / "arrow" / "array.py"
    if not array_path.exists():
        return methods

    source = array_path.read_text()

    # Find _reduce_pyarrow method and extract the name mapping
    reduce_match = re.search(
        r"def _reduce_pyarrow\([^)]+\).*?(?=\n    def |\nclass |\Z)",
        source,
        re.DOTALL,
    )
    if not reduce_match:
        return methods

    reduce_source = reduce_match.group(0)

    # Extract the name mapping from source: {"median": "quantile", ...}
    name_map_match = re.search(
        r"pyarrow_name\s*=\s*\{([^}]+)\}\.get\(name,\s*name\)",
        reduce_source,
    )
    name_mapping: dict[str, str] = {}
    if name_map_match:
        # Parse the mapping from source dynamically
        map_content = name_map_match.group(1)
        for match in re.finditer(r'"(\w+)":\s*"(\w+)"', map_content):
            name_mapping[match.group(1)] = match.group(2)

    # Extract aggregation method names from the docstring dynamically
    # Pattern: "Name of the function, supported values are: { X, Y, Z }."
    docstring_match = re.search(
        r"supported values are:\s*\{\s*([^}]+)\s*\}",
        reduce_source,
    )
    aggregation_names: list[str] = []
    if docstring_match:
        # Parse comma-separated method names from docstring
        names_str = docstring_match.group(1)
        aggregation_names = [n.strip() for n in names_str.split(",") if n.strip()]

    for name in aggregation_names:
        info = _analyze_aggregation_method(name, reduce_source, name_mapping)
        if info:
            methods.append(info)

    return methods


def _analyze_aggregation_method(
    name: str, reduce_source: str, name_mapping: dict[str, str]
) -> MethodInfo | None:
    """Analyze a single aggregation method programmatically from source patterns."""
    # Get the PyArrow function name from the extracted mapping
    pyarrow_name = name_mapping.get(name, name)

    # Determine arrow function and fallback type based on source patterns
    arrow_func = None
    fallback_type = FallbackType.ARROW_NATIVE
    notes_parts = ["Series-level only"]

    # Check if this method has a GENERAL custom def pyarrow_meth block in source
    # Pattern: if name == "X": followed by def pyarrow_meth (no extra conditions)
    # We exclude patterns like "if name == 'sum' and pa.types.is_string" which are
    # type-specific overrides, not general implementations
    general_custom_pattern = rf'if name == "{name}"\s*:[\s\S]*?def pyarrow_meth'
    has_custom_impl = re.search(general_custom_pattern, reduce_source)

    # Check if method uses type conversion for any/all pattern
    # This is specifically the pattern where any/all converts non-bool to bool
    type_conversion_pattern = r'if name in \["any", "all"\].*?not_equal'
    uses_type_conversion = name in ("any", "all") and re.search(
        type_conversion_pattern, reduce_source, re.DOTALL
    )

    # Check for getattr(pc, ..., None) pattern - indicates runtime availability check
    uses_getattr_check = "getattr(pc," in reduce_source and "None)" in reduce_source

    if has_custom_impl:
        # Method has custom implementation - extract pc.X calls from the if block
        # Find the complete if block for this method
        block_end = r"(?=\n        elif|\n        else:|\n    def )"
        block_pattern = rf'if name == "{name}"[^:]*:.*?{block_end}'
        block_match = re.search(block_pattern, reduce_source, re.DOTALL)
        if block_match:
            custom_block = block_match.group(0)
            # Look for pc.X calls in the custom block
            pc_calls = re.findall(r"pc\.(\w+)", custom_block)
            if pc_calls:
                unique_calls = list(dict.fromkeys(pc_calls))
                arrow_func = " / ".join(f"pc.{c}" for c in unique_calls)
                notes_parts.append(f"computed from {', '.join(unique_calls)}")
            else:
                arrow_func = f"pc.{pyarrow_name}"
        else:
            arrow_func = f"pc.{pyarrow_name}"
    elif uses_type_conversion:
        # Method converts types before applying (like any/all for non-bool)
        arrow_func = f"pc.{pyarrow_name}"
        fallback_type = FallbackType.CONDITIONAL
        notes_parts.append("converts non-bool to bool first")
    elif uses_getattr_check:
        # Method uses runtime check - may be version-gated or not implemented
        # Check test file for version requirements
        test_path = REPO_ROOT / "pandas" / "tests" / "extension" / "test_arrow.py"
        if test_path.exists():
            test_source = test_path.read_text()
            # Look for version guard patterns
            version_pattern = rf"pa_version_under(\d+)p(\d+).*{name}"
            version_match = re.search(version_pattern, test_source)
            if version_match:
                major = version_match.group(1)
                minor = version_match.group(2)
                arrow_func = f"pc.{pyarrow_name}"
                fallback_type = FallbackType.VERSION_GATED
                notes_parts.append(f"requires PyArrow >= {major}.{minor}")
            else:
                # No version guard - check if PyArrow has the function
                try:
                    import pyarrow.compute as pc

                    if hasattr(pc, pyarrow_name):
                        arrow_func = f"pc.{pyarrow_name}"
                    else:
                        arrow_func = None
                        fallback_type = FallbackType.NOT_IMPLEMENTED
                        notes_parts = ["Not supported for Arrow arrays"]
                except ImportError:
                    arrow_func = None
                    fallback_type = FallbackType.NOT_IMPLEMENTED
                    notes_parts = ["Not supported for Arrow arrays"]
        else:
            arrow_func = f"pc.{pyarrow_name}"
    else:
        # Standard case - uses pc.{pyarrow_name}
        arrow_func = f"pc.{pyarrow_name}"

    # Check for type-specific handling (data preparation, not fallback behavior)
    type_handling_pattern = rf'name in \[.*?"{name}".*?\].*?(?:is_duration|is_temporal)'
    if re.search(type_handling_pattern, reduce_source, re.DOTALL):
        notes_parts.append("special handling for temporal types")

    notes = "; ".join(notes_parts) if notes_parts else None

    return MethodInfo(
        name=name,
        public_name=name,
        category=MethodCategory.AGGREGATION,
        fallback_type=fallback_type,
        arrow_function=arrow_func,
        notes=notes,
        performance_notes=_get_performance_note(fallback_type),
    )


# =============================================================================
# List Accessor Analysis - AST-based
# =============================================================================

# Public name mappings for list accessor methods
LIST_METHOD_NAME_MAP = {
    "len": "len()",
    "__getitem__": "[index] / [slice]",
    "flatten": "flatten()",
}


def analyze_list_accessor() -> list[MethodInfo]:
    """Analyze ListAccessor for Arrow list methods using AST."""
    methods = []

    accessor_path = REPO_ROOT / "pandas" / "core" / "arrays" / "arrow" / "accessors.py"
    if not accessor_path.exists():
        return methods

    source = accessor_path.read_text()
    tree = ast.parse(source)

    # Find the ListAccessor class
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == "ListAccessor":
            for item in node.body:
                if isinstance(item, ast.FunctionDef):
                    method_name = item.name
                    is_public = not method_name.startswith("_")
                    if method_name in LIST_METHOD_NAME_MAP or is_public:
                        method_source = ast.get_source_segment(source, item)
                        if method_source:
                            public_name = LIST_METHOD_NAME_MAP.get(
                                method_name, f"{method_name}()"
                            )
                            info = _analyze_accessor_method(
                                method_name,
                                public_name,
                                method_source,
                                MethodCategory.LIST,
                            )
                            if info:
                                methods.append(info)

    return methods


# =============================================================================
# Struct Accessor Analysis - AST-based
# =============================================================================

# Public name mappings for struct accessor methods
STRUCT_METHOD_NAME_MAP = {
    "dtypes": "dtypes (property)",
    "field": "field(name_or_index)",
    "explode": "explode()",
}


def analyze_struct_accessor() -> list[MethodInfo]:
    """Analyze StructAccessor for Arrow struct methods using AST."""
    methods = []

    accessor_path = REPO_ROOT / "pandas" / "core" / "arrays" / "arrow" / "accessors.py"
    if not accessor_path.exists():
        return methods

    source = accessor_path.read_text()
    tree = ast.parse(source)

    # Find the StructAccessor class
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == "StructAccessor":
            for item in node.body:
                if isinstance(item, ast.FunctionDef):
                    method_name = item.name
                    is_public = not method_name.startswith("_")
                    if method_name in STRUCT_METHOD_NAME_MAP or is_public:
                        method_source = ast.get_source_segment(source, item)
                        if method_source:
                            public_name = STRUCT_METHOD_NAME_MAP.get(
                                method_name, f"{method_name}()"
                            )
                            info = _analyze_accessor_method(
                                method_name,
                                public_name,
                                method_source,
                                MethodCategory.STRUCT,
                            )
                            if info:
                                methods.append(info)

    return methods


def _analyze_accessor_method(
    method_name: str, public_name: str, source: str, category: MethodCategory
) -> MethodInfo | None:
    """Analyze an accessor method's source to determine Arrow support."""
    # Skip internal/validation methods
    if method_name.startswith("_") and method_name not in ("__getitem__", "__iter__"):
        return None
    if method_name == "__iter__":
        return None  # Not a real data method

    arrow_func = _extract_arrow_function_for_accessor(source)
    notes = _extract_accessor_notes(source, method_name)

    # Determine fallback type
    if "pc." in source or "pa.compute" in source:
        fallback_type = FallbackType.ARROW_NATIVE
        # Check for conditional patterns
        if "isinstance(key, int)" in source and "isinstance(key, slice)" in source:
            fallback_type = FallbackType.CONDITIONAL
            if not notes:
                notes = "Different Arrow functions for int vs slice"
    elif "._pa_array" in source:
        fallback_type = FallbackType.ARROW_NATIVE
    else:
        # Accessors are generally Arrow-native
        fallback_type = FallbackType.ARROW_NATIVE

    return MethodInfo(
        name=method_name,
        public_name=public_name,
        category=category,
        fallback_type=fallback_type,
        arrow_function=arrow_func,
        notes=notes,
        performance_notes=_get_performance_note(fallback_type),
        github_issues=_extract_github_issues(source),
    )


def _extract_arrow_function_for_accessor(source: str) -> str | None:
    """Extract the PyArrow function used for an accessor method programmatically."""
    # Extract pc.X function calls from source
    pc_calls = re.findall(r"pc\.(\w+)\s*\(", source)
    if pc_calls:
        unique_calls = list(dict.fromkeys(pc_calls))
        # Filter out utility functions
        utility_funcs = {"is_null", "if_else", "invert"}
        main_funcs = [c for c in unique_calls if c not in utility_funcs]

        if main_funcs:
            if len(main_funcs) == 1:
                return f"pc.{main_funcs[0]}"
            # Multiple functions - return them joined
            return " / ".join(f"pc.{f}" for f in main_funcs[:2])
        elif unique_calls:
            return f"pc.{unique_calls[0]}"

    # Check for pa.compute.X calls - prefer the more specific function
    pa_compute_calls = re.findall(r"pa\.compute\.(\w+)\s*\(", source)
    if pa_compute_calls:
        unique_pa_calls = list(dict.fromkeys(pa_compute_calls))
        # Filter out utility/helper functions to get the main one
        utility_funcs = {"list_value_length"}
        main_funcs = [c for c in unique_pa_calls if c not in utility_funcs]
        if main_funcs:
            return f"pa.compute.{main_funcs[0]}"
        return f"pa.compute.{unique_pa_calls[0]}"

    return None


def _extract_accessor_notes(source: str, method_name: str) -> str | None:
    """Extract notes from accessor method source programmatically."""
    notes = []

    # Check for TODO comments indicating limitations
    if "# TODO" in source:
        if "negative" in source.lower():
            notes.append("Negative indices not supported")

    # Check for nested/list-like support
    if "is_list_like" in source:
        notes.append("Supports nested field access via list")

    # Check for DataFrame return
    if "concat" in source and "DataFrame" in source:
        notes.append("Returns all fields as DataFrame")

    # Check for property-style methods (dtypes only - returns type info)
    if method_name == "dtypes":
        notes.append("Returns dtype info (property)")

    # Check for isinstance branching
    if "isinstance(key, int)" in source and "isinstance(key, slice)" in source:
        notes.append("Different handling for int vs slice")

    # Check for GH issue references
    gh_matches = re.findall(r"GH#?(\d+)", source)
    if gh_matches:
        notes.append(f"Related: GH#{gh_matches[0]}")

    return "; ".join(notes) if notes else None


# =============================================================================
# Array Methods Analysis (General Series/Array operations) - AST-based
# =============================================================================

# Methods to analyze from ArrowExtensionArray with their public names
ARRAY_METHOD_NAME_MAP = {
    "argsort": "argsort()",
    "dropna": "dropna()",
    "fillna": "fillna(value)",
    "round": "round(decimals)",
    "unique": "unique()",
    "value_counts": "value_counts()",
    "duplicated": "duplicated()",
    "interpolate": "interpolate()",
    "isna": "isna() / isnull()",
    "_accumulate": "cumsum/cumprod/cummin/cummax",
    "_reduce": "sum/mean/min/max/etc (internal)",
    "_pad_or_backfill": "ffill/bfill",
    "factorize": "factorize()",
    "searchsorted": "searchsorted(value)",
    "take": "take(indices)",
    "copy": "copy()",
    "view": "view(dtype)",
    "_cmp_method": "comparison operators",
    "_arith_method": "arithmetic operators",
}


def analyze_array_methods() -> list[MethodInfo]:
    """
    Analyze ArrowExtensionArray for general array method support using AST.

    Reads the source code of array.py and analyzes each method to determine
    whether it uses PyArrow compute functions or falls back to NumPy/Python.
    """
    methods = []

    array_path = REPO_ROOT / "pandas" / "core" / "arrays" / "arrow" / "array.py"
    if not array_path.exists():
        return methods

    source = array_path.read_text()
    tree = ast.parse(source)

    # Find the ArrowExtensionArray class
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == "ArrowExtensionArray":
            for item in node.body:
                if isinstance(item, ast.FunctionDef):
                    method_name = item.name
                    if method_name in ARRAY_METHOD_NAME_MAP:
                        method_source = ast.get_source_segment(source, item)
                        if method_source:
                            info = _analyze_array_method_source(
                                method_name,
                                ARRAY_METHOD_NAME_MAP[method_name],
                                method_source,
                            )
                            if info:
                                methods.append(info)

    # Add cumulative methods (they go through _accumulate)
    methods.extend(_get_cumulative_methods(source))

    # Add null handling methods derived from _pad_or_backfill
    methods.extend(_get_fill_methods(source))

    return methods


def _analyze_array_method_source(
    method_name: str, public_name: str, source: str
) -> MethodInfo | None:
    """Analyze an array method's source to determine fallback behavior."""
    # Skip internal dispatch methods - they're analyzed separately
    if method_name in (
        "_accumulate",
        "_reduce",
        "_pad_or_backfill",
        "_cmp_method",
        "_arith_method",
    ):
        return None

    arrow_func = _extract_arrow_function_for_array(source)
    fallback_type = _determine_fallback_type_from_source(source, arrow_func)
    conditions = _extract_array_fallback_conditions(source)
    notes = _extract_array_method_notes(source, method_name, fallback_type)

    return MethodInfo(
        name=method_name,
        public_name=public_name,
        category=MethodCategory.ARRAY_METHODS,
        fallback_type=fallback_type,
        arrow_function=arrow_func,
        fallback_conditions=conditions,
        notes=notes,
        performance_notes=_get_performance_note(fallback_type),
        github_issues=_extract_github_issues(source),
    )


def _determine_fallback_type_from_source(
    source: str, arrow_func: str | None
) -> FallbackType:
    """Determine the fallback type based on source code patterns."""
    # Check for conditional fallback patterns
    has_super_call = "super()." in source or "super(" in source
    has_to_numpy = "to_numpy(" in source or ".to_numpy()" in source
    has_limit_check = "if limit" in source or "limit is not None" in source

    # Check for PyArrow native patterns (not just pc.*)
    # Note: ._pa_array.type access doesn't mean Arrow-native - it's just type checking
    has_pa_compute = "pc." in source
    has_pa_array_result_method = (
        ".is_null()" in source and "return" in source.split(".is_null()")[1][:50]
    )
    has_chunked_array = "ChunkedArray" in source or "pa.chunked_array" in source

    # If explicitly marked as no Arrow function and uses to_numpy, it's object fallback
    if arrow_func is None and has_to_numpy:
        return FallbackType.OBJECT_FALLBACK

    if arrow_func:
        if has_super_call or has_limit_check:
            return FallbackType.CONDITIONAL
        return FallbackType.ARROW_NATIVE
    elif has_pa_compute or has_pa_array_result_method:
        if has_super_call or has_limit_check:
            return FallbackType.CONDITIONAL
        return FallbackType.ARROW_NATIVE
    elif has_to_numpy:
        return FallbackType.OBJECT_FALLBACK
    elif has_super_call:
        return FallbackType.CONDITIONAL
    elif has_chunked_array:
        return FallbackType.ARROW_NATIVE
    else:
        return FallbackType.ARROW_NATIVE


def _extract_array_fallback_conditions(source: str) -> list[str]:
    """Extract conditions that trigger fallback for array methods."""
    conditions = []

    if "if limit is not None" in source or "limit is not None" in source:
        conditions.append("limit parameter specified")

    if "limit_area" in source:
        conditions.append("limit_area parameter specified")

    if "is_string_dtype" in source:
        conditions.append("string dtype (uses different path)")

    if "pa.types.is_temporal" in source:
        conditions.append("temporal types have special handling")

    if "pa.types.is_duration" in source:
        conditions.append("duration types have special handling")

    if "super()." in source:
        conditions.append("falls back to base implementation for some cases")

    return conditions


def _extract_array_method_notes(
    source: str, method_name: str, fallback_type: FallbackType | None = None
) -> str | None:
    """Extract relevant notes from method source."""
    notes = []

    # Add fallback flavor description for object fallback methods
    if fallback_type == FallbackType.OBJECT_FALLBACK:
        flavor = _determine_fallback_flavor(source)
        notes.append(_get_fallback_flavor_description(flavor))
    elif fallback_type == FallbackType.ELEMENTWISE:
        notes.append("element-wise loop")

    if "GH#" in source or "GH " in source:
        # Extract issue references
        issues = _extract_github_issues(source)
        if issues:
            notes.append(f"GH{issues[0]}")

    # Method-specific notes based on code patterns
    if method_name == "isna":
        if "null_count == 0" in source:
            notes.append("fast null path")

    if method_name == "duplicated":
        if "to_numpy" in source and fallback_type != FallbackType.OBJECT_FALLBACK:
            notes.append("via numpy")

    if method_name == "searchsorted":
        if "to_numpy" in source and fallback_type != FallbackType.OBJECT_FALLBACK:
            notes.append("via numpy")

    return "; ".join(notes) if notes else None


def _extract_arrow_function_for_array(source: str) -> str | None:
    """Extract the PyArrow function used for an array method programmatically."""
    # Check for fallback patterns first - if present, this method doesn't use Arrow
    fallback_indicators = [
        "to_numpy(",
        ".to_numpy()",
    ]
    # Only return None if the method ONLY uses fallback (no Arrow-native patterns)
    has_fallback = any(ind in source for ind in fallback_indicators)
    has_arrow_path = _has_arrow_native_pattern(source)

    if has_fallback and not has_arrow_path:
        return None

    # Extract pc.X function calls from source
    pc_calls = re.findall(r"pc\.(\w+)\s*\(", source)
    if pc_calls:
        unique_calls = list(dict.fromkeys(pc_calls))
        # Filter out utility/type-checking functions
        utility_funcs = {"is_null", "is_valid", "if_else", "invert"}
        main_funcs = [c for c in unique_calls if c not in utility_funcs]

        if main_funcs:
            if len(main_funcs) == 1:
                return f"pc.{main_funcs[0]}"
            # Multiple functions - return the most significant ones
            return " / ".join(f"pc.{f}" for f in main_funcs[:2])
        elif unique_calls:
            return f"pc.{unique_calls[0]}"

    # Check for ChunkedArray method calls
    chunked_methods = re.findall(r"\.(\w+)\s*\(", source)
    arrow_chunked_methods = {
        "is_null",
        "is_valid",
        "value_counts",
        "combine_chunks",
        "view",
        "cast",
        "to_pylist",
    }
    for method in chunked_methods:
        if method in arrow_chunked_methods:
            return f"ChunkedArray.{method}"

    return None


def _get_cumulative_methods(source: str) -> list[MethodInfo]:
    """Generate MethodInfo for cumulative operations based on _accumulate."""
    methods = []

    # Extract the pyarrow_name mapping from source
    # Pattern: {"cummax": "cumulative_max", "cummin": "cumulative_min", ...}
    mapping_match = re.search(
        r"pyarrow_name\s*=\s*\{([^}]+)\}\.get\(name,\s*name\)",
        source,
    )

    cumulative_ops = {}
    if mapping_match:
        map_content = mapping_match.group(1)
        # Parse entries like "cummax": "cumulative_max"
        for match in re.finditer(r'"(\w+)":\s*"(\w+)"', map_content):
            op_name = match.group(1)
            pyarrow_func = match.group(2)
            cumulative_ops[op_name] = (f"{op_name}()", f"pc.{pyarrow_func}")
    else:
        # Fallback to hardcoded if pattern not found
        cumulative_ops = {
            "cumsum": ("cumsum()", "pc.cumulative_sum_checked"),
            "cumprod": ("cumprod()", "pc.cumulative_prod_checked"),
            "cummin": ("cummin()", "pc.cumulative_min"),
            "cummax": ("cummax()", "pc.cumulative_max"),
        }

    # Check if string dtype has special handling
    has_string_fallback = "is_string_dtype" in source and "_str_accumulate" in source

    for name, (public_name, arrow_func) in cumulative_ops.items():
        notes = None
        if has_string_fallback:
            if name == "cumprod":
                notes = "Not supported for string dtype"
            else:
                notes = "String dtype uses NumPy fallback"

        methods.append(
            MethodInfo(
                name=name,
                public_name=public_name,
                category=MethodCategory.ARRAY_METHODS,
                fallback_type=FallbackType.CONDITIONAL,
                arrow_function=arrow_func,
                notes=notes,
                performance_notes=_get_performance_note(FallbackType.CONDITIONAL),
            )
        )

    return methods


def _get_fill_methods(source: str) -> list[MethodInfo]:
    """Generate MethodInfo for ffill/bfill based on _pad_or_backfill."""
    methods = []

    # Check what conditions trigger fallback by parsing the source
    # Pattern: if limit is None and limit_area is None: (use Arrow)
    has_limit_check = "limit is None" in source or "limit is not None" in source
    has_limit_area_check = (
        "limit_area is None" in source or "limit_area is not None" in source
    )

    conditions = []
    if has_limit_check:
        conditions.append("limit parameter specified")
    if has_limit_area_check:
        conditions.append("limit_area parameter specified")

    notes = f"Falls back if {' or '.join(conditions)}" if conditions else None

    # Extract the PyArrow functions used from source
    # Pattern: pc.fill_null_forward, pc.fill_null_backward
    fill_ops = {}

    # Look for pad/ffill -> fill_null_forward
    if "pc.fill_null_forward" in source:
        fill_ops["ffill"] = ("ffill() / pad()", "pc.fill_null_forward")

    # Look for backfill/bfill -> fill_null_backward
    if "pc.fill_null_backward" in source:
        fill_ops["bfill"] = ("bfill() / backfill()", "pc.fill_null_backward")

    # Fallback if pattern not found
    if not fill_ops:
        fill_ops = {
            "ffill": ("ffill() / pad()", "pc.fill_null_forward"),
            "bfill": ("bfill() / backfill()", "pc.fill_null_backward"),
        }

    for name, (public_name, arrow_func) in fill_ops.items():
        methods.append(
            MethodInfo(
                name=name,
                public_name=public_name,
                category=MethodCategory.ARRAY_METHODS,
                fallback_type=FallbackType.CONDITIONAL
                if conditions
                else FallbackType.ARROW_NATIVE,
                arrow_function=arrow_func,
                fallback_conditions=conditions,
                notes=notes,
                performance_notes=_get_performance_note(
                    FallbackType.CONDITIONAL
                    if conditions
                    else FallbackType.ARROW_NATIVE
                ),
            )
        )

    return methods


# =============================================================================
# Helper Functions
# =============================================================================


def _extract_arrow_function(source: str) -> str | None:
    """Extract the PyArrow compute function name from source."""
    matches = re.findall(r"pc\.(\w+)", source)
    if matches:
        significant = [
            m for m in matches if not m.startswith("is_") and m not in ("if_else",)
        ]
        return significant[0] if significant else matches[0]
    return None


def _extract_github_issues(source: str) -> list[str]:
    """Extract GitHub issue references from source."""
    matches = re.findall(r"GH#?(\d+)", source)
    return [f"#{m}" for m in matches]


class FallbackFlavor(Enum):
    """
    More specific classification for OBJECT_FALLBACK methods.

    These represent different fallback mechanisms with distinct performance
    characteristics:

    OBJECT_STRING_MIXIN:
        Uses super()._str_*() to delegate to ObjectStringArrayMixin methods.
        The ObjectStringArrayMixin uses lib.map_infer_mask() which iterates
        over elements and applies Python string methods.
        Fallback method: ObjectStringArrayMixin._str_map() in object_array.py

    NUMPY_MATERIALIZATION:
        Converts Arrow array to NumPy via to_numpy(), to_pylist(),
        _to_datetimearray(), or _to_timedeltaarray() before processing.
        Fallback methods: ChunkedArray.to_numpy(), ChunkedArray.to_pylist()

    PYTHON_LOOP:
        Uses _apply_elementwise(func) which iterates over each chunk,
        calls to_numpy() on it, then applies a Python function to each element.
        Fallback method: ArrowExtensionArray._apply_elementwise() in array.py
    """

    # Delegates via super() to ObjectStringArrayMixin._str_map()
    OBJECT_STRING_MIXIN = "object_string_mixin"
    # Calls to_numpy()/to_pylist()/_to_datetimearray()/_to_timedeltaarray()
    NUMPY_MATERIALIZATION = "numpy_materialization"
    # Uses _apply_elementwise(func) for per-element Python function calls
    PYTHON_LOOP = "python_loop"
    # Could not determine specific flavor
    UNKNOWN = "unknown"


def _determine_fallback_flavor(source: str) -> FallbackFlavor:
    """
    Determine the specific flavor of object fallback from source code patterns.

    Detection patterns:
    - PYTHON_LOOP: Uses _apply_elementwise(func) which applies a Python
      callable to each element. Found in: _str_partition, _str_casefold,
      _str_encode, _str_findall, _str_normalize, _str_wrap, _str_translate
    - NUMPY_MATERIALIZATION: Converts to NumPy/Python via to_numpy(),
      to_pylist(), _to_datetimearray(), _to_timedeltaarray(). Found in:
      datetime component accessors, value_counts, _cmp_method
    - OBJECT_STRING_MIXIN: Calls super()._str_*() which delegates to
      ObjectStringArrayMixin. Found in: _str_contains (with flags),
      _str_match (with flags), _str_replace (with callable/flags),
      _str_count (with flags), _str_repeat (with array-like)
    """
    has_super_call = "super()." in source or "super(" in source
    has_to_numpy = "to_numpy(" in source or ".to_numpy()" in source
    has_elementwise = "_apply_elementwise" in source

    # Check for conversion patterns (datetime/timedelta specific)
    conversion_patterns = [
        "_to_timedeltaarray()",
        "_to_datetimearray()",
        "to_pylist()",
    ]
    has_conversion = any(pattern in source for pattern in conversion_patterns)

    # Priority order: elementwise is most specific, then numpy/conversion,
    # then super()
    if has_elementwise:
        return FallbackFlavor.PYTHON_LOOP
    elif has_to_numpy or has_conversion:
        return FallbackFlavor.NUMPY_MATERIALIZATION
    elif has_super_call:
        return FallbackFlavor.OBJECT_STRING_MIXIN
    else:
        return FallbackFlavor.UNKNOWN


def _get_fallback_flavor_description(flavor: FallbackFlavor) -> str:
    """
    Get a human-readable description for a fallback flavor.

    These descriptions indicate the actual fallback mechanism used:
    - OBJECT_STRING_MIXIN: Delegates to ObjectStringArrayMixin._str_map()
      which uses lib.map_infer_mask() for element-wise Python string ops
    - NUMPY_MATERIALIZATION: Calls to_numpy()/to_pylist() to materialize
      the Arrow array before processing with NumPy/Python
    - PYTHON_LOOP: Uses _apply_elementwise(func) to apply a Python callable
      per element while maintaining chunked structure
    """
    descriptions = {
        FallbackFlavor.OBJECT_STRING_MIXIN: "via ObjectStringArrayMixin",
        FallbackFlavor.NUMPY_MATERIALIZATION: "via to_numpy()",
        FallbackFlavor.PYTHON_LOOP: "via _apply_elementwise()",
        FallbackFlavor.UNKNOWN: "Object fallback",
    }
    return descriptions.get(flavor, "Object fallback")


def _get_performance_note(fallback_type: FallbackType) -> str:
    """Get performance note for a method based on its fallback type."""
    notes = {
        FallbackType.ARROW_NATIVE: (
            "Optimized: Uses PyArrow compute kernel (vectorized)"
        ),
        FallbackType.CONDITIONAL: (
            "Mixed: Uses Arrow when possible, falls back for edge cases"
        ),
        FallbackType.ELEMENTWISE: "Slow: Element-by-element Python iteration",
        FallbackType.OBJECT_FALLBACK: (
            "Slow: Converts to object dtype, uses Python iteration"
        ),
        FallbackType.VERSION_GATED: (
            "Version-dependent: Newer PyArrow versions use optimized kernel"
        ),
        FallbackType.NOT_IMPLEMENTED: (
            "Not supported: Raises NotImplementedError or TypeError"
        ),
    }
    return notes.get(fallback_type, "")


# =============================================================================
# Build Tables
# =============================================================================


def build_string_method_table() -> dict[str, MethodInfo]:
    """Build complete table of string methods with their Arrow support status."""
    methods = {}

    # 1. First add methods from ArrowStringArrayMixin (_arrow_string_mixins.py)
    for info in analyze_arrow_string_mixin():
        methods[info.name] = info

    # 2. Add methods from ArrowStringArray (string_arrow.py)
    #    Track which methods are actually implemented in ArrowStringArray
    arrow_string_array_methods = set()
    for info in analyze_arrow_string_array():
        arrow_string_array_methods.add(info.name)
        if info.name in methods:
            existing = methods[info.name]
            existing.fallback_conditions.extend(info.fallback_conditions)
            # Dedupe while preserving order for stable output
            existing.fallback_conditions = _dedupe_preserve_order(
                existing.fallback_conditions
            )
            existing.fallback_type = FallbackType.CONDITIONAL
        else:
            methods[info.name] = info

    # 3. Add string methods from ArrowExtensionArray (array.py)
    #    These are typically methods using _apply_elementwise()
    #    IMPORTANT: Methods only in ArrowExtensionArray (not in ArrowStringArray)
    #    should be classified as OBJECT_FALLBACK for dtype="string[pyarrow]" users
    for info in analyze_arrow_extension_array_string_methods():
        if info.name in methods:
            existing = methods[info.name]
            # Merge conditions and keep the more specific fallback type
            existing.fallback_conditions.extend(info.fallback_conditions)
            existing.fallback_conditions = _dedupe_preserve_order(
                existing.fallback_conditions
            )
            # Update notes if the new info has more specific notes
            if info.notes and (
                not existing.notes or "Object fallback" in existing.notes
            ):
                existing.notes = info.notes
        else:
            # Method only exists in ArrowExtensionArray, not ArrowStringArray
            # For dtype="string[pyarrow]", this falls back to ObjectStringArrayMixin
            if info.name not in arrow_string_array_methods:
                info.fallback_type = FallbackType.OBJECT_FALLBACK
                info.notes = (
                    "Falls back for dtype='string[pyarrow]'; "
                    "ArrowDtype uses _apply_elementwise()"
                )
                info.performance_notes = _get_performance_note(
                    FallbackType.OBJECT_FALLBACK
                )
            methods[info.name] = info

    # 4. Finally add methods from ObjectStringArrayMixin (object_array.py)
    #    These are pure Python fallback methods
    object_methods = analyze_object_string_mixin()
    for info in object_methods:
        if info.name not in methods:
            methods[info.name] = info

    return methods


def build_all_methods_table() -> dict[str, dict[str, MethodInfo]]:
    """Build complete table of all methods organized by category."""
    return {
        "string": build_string_method_table(),
        "arithmetic": {m.name: m for m in analyze_arithmetic_operations()},
        "datetime": {m.name: m for m in analyze_datetime_methods()},
        "aggregation": {m.name: m for m in analyze_aggregation_methods()},
        "array_methods": {m.name: m for m in analyze_array_methods()},
        "list": {m.name: m for m in analyze_list_accessor()},
        "struct": {m.name: m for m in analyze_struct_accessor()},
    }


# =============================================================================
# Output Formatters
# =============================================================================


def _generate_summary_stats(all_methods: dict[str, dict[str, MethodInfo]]) -> list[str]:
    """Generate summary statistics section for the RST output."""
    # Calculate stats per category
    category_stats: dict[str, dict[str, int]] = {}
    total_stats: dict[str, int] = {
        "arrow": 0,
        "conditional": 0,
        "version_gated": 0,
        "elementwise": 0,
        "object": 0,
        "not_implemented": 0,
        "total": 0,
    }

    category_display_names = {
        "string": "String methods",
        "arithmetic": "Arithmetic & comparison",
        "datetime": "Datetime methods",
        "aggregation": "Aggregation methods",
        "array_methods": "Array methods",
        "list": "List accessor",
        "struct": "Struct accessor",
    }

    for cat_key, methods in all_methods.items():
        if not methods:
            continue

        stats = {
            "arrow": 0,
            "conditional": 0,
            "version_gated": 0,
            "elementwise": 0,
            "object": 0,
            "not_implemented": 0,
            "total": 0,
        }

        for method_info in methods.values():
            stats["total"] += 1
            total_stats["total"] += 1

            if method_info.fallback_type == FallbackType.ARROW_NATIVE:
                stats["arrow"] += 1
                total_stats["arrow"] += 1
            elif method_info.fallback_type == FallbackType.CONDITIONAL:
                stats["conditional"] += 1
                total_stats["conditional"] += 1
            elif method_info.fallback_type == FallbackType.VERSION_GATED:
                stats["version_gated"] += 1
                total_stats["version_gated"] += 1
            elif method_info.fallback_type == FallbackType.ELEMENTWISE:
                stats["elementwise"] += 1
                total_stats["elementwise"] += 1
            elif method_info.fallback_type == FallbackType.OBJECT_FALLBACK:
                stats["object"] += 1
                total_stats["object"] += 1
            elif method_info.fallback_type == FallbackType.NOT_IMPLEMENTED:
                stats["not_implemented"] += 1
                total_stats["not_implemented"] += 1

        category_stats[cat_key] = stats

    # Build summary section
    # Note: The "Python fallback" column combines elementwise + object fallback types
    lines = [
        "",
        "Summary",
        "=======",
        "",
        ".. list-table::",
        "   :widths: 30 12 12 12 12 12 10",
        "   :header-rows: 1",
        "",
        "   * - Category",
        "     - |arrow|",
        "     - |conditional|",
        "     - |version|",
        "     - |pyfallback|",
        "     - |notimpl|",
        "     - Total",
    ]

    # Add rows for each category
    for cat_key in [
        "string",
        "arithmetic",
        "datetime",
        "aggregation",
        "array_methods",
        "list",
        "struct",
    ]:
        if cat_key not in category_stats:
            continue
        stats = category_stats[cat_key]
        display_name = category_display_names.get(cat_key, cat_key)
        lines.append(f"   * - {display_name}")
        lines.append(f"     - {stats['arrow']}")
        lines.append(f"     - {stats['conditional']}")
        lines.append(f"     - {stats['version_gated']}")
        # Combine elementwise + object for display
        fallback_count = stats["elementwise"] + stats["object"]
        lines.append(f"     - {fallback_count}")
        lines.append(f"     - {stats['not_implemented']}")
        lines.append(f"     - {stats['total']}")

    # Add total row
    lines.append("   * - **Total**")
    lines.append(f"     - **{total_stats['arrow']}**")
    lines.append(f"     - **{total_stats['conditional']}**")
    lines.append(f"     - **{total_stats['version_gated']}**")
    fallback_total = total_stats["elementwise"] + total_stats["object"]
    lines.append(f"     - **{fallback_total}**")
    lines.append(f"     - **{total_stats['not_implemented']}**")
    lines.append(f"     - **{total_stats['total']}**")

    # Calculate and display percentages
    total = total_stats["total"]
    if total > 0:
        arrow_pct = (total_stats["arrow"] / total) * 100
        conditional_pct = (total_stats["conditional"] / total) * 100
        native_pct = arrow_pct + conditional_pct
        arrow_count = total_stats["arrow"]
        cond_count = total_stats["conditional"]
        lines.extend(
            [
                "",
                f"**{arrow_count}** methods ({arrow_pct:.0f}%) use native Arrow.",
                f"**{cond_count}** methods ({conditional_pct:.0f}%) use Arrow "
                "with conditional fallbacks.",
                f"Overall, **{native_pct:.0f}%** have Arrow-native code paths.",
                "",
            ]
        )

    return lines


def format_rst_table(
    all_methods: dict[str, dict[str, MethodInfo]], pyarrow_versions: dict[str, str]
) -> str:
    """Format all methods as RST table for pandas documentation."""
    lines = [
        ".. _arrow-fallbacks:",
        "",
        "{{ header }}",
        "",
        "******************************",
        "Arrow Method Support Reference",
        "******************************",
        "",
        "This document provides a comprehensive reference of which pandas methods",
        "use native PyArrow compute functions vs. falling back to Python/NumPy",
        "implementations when using Arrow-backed arrays.",
        "",
        (
            "Minimum PyArrow version: "
            f"**{pyarrow_versions.get('PYARROW_MIN_VERSION', '13.0.0')}**"
        ),
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
        "     - Uses native PyArrow compute kernel",
        "   * - |conditional|",
        "     - Uses Arrow by default, falls back to Python for certain parameters",
        "   * - |version|",
        "     - Behavior depends on PyArrow version",
        "   * - |pyfallback|",
        "     - Python fallback (combines |elementwise| + |object| in summary table)",
        "   * - |elementwise|",
        "     - Processes elements one-by-one in Python (not vectorized)",
        "   * - |object|",
        "     - Converts to Python objects, then processes element-by-element",
        "   * - |notimpl|",
        "     - Not implemented for Arrow arrays",
        "",
        # Define substitutions once - plain text works for all output formats
        # Note: Sphinx doesn't support conditional substitution definitions well
        ".. |arrow| replace:: Arrow",
        ".. |conditional| replace:: Conditional",
        ".. |version| replace:: Version-gated",
        ".. |pyfallback| replace:: Py fallback",
        ".. |elementwise| replace:: Element-wise",
        ".. |object| replace:: Object fallback",
        ".. |notimpl| replace:: Not implemented",
        "",
    ]

    # Add Summary section
    lines.extend(_generate_summary_stats(all_methods))

    # Add String Methods section
    if all_methods.get("string"):
        lines.extend(
            _format_rst_category_section(
                "String Methods (Series.str.*)", "string", all_methods["string"], "str"
            )
        )

    # Add Arithmetic Operations section
    if all_methods.get("arithmetic"):
        lines.extend(
            _format_rst_category_section(
                "Arithmetic & Comparison Operations",
                "arithmetic",
                all_methods["arithmetic"],
                "",
            )
        )

    # Add Datetime Methods section
    if all_methods.get("datetime"):
        lines.extend(
            _format_rst_category_section(
                "Datetime Methods (Series.dt.*)",
                "datetime",
                all_methods["datetime"],
                "dt",
            )
        )

    # Add Aggregation Methods section
    if all_methods.get("aggregation"):
        lines.extend(
            _format_rst_category_section(
                "Series Aggregation Methods (Series.sum(), etc.)",
                "aggregation",
                all_methods["aggregation"],
                "",
            )
        )
        # Add note about groupby
        lines.extend(
            [
                "",
                ".. note::",
                "",
                "   The aggregation methods above apply to **Series-level** operations",
                "   only (e.g., ``ser.sum()``). As of this pandas version, **GroupBy**",
                "   **aggregations** (``df.groupby(...).sum()``) use the standard",
                "   pandas implementation rather than Arrow compute kernels. This may",
                "   change in future versions.",
                "",
            ]
        )

    # Add Array Methods section (cumulative, null handling, etc.)
    if all_methods.get("array_methods"):
        lines.extend(
            _format_rst_category_section(
                "General Array Methods",
                "array_methods",
                all_methods["array_methods"],
                "",
            )
        )

    # Add List Accessor section
    if all_methods.get("list"):
        lines.extend(
            _format_rst_category_section(
                "List Accessor Methods (Series.list.*)",
                "list",
                all_methods["list"],
                "list",
            )
        )

    # Add Struct Accessor section
    if all_methods.get("struct"):
        lines.extend(
            _format_rst_category_section(
                "Struct Accessor Methods (Series.struct.*)",
                "struct",
                all_methods["struct"],
                "struct",
            )
        )

    return "\n".join(lines) + "\n"


def _format_rst_category_section(
    title: str, category: str, methods: dict[str, MethodInfo], accessor: str
) -> list[str]:
    """Format a category section for RST output."""
    lines = [
        "",
        title,
        "=" * len(title),
        "",
        ".. list-table::",
        "   :widths: 20 12 25 43",
        "   :header-rows: 1",
        "",
        "   * - Method",
        "     - Status",
        "     - Arrow Function",
        "     - Notes",
    ]

    sorted_methods = sorted(methods.values(), key=lambda m: m.public_name)

    for info in sorted_methods:
        status = _get_rst_status(info.fallback_type)
        arrow_func = f"``{info.arrow_function}``" if info.arrow_function else "-"

        # Format notes
        notes_parts = []
        if info.min_pyarrow_version:
            notes_parts.append(f"Requires PyArrow >= {info.min_pyarrow_version}")
        if info.fallback_conditions:
            # Dedupe and sort for stable, readable output
            conditions = sorted(set(info.fallback_conditions))
            cond_str = "; ".join(conditions[:2])
            if len(conditions) > 2:
                cond_str += f" (+{len(conditions) - 2} more)"
            notes_parts.append(f"Falls back: {cond_str}")
        if info.notes:
            notes_parts.append(info.notes)

        notes = " | ".join(notes_parts) if notes_parts else "-"

        # Format method name
        if accessor:
            method_display = f"``{accessor}.{info.public_name}``"
        else:
            method_display = f"``{info.public_name}``"

        lines.append(f"   * - {method_display}")
        lines.append(f"     - {status}")
        lines.append(f"     - {arrow_func}")
        lines.append(f"     - {notes}")

    return lines


def _get_rst_status(fallback_type: FallbackType) -> str:
    """Get RST status indicator for fallback type."""
    status_map = {
        FallbackType.ARROW_NATIVE: "|arrow|",
        FallbackType.CONDITIONAL: "|conditional|",
        FallbackType.VERSION_GATED: "|version|",
        FallbackType.ELEMENTWISE: "|elementwise|",
        FallbackType.OBJECT_FALLBACK: "|object|",
        FallbackType.NOT_IMPLEMENTED: "|notimpl|",
    }
    return status_map.get(fallback_type, "Unknown")


def format_json(
    all_methods: dict[str, dict[str, MethodInfo]], pyarrow_versions: dict[str, str]
) -> str:
    """Format all methods as JSON for programmatic use."""
    # Calculate summary statistics
    summary: dict[str, dict[str, int]] = {}
    totals = {
        "arrow": 0,
        "conditional": 0,
        "version_gated": 0,
        "elementwise": 0,
        "object": 0,
        "not_implemented": 0,
        "total": 0,
    }

    for category, methods in all_methods.items():
        if not methods:
            continue
        cat_stats = {
            "arrow": 0,
            "conditional": 0,
            "version_gated": 0,
            "elementwise": 0,
            "object": 0,
            "not_implemented": 0,
            "total": 0,
        }
        for info in methods.values():
            cat_stats["total"] += 1
            totals["total"] += 1
            status = info.fallback_type.value
            if status in cat_stats:
                cat_stats[status] += 1
                totals[status] += 1
        summary[category] = cat_stats

    # Calculate percentages
    total = totals["total"]
    arrow_pct = (totals["arrow"] / total * 100) if total else 0
    conditional_pct = (totals["conditional"] / total * 100) if total else 0

    data: dict = {
        "pyarrow_min_version": pyarrow_versions.get("PYARROW_MIN_VERSION", "13.0.0"),
        "summary": {
            "by_category": summary,
            "totals": totals,
            "percentages": {
                "arrow_native": round(arrow_pct, 1),
                "conditional": round(conditional_pct, 1),
                "arrow_total": round(arrow_pct + conditional_pct, 1),
            },
        },
        "categories": {},
    }

    for category, methods in all_methods.items():
        data["categories"][category] = {}
        for name, info in sorted(methods.items()):
            data["categories"][category][info.public_name] = {
                "internal_name": info.name,
                "status": info.fallback_type.value,
                "arrow_function": info.arrow_function,
                "min_pyarrow_version": info.min_pyarrow_version,
                "fallback_conditions": info.fallback_conditions,
                "performance_notes": info.performance_notes,
                "github_issues": info.github_issues,
                "notes": info.notes,
            }

    return json.dumps(data, indent=2)


def format_markdown_table(
    all_methods: dict[str, dict[str, MethodInfo]], pyarrow_versions: dict[str, str]
) -> str:
    """Format all methods as Markdown tables."""
    min_version = pyarrow_versions.get("PYARROW_MIN_VERSION", "13.0.0")
    lines = [
        "# Arrow Method Support Reference",
        "",
        f"Minimum PyArrow version: **{min_version}**",
        "",
    ]

    category_titles = {
        "string": "String Methods (Series.str.*)",
        "arithmetic": "Arithmetic & Comparison Operations",
        "datetime": "Datetime Methods (Series.dt.*)",
        "aggregation": "Series Aggregation Methods (Series.sum(), etc.)",
        "array_methods": "General Array Methods",
        "list": "List Accessor Methods (Series.list.*)",
        "struct": "Struct Accessor Methods (Series.struct.*)",
    }

    for category, methods in all_methods.items():
        if not methods:
            continue

        lines.extend(
            [
                f"## {category_titles.get(category, category)}",
                "",
                "| Method | Status | Arrow Function | Notes |",
                "|--------|--------|----------------|-------|",
            ]
        )

        sorted_methods = sorted(methods.values(), key=lambda m: m.public_name)

        for info in sorted_methods:
            status = info.fallback_type.value.replace("_", " ").title()
            arrow_func = f"`{info.arrow_function}`" if info.arrow_function else "-"

            notes_parts = []
            if info.min_pyarrow_version:
                notes_parts.append(f"PA >= {info.min_pyarrow_version}")
            if info.notes:
                notes_parts.append(info.notes)
            notes = "; ".join(notes_parts) if notes_parts else "-"

            lines.append(
                f"| `{info.public_name}` | {status} | {arrow_func} | {notes} |"
            )

        lines.append("")

    return "\n".join(lines)


# =============================================================================
# Main
# =============================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Generate comprehensive Arrow fallback documentation"
    )
    parser.add_argument(
        "--format",
        choices=["rst", "md", "json"],
        default="rst",
        help="Output format (default: rst)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output file path (default: stdout)",
    )
    parser.add_argument(
        "--category",
        choices=[
            "all",
            "string",
            "arithmetic",
            "datetime",
            "aggregation",
            "array_methods",
            "list",
            "struct",
        ],
        default="all",
        help="Which category to include (default: all)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Check if generated docs match existing file (exit 1 if not)",
    )
    args = parser.parse_args()

    # Build method tables
    category_builders = {
        "string": lambda: {"string": build_string_method_table()},
        "arithmetic": lambda: {
            "arithmetic": {m.name: m for m in analyze_arithmetic_operations()}
        },
        "datetime": lambda: {
            "datetime": {m.name: m for m in analyze_datetime_methods()}
        },
        "aggregation": lambda: {
            "aggregation": {m.name: m for m in analyze_aggregation_methods()}
        },
        "array_methods": lambda: {
            "array_methods": {m.name: m for m in analyze_array_methods()}
        },
        "list": lambda: {"list": {m.name: m for m in analyze_list_accessor()}},
        "struct": lambda: {"struct": {m.name: m for m in analyze_struct_accessor()}},
    }

    if args.category == "all":
        all_methods = build_all_methods_table()
    else:
        all_methods = category_builders.get(args.category, build_all_methods_table)()

    pyarrow_versions = get_pyarrow_versions()

    # Format output
    if args.format == "rst":
        output = format_rst_table(all_methods, pyarrow_versions)
    elif args.format == "md":
        output = format_markdown_table(all_methods, pyarrow_versions)
    else:
        output = format_json(all_methods, pyarrow_versions)

    # Check mode: compare with existing file
    if args.check:
        target = (
            REPO_ROOT / "doc" / "source" / "user_guide" / "arrow_string_fallbacks.rst"
        )
        relative_target = "doc/source/user_guide/arrow_string_fallbacks.rst"
        if not target.exists():
            print(f"ERROR: {relative_target} does not exist.", file=sys.stderr)
            print(
                "Run: python scripts/generate_arrow_fallback_table.py > "
                f"{relative_target}",
                file=sys.stderr,
            )
            sys.exit(1)
        current = target.read_text()
        if current != output:
            print(f"ERROR: {relative_target} is out of date.", file=sys.stderr)
            print(
                "Run: python scripts/generate_arrow_fallback_table.py > "
                f"{relative_target}",
                file=sys.stderr,
            )
            sys.exit(1)
        print("Arrow fallback documentation is up to date.")
        sys.exit(0)

    # Write output
    if args.output:
        args.output.write_text(output)
        print(f"Written to {args.output}")
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
