"""Compatibility utilities for cross-platform functionality."""

from __future__ import annotations


def normalize_path_for_assertion(path: str) -> str:
    """Normalize path separators for cross-platform assertions.

    On Windows, this converts backslashes to forward slashes to ensure
    path comparisons work correctly. On other platforms, returns the path unchanged.
    The length of the string is not changed by this operation.

    Args:
        path: The path string to normalize

    Returns:
        The path with normalized separators
    """
    return path.replace("\\", "/")


def strip_path_suffix(
    full_path: str, suffix_path: str, error_msg: str | None = None
) -> str:
    """Strip a suffix from a path, with cross-platform path separator handling.

    This function first normalizes path separators for Windows compatibility,
    then asserts that the full path ends with the suffix, and finally returns
    the path with the suffix removed. This is the common pattern used for
    computing parent directories from git output.

    Args:
        full_path: The full path string
        suffix_path: The suffix path to strip from the end
        error_msg: Optional custom error message for the assertion

    Returns:
        The prefix path with the suffix removed

    Raises:
        AssertionError: If the full path doesn't end with the suffix
    """
    normalized_full = normalize_path_for_assertion(full_path)

    if error_msg:
        assert normalized_full.endswith(suffix_path), error_msg
    else:
        assert normalized_full.endswith(suffix_path), (
            f"Path assertion failed: {full_path!r} does not end with {suffix_path!r}"
        )

    return full_path[: -len(suffix_path)]


# Legacy aliases for backward compatibility during transition
def assert_path_endswith(
    full_path: str, suffix_path: str, error_msg: str | None = None
) -> None:
    """Legacy alias - use strip_path_suffix instead."""
    strip_path_suffix(full_path, suffix_path, error_msg)


def compute_path_prefix(full_path: str, suffix_path: str) -> str:
    """Legacy alias - use strip_path_suffix instead."""
    return strip_path_suffix(full_path, suffix_path)
