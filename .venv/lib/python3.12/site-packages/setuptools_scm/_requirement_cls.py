from __future__ import annotations

__all__ = ["Requirement", "extract_package_name"]

try:
    from packaging.requirements import Requirement
    from packaging.utils import canonicalize_name
except ImportError:
    from setuptools.extern.packaging.requirements import (  # type: ignore[import-not-found,no-redef]
        Requirement as Requirement,
    )
    from setuptools.extern.packaging.utils import (  # type: ignore[import-not-found,no-redef]
        canonicalize_name as canonicalize_name,
    )

from . import _log

log = _log.log.getChild("requirement_cls")


def extract_package_name(requirement_string: str) -> str:
    """Extract the canonical package name from a requirement string.

    This function uses packaging.requirements.Requirement to properly parse
    the requirement and extract the package name, handling all edge cases
    that the custom regex-based approach might miss.

    Args:
        requirement_string: The requirement string to parse

    Returns:
        The package name as a string
    """
    return canonicalize_name(Requirement(requirement_string).name)
