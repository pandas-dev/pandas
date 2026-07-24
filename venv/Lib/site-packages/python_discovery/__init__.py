"""Self-contained Python interpreter discovery."""

from __future__ import annotations

from importlib.metadata import version

from ._cache import ContentStore, DiskCache, PyInfoCache
from ._discovery import get_interpreter, iter_interpreters
from ._py_info import KNOWN_ARCHITECTURES, PythonInfo, normalize_isa
from ._py_spec import KNOWN_IMPLEMENTATIONS, PythonSpec
from ._specifier import SimpleSpecifier, SimpleSpecifierSet, SimpleVersion

__version__ = version("python-discovery")

__all__ = [
    "KNOWN_ARCHITECTURES",
    "KNOWN_IMPLEMENTATIONS",
    "ContentStore",
    "DiskCache",
    "PyInfoCache",
    "PythonInfo",
    "PythonSpec",
    "SimpleSpecifier",
    "SimpleSpecifierSet",
    "SimpleVersion",
    "__version__",
    "get_interpreter",
    "iter_interpreters",
    "normalize_isa",
]
