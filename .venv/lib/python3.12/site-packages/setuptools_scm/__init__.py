"""
:copyright: 2010-2023 by Ronny Pfannschmidt
:license: MIT
"""

from __future__ import annotations

from ._config import DEFAULT_LOCAL_SCHEME
from ._config import DEFAULT_VERSION_SCHEME
from ._config import Configuration
from ._get_version_impl import _get_version
from ._get_version_impl import get_version
from ._integration.dump_version import dump_version  # soft deprecated
from ._version_cls import NonNormalizedVersion
from ._version_cls import Version
from .version import ScmVersion

# Public API
__all__ = [
    "DEFAULT_LOCAL_SCHEME",
    "DEFAULT_VERSION_SCHEME",
    "Configuration",
    "NonNormalizedVersion",
    "ScmVersion",
    "Version",
    "_get_version",
    "dump_version",
    # soft deprecated imports, left for backward compatibility
    "get_version",
]
