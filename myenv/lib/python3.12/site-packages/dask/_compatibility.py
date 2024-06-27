from __future__ import annotations

import sys
import warnings

if sys.version_info >= (3, 12):
    import importlib.metadata as importlib_metadata
else:
    import importlib_metadata
from packaging.version import Version

PY_VERSION = Version(".".join(map(str, sys.version_info[:3])))

EMSCRIPTEN = sys.platform == "emscripten"


def entry_points(group=None):
    warnings.warn(
        "`dask._compatibility.entry_points` has been replaced by `importlib_metadata.entry_points` and will be removed "
        "in a future version. Please use `importlib_metadata.entry_points` instead.",
        DeprecationWarning,
        stacklevel=2,
    )
    return importlib_metadata.entry_points(group=group)
