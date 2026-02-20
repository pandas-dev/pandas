from __future__ import annotations

import os
from typing import Final

# Earliest fully supported Python 3.x version. Used as the default Python
# version in tests. Mypy wheels should be built starting with this version,
# and CI tests should be run on this version (and later versions).
PYTHON3_VERSION: Final = (3, 9)

# Earliest Python 3.x version supported via --python-version 3.x. To run
# mypy, at least version PYTHON3_VERSION is needed.
PYTHON3_VERSION_MIN: Final = (3, 9)  # Keep in sync with typeshed's python support

CACHE_DIR: Final = ".mypy_cache"

CONFIG_NAMES: Final = ["mypy.ini", ".mypy.ini"]
SHARED_CONFIG_NAMES: Final = ["pyproject.toml", "setup.cfg"]

USER_CONFIG_FILES: list[str] = ["~/.config/mypy/config", "~/.mypy.ini"]
if os.environ.get("XDG_CONFIG_HOME"):
    USER_CONFIG_FILES.insert(0, os.path.join(os.environ["XDG_CONFIG_HOME"], "mypy/config"))
USER_CONFIG_FILES = [os.path.expanduser(f) for f in USER_CONFIG_FILES]

# This must include all reporters defined in mypy.report. This is defined here
# to make reporter names available without importing mypy.report -- this speeds
# up startup.
REPORTER_NAMES: Final = [
    "linecount",
    "any-exprs",
    "linecoverage",
    "memory-xml",
    "cobertura-xml",
    "xml",
    "xslt-html",
    "xslt-txt",
    "html",
    "txt",
    "lineprecision",
]

# Threshold after which we sometimes filter out most errors to avoid very
# verbose output. The default is to show all errors.
MANY_ERRORS_THRESHOLD: Final = -1
