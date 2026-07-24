from __future__ import annotations

import os
import sys
from typing import Final

# Earliest fully supported Python 3.x version. Used as the default Python
# version in tests. Mypy wheels should be built starting with this version,
# and CI tests should be run on this version (and later versions).
PYTHON3_VERSION: Final = (3, 10)

# Earliest Python 3.x version supported via --python-version 3.x. To run
# mypy, at least version PYTHON3_VERSION is needed.
PYTHON3_VERSION_MIN: Final = (3, 10)  # Keep in sync with supported target versions

CACHE_DIR: Final = ".mypy_cache"
SQLITE_NUM_SHARDS: Final = 16

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

RECURSION_LIMIT: Final = 2**14

# It looks like Windows is slow with processes, causing test flakiness even
# with our generous timeouts, so we set them higher.
WORKER_START_INTERVAL: Final = 0.01 if sys.platform != "win32" else 0.03
WORKER_START_TIMEOUT: Final = 3 if sys.platform != "win32" else 10
WORKER_SHUTDOWN_TIMEOUT: Final = 1 if sys.platform != "win32" else 3

WORKER_CONNECTION_TIMEOUT: Final = 10
WORKER_IDLE_TIMEOUT: Final = 600
WORKER_DONE_TIMEOUT: Final = 600
