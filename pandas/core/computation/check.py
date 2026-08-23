from __future__ import annotations

import warnings

from pandas.compat._optional import import_optional_dependency
from pandas.util._exceptions import find_stack_level

from pandas.util.version import Version

# GH#65408 these numexpr releases ship binaries whose multithreaded evaluation
#  path can intermittently and silently return wrong results (e.g. an all-zero
#  array) for large inputs, see https://github.com/pydata/numexpr/issues/566.
#  Only some builds of an affected release are broken, but we have no way to
#  tell them apart, so we refuse the release as a whole.
BLOCKED_NUMEXPR_VERSIONS = frozenset({"2.14.1"})

ne = import_optional_dependency("numexpr", errors="warn")

if ne is not None and Version(ne.__version__).base_version in BLOCKED_NUMEXPR_VERSIONS:
    warnings.warn(
        f"numexpr version '{ne.__version__}' can silently return incorrect "
        "results and will not be used by pandas. Install numexpr 2.14.2 or "
        "newer to re-enable numexpr acceleration.",
        UserWarning,
        stacklevel=find_stack_level(),
    )
    ne = None

NUMEXPR_INSTALLED = ne is not None

__all__ = ["BLOCKED_NUMEXPR_VERSIONS", "NUMEXPR_INSTALLED"]
