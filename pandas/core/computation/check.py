from __future__ import annotations

import warnings

from pandas.compat._optional import import_optional_dependency
from pandas.util._exceptions import find_stack_level

from pandas.util.version import Version

# https://github.com/pandas-dev/pandas/issues/65408
# Avoid bad numexpr binaries producing incorrect results
# https://github.com/pydata/numexpr/issues/557
# https://github.com/pydata/numexpr/issues/566
BLOCKED_NUMEXPR_VERSIONS = frozenset({"2.14.1"})

ne = import_optional_dependency("numexpr", errors="warn")

# the installed version if pandas refuses to use it, otherwise None
NUMEXPR_BLOCKED_VERSION: str | None = None

if ne is not None and Version(ne.__version__).base_version in BLOCKED_NUMEXPR_VERSIONS:
    NUMEXPR_BLOCKED_VERSION = ne.__version__
    ne = None

NUMEXPR_INSTALLED = ne is not None

_warned_blocked = False


def warn_numexpr_blocked() -> None:
    """
    Warn at most once that the installed numexpr is unusable.

    Called where numexpr would otherwise have been used, so that users who never
    reach such a code path are not warned (GH#66956).
    """
    global _warned_blocked

    if NUMEXPR_BLOCKED_VERSION is None or _warned_blocked:
        return

    _warned_blocked = True
    warnings.warn(
        f"numexpr version '{NUMEXPR_BLOCKED_VERSION}' can silently return incorrect "
        "results and will not be used by pandas. Install numexpr 2.14.2 or "
        "newer to re-enable numexpr acceleration.",
        UserWarning,
        stacklevel=find_stack_level(),
    )


__all__ = [
    "BLOCKED_NUMEXPR_VERSIONS",
    "NUMEXPR_BLOCKED_VERSION",
    "NUMEXPR_INSTALLED",
    "warn_numexpr_blocked",
]
