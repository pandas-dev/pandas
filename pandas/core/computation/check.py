from __future__ import annotations

import warnings

from pandas.compat._optional import (
    VERSIONS,
    get_version,
)
from pandas.util._exceptions import find_stack_level

from pandas.util.version import Version

# the minimum also screens out builds that can silently return incorrect
#  results, GH#62468, GH#65408
NUMEXPR_MIN_VERSION = VERSIONS["numexpr"]

# not import_optional_dependency, which either warns at import (GH#66956) or
#  hides that a too-old numexpr is installed
try:
    import numexpr as ne
except ImportError:
    ne = None

# the installed version if pandas refuses to use it, otherwise None
NUMEXPR_BLOCKED_VERSION: str | None = None

if ne is not None and Version(get_version(ne)) < Version(NUMEXPR_MIN_VERSION):
    NUMEXPR_BLOCKED_VERSION = get_version(ne)
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

    warnings.warn(
        f"numexpr version '{NUMEXPR_BLOCKED_VERSION}' is older than the minimum "
        f"supported version '{NUMEXPR_MIN_VERSION}' and will not be used by "
        "pandas. Upgrade numexpr to re-enable numexpr acceleration.",
        UserWarning,
        stacklevel=find_stack_level(),
    )
    # set after warning, so that a caller suppressing warnings does not use up
    #  the one-shot
    _warned_blocked = True


__all__ = [
    "NUMEXPR_BLOCKED_VERSION",
    "NUMEXPR_INSTALLED",
    "NUMEXPR_MIN_VERSION",
    "warn_numexpr_blocked",
]
