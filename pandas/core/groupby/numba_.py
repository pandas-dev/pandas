from __future__ import annotations

from pandas._core.groupby import numba_
from pandas.core.common import _depr_core

_depr_core()

_globals = globals()

for item in numba_.__dir__():
    _globals[item] = getattr(numba_, item)
