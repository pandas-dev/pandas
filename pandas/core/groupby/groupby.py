from __future__ import annotations

from pandas._core.groupby import groupby
from pandas.core.common import _depr_core

_depr_core()

_globals = globals()

for item in groupby.__dir__():
    _globals[item] = getattr(groupby, item)
