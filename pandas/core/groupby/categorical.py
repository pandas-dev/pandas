from __future__ import annotations

from pandas._core.groupby import categorical
from pandas.core.common import _depr_core

_depr_core()

_globals = globals()

for item in categorical.__dir__():
    _globals[item] = getattr(categorical, item)
