from types import ModuleType
from typing import Generic
from typing_extensions import TypeVar

import numpy as np

from ._base import FixedRule as FixedRule
from scipy.special import roots_legendre as roots_legendre

_XPT_co = TypeVar("_XPT_co", default=ModuleType, covariant=True)

###

class GaussLegendreQuadrature(FixedRule[_XPT_co, np.float64], Generic[_XPT_co]):  # undocumented
    npoints: int
    def __init__(self, /, npoints: int, xp: _XPT_co | None = None) -> None: ...
