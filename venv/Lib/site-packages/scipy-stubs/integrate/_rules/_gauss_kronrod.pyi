from types import ModuleType
from typing import Generic
from typing_extensions import TypeVar

import numpy as np

from ._base import NestedFixedRule
from ._gauss_legendre import GaussLegendreQuadrature

_XPT_co = TypeVar("_XPT_co", default=ModuleType, covariant=True)

###

class GaussKronrodQuadrature(NestedFixedRule[_XPT_co, np.float64], Generic[_XPT_co]):  # undocumented
    npoints: int
    gauss: GaussLegendreQuadrature[_XPT_co]
    def __init__(self, /, npoints: int, xp: _XPT_co | None = None) -> None: ...
