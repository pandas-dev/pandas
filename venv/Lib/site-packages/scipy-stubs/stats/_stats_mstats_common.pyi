from typing import Any, Generic, Self, override
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from . import distributions as distributions
from ._typing import BunchMixin

__all__ = ["_find_repeats"]

_ResultT_co = TypeVar("_ResultT_co", bound=np.float64 | onp.ArrayND[np.float64], default=np.float64 | Any, covariant=True)

###

class SiegelslopesResult(BunchMixin[tuple[_ResultT_co, _ResultT_co]], tuple[_ResultT_co, _ResultT_co], Generic[_ResultT_co]):
    def __new__(_cls, slope: _ResultT_co, intercept: _ResultT_co) -> Self: ...
    def __init__(self, /, slope: _ResultT_co, intercept: _ResultT_co) -> None: ...
    @property
    def slope(self, /) -> _ResultT_co: ...
    @property
    def intercept(self, /) -> _ResultT_co: ...

class TheilslopesResult(
    BunchMixin[tuple[_ResultT_co, _ResultT_co, _ResultT_co, _ResultT_co]],
    tuple[_ResultT_co, _ResultT_co, _ResultT_co, _ResultT_co],
    Generic[_ResultT_co],
):
    @override
    def __new__(_cls, slope: _ResultT_co, intercept: _ResultT_co, low_slope: _ResultT_co, high_slope: _ResultT_co) -> Self: ...  # pyrefly:ignore[bad-override]
    @override
    def __init__(  # pyrefly:ignore[bad-override]
        self, /, slope: _ResultT_co, intercept: _ResultT_co, low_slope: _ResultT_co, high_slope: _ResultT_co
    ) -> None: ...

    #
    @property
    def slope(self, /) -> _ResultT_co: ...
    @property
    def intercept(self, /) -> _ResultT_co: ...
    @property
    def low_slope(self, /) -> _ResultT_co: ...
    @property
    def high_slope(self, /) -> _ResultT_co: ...

###

def _find_repeats(arr: onp.ArrayND[npc.number]) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.intp]]: ...  # undocumented
