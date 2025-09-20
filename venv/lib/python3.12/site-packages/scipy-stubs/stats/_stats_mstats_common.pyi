from typing import Any, Generic, Literal, Self, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from . import distributions as distributions
from ._typing import BunchMixin, NanPolicy

__all__ = ["_find_repeats", "siegelslopes", "theilslopes"]

_ResultT_co = TypeVar("_ResultT_co", bound=np.float64 | onp.ArrayND[np.float64], default=np.float64 | Any, covariant=True)

_Method: TypeAlias = Literal["hierarchical", "separate"]

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
    def __new__(_cls, slope: _ResultT_co, intercept: _ResultT_co, low_slope: _ResultT_co, high_slope: _ResultT_co) -> Self: ...
    def __init__(
        self, /, slope: _ResultT_co, intercept: _ResultT_co, low_slope: _ResultT_co, high_slope: _ResultT_co
    ) -> None: ...
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

#
@overload
def siegelslopes(
    y: onp.ToFloatND,
    x: onp.ToFloatND | None = None,
    method: _Method = "hierarchical",
    *,
    axis: None = None,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> SiegelslopesResult[np.float64]: ...
@overload
def siegelslopes(
    y: onp.ToFloatStrict1D,
    x: onp.ToFloatStrict1D | None = None,
    method: _Method = "hierarchical",
    *,
    axis: int | None = None,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> SiegelslopesResult[np.float64]: ...
@overload
def siegelslopes(
    y: onp.ToFloatStrict2D,
    x: onp.ToFloatStrict2D | None = None,
    method: _Method = "hierarchical",
    *,
    axis: int,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> SiegelslopesResult[onp.Array1D[np.float64]]: ...
@overload
def siegelslopes(
    y: onp.ToFloatStrict3D,
    x: onp.ToFloatStrict3D | None = None,
    method: _Method = "hierarchical",
    *,
    axis: int,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> SiegelslopesResult[onp.Array2D[np.float64]]: ...
@overload
def siegelslopes(
    y: onp.ToFloatND,
    x: onp.ToFloatND | None = None,
    method: _Method = "hierarchical",
    *,
    axis: int | None = None,
    keepdims: Literal[True],
    nan_policy: NanPolicy = "propagate",
) -> SiegelslopesResult[onp.ArrayND[np.float64]]: ...
@overload
def siegelslopes(
    y: onp.ToFloatND,
    x: onp.ToFloatND | None = None,
    method: _Method = "hierarchical",
    *,
    axis: int | None = None,
    keepdims: bool = False,
    nan_policy: NanPolicy = "propagate",
) -> SiegelslopesResult[np.float64 | Any]: ...

#
@overload
def theilslopes(
    y: onp.ToFloatND,
    x: onp.ToFloatND | None = None,
    alpha: onp.ToJustFloat = 0.95,
    method: _Method = "separate",
    *,
    axis: None = None,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> TheilslopesResult[np.float64]: ...
@overload
def theilslopes(
    y: onp.ToFloatStrict1D,
    x: onp.ToFloatStrict1D | None = None,
    alpha: onp.ToJustFloat = 0.95,
    method: _Method = "separate",
    *,
    axis: int | None = None,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> TheilslopesResult[np.float64]: ...
@overload
def theilslopes(
    y: onp.ToFloatStrict2D,
    x: onp.ToFloatStrict2D | None = None,
    alpha: onp.ToJustFloat = 0.95,
    method: _Method = "separate",
    *,
    axis: int,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> TheilslopesResult[onp.Array1D[np.float64]]: ...
@overload
def theilslopes(
    y: onp.ToFloatStrict3D,
    x: onp.ToFloatStrict3D | None = None,
    alpha: onp.ToJustFloat = 0.95,
    method: _Method = "separate",
    *,
    axis: int,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> TheilslopesResult[onp.Array2D[np.float64]]: ...
@overload
def theilslopes(
    y: onp.ToFloatND,
    x: onp.ToFloatND | None = None,
    alpha: onp.ToJustFloat = 0.95,
    method: _Method = "separate",
    *,
    axis: int | None = None,
    keepdims: Literal[True],
    nan_policy: NanPolicy = "propagate",
) -> TheilslopesResult[onp.ArrayND[np.float64]]: ...
@overload
def theilslopes(
    y: onp.ToFloatND,
    x: onp.ToFloatND | None = None,
    alpha: onp.ToJustFloat = 0.95,
    method: _Method = "separate",
    *,
    axis: int | None = None,
    keepdims: bool = False,
    nan_policy: NanPolicy = "propagate",
) -> TheilslopesResult[np.float64 | Any]: ...
