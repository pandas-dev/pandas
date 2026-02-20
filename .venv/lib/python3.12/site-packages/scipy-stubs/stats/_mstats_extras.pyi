from typing import TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp

__all__ = [
    "compare_medians_ms",
    "hdmedian",
    "hdquantiles",
    "hdquantiles_sd",
    "idealfourths",
    "median_cihs",
    "mjci",
    "mquantiles_cimj",
    "rsh",
    "trimmed_mean_ci",
]

_T = TypeVar("_T")
_Tuple2: TypeAlias = tuple[_T, _T]
_FloatND: TypeAlias = onp.ArrayND[np.float64]

_ToProb: TypeAlias = onp.ToFloat | onp.ToFloatND
_ToAxis: TypeAlias = op.CanIndex | None

###

@overload
def hdquantiles(
    data: onp.ToFloat1D, prob: _ToProb = (0.25, 0.5, 0.75), axis: _ToAxis = None, var: onp.ToFalse = False
) -> onp.MArray1D[np.float64]: ...
@overload
def hdquantiles(data: onp.ToFloat1D, prob: _ToProb, axis: _ToAxis, var: onp.ToTrue) -> onp.MArray2D[np.float64]: ...
@overload
def hdquantiles(
    data: onp.ToFloat1D, prob: _ToProb = (0.25, 0.5, 0.75), axis: _ToAxis = None, *, var: onp.ToTrue
) -> onp.MArray2D[np.float64]: ...
@overload
def hdquantiles(
    data: onp.ToFloatND, prob: _ToProb = (0.25, 0.5, 0.75), axis: _ToAxis = None, var: bool = False
) -> onp.MArray[np.float64]: ...

#
@overload
def hdmedian(data: onp.ToFloatND, axis: _ToAxis = -1, var: onp.ToFalse = False) -> onp.MArray[np.float64]: ...
@overload
def hdmedian(data: onp.ToFloatND, axis: _ToAxis, var: onp.ToTrue) -> onp.MArray[np.float64]: ...
@overload
def hdmedian(data: onp.ToFloatND, axis: _ToAxis = -1, *, var: onp.ToTrue) -> onp.MArray[np.float64]: ...

#
def hdquantiles_sd(data: onp.ToFloatND, prob: _ToProb = (0.25, 0.5, 0.75), axis: _ToAxis = None) -> onp.MArray[np.float64]: ...

#
def trimmed_mean_ci(
    data: onp.ToFloatND,
    limits: _Tuple2[onp.ToFloat] | None = (0.2, 0.2),
    inclusive: _Tuple2[op.CanBool] = (True, True),
    alpha: onp.ToJustFloat = 0.05,
    axis: _ToAxis = None,
) -> _FloatND: ...

#
def mjci(data: onp.ToFloatND, prob: _ToProb = (0.25, 0.5, 0.75), axis: _ToAxis = None) -> _FloatND: ...

#
def mquantiles_cimj(
    data: onp.ToFloatND, prob: _ToProb = (0.25, 0.5, 0.75), alpha: onp.ToJustFloat = 0.05, axis: _ToAxis = None
) -> _Tuple2[_FloatND]: ...

#
@overload
def median_cihs(data: onp.ToFloatND, alpha: onp.ToJustFloat = 0.05, axis: None = None) -> _Tuple2[np.float64]: ...
@overload
def median_cihs(data: onp.ToFloatND, alpha: onp.ToJustFloat, axis: op.CanIndex) -> _Tuple2[np.float64 | _FloatND]: ...
@overload
def median_cihs(data: onp.ToFloatND, alpha: onp.ToJustFloat = 0.05, *, axis: op.CanIndex) -> _Tuple2[np.float64 | _FloatND]: ...

#
@overload
def compare_medians_ms(group_1: onp.ToFloatND, group_2: onp.ToFloatND, axis: None = None) -> np.float64: ...
@overload
def compare_medians_ms(group_1: onp.ToFloatND, group_2: onp.ToFloatND, axis: op.CanIndex) -> _FloatND: ...

#
@overload
def idealfourths(data: onp.ToFloatND, axis: None = None) -> list[np.float64]: ...
@overload
def idealfourths(data: onp.ToFloatND, axis: op.CanIndex) -> onp.MArray[np.float64]: ...

#
def rsh(data: onp.ToFloatND, points: onp.ToFloatND | None = None) -> np.float64: ...
