from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, Concatenate, Final, Generic, Literal, NamedTuple, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp

from ._common import ConfidenceInterval
from ._stats_py import SignificanceResult
from ._typing import Alternative, NanPolicy

__all__ = [
    "barnard_exact",
    "boschloo_exact",
    "cramervonmises",
    "cramervonmises_2samp",
    "epps_singleton_2samp",
    "poisson_means_test",
    "somersd",
    "tukey_hsd",
]

_Float2D: TypeAlias = onp.Array2D[np.float64]
_FloatND: TypeAlias = onp.ArrayND[np.float64]
_FloatOrND: TypeAlias = float | _FloatND
_FloatOrNDT = TypeVar("_FloatOrNDT", bound=_FloatOrND, default=Any)

_ToCDF: TypeAlias = str | Callable[Concatenate[float, ...], float | np.float32]
_ToCDFArgs: TypeAlias = tuple[onp.ToFloat, ...]
_CV2Method: TypeAlias = Literal["auto", "asymptotic", "exact"]

###

class Epps_Singleton_2sampResult(NamedTuple, Generic[_FloatOrNDT]):
    statistic: _FloatOrNDT  # readonly
    pvalue: _FloatOrNDT  # readonly

@overload
def epps_singleton_2samp(
    x: onp.ToFloatStrict1D,
    y: onp.ToFloatStrict1D,
    t: onp.ToFloatStrict1D = (0.4, 0.8),
    *,
    axis: Literal[0, -1] | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: onp.ToFalse = False,
) -> Epps_Singleton_2sampResult[float]: ...
@overload
def epps_singleton_2samp(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    t: onp.ToFloatND = (0.4, 0.8),
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: onp.ToFalse = False,
) -> Epps_Singleton_2sampResult[float]: ...
@overload
def epps_singleton_2samp(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    t: onp.ToFloatND = (0.4, 0.8),
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: onp.ToTrue,
) -> Epps_Singleton_2sampResult[_FloatND]: ...
@overload
def epps_singleton_2samp(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    t: onp.ToFloatND = (0.4, 0.8),
    *,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> Epps_Singleton_2sampResult: ...

class CramerVonMisesResult(Generic[_FloatOrNDT]):
    statistic: _FloatOrNDT  # readonly
    pvalue: _FloatOrNDT  # readonly
    def __init__(self, /, statistic: _FloatOrNDT, pvalue: _FloatOrNDT) -> None: ...

@overload
def cramervonmises(
    rvs: onp.ToFloatStrict1D,
    cdf: _ToCDF,
    args: _ToCDFArgs = (),
    *,
    axis: Literal[0, -1] | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: onp.ToFalse = False,
) -> CramerVonMisesResult[float]: ...
@overload
def cramervonmises(
    rvs: onp.ToFloatND,
    cdf: _ToCDF,
    args: _ToCDFArgs = (),
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: onp.ToFalse = False,
) -> CramerVonMisesResult[float]: ...
@overload
def cramervonmises(
    rvs: onp.ToFloatND,
    cdf: _ToCDF,
    args: _ToCDFArgs = (),
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: onp.ToTrue,
) -> CramerVonMisesResult[_FloatND]: ...
@overload
def cramervonmises(
    rvs: onp.ToFloatND,
    cdf: _ToCDF,
    args: _ToCDFArgs = (),
    *,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> CramerVonMisesResult: ...

#
@overload
def cramervonmises_2samp(
    x: onp.ToFloatStrict1D,
    y: onp.ToFloatStrict1D,
    method: _CV2Method = "auto",
    *,
    axis: Literal[0, -1] | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: onp.ToFalse = False,
) -> CramerVonMisesResult[float]: ...
@overload
def cramervonmises_2samp(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    method: _CV2Method = "auto",
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: onp.ToFalse = False,
) -> CramerVonMisesResult[float]: ...
@overload
def cramervonmises_2samp(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    method: _CV2Method = "auto",
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: onp.ToTrue,
) -> CramerVonMisesResult[_FloatND]: ...
@overload
def cramervonmises_2samp(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    method: _CV2Method = "auto",
    *,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> CramerVonMisesResult: ...

#
def poisson_means_test(
    k1: int, n1: float, k2: int, n2: float, *, diff: float = 0, alternative: Alternative = "two-sided"
) -> SignificanceResult[np.float64]: ...

#
@dataclass
class SomersDResult:
    statistic: Final[float]
    pvalue: Final[float]
    table: Final[_Float2D]

def somersd(
    x: onp.ToFloat1D | onp.ToFloat2D, y: onp.ToFloat1D | None = None, alternative: Alternative = "two-sided"
) -> SomersDResult: ...

#
@dataclass
class BarnardExactResult:
    statistic: Final[float]
    pvalue: Final[float]

def barnard_exact(
    table: onp.ToInt2D, alternative: Alternative = "two-sided", pooled: bool = True, n: op.JustInt = 32
) -> BarnardExactResult: ...

#
@dataclass
class BoschlooExactResult:
    statistic: Final[float]
    pvalue: Final[float]

def boschloo_exact(table: onp.ToInt2D, alternative: Alternative = "two-sided", n: op.JustInt = 32) -> BoschlooExactResult: ...

#
class TukeyHSDResult:
    statistic: Final[_Float2D]
    pvalue: Final[_Float2D]
    _ntreatments: Final[int]
    _df: Final[int]
    _stand_err: Final[float]
    def __init__(self, /, statistic: _Float2D, pvalue: _Float2D, _ntreatments: int, _df: int, _stand_err: float) -> None: ...

    #
    _ci: ConfidenceInterval | None
    _ci_cl: float | None
    def confidence_interval(self, /, confidence_level: op.JustFloat | np.float64 = 0.95) -> ConfidenceInterval: ...

def tukey_hsd(arg0: onp.ToFloatND, arg1: onp.ToFloatND, /, *args: onp.ToFloatND, equal_var: bool = True) -> TukeyHSDResult: ...
