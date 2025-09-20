from collections.abc import Callable, Sequence
from dataclasses import dataclass
from types import ModuleType
from typing import Any, Generic, Literal as L, Never, Protocol, Self, TypeAlias, overload, type_check_only
from typing_extensions import NamedTuple, TypeVar, deprecated

import numpy as np
import numpy.typing as npt
import numpy_typing_compat as nptc
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._resampling import BootstrapMethod, ResamplingMethod
from ._stats_mstats_common import siegelslopes, theilslopes
from ._typing import Alternative, BaseBunch, BunchMixin, NanPolicy, PowerDivergenceStatistic

__all__ = [
    "alexandergovern",
    "brunnermunzel",
    "chisquare",
    "combine_pvalues",
    "cumfreq",
    "describe",
    "energy_distance",
    "expectile",
    "f_oneway",
    "find_repeats",
    "fisher_exact",
    "friedmanchisquare",
    "gmean",
    "gstd",
    "gzscore",
    "hmean",
    "iqr",
    "jarque_bera",
    "kendalltau",
    "kruskal",
    "ks_1samp",
    "ks_2samp",
    "kstest",
    "kurtosis",
    "kurtosistest",
    "linregress",
    "lmoment",
    "median_abs_deviation",
    "mode",
    "moment",
    "normaltest",
    "obrientransform",
    "pearsonr",
    "percentileofscore",
    "pmean",
    "pointbiserialr",
    "power_divergence",
    "quantile_test",
    "rankdata",
    "ranksums",
    "relfreq",
    "scoreatpercentile",
    "sem",
    "siegelslopes",
    "sigmaclip",
    "skew",
    "skewtest",
    "spearmanr",
    "theilslopes",
    "tiecorrect",
    "tmax",
    "tmean",
    "tmin",
    "trim1",
    "trim_mean",
    "trimboth",
    "tsem",
    "tstd",
    "ttest_1samp",
    "ttest_ind",
    "ttest_ind_from_stats",
    "ttest_rel",
    "tvar",
    "wasserstein_distance",
    "wasserstein_distance_nd",
    "weightedtau",
    "zmap",
    "zscore",
]

###

_SCT = TypeVar("_SCT", bound=np.generic)

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])
_InexactT = TypeVar("_InexactT", bound=npc.inexact)
_FloatT = TypeVar("_FloatT", bound=npc.floating, default=npc.floating)
_RealT = TypeVar("_RealT", bound=_Real0D, default=_Real0D)
_RealT_co = TypeVar("_RealT_co", bound=_Real0D, default=_Real0D, covariant=True)

_IntOrArrayT_co = TypeVar("_IntOrArrayT_co", bound=_ScalarOrND[np.intp], default=_ScalarOrND[np.intp], covariant=True)
_FloatOrArrayT = TypeVar("_FloatOrArrayT", bound=_ScalarOrND[npc.floating])
_FloatOrArrayT_co = TypeVar(
    "_FloatOrArrayT_co", bound=float | _ScalarOrND[npc.floating], default=float | onp.ArrayND[np.float64], covariant=True
)
_FloatOrArrayT2_co = TypeVar(
    "_FloatOrArrayT2_co", bound=float | _ScalarOrND[npc.floating], default=float | onp.ArrayND[np.float64], covariant=True
)
_F64OrArrayT_co = TypeVar(
    "_F64OrArrayT_co", bound=np.float64 | onp.ArrayND[np.float64], default=np.float64 | onp.ArrayND[np.float64], covariant=True
)
_RealOrArrayT_co = TypeVar("_RealOrArrayT_co", bound=_ScalarOrND[_Real0D], default=_ScalarOrND[Any], covariant=True)

_Real0D: TypeAlias = npc.integer | npc.floating

_ScalarOrND: TypeAlias = _SCT | onp.ArrayND[_SCT]
_FloatOrND: TypeAlias = _ScalarOrND[_FloatT]
_RealOrND: TypeAlias = _ScalarOrND[_RealT]

_InterpolationMethod: TypeAlias = L["linear", "lower", "higher", "nearest", "midpoint"]
_TrimTail: TypeAlias = L["left", "right"]
_KendallTauMethod: TypeAlias = L["auto", "asymptotic", "exact"]
_KendallTauVariant: TypeAlias = L["b", "c"]
_KS1TestMethod: TypeAlias = L[_KS2TestMethod, "approx"]
_KS2TestMethod: TypeAlias = L["auto", "exact", "asymp"]
_CombinePValuesMethod: TypeAlias = L["fisher", "pearson", "tippett", "stouffer", "mudholkar_george"]
_RankMethod: TypeAlias = L["average", "min", "max", "dense", "ordinal"]

_LMomentOrder: TypeAlias = L[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16] | npc.integer
_LMomentOrder1D: TypeAlias = Sequence[_LMomentOrder] | onp.CanArrayND[npc.integer]
_RealLimits: TypeAlias = tuple[float | _Real0D, float | _Real0D]
_Weigher: TypeAlias = Callable[[int], float | _Real0D]

@type_check_only
class _RVSCallable(Protocol):
    def __call__(self, /, *, size: int | tuple[int, ...]) -> onp.ArrayND[npc.floating]: ...

@type_check_only
class _MADCenterFunc(Protocol):
    def __call__(self, x: onp.Array1D[np.float64], /, *, axis: int | None) -> onp.ToFloat: ...

@type_check_only
class _TestResultTuple(NamedTuple, Generic[_FloatOrArrayT_co]):
    statistic: _FloatOrArrayT_co
    pvalue: _FloatOrArrayT_co

@type_check_only
class _TestResultBunch(BaseBunch[_FloatOrArrayT_co, _FloatOrArrayT2_co], Generic[_FloatOrArrayT_co, _FloatOrArrayT2_co]):
    @property
    def statistic(self, /) -> _FloatOrArrayT_co: ...
    @property
    def pvalue(self, /) -> _FloatOrArrayT2_co: ...
    def __new__(_cls, statistic: _FloatOrArrayT_co, pvalue: _FloatOrArrayT2_co) -> Self: ...
    def __init__(self, /, statistic: _FloatOrArrayT_co, pvalue: _FloatOrArrayT2_co) -> None: ...

###

class SkewtestResult(_TestResultTuple[_FloatOrArrayT_co], Generic[_FloatOrArrayT_co]): ...
class KurtosistestResult(_TestResultTuple[_FloatOrArrayT_co], Generic[_FloatOrArrayT_co]): ...
class NormaltestResult(_TestResultTuple[_FloatOrArrayT_co], Generic[_FloatOrArrayT_co]): ...
class Ttest_indResult(_TestResultTuple[_FloatOrArrayT_co], Generic[_FloatOrArrayT_co]): ...
class Power_divergenceResult(_TestResultTuple[_FloatOrArrayT_co], Generic[_FloatOrArrayT_co]): ...
class RanksumsResult(_TestResultTuple[_FloatOrArrayT_co], Generic[_FloatOrArrayT_co]): ...
class KruskalResult(_TestResultTuple[_FloatOrArrayT_co], Generic[_FloatOrArrayT_co]): ...
class FriedmanchisquareResult(_TestResultTuple[_FloatOrArrayT_co], Generic[_FloatOrArrayT_co]): ...
class BrunnerMunzelResult(_TestResultTuple[_FloatOrArrayT_co], Generic[_FloatOrArrayT_co]): ...
class F_onewayResult(_TestResultTuple[_FloatOrArrayT_co], Generic[_FloatOrArrayT_co]): ...

class ConfidenceInterval(NamedTuple, Generic[_FloatOrArrayT_co]):
    low: _FloatOrArrayT_co
    high: _FloatOrArrayT_co

class DescribeResult(NamedTuple, Generic[_RealOrArrayT_co, _FloatOrArrayT_co]):
    nobs: int
    minmax: tuple[_RealOrArrayT_co, _RealOrArrayT_co]
    mean: _FloatOrArrayT_co
    variance: _FloatOrArrayT_co
    skewness: _FloatOrArrayT_co
    kurtosis: _FloatOrArrayT_co

class ModeResult(NamedTuple, Generic[_RealOrArrayT_co, _IntOrArrayT_co]):
    mode: _RealOrArrayT_co
    count: _IntOrArrayT_co  # type: ignore[assignment]  # pyright: ignore[reportIncompatibleMethodOverride]

class HistogramResult(NamedTuple):
    count: onp.Array1D[np.float64]  # type: ignore[assignment]  # pyright: ignore[reportIncompatibleMethodOverride]
    lowerlimit: L[0] | npc.floating
    binsize: onp.Array1D[np.float64]
    extrapoints: int

class CumfreqResult(NamedTuple):
    cumcount: onp.Array1D[np.float64]
    lowerlimit: L[0] | npc.floating
    binsize: onp.Array1D[np.float64]
    extrapoints: int

class RelfreqResult(NamedTuple):
    frequency: onp.Array1D[np.float64]
    lowerlimit: L[0] | npc.floating
    binsize: onp.Array1D[np.float64]
    extrapoints: int

class SigmaclipResult(NamedTuple, Generic[_RealT_co, _FloatOrArrayT_co]):
    clipped: onp.Array1D[_RealT_co]
    lower: _FloatOrArrayT_co
    upper: _FloatOrArrayT_co

class RepeatedResults(NamedTuple):
    values: onp.Array1D[np.float64]
    counts: onp.Array1D[np.intp]

@dataclass
class AlexanderGovernResult:
    statistic: float
    pvalue: float

@dataclass
class QuantileTestResult:
    statistic: float
    statistic_type: int
    pvalue: float
    _alternative: list[str]
    _x: onp.ArrayND[_Real0D]
    _p: float
    def confidence_interval(self, /, confidence_level: float = 0.95) -> float: ...

class SignificanceResult(_TestResultBunch[_FloatOrArrayT_co, _FloatOrArrayT_co], Generic[_FloatOrArrayT_co]): ...
class PearsonRResultBase(_TestResultBunch[_FloatOrArrayT_co, _F64OrArrayT_co], Generic[_FloatOrArrayT_co, _F64OrArrayT_co]): ...

class PearsonRResult(PearsonRResultBase[_FloatOrArrayT_co, _F64OrArrayT_co], Generic[_FloatOrArrayT_co, _F64OrArrayT_co]):
    _alternative: Alternative
    _n: int
    _x: onp.ArrayND[_Real0D]
    _y: onp.ArrayND[_Real0D]
    _axis: int
    correlation: _FloatOrArrayT_co  # alias for `statistic`

    def __init__(  # pyright: ignore[reportInconsistentConstructor]
        self,
        /,
        statistic: _FloatOrArrayT_co,
        pvalue: _F64OrArrayT_co,
        alternative: Alternative,
        n: int,
        x: onp.ArrayND[_Real0D],
        y: onp.ArrayND[_Real0D],
        axis: int,
    ) -> None: ...
    def confidence_interval(
        self, /, confidence_level: float = 0.95, method: BootstrapMethod | None = None
    ) -> ConfidenceInterval[_FloatOrArrayT_co]: ...

class TtestResultBase(_TestResultBunch[_FloatOrArrayT_co, _FloatOrArrayT_co], Generic[_FloatOrArrayT_co]):
    @property
    def df(self, /) -> _FloatOrArrayT_co: ...
    def __new__(_cls, statistic: _FloatOrArrayT_co, pvalue: _FloatOrArrayT_co, *, df: _FloatOrArrayT_co) -> Self: ...
    def __init__(self, /, statistic: _FloatOrArrayT_co, pvalue: _FloatOrArrayT_co, *, df: _FloatOrArrayT_co) -> None: ...

class TtestResult(TtestResultBase[_FloatOrArrayT_co], Generic[_FloatOrArrayT_co]):
    _alternative: Alternative
    _standard_error: _FloatOrArrayT_co
    _estimate: _FloatOrArrayT_co
    _statistic_np: _FloatOrArrayT_co
    _dtype: np.dtype[npc.floating]
    _xp: ModuleType

    def __init__(  # pyright: ignore[reportInconsistentConstructor]
        self,
        /,
        statistic: _FloatOrArrayT_co,
        pvalue: _FloatOrArrayT_co,
        df: _FloatOrArrayT_co,
        alternative: Alternative,
        standard_error: _FloatOrArrayT_co,
        estimate: _FloatOrArrayT_co,
        statistic_np: _FloatOrArrayT_co | None = None,
        xp: ModuleType | None = None,
    ) -> None: ...
    def confidence_interval(self, /, confidence_level: float = 0.95) -> ConfidenceInterval[_FloatOrArrayT_co]: ...

class KstestResult(_TestResultBunch[np.float64, np.float64]):
    @property
    def statistic_location(self, /) -> np.float64: ...
    @property
    def statistic_sign(self, /) -> np.int8: ...
    def __new__(
        _cls, statistic: np.float64, pvalue: np.float64, *, statistic_location: np.float64, statistic_sign: np.int8
    ) -> Self: ...
    def __init__(
        self, /, statistic: np.float64, pvalue: np.float64, *, statistic_location: np.float64, statistic_sign: np.int8
    ) -> None: ...

Ks_2sampResult = KstestResult

class LinregressResult(
    BunchMixin[
        tuple[_FloatOrArrayT_co, _FloatOrArrayT_co, _FloatOrArrayT_co, _FloatOrArrayT_co, _FloatOrArrayT_co, _FloatOrArrayT_co]
    ],
    tuple[_FloatOrArrayT_co, _FloatOrArrayT_co, _FloatOrArrayT_co, _FloatOrArrayT_co, _FloatOrArrayT_co, _FloatOrArrayT_co],
    Generic[_FloatOrArrayT_co],
):
    def __new__(
        _cls,
        slope: _FloatOrArrayT_co,
        intercept: _FloatOrArrayT_co,
        rvalue: _FloatOrArrayT_co,
        pvalue: _FloatOrArrayT_co,
        stderr: _FloatOrArrayT_co,
        *,
        intercept_stderr: _FloatOrArrayT_co,
    ) -> Self: ...
    def __init__(
        self,
        /,
        slope: _FloatOrArrayT_co,
        intercept: _FloatOrArrayT_co,
        rvalue: _FloatOrArrayT_co,
        pvalue: _FloatOrArrayT_co,
        stderr: _FloatOrArrayT_co,
        *,
        intercept_stderr: _FloatOrArrayT_co,
    ) -> None: ...
    @property
    def slope(self, /) -> _FloatOrArrayT_co: ...
    @property
    def intercept(self, /) -> _FloatOrArrayT_co: ...
    @property
    def rvalue(self, /) -> _FloatOrArrayT_co: ...
    @property
    def pvalue(self, /) -> _FloatOrArrayT_co: ...
    @property
    def stderr(self, /) -> _FloatOrArrayT_co: ...
    @property
    def intercept_stderr(self, /) -> _FloatOrArrayT_co: ...

# TODO(jorenham): improve
def gmean(
    a: onp.ToFloatND,
    axis: int | None = 0,
    dtype: npt.DTypeLike | None = None,
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> _RealOrND: ...

# TODO(jorenham): improve
def hmean(
    a: onp.ToFloatND,
    axis: int | None = 0,
    dtype: npt.DTypeLike | None = None,
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> _RealOrND: ...

# TODO(jorenham): improve
def pmean(
    a: onp.ToFloatND,
    p: float | _Real0D,
    *,
    axis: int | None = 0,
    dtype: npt.DTypeLike | None = None,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> _RealOrND: ...

# NOTE: The two mypy `overload-overlap` errors are false positive
@overload  # int {0,1}d, keepdims=False (default)
def mode(
    a: int | Sequence[int], axis: int | None = 0, nan_policy: NanPolicy = "propagate", keepdims: L[False] = False
) -> ModeResult[np.int_, np.intp]: ...
@overload  # int ?d, axis=None, keepdims=False (default)
def mode(
    a: int | onp.SequenceND[int], axis: None, nan_policy: NanPolicy = "propagate", keepdims: L[False] = False
) -> ModeResult[np.int_, np.intp]: ...
@overload  # int ?d, keepdims=True (keyword)
def mode(
    a: int | onp.SequenceND[int], axis: int | None = 0, nan_policy: NanPolicy = "propagate", *, keepdims: L[True]
) -> ModeResult[onp.ArrayND[np.int_], onp.ArrayND[np.intp]]: ...
@overload  # int >1d, axis: int (default)
def mode(  # type: ignore[overload-overlap]
    a: Sequence[onp.SequenceND[int]], axis: int = 0, nan_policy: NanPolicy = "propagate", keepdims: bool = False
) -> ModeResult[onp.ArrayND[np.int_], onp.ArrayND[np.intp]]: ...
@overload  # float {0,1}d, keepdims=False (default)
def mode(
    a: op.JustFloat | Sequence[op.JustFloat],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> ModeResult[np.float64, np.intp]: ...
@overload  # float ?d, axis=None, keepdims=False (default)
def mode(
    a: op.JustFloat | onp.SequenceND[op.JustFloat], axis: None, nan_policy: NanPolicy = "propagate", keepdims: L[False] = False
) -> ModeResult[np.float64, np.intp]: ...
@overload  # float ?d, keepdims=True (keyword)
def mode(
    a: op.JustFloat | onp.SequenceND[op.JustFloat],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> ModeResult[onp.ArrayND[np.float64], onp.ArrayND[np.intp]]: ...
@overload  # float >1d, axis: int (default)
def mode(  # type: ignore[overload-overlap]
    a: Sequence[onp.SequenceND[op.JustFloat]], axis: int = 0, nan_policy: NanPolicy = "propagate", keepdims: bool = False
) -> ModeResult[onp.ArrayND[np.float64], onp.ArrayND[np.intp]]: ...
@overload  # T@real {0,1}d, keepdims=False (default)
def mode(
    a: _RealT | onp.ToArrayStrict1D[Never, _RealT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> ModeResult[_RealT, np.intp]: ...
@overload  # T@real ?d, axis=None, keepdims=False (default)
def mode(
    a: _RealT | onp.ToArrayND[Never, _RealT], axis: None, nan_policy: NanPolicy = "propagate", keepdims: L[False] = False
) -> ModeResult[_RealT, np.intp]: ...
@overload  # T@real ?d, keepdims=True (keyword)
def mode(
    a: _RealT | onp.ToArrayND[Never, _RealT], axis: int | None = 0, nan_policy: NanPolicy = "propagate", *, keepdims: L[True]
) -> ModeResult[onp.ArrayND[_RealT], onp.ArrayND[np.intp]]: ...
@overload  # T@real >1d, axis: int (default)
def mode(
    a: onp.CanArray[onp.AtLeast2D, np.dtype[_RealT]], axis: int = 0, nan_policy: NanPolicy = "propagate", keepdims: bool = False
) -> ModeResult[onp.ArrayND[_RealT], onp.ArrayND[np.intp]]: ...
@overload  # real ?d, axis=None, keepdims=False (default)
def mode(
    a: onp.ToFloat | onp.ToFloatND, axis: None, nan_policy: NanPolicy = "propagate", keepdims: L[False] = False
) -> ModeResult[np.float64 | Any, np.intp]: ...
@overload  # real ?d, keepdims=True (keyword)
def mode(
    a: onp.ToFloat | onp.ToFloatND, axis: int | None = 0, nan_policy: NanPolicy = "propagate", *, keepdims: L[True]
) -> ModeResult[onp.ArrayND[np.float64 | Any], onp.ArrayND[np.intp]]: ...
@overload  # real ?d
def mode(
    a: onp.ToFloat | onp.ToFloatND, axis: int | None = 0, nan_policy: NanPolicy = "propagate", keepdims: bool = False
) -> ModeResult: ...

# TODO(jorenham): improve
def tmean(
    a: onp.ToFloatND,
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> _FloatOrND: ...

# TODO(jorenham): improve
def tvar(
    a: onp.ToFloatND,
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int | None = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> _FloatOrND: ...

# TODO(jorenham): improve
def tmin(
    a: onp.ToFloatND,
    lowerlimit: float | _Real0D | None = None,
    axis: int | None = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: bool = False,
) -> _RealOrND: ...

# TODO(jorenham): improve
def tmax(
    a: onp.ToFloatND,
    upperlimit: float | _Real0D | None = None,
    axis: int | None = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: bool = False,
) -> _RealOrND: ...

# TODO(jorenham): improve
def tstd(
    a: onp.ToFloatND,
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int | None = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> _FloatOrND: ...

# TODO(jorenham): improve
def tsem(
    a: onp.ToFloatND,
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int | None = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> _FloatOrND: ...

#
@overload
def gstd(
    a: onp.ToFloatND, axis: None, ddof: int = 1, *, keepdims: L[False] = False, nan_policy: NanPolicy = "propagate"
) -> np.float64: ...
@overload
def gstd(
    a: onp.ToFloatStrict1D,
    axis: int | None = 0,
    ddof: int = 1,
    *,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> np.float64: ...
@overload
def gstd(
    a: onp.ToFloatND, axis: int | None = 0, ddof: int = 1, *, keepdims: L[True], nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[np.float64]: ...
@overload
def gstd(
    a: onp.ToFloatND, axis: int | None = 0, ddof: int = 1, *, keepdims: bool = False, nan_policy: NanPolicy = "propagate"
) -> np.float64 | onp.ArrayND[np.float64]: ...

# TODO(jorenham): improve
def moment(
    a: onp.ToFloatND,
    order: int = 1,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | npc.floating | None = None,
    keepdims: bool = False,
) -> _FloatOrND: ...

# TODO(jorenham): improve
def skew(
    a: onp.ToFloatND, axis: int | None = 0, bias: bool = True, nan_policy: NanPolicy = "propagate", *, keepdims: bool = False
) -> _FloatOrND: ...

# TODO(jorenham): improve
def kurtosis(
    a: onp.ToFloatND,
    axis: int | None = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: bool = False,
) -> _FloatOrND: ...

#
def describe(
    a: onp.ToFloatND, axis: int | None = 0, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult: ...

# TODO(jorenham): improve
def skewtest(
    a: onp.ToFloatND,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: bool = False,
) -> SkewtestResult: ...

# TODO(jorenham): improve
def kurtosistest(
    a: onp.ToFloatND,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: bool = False,
) -> KurtosistestResult: ...

# TODO(jorenham): improve
def normaltest(
    a: onp.ToFloatND, axis: int | None = 0, nan_policy: NanPolicy = "propagate", *, keepdims: bool = False
) -> NormaltestResult: ...

# TODO(jorenham): improve
def jarque_bera(
    x: onp.ToFloatND, *, axis: int | None = None, nan_policy: NanPolicy = "propagate", keepdims: bool = False
) -> SignificanceResult: ...

# TODO(jorenham): improve
def scoreatpercentile(
    a: onp.ToFloat1D,
    per: onp.ToFloat | onp.ToFloatND,
    limit: _RealLimits | tuple[()] = (),
    interpolation_method: L["fraction", "lower", "higher"] = "fraction",
    axis: int | None = None,
) -> _FloatOrND: ...

#
def percentileofscore(
    a: onp.ToFloat1D,
    score: onp.ToFloat | onp.ToFloatND,
    kind: L["rank", "weak", "strict", "mean"] = "rank",
    nan_policy: NanPolicy = "propagate",
) -> np.float64: ...

#
def cumfreq(
    a: onp.ToFloatND, numbins: int = 10, defaultreallimits: _RealLimits | None = None, weights: onp.ToFloatND | None = None
) -> CumfreqResult: ...
def relfreq(
    a: onp.ToFloatND, numbins: int = 10, defaultreallimits: _RealLimits | None = None, weights: onp.ToFloatND | None = None
) -> RelfreqResult: ...

#
def obrientransform(*samples: onp.ToFloatND) -> onp.Array2D[npc.floating] | onp.Array1D[np.object_]: ...

#
@overload  # 1d ~inexact64 | +integer, keepdims=False (default)
def sem(
    a: onp.ToArrayStrict1D[complex, npc.inexact64 | npc.integer | np.bool_],
    axis: L[0, -1] | None = 0,
    ddof: int = 1,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # >1d ~inexact64 | +integer, axis: int (default)
def sem(
    a: onp.CanArray[onp.AtLeast2D, np.dtype[npc.inexact64 | npc.integer | np.bool_]] | Sequence[onp.SequenceND[complex]],
    axis: int = 0,
    ddof: int = 1,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d ~inexact64 | +integer, axis=None, keepdims=False (default)
def sem(
    a: onp.ToArrayND[complex, npc.inexact64 | npc.integer | np.bool_],
    axis: None,
    ddof: int = 1,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # ?d ~inexact64 | +integer, keepdims=True
def sem(
    a: onp.ToArrayND[complex, npc.inexact64 | npc.integer | np.bool_],
    axis: int | None = 0,
    ddof: int = 1,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # 1d +complex, keepdims=False (default)
def sem(
    a: onp.ToComplexStrict1D,
    axis: L[0, -1] | None = 0,
    ddof: int = 1,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> npc.floating: ...
@overload  # >1d +complex, axis: int (default)
def sem(
    a: onp.CanArray[onp.AtLeast2D, np.dtype[npc.number]],
    axis: int = 0,
    ddof: int = 1,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: bool = False,
) -> onp.ArrayND[npc.floating]: ...
@overload  # ?d +complex, axis=None, keepdims=False (default)
def sem(
    a: onp.ToComplexND, axis: None, ddof: int = 1, nan_policy: NanPolicy = "propagate", *, keepdims: L[False] = False
) -> npc.floating: ...
@overload  # ?d +complex, keepdims=True
def sem(
    a: onp.ToComplexND, axis: int | None = 0, ddof: int = 1, nan_policy: NanPolicy = "propagate", *, keepdims: L[True]
) -> onp.ArrayND[npc.floating]: ...
@overload  # ?d +complex
def sem(
    a: onp.ToComplexND, axis: int | None = 0, ddof: int = 1, nan_policy: NanPolicy = "propagate", *, keepdims: bool = False
) -> _FloatOrND: ...

# NOTE: keep in sync with `gzscore`
@overload  # +integer, known shape
def zscore(
    a: nptc.CanArray[_ShapeT, np.dtype[npc.integer | np.bool_]],
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64, _ShapeT]: ...
@overload  # known inexact dtype, known shape
def zscore(
    a: nptc.CanArray[_ShapeT, np.dtype[_InexactT]], axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[_InexactT, _ShapeT]: ...
@overload  # float 1d
def zscore(
    a: Sequence[float], axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array1D[np.float64]: ...
@overload  # float 2d
def zscore(
    a: Sequence[Sequence[float]], axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array2D[np.float64]: ...
@overload  # float 3d
def zscore(
    a: Sequence[Sequence[Sequence[float]]], axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array3D[np.float64]: ...
@overload  # complex 1d
def zscore(
    a: Sequence[op.JustComplex], axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array1D[np.complex128]: ...
@overload  # complex 2d
def zscore(
    a: Sequence[Sequence[op.JustComplex]], axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array2D[np.complex128]: ...
@overload  # complex 3d
def zscore(
    a: Sequence[Sequence[Sequence[op.JustComplex]]], axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array3D[np.complex128]: ...
@overload  # floating fallback
def zscore(  # the weird shape-type is a workaround for a bug in pyright's overlapping overload detection on numpy<2.1
    a: onp.ToFloatND, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[npc.floating, tuple[int] | tuple[Any, ...]]: ...
@overload  # complex fallback
def zscore(
    a: onp.ToJustComplexND, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[npc.complexfloating]: ...

# NOTE: keep in sync with `zscore`
@overload  # +integer, known shape
def gzscore(
    a: nptc.CanArray[_ShapeT, np.dtype[npc.integer | np.bool_]],
    *,
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64, _ShapeT]: ...
@overload  # known inexact dtype, known shape
def gzscore(
    a: nptc.CanArray[_ShapeT, np.dtype[_InexactT]], *, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[_InexactT, _ShapeT]: ...
@overload  # float 1d
def gzscore(
    a: Sequence[float], *, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array1D[np.float64]: ...
@overload  # float 2d
def gzscore(
    a: Sequence[Sequence[float]], *, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array2D[np.float64]: ...
@overload  # float 3d
def gzscore(
    a: Sequence[Sequence[Sequence[float]]], *, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array3D[np.float64]: ...
@overload  # complex 1d
def gzscore(
    a: Sequence[op.JustComplex], *, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array1D[np.complex128]: ...
@overload  # complex 2d
def gzscore(
    a: Sequence[Sequence[op.JustComplex]], *, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array2D[np.complex128]: ...
@overload  # complex 3d
def gzscore(
    a: Sequence[Sequence[Sequence[op.JustComplex]]], *, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array3D[np.complex128]: ...
@overload  # floating fallback
def gzscore(  # the weird shape-type is a workaround for a bug in pyright's overlapping overload detection on numpy<2.1
    a: onp.ToFloatND, *, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[npc.floating, tuple[int] | tuple[Any, ...]]: ...
@overload  # complex fallback
def gzscore(
    a: onp.ToJustComplexND, *, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[npc.complexfloating]: ...

# TODO(jorenham): improve like zscore
@overload  # (real vector-like, real vector-like) -> floating vector
def zmap(
    scores: onp.ToFloat1D, compare: onp.ToFloat1D, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array1D[npc.floating]: ...
@overload  # (real array-like, real array-like) -> floating array
def zmap(
    scores: onp.ToFloatND, compare: onp.ToFloatND, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[npc.floating]: ...
@overload  # (just complex vector-like, complex vector-like) -> floating vector
def zmap(
    scores: onp.ToJustComplex1D,
    compare: onp.ToComplex1D,
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[npc.complexfloating]: ...
@overload  # (complex vector-like, just complex vector-like) -> floating vector
def zmap(
    scores: onp.ToComplex1D,
    compare: onp.ToJustComplex1D,
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[npc.complexfloating]: ...
@overload  # (just complex array-like, complex array-like) -> floating array
def zmap(
    scores: onp.ToJustComplexND,
    compare: onp.ToComplexND,
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[npc.complexfloating]: ...
@overload  # (complex array-like, just complex array-like) -> floating array
def zmap(
    scores: onp.ToComplexND,
    compare: onp.ToJustComplexND,
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[npc.complexfloating]: ...

# TODO(jorenham): improve
def iqr(
    x: onp.ToFloatND,
    axis: int | Sequence[int] | None = None,
    rng: tuple[float, float] = (25, 75),
    scale: L["normal"] | onp.ToFloat | onp.ToFloatND = 1.0,
    nan_policy: NanPolicy = "propagate",
    interpolation: _InterpolationMethod = "linear",
    keepdims: bool = False,
) -> _FloatOrND: ...

# TODO(jorenham): improve
def median_abs_deviation(
    x: onp.ToFloatND,
    axis: int | None = 0,
    center: np.ufunc | _MADCenterFunc = ...,
    scale: L["normal"] | onp.ToFloat = 1.0,
    nan_policy: NanPolicy = "propagate",
) -> _FloatOrND: ...

#
def sigmaclip(a: onp.ToFloatND, low: float = 4.0, high: float = 4.0) -> SigmaclipResult: ...

# TODO(jorenham): improve
def trimboth(a: onp.ToFloatND, proportiontocut: float, axis: int | None = 0) -> onp.ArrayND[_Real0D]: ...

# TODO(jorenham): improve
def trim1(a: onp.ToFloatND, proportiontocut: float, tail: _TrimTail = "right", axis: int | None = 0) -> onp.ArrayND[_Real0D]: ...

# TODO(jorenham): improve
def trim_mean(a: onp.ToFloatND, proportiontocut: float, axis: int | None = 0) -> _FloatOrND: ...

# TODO(jorenham): improve
def f_oneway(
    *samples: onp.ToFloatND,
    axis: int | None = 0,
    keepdims: bool = False,
    nan_policy: NanPolicy = "propagate",
    equal_var: bool = True,
) -> F_onewayResult: ...

# TODO(jorenham): improve
def alexandergovern(
    *samples: onp.ToFloatND, axis: int | None = 0, keepdims: bool = False, nan_policy: NanPolicy = "propagate"
) -> AlexanderGovernResult: ...

#
@overload  # 1d +integer | ~float64, +floating
def pearsonr(
    x: onp.ToJustFloat64Strict1D | onp.ToIntStrict1D,
    y: onp.ToFloatStrict1D,
    *,
    axis: L[0, -1] | None = 0,
    alternative: Alternative = "two-sided",
    method: ResamplingMethod | None = None,
) -> PearsonRResult[np.float64, np.float64]: ...
@overload  # 1d +floating, +integer | ~float64
def pearsonr(
    x: onp.ToFloatStrict1D,
    y: onp.ToJustFloat64Strict1D | onp.ToIntStrict1D,
    *,
    axis: L[0, -1] | None = 0,
    alternative: Alternative = "two-sided",
    method: ResamplingMethod | None = None,
) -> PearsonRResult[np.float64, np.float64]: ...
@overload  # 1d +floating, +floating
def pearsonr(
    x: onp.ToFloatStrict1D,
    y: onp.ToFloatStrict1D,
    *,
    axis: L[0, -1] | None = 0,
    alternative: Alternative = "two-sided",
    method: ResamplingMethod | None = None,
) -> PearsonRResult[npc.floating, np.float64]: ...
@overload  # ?d +integer | ~float64, +floating, axis=None
def pearsonr(
    x: onp.ToJustFloat64_ND | onp.ToIntND,
    y: onp.ToFloatND,
    *,
    axis: None,
    alternative: Alternative = "two-sided",
    method: ResamplingMethod | None = None,
) -> PearsonRResult[np.float64, np.float64]: ...
@overload  # ?d +floating, +integer | ~float64, axis=None
def pearsonr(
    x: onp.ToFloatND,
    y: onp.ToJustFloat64_ND | onp.ToIntND,
    *,
    axis: None,
    alternative: Alternative = "two-sided",
    method: ResamplingMethod | None = None,
) -> PearsonRResult[np.float64, np.float64]: ...
@overload  # ?d +floating, +floating, axis=None
def pearsonr(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    *,
    axis: None,
    alternative: Alternative = "two-sided",
    method: ResamplingMethod | None = None,
) -> PearsonRResult[npc.floating, np.float64]: ...
@overload  # >=2d +integer | ~float64, +floating
def pearsonr(
    x: onp.CanArray[onp.AtLeast2D, np.dtype[npc.integer | npc.floating]] | Sequence[onp.SequenceND[float]],
    y: onp.CanArray[onp.AtLeast2D, np.dtype[npc.integer | np.float64]] | Sequence[onp.SequenceND[float]],
    *,
    axis: int = 0,
    alternative: Alternative = "two-sided",
    method: ResamplingMethod | None = None,
) -> PearsonRResult[onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
@overload  # >=2d +floating, +integer | ~float64
def pearsonr(
    x: onp.CanArray[onp.AtLeast2D, np.dtype[npc.integer | np.float64]] | Sequence[onp.SequenceND[float]],
    y: onp.CanArray[onp.AtLeast2D, np.dtype[npc.integer | npc.floating]] | Sequence[onp.SequenceND[float]],
    *,
    axis: int = 0,
    alternative: Alternative = "two-sided",
    method: ResamplingMethod | None = None,
) -> PearsonRResult[onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
@overload  # >=2d +floating, +floating
def pearsonr(
    x: onp.CanArray[onp.AtLeast2D, np.dtype[npc.integer | npc.floating]] | Sequence[onp.SequenceND[float]],
    y: onp.CanArray[onp.AtLeast2D, np.dtype[npc.integer | npc.floating]] | Sequence[onp.SequenceND[float]],
    *,
    axis: int = 0,
    alternative: Alternative = "two-sided",
    method: ResamplingMethod | None = None,
) -> PearsonRResult[onp.ArrayND[npc.floating], onp.ArrayND[np.float64]]: ...
@overload  # fallback
def pearsonr(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    *,
    axis: int | None = 0,
    alternative: Alternative = "two-sided",
    method: ResamplingMethod | None = None,
) -> PearsonRResult[npc.floating, np.float64] | PearsonRResult[onp.ArrayND[npc.floating], onp.ArrayND[np.float64]]: ...

#
def fisher_exact(
    table: onp.ArrayND[_Real0D], alternative: Alternative | None = None, *, method: ResamplingMethod | None = None
) -> SignificanceResult[float]: ...

# TODO(jorenham): improve like `pearsonr` (but always return `f64`)
def spearmanr(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None = None,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
) -> SignificanceResult: ...

# TODO(jorenham): improve like `pearsonr` (but return `SignificanceResult`, not `PearsonRResult`)
@overload
def pointbiserialr(
    x: onp.ToBoolND, y: onp.ToFloatND, *, axis: None, nan_policy: NanPolicy = "propagate", keepdims: L[False] = False
) -> SignificanceResult[np.float64]: ...
@overload
def pointbiserialr(
    x: onp.ToBoolStrict1D,
    y: onp.ToFloatStrict1D,
    *,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[np.float64]: ...
@overload
def pointbiserialr(
    x: onp.ToBoolStrict2D,
    y: onp.ToFloatStrict2D,
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[onp.Array1D[np.float64]]: ...
@overload
def pointbiserialr(
    x: onp.ToBoolStrict3D,
    y: onp.ToFloatStrict3D,
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[onp.Array2D[np.float64]]: ...
@overload
def pointbiserialr(
    x: onp.ToBoolND, y: onp.ToFloatND, *, axis: int | None = 0, nan_policy: NanPolicy = "propagate", keepdims: L[True]
) -> SignificanceResult[onp.ArrayND[np.float64]]: ...
@overload
def pointbiserialr(
    x: onp.ToBoolND, y: onp.ToFloatND, *, axis: int | None = 0, nan_policy: NanPolicy = "propagate", keepdims: bool = False
) -> SignificanceResult[np.float64 | Any]: ...

#
@overload
def kendalltau(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    *,
    axis: None,
    keepdims: L[False] = False,
    method: _KendallTauMethod = "auto",
    variant: _KendallTauVariant = "b",
    alternative: Alternative = "two-sided",
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[np.float64]: ...
@overload
def kendalltau(
    x: onp.ToFloatStrict1D,
    y: onp.ToFloatStrict1D,
    *,
    axis: int | None = 0,
    keepdims: L[False] = False,
    method: _KendallTauMethod = "auto",
    variant: _KendallTauVariant = "b",
    alternative: Alternative = "two-sided",
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[np.float64]: ...
@overload
def kendalltau(
    x: onp.ToFloatStrict2D,
    y: onp.ToFloatStrict2D,
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    method: _KendallTauMethod = "auto",
    variant: _KendallTauVariant = "b",
    alternative: Alternative = "two-sided",
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[onp.Array1D[np.float64]]: ...
@overload
def kendalltau(
    x: onp.ToFloatStrict3D,
    y: onp.ToFloatStrict3D,
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    method: _KendallTauMethod = "auto",
    variant: _KendallTauVariant = "b",
    alternative: Alternative = "two-sided",
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[onp.Array2D[np.float64]]: ...
@overload
def kendalltau(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    *,
    axis: int | None = 0,
    keepdims: L[True],
    method: _KendallTauMethod = "auto",
    variant: _KendallTauVariant = "b",
    alternative: Alternative = "two-sided",
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[onp.ArrayND[np.float64]]: ...
@overload
def kendalltau(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    *,
    axis: int | None = 0,
    keepdims: bool = False,
    method: _KendallTauMethod = "auto",
    variant: _KendallTauVariant = "b",
    alternative: Alternative = "two-sided",
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[np.float64 | Any]: ...

#
@overload
def weightedtau(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    rank: onp.ToInt | onp.ToIntND = True,
    weigher: _Weigher | None = None,
    additive: bool = True,
    *,
    axis: None,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[np.float64]: ...
@overload
def weightedtau(
    x: onp.ToFloatStrict1D,
    y: onp.ToFloatStrict1D,
    rank: onp.ToInt | onp.ToIntND = True,
    weigher: _Weigher | None = None,
    additive: bool = True,
    *,
    axis: int | None = 0,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[np.float64]: ...
@overload
def weightedtau(
    x: onp.ToFloatStrict2D,
    y: onp.ToFloatStrict2D,
    rank: onp.ToInt | onp.ToIntND = True,
    weigher: _Weigher | None = None,
    additive: bool = True,
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[onp.Array1D[np.float64]]: ...
@overload
def weightedtau(
    x: onp.ToFloatStrict3D,
    y: onp.ToFloatStrict3D,
    rank: onp.ToInt | onp.ToIntND = True,
    weigher: _Weigher | None = None,
    additive: bool = True,
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[onp.Array2D[np.float64]]: ...
@overload
def weightedtau(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    rank: onp.ToInt | onp.ToIntND = True,
    weigher: _Weigher | None = None,
    additive: bool = True,
    *,
    axis: int | None = 0,
    keepdims: L[True],
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[onp.ArrayND[np.float64]]: ...
@overload
def weightedtau(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    rank: onp.ToInt | onp.ToIntND = True,
    weigher: _Weigher | None = None,
    additive: bool = True,
    *,
    axis: int | None = 0,
    keepdims: bool = False,
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[np.float64 | Any]: ...

#
def pack_TtestResult(
    statistic: _FloatOrArrayT,
    pvalue: _FloatOrArrayT,
    df: _FloatOrArrayT,
    alternative: Alternative,
    standard_error: _FloatOrArrayT,
    estimate: _FloatOrArrayT,
) -> TtestResult[_FloatOrArrayT]: ...  # undocumented

#
def unpack_TtestResult(
    res: TtestResult[_FloatOrArrayT], _: int
) -> tuple[
    _FloatOrArrayT,  # statistic
    _FloatOrArrayT,  # pvalue
    _FloatOrArrayT,  # df
    Alternative,  # _alternative
    _FloatOrArrayT,  # _standard_error
    _FloatOrArrayT,  # _estimate
]: ...  # undocumented

# TODO(jorenham): improve
def ttest_1samp(
    a: onp.ToFloatND,
    popmean: onp.ToFloat | onp.ToFloatND,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: bool = False,
) -> TtestResult: ...

# TODO(jorenham): improve
def ttest_ind_from_stats(
    mean1: onp.ToFloat | onp.ToFloatND,
    std1: onp.ToFloat | onp.ToFloatND,
    nobs1: onp.ToInt | onp.ToIntND,
    mean2: onp.ToFloat | onp.ToFloatND,
    std2: onp.ToFloat | onp.ToFloatND,
    nobs2: onp.ToInt | onp.ToIntND,
    equal_var: bool = True,
    alternative: Alternative = "two-sided",
) -> Ttest_indResult: ...

#
@overload
def ttest_ind(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    *,
    axis: int | None = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    permutations: None = None,
    random_state: None = None,
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: bool = False,
) -> TtestResult: ...
@overload
@deprecated(
    "Argument `random_state` is deprecated, and will be removed in SciPy 1.17. "
    "Use `method to perform a permutation test."
)  # fmt: skip
def ttest_ind(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    *,
    axis: int | None = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    permutations: None = None,
    random_state: onp.random.ToRNG | None,
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: bool = False,
) -> TtestResult: ...
@overload
@deprecated(
    "Argument `permutations` is deprecated, and will be removed in SciPy 1.17. "
    "Use method` to perform a permutation test."
)  # fmt: skip
def ttest_ind(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    *,
    axis: int | None = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    permutations: onp.ToFloat,
    random_state: None = None,
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: bool = False,
) -> TtestResult: ...
@overload
@deprecated(
    "Arguments {'random_state', 'permutations'} are deprecated, and will be removed in SciPy 1.17. "
    "Use `method` to perform a permutation test."
)
def ttest_ind(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    *,
    axis: int | None = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    permutations: onp.ToFloat,
    random_state: onp.random.ToRNG | None,
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: bool = False,
) -> TtestResult: ...

# TODO(jorenham): improve
def ttest_rel(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: bool = False,
) -> TtestResult: ...

#
@overload
def power_divergence(
    f_obs: onp.ToFloatStrict1D,
    f_exp: onp.ToFloatStrict1D | None = None,
    ddof: int = 0,
    axis: int | None = 0,
    lambda_: PowerDivergenceStatistic | float | None = None,
    *,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> Power_divergenceResult[np.float64]: ...
@overload
def power_divergence(
    f_obs: onp.ToFloatND,
    f_exp: onp.ToFloatND | None,
    ddof: int,
    axis: None,
    lambda_: PowerDivergenceStatistic | float | None = None,
    *,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> Power_divergenceResult[np.float64]: ...
@overload
def power_divergence(
    f_obs: onp.ToFloatND,
    f_exp: onp.ToFloatND | None = None,
    ddof: int = 0,
    *,
    axis: None,
    lambda_: PowerDivergenceStatistic | float | None = None,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> Power_divergenceResult[np.float64]: ...
@overload
def power_divergence(
    f_obs: onp.ToFloatND,
    f_exp: onp.ToFloatND | None = None,
    ddof: int = 0,
    axis: int | None = 0,
    lambda_: PowerDivergenceStatistic | float | None = None,
    *,
    keepdims: L[True],
    nan_policy: NanPolicy = "propagate",
) -> Power_divergenceResult[onp.ArrayND[np.float64]]: ...
@overload
def power_divergence(
    f_obs: onp.ToFloatND,
    f_exp: onp.ToFloatND | None = None,
    ddof: int = 0,
    axis: int | None = 0,
    lambda_: PowerDivergenceStatistic | float | None = None,
    *,
    keepdims: bool = False,
    nan_policy: NanPolicy = "propagate",
) -> Power_divergenceResult[np.float64 | Any]: ...

#
@overload
def chisquare(
    f_obs: onp.ToFloatStrict1D,
    f_exp: onp.ToFloatStrict1D | None = None,
    ddof: int = 0,
    axis: int | None = 0,
    *,
    sum_check: bool = True,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> Power_divergenceResult[np.float64]: ...
@overload
def chisquare(
    f_obs: onp.ToFloatND,
    f_exp: onp.ToFloatND | None,
    ddof: int,
    axis: None,
    *,
    sum_check: bool = True,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> Power_divergenceResult[np.float64]: ...
@overload
def chisquare(
    f_obs: onp.ToFloatND,
    f_exp: onp.ToFloatND | None = None,
    ddof: int = 0,
    *,
    axis: None,
    sum_check: bool = True,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> Power_divergenceResult[np.float64]: ...
@overload
def chisquare(
    f_obs: onp.ToFloatND,
    f_exp: onp.ToFloatND | None = None,
    ddof: int = 0,
    axis: int | None = 0,
    *,
    sum_check: bool = True,
    keepdims: L[True],
    nan_policy: NanPolicy = "propagate",
) -> Power_divergenceResult[onp.ArrayND[np.float64]]: ...
@overload
def chisquare(
    f_obs: onp.ToFloatND,
    f_exp: onp.ToFloatND | None = None,
    ddof: int = 0,
    axis: int | None = 0,
    *,
    sum_check: bool = True,
    keepdims: bool = False,
    nan_policy: NanPolicy = "propagate",
) -> Power_divergenceResult[np.float64 | Any]: ...

# TODO(jorenham): improve
def ks_1samp(
    x: onp.ToFloatND,
    cdf: Callable[[float], float | _Real0D],
    args: tuple[object, ...] = (),
    alternative: Alternative = "two-sided",
    method: _KS1TestMethod = "auto",
    *,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> KstestResult: ...

# TODO(jorenham): improve
def ks_2samp(
    data1: onp.ToFloatND,
    data2: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    method: _KS2TestMethod = "auto",
    *,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> KstestResult: ...

# TODO(jorenham): improve
def kstest(
    rvs: str | onp.ToFloatND | _RVSCallable,
    cdf: str | onp.ToFloatND | Callable[[float], float | npc.floating],
    args: tuple[object, ...] = (),
    N: int = 20,
    alternative: Alternative = "two-sided",
    method: _KS1TestMethod = "auto",
    *,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> KstestResult: ...

#
def tiecorrect(rankvals: onp.ToIntND) -> float: ...

# TODO(jorenham): improve
def ranksums(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    *,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> RanksumsResult: ...

# TODO(jorenham): improve
def kruskal(
    *samples: onp.ToFloatND, nan_policy: NanPolicy = "propagate", axis: int | None = 0, keepdims: bool = False
) -> KruskalResult: ...

# TODO(jorenham): improve
def friedmanchisquare(
    *samples: onp.ToFloatND, axis: int | None = 0, nan_policy: NanPolicy = "propagate", keepdims: bool = False
) -> FriedmanchisquareResult: ...

# TODO(jorenham): improve
def brunnermunzel(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    distribution: L["t", "normal"] = "t",
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: bool = False,
    axis: int | None = 0,
) -> BrunnerMunzelResult: ...

# TODO(jorenham): improve
def combine_pvalues(
    pvalues: onp.ToFloatND,
    method: _CombinePValuesMethod = "fisher",
    weights: onp.ToFloatND | None = None,
    *,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> SignificanceResult: ...

#
def quantile_test_iv(
    x: onp.ToFloatND, q: float | _Real0D, p: float | npc.floating, alternative: Alternative
) -> tuple[onp.ArrayND[_Real0D], _Real0D, npc.floating, Alternative]: ...  # undocumented

#
def quantile_test(
    x: onp.ToFloatND, *, q: float | _Real0D = 0, p: float | npc.floating = 0.5, alternative: Alternative = "two-sided"
) -> QuantileTestResult: ...

#
def wasserstein_distance_nd(
    u_values: onp.ToFloatND,
    v_values: onp.ToFloatND,
    u_weights: onp.ToFloatND | None = None,
    v_weights: onp.ToFloatND | None = None,
) -> np.float64: ...
def wasserstein_distance(
    u_values: onp.ToFloatND,
    v_values: onp.ToFloatND,
    u_weights: onp.ToFloatND | None = None,
    v_weights: onp.ToFloatND | None = None,
) -> np.float64: ...
def energy_distance(
    u_values: onp.ToFloatND,
    v_values: onp.ToFloatND,
    u_weights: onp.ToFloatND | None = None,
    v_weights: onp.ToFloatND | None = None,
) -> np.float64: ...

# TODO(jorenham): improve
def rankdata(
    a: onp.ToFloatND, method: _RankMethod = "average", *, axis: int | None = None, nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[_Real0D]: ...

#
def expectile(a: onp.ToFloatND, alpha: float = 0.5, *, weights: onp.ToFloatND | None = None) -> np.float64: ...

#
@overload
def linregress(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    *,
    axis: None,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> LinregressResult[np.float64]: ...
@overload
def linregress(
    x: onp.ToFloatStrict1D,
    y: onp.ToFloatStrict1D,
    alternative: Alternative = "two-sided",
    *,
    axis: int | None = 0,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> LinregressResult[np.float64]: ...
@overload
def linregress(
    x: onp.ToFloatStrict2D,
    y: onp.ToFloatStrict2D,
    alternative: Alternative = "two-sided",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> LinregressResult[onp.Array1D[np.float64]]: ...
@overload
def linregress(
    x: onp.ToFloatStrict3D,
    y: onp.ToFloatStrict3D,
    alternative: Alternative = "two-sided",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> LinregressResult[onp.Array2D[np.float64]]: ...
@overload
def linregress(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    *,
    axis: int | None = 0,
    keepdims: L[True],
    nan_policy: NanPolicy = "propagate",
) -> LinregressResult[onp.ArrayND[np.float64]]: ...
@overload
def linregress(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    *,
    axis: int | None = 0,
    keepdims: bool = False,
    nan_policy: NanPolicy = "propagate",
) -> LinregressResult[np.float64 | Any]: ...

#
@deprecated(
    "`scipy.stats.find_repeats` is deprecated as of SciPy 1.15.0 and will be removed in SciPy 1.17.0. "
    "Please use `numpy.unique`/`numpy.unique_counts` instead."
)
def find_repeats(arr: onp.ToFloatND) -> RepeatedResults: ...

# NOTE: `lmoment` is currently numerically unstable after `order > 16`.
# See https://github.com/jorenham/Lmo/ for a more stable implementation that additionally supports generalized trimmed TL-moments,
# multivariate L- and TL-comoments, theoretical L- and TL-moments or `scipy.stats` distributions, and much more ;)

# TODO(jorenham): improve dtype
@overload  # sample: 1-d, order: 0-d, keepdims: falsy
def lmoment(
    sample: onp.ToFloatStrict1D,
    order: _LMomentOrder,
    *,
    axis: L[0, -1] | None = 0,
    keepdims: onp.ToFalse = False,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> np.float32 | np.float64: ...
@overload  # sample: 1-d, order: 0-d, keepdims: truthy
def lmoment(
    sample: onp.ToFloatStrict1D,
    order: _LMomentOrder,
    *,
    axis: L[0, -1] | None = 0,
    keepdims: onp.ToTrue,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float32 | np.float64]: ...
@overload  # sample: 1-d, order: 1-d, keepdims: falsy
def lmoment(
    sample: onp.ToFloatStrict1D,
    order: _LMomentOrder1D | None = None,
    *,
    axis: L[0, -1] | None = 0,
    keepdims: onp.ToFalse = False,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float32 | np.float64]: ...
@overload  # sample: 1-d, order: 1-d, keepdims: truthy
def lmoment(
    sample: onp.ToFloatStrict1D,
    order: _LMomentOrder1D | None = None,
    *,
    axis: L[0, -1] | None = 0,
    keepdims: onp.ToTrue,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array2D[np.float32 | np.float64]: ...
@overload  # sample: 2-d, order: 0-d, keepdims: falsy
def lmoment(
    sample: onp.ToFloatStrict2D,
    order: _LMomentOrder,
    *,
    axis: int = 0,
    keepdims: onp.ToFalse = False,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float32 | np.float64]: ...
@overload  # sample: 2-d, order: 0-d, keepdims: truthy
def lmoment(
    sample: onp.ToFloatStrict2D,
    order: _LMomentOrder,
    *,
    axis: int | None = 0,
    keepdims: onp.ToTrue,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array2D[np.float32 | np.float64]: ...
@overload  # sample: 2-d, order: 1-d, keepdims: falsy
def lmoment(
    sample: onp.ToFloatStrict2D,
    order: _LMomentOrder1D | None = None,
    *,
    axis: int = 0,
    keepdims: onp.ToFalse = False,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array2D[np.float32 | np.float64]: ...
@overload  # sample: 2-d, order: 1-d, keepdims: truthy
def lmoment(
    sample: onp.ToFloatStrict2D,
    order: _LMomentOrder1D | None = None,
    *,
    axis: int | None = 0,
    keepdims: onp.ToTrue,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array3D[np.float32 | np.float64]: ...
@overload  # sample: 3-d, order: 0-d, keepdims: falsy
def lmoment(
    sample: onp.ToFloatStrict3D,
    order: _LMomentOrder,
    *,
    axis: int = 0,
    keepdims: onp.ToFalse = False,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array2D[np.float32 | np.float64]: ...
@overload  # sample: 3-d, order: 0-d, keepdims: truthy
def lmoment(
    sample: onp.ToFloatStrict3D,
    order: _LMomentOrder,
    *,
    axis: int | None = 0,
    keepdims: onp.ToTrue,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array3D[np.float32 | np.float64]: ...
@overload  # sample: 2-d, order: 1-d, keepdims: falsy
def lmoment(
    sample: onp.ToFloatStrict3D,
    order: _LMomentOrder1D | None = None,
    *,
    axis: int = 0,
    keepdims: onp.ToFalse = False,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array3D[np.float32 | np.float64]: ...
@overload  # sample: 3-d, order: 1-d, keepdims: truthy
def lmoment(
    sample: onp.ToFloatStrict3D,
    order: _LMomentOrder1D | None = None,
    *,
    axis: int | None = 0,
    keepdims: onp.ToTrue,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array[tuple[int, int, int, int], np.float32 | np.float64]: ...
@overload  # sample: N-d, order: 0-d, keepdims: falsy, axis: None
def lmoment(
    sample: onp.ToFloatND,
    order: _LMomentOrder,
    *,
    axis: None,
    keepdims: onp.ToFalse = False,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> np.float32 | np.float64: ...
@overload  # sample: N-d, order: 1-d, keepdims: falsy, axis: None
def lmoment(
    sample: onp.ToFloatND,
    order: _LMomentOrder1D | None = None,
    *,
    axis: None,
    keepdims: onp.ToFalse = False,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float32 | np.float64]: ...
@overload  # sample: N-d, keepdims: truthy
def lmoment(
    sample: onp.ToFloatND,
    order: _LMomentOrder | _LMomentOrder1D | None = None,
    *,
    axis: int | None = 0,
    keepdims: onp.ToTrue,
    sorted: op.CanBool = False,
    standardize: op.CanBool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float32 | np.float64]: ...
