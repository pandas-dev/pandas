from collections.abc import Callable, Sequence
from dataclasses import dataclass
from types import ModuleType
from typing import Any, Generic, Literal as L, Never, Protocol, Self, TypeAlias, overload, type_check_only
from typing_extensions import NamedTuple, TypeVar

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
    "_FloatOrArrayT_co",
    bound=float | npc.floating | onp.ArrayND[npc.floating, Any],
    default=float | onp.ArrayND[np.float64],
    covariant=True,
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

_RealLimits: TypeAlias = tuple[float | _Real0D, float | _Real0D]
_Weigher: TypeAlias = Callable[[int], float | _Real0D]

_JustAnyShape: TypeAlias = tuple[Never, Never, Never, Never]  # workaround for https://github.com/microsoft/pyright/issues/10232
_AsFloat64_1D: TypeAlias = onp.ToArrayStrict1D[float, npc.floating64 | npc.integer]
_AsFloat64_2D: TypeAlias = onp.ToArrayStrict2D[float, npc.floating64 | npc.integer]
_AsFloat64_ND: TypeAlias = onp.ToArrayND[float, npc.floating64 | npc.integer]
_AsFloat32_1D: TypeAlias = onp.ToArrayStrict1D[np.float32, np.float32 | np.float16]
_AsFloat32_2D: TypeAlias = onp.ToArrayStrict2D[np.float32, np.float32 | np.float16]
_AsFloat32_ND: TypeAlias = onp.ToArrayND[Never, np.float32 | np.float16]

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

#
@overload  # ?d ~f64, order: 0d
def moment(
    a: onp.ArrayND[npc.floating64 | npc.integer | np.bool_, _JustAnyShape],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # ?d ~T: floating, order: 0d
def moment(
    a: onp.ArrayND[_FloatT, _JustAnyShape],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> _FloatT | onp.ArrayND[_FloatT]: ...
@overload  # 1d ~f64, order: 0d
def moment(
    a: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool_],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # 1d ~T: floating, order: 0d
def moment(
    a: onp.ToArrayStrict1D[_FloatT, _FloatT],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> _FloatT: ...
@overload  # 2d ~f64, order: 0d
def moment(
    a: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool_],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~T: floating, order: 0d
def moment(
    a: onp.ToArrayStrict2D[_FloatT, _FloatT],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.Array1D[_FloatT]: ...
@overload  # 3d ~f64, order: 0d
def moment(
    a: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool_],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.Array2D[np.float64]: ...
@overload  # 3d ~T: floating, order: 0d
def moment(
    a: onp.ToArrayStrict3D[_FloatT, _FloatT],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.Array2D[_FloatT]: ...
@overload  # nd ~f64, order: 0d
def moment(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64] | Any: ...
@overload  # nd ~f64, order: 0d, axis=None  (positional)
def moment(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    order: int,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # nd ~f64, order: 0d, axis=None  (keyword)
def moment(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    order: int = 1,
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    center: float | None = None,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # nd ~f64, order: nd
def moment(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    order: onp.ToIntND,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64]: ...
@overload  # nd ~f64, keepdims=True
def moment(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    order: int | onp.ToIntND = 1,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # nd ~T: floating, order: 0d
def moment(
    a: onp.ToArrayND[_FloatT, _FloatT],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.ArrayND[_FloatT] | Any: ...
@overload  # nd ~T: floating, order: 0d, axis=None  (positional)
def moment(
    a: onp.ToArrayND[_FloatT, _FloatT],
    order: int,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> _FloatT: ...
@overload  # nd ~T: floating, order: 0d, axis=None  (keyword)
def moment(
    a: onp.ToArrayND[_FloatT, _FloatT],
    order: int = 1,
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    center: float | None = None,
    keepdims: L[False] = False,
) -> _FloatT: ...
@overload  # nd ~T: floating, order: nd
def moment(
    a: onp.ToArrayND[_FloatT, _FloatT],
    order: onp.ToIntND,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.ArrayND[_FloatT]: ...
@overload  # nd ~T: floating, keepdims=True
def moment(
    a: onp.ToArrayND[_FloatT, _FloatT],
    order: int | onp.ToIntND = 1,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[True],
) -> onp.ArrayND[_FloatT]: ...
@overload  # nd +floating, order: 0d
def moment(
    a: onp.ToFloatND,
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> npc.floating | onp.ArrayND[npc.floating]: ...
@overload  # nd +floating, order: 0d, axis=None  (positional)
def moment(
    a: onp.ToFloatND,
    order: int,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> npc.floating: ...
@overload  # nd +floating, order: 0d, axis=None  (keyword)
def moment(
    a: onp.ToFloatND,
    order: int = 1,
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    center: float | None = None,
    keepdims: L[False] = False,
) -> npc.floating: ...
@overload  # nd +floating, order: nd
def moment(
    a: onp.ToFloatND,
    order: onp.ToIntND,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.ArrayND[npc.floating]: ...
@overload  # nd +floating, keepdims=True
def moment(
    a: onp.ToFloatND,
    order: int | onp.ToIntND = 1,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[True],
) -> onp.ArrayND[npc.floating]: ...

# keep in sync with kurtosis
@overload  # ?d ~f64
def skew(
    a: onp.ArrayND[npc.floating64 | npc.integer | np.bool_, _JustAnyShape],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # ?d ~T
def skew(
    a: onp.ArrayND[_FloatT, _JustAnyShape],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> _FloatT | onp.ArrayND[_FloatT]: ...
@overload  # 1d ~f64
def skew(
    a: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool_],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # 1d ~T
def skew(
    a: onp.ToArrayStrict1D[_FloatT, _FloatT],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> _FloatT: ...
@overload  # 2d ~f64
def skew(
    a: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool_],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~T
def skew(
    a: onp.ToArrayStrict2D[_FloatT, _FloatT],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array1D[_FloatT]: ...
@overload  # 3d ~f64
def skew(
    a: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool_],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array2D[np.float64]: ...
@overload  # 3d ~T
def skew(
    a: onp.ToArrayStrict3D[_FloatT, _FloatT],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array2D[_FloatT]: ...
@overload  # nd ~f64
def skew(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64] | Any: ...
@overload  # nd ~f64, axis=None
def skew(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    axis: None,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # nd ~f64, keepdims=True
def skew(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    axis: int | None = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # nd ~T
def skew(
    a: onp.ToArrayND[_FloatT, _FloatT],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.ArrayND[_FloatT] | Any: ...
@overload  # nd ~T, axis=None
def skew(
    a: onp.ToArrayND[_FloatT, _FloatT],
    axis: None,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> _FloatT: ...
@overload  # nd ~T, keepdims=True
def skew(
    a: onp.ToArrayND[_FloatT, _FloatT],
    axis: int | None = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[_FloatT]: ...
@overload  # nd +floating
def skew(
    a: onp.ToFloatND, axis: int = 0, bias: bool = True, nan_policy: NanPolicy = "propagate", *, keepdims: L[False] = False
) -> onp.ArrayND[npc.floating] | Any: ...
@overload  # nd +floating, axis=None
def skew(
    a: onp.ToFloatND, axis: None, bias: bool = True, nan_policy: NanPolicy = "propagate", *, keepdims: L[False] = False
) -> npc.floating: ...
@overload  # nd +floating, keepdims=True
def skew(
    a: onp.ToFloatND, axis: int | None = 0, bias: bool = True, nan_policy: NanPolicy = "propagate", *, keepdims: L[True]
) -> onp.ArrayND[npc.floating]: ...

# keep in sync with skew
@overload  # ?d ~f64
def kurtosis(
    a: onp.ArrayND[npc.floating64 | npc.integer | np.bool_, _JustAnyShape],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # ?d ~T
def kurtosis(
    a: onp.ArrayND[_FloatT, _JustAnyShape],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> _FloatT | onp.ArrayND[_FloatT]: ...
@overload  # 1d ~f64
def kurtosis(
    a: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool_],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # 1d ~T
def kurtosis(
    a: onp.ToArrayStrict1D[_FloatT, _FloatT],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> _FloatT: ...
@overload  # 2d ~f64
def kurtosis(
    a: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool_],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~T
def kurtosis(
    a: onp.ToArrayStrict2D[_FloatT, _FloatT],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array1D[_FloatT]: ...
@overload  # 3d ~f64
def kurtosis(
    a: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool_],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array2D[np.float64]: ...
@overload  # 3d ~T
def kurtosis(
    a: onp.ToArrayStrict3D[_FloatT, _FloatT],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array2D[_FloatT]: ...
@overload  # nd ~f64
def kurtosis(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64] | Any: ...
@overload  # nd ~f64, axis=None
def kurtosis(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    axis: None,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # nd ~f64, keepdims=True
def kurtosis(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    axis: int | None = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # nd ~T
def kurtosis(
    a: onp.ToArrayND[_FloatT, _FloatT],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.ArrayND[_FloatT] | Any: ...
@overload  # nd ~T, axis=None
def kurtosis(
    a: onp.ToArrayND[_FloatT, _FloatT],
    axis: None,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> _FloatT: ...
@overload  # nd ~T, keepdims=True
def kurtosis(
    a: onp.ToArrayND[_FloatT, _FloatT],
    axis: int | None = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[_FloatT]: ...
@overload  # nd +floating
def kurtosis(
    a: onp.ToFloatND,
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.ArrayND[npc.floating] | Any: ...
@overload  # nd +floating, axis=None
def kurtosis(
    a: onp.ToFloatND,
    axis: None,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> npc.floating: ...
@overload  # nd +floating, keepdims=True
def kurtosis(
    a: onp.ToFloatND,
    axis: int | None = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[npc.floating]: ...

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

#
@overload
def median_abs_deviation(
    x: onp.ToFloatStrict1D,
    axis: int = 0,
    center: np.ufunc | _MADCenterFunc | None = None,
    scale: L["normal"] | float = 1.0,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload
def median_abs_deviation(
    x: onp.ToFloatND,
    axis: None,
    center: np.ufunc | _MADCenterFunc | None = None,
    scale: L["normal"] | float = 1.0,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64]: ...
@overload
def median_abs_deviation(
    x: onp.ToFloatND,
    axis: int = 0,
    center: np.ufunc | _MADCenterFunc | None = None,
    scale: L["normal"] | float = 1.0,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64] | Any: ...
@overload
def median_abs_deviation(
    x: onp.ToFloatND,
    axis: int | None = 0,
    center: np.ufunc | _MADCenterFunc | None = None,
    scale: L["normal"] | float = 1.0,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...

#
def sigmaclip(a: onp.ToFloatND, low: float = 4.0, high: float = 4.0) -> SigmaclipResult: ...

# TODO(jorenham): improve
def trimboth(a: onp.ToFloatND, proportiontocut: float, axis: int | None = 0) -> onp.ArrayND[_Real0D]: ...

# TODO(jorenham): improve
def trim1(a: onp.ToFloatND, proportiontocut: float, tail: _TrimTail = "right", axis: int | None = 0) -> onp.ArrayND[_Real0D]: ...

#
@overload
def trim_mean(
    a: onp.ToFloatStrict1D,
    proportiontocut: float,
    axis: int | None = 0,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload
def trim_mean(
    a: onp.ToFloatND, proportiontocut: float, axis: None, *, nan_policy: NanPolicy = "propagate", keepdims: L[False] = False
) -> np.float64: ...
@overload
def trim_mean(
    a: onp.ToFloatND, proportiontocut: float, axis: int = 0, *, nan_policy: NanPolicy = "propagate", keepdims: L[False] = False
) -> _FloatOrND: ...
@overload
def trim_mean(
    a: onp.ToFloatND, proportiontocut: float, axis: int = 0, *, nan_policy: NanPolicy = "propagate", keepdims: L[True]
) -> onp.ArrayND[np.float64]: ...

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
@overload  # ?d, ?d
def spearmanr(
    a: onp.ArrayND[npc.floating | npc.integer | np.bool_, _JustAnyShape],
    b: onp.ArrayND[npc.floating | npc.integer | np.bool_, _JustAnyShape],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
) -> SignificanceResult[np.float64 | onp.Array2D[np.float64]]: ...
@overload  # 1d, 1d
def spearmanr(
    a: onp.ToFloatStrict1D,
    b: onp.ToFloatStrict1D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
) -> SignificanceResult[np.float64]: ...
@overload  # 2d, 2d
def spearmanr(
    a: onp.ToFloatStrict2D,
    b: onp.ToFloatStrict2D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
) -> SignificanceResult[onp.Array2D[np.float64]]: ...
@overload  # axis=None
def spearmanr(
    a: onp.ToFloatND, b: onp.ToFloatND, axis: None, nan_policy: NanPolicy = "propagate", alternative: Alternative = "two-sided"
) -> SignificanceResult[np.float64]: ...
@overload  # 2d, None
def spearmanr(
    a: onp.ToFloat2D, b: None = None, axis: int = 0, nan_policy: NanPolicy = "propagate", alternative: Alternative = "two-sided"
) -> SignificanceResult[np.float64 | onp.Array2D[np.float64]]: ...

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
@overload  # nd, axis=None (default)
def kendalltau(
    x: onp.ToComplexND,
    y: onp.ToComplexND,
    *,
    axis: None = None,
    keepdims: L[False] = False,
    method: _KendallTauMethod = "auto",
    variant: _KendallTauVariant = "b",
    alternative: Alternative = "two-sided",
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[np.float64]: ...
@overload  # ?d, axis: int
def kendalltau(
    x: onp.ArrayND[npc.number | np.bool_, _JustAnyShape],
    y: onp.ArrayND[npc.number | np.bool_, _JustAnyShape],
    *,
    axis: int,
    keepdims: L[False] = False,
    method: _KendallTauMethod = "auto",
    variant: _KendallTauVariant = "b",
    alternative: Alternative = "two-sided",
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[onp.ArrayND[np.float64] | Any]: ...
@overload  # 1d, axis: int
def kendalltau(
    x: onp.ToComplexStrict1D,
    y: onp.ToComplexStrict1D,
    *,
    axis: int,
    keepdims: L[False] = False,
    method: _KendallTauMethod = "auto",
    variant: _KendallTauVariant = "b",
    alternative: Alternative = "two-sided",
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[np.float64]: ...
@overload  # 2d, axis: int
def kendalltau(
    x: onp.ToComplexStrict2D,
    y: onp.ToComplexStrict2D,
    *,
    axis: int,
    keepdims: L[False] = False,
    method: _KendallTauMethod = "auto",
    variant: _KendallTauVariant = "b",
    alternative: Alternative = "two-sided",
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[onp.Array1D[np.float64]]: ...
@overload  # 3d, axis: int
def kendalltau(
    x: onp.ToComplexStrict3D,
    y: onp.ToComplexStrict3D,
    *,
    axis: int,
    keepdims: L[False] = False,
    method: _KendallTauMethod = "auto",
    variant: _KendallTauVariant = "b",
    alternative: Alternative = "two-sided",
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[onp.Array2D[np.float64]]: ...
@overload  # nd, axis: int
def kendalltau(
    x: onp.ToComplexND,
    y: onp.ToComplexND,
    *,
    axis: int,
    keepdims: bool = False,
    method: _KendallTauMethod = "auto",
    variant: _KendallTauVariant = "b",
    alternative: Alternative = "two-sided",
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[onp.ArrayND[np.float64] | Any]: ...
@overload  # ?d, keepdims=True
def kendalltau(
    x: onp.ToComplexND,
    y: onp.ToComplexND,
    *,
    axis: int | None = None,
    keepdims: L[True],
    method: _KendallTauMethod = "auto",
    variant: _KendallTauVariant = "b",
    alternative: Alternative = "two-sided",
    nan_policy: NanPolicy = "propagate",
) -> SignificanceResult[onp.ArrayND[np.float64]]: ...

#
@overload
def weightedtau(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    rank: onp.ToInt | onp.ToIntND = True,
    weigher: _Weigher | None = None,
    additive: bool = True,
    *,
    axis: None = None,
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
    axis: int | None = None,
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
    axis: int,
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
    axis: int,
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
    axis: int | None = None,
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
    axis: int | None = None,
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

_AnyFloatSub64T = TypeVar("_AnyFloatSub64T", bound=np.float32 | np.float16)

# keep in sync with `ttest_rel`
@overload  # ?d ~float64
def ttest_ind(
    a: onp.ArrayND[npc.floating64 | npc.integer | np.bool_, _JustAnyShape],
    b: onp.ArrayND[npc.floating64 | npc.integer | np.bool_, _JustAnyShape],
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[np.float64 | Any]: ...
@overload  # ?d ~T
def ttest_ind(
    a: onp.ArrayND[_AnyFloatSub64T, _JustAnyShape],
    b: onp.ArrayND[_AnyFloatSub64T, _JustAnyShape],
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[_AnyFloatSub64T | Any]: ...
@overload  # 1d ~f64
def ttest_ind(
    a: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool_],
    b: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool_],
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[np.float64]: ...
@overload  # 1d ~T
def ttest_ind(
    a: onp.ToArrayStrict1D[_AnyFloatSub64T, _AnyFloatSub64T],
    b: onp.ToArrayStrict1D[_AnyFloatSub64T, _AnyFloatSub64T],
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[_AnyFloatSub64T]: ...
@overload  # 1d +floating
def ttest_ind(
    a: onp.ToFloatStrict1D,
    b: onp.ToFloatStrict1D,
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[np.float64 | Any]: ...
@overload  # 2d ~f64
def ttest_ind(
    a: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool_],
    b: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool_],
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array1D[np.float64]]: ...
@overload  # 2d ~T
def ttest_ind(
    a: onp.ToArrayStrict2D[_AnyFloatSub64T, _AnyFloatSub64T],
    b: onp.ToArrayStrict2D[_AnyFloatSub64T, _AnyFloatSub64T],
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array1D[_AnyFloatSub64T]]: ...
@overload  # 2d +floating
def ttest_ind(
    a: onp.ToFloatStrict2D,
    b: onp.ToFloatStrict2D,
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array1D[np.float64 | Any]]: ...
@overload  # 3d ~f64
def ttest_ind(
    a: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool_],
    b: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool_],
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array2D[np.float64]]: ...
@overload  # 3d ~T
def ttest_ind(
    a: onp.ToArrayStrict3D[_AnyFloatSub64T, _AnyFloatSub64T],
    b: onp.ToArrayStrict3D[_AnyFloatSub64T, _AnyFloatSub64T],
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array2D[_AnyFloatSub64T]]: ...
@overload  # 3d +floating
def ttest_ind(
    a: onp.ToFloatStrict3D,
    b: onp.ToFloatStrict3D,
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array2D[np.float64 | Any]]: ...
@overload  # nd ~f64, axis=None
def ttest_ind(  # type: ignore[overload-overlap]
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    b: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    *,
    axis: None,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[np.float64]: ...
@overload  # nd ~f64, keepdims=True
def ttest_ind(  # type: ignore[overload-overlap]
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    b: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    *,
    axis: int | None = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[True],
) -> TtestResult[onp.ArrayND[np.float64]]: ...
@overload  # nd ~T, axis=None
def ttest_ind(
    a: onp.ToArrayND[_AnyFloatSub64T, _AnyFloatSub64T],
    b: onp.ToArrayND[_AnyFloatSub64T, _AnyFloatSub64T],
    *,
    axis: None,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[_AnyFloatSub64T]: ...
@overload  # nd ~T, keepdims=True
def ttest_ind(
    a: onp.ToArrayND[_AnyFloatSub64T, _AnyFloatSub64T],
    b: onp.ToArrayND[_AnyFloatSub64T, _AnyFloatSub64T],
    *,
    axis: int | None = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[True],
) -> TtestResult[onp.ArrayND[_AnyFloatSub64T]]: ...
@overload  # nd +floating, axis=None
def ttest_ind(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    *,
    axis: None,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[np.float64 | Any]: ...
@overload  # nd +floating, keepdims=True
def ttest_ind(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    *,
    axis: int | None = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[True],
) -> TtestResult[onp.ArrayND[np.float64 | Any]]: ...
@overload  # nd +floating
def ttest_ind(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[onp.ArrayND[np.float64 | Any] | np.float64 | Any]: ...

# keep in sync with `ttest_ind`
@overload  # ?d ~float64
def ttest_rel(
    a: onp.ArrayND[npc.floating64 | npc.integer | np.bool_, _JustAnyShape],
    b: onp.ArrayND[npc.floating64 | npc.integer | np.bool_, _JustAnyShape],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[np.float64 | Any]: ...
@overload  # ?d ~T
def ttest_rel(
    a: onp.ArrayND[_AnyFloatSub64T, _JustAnyShape],
    b: onp.ArrayND[_AnyFloatSub64T, _JustAnyShape],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[_AnyFloatSub64T | Any]: ...
@overload  # 1d ~f64
def ttest_rel(
    a: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool_],
    b: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool_],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[np.float64]: ...
@overload  # 1d ~T
def ttest_rel(
    a: onp.ToArrayStrict1D[_AnyFloatSub64T, _AnyFloatSub64T],
    b: onp.ToArrayStrict1D[_AnyFloatSub64T, _AnyFloatSub64T],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[_AnyFloatSub64T]: ...
@overload  # 1d +floating
def ttest_rel(
    a: onp.ToFloatStrict1D,
    b: onp.ToFloatStrict1D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[np.float64 | Any]: ...
@overload  # 2d ~f64
def ttest_rel(
    a: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool_],
    b: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool_],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array1D[np.float64]]: ...
@overload  # 2d ~T
def ttest_rel(
    a: onp.ToArrayStrict2D[_AnyFloatSub64T, _AnyFloatSub64T],
    b: onp.ToArrayStrict2D[_AnyFloatSub64T, _AnyFloatSub64T],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array1D[_AnyFloatSub64T]]: ...
@overload  # 2d +floating
def ttest_rel(
    a: onp.ToFloatStrict2D,
    b: onp.ToFloatStrict2D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array1D[np.float64 | Any]]: ...
@overload  # 3d ~f64
def ttest_rel(
    a: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool_],
    b: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool_],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array2D[np.float64]]: ...
@overload  # 3d ~T
def ttest_rel(
    a: onp.ToArrayStrict3D[_AnyFloatSub64T, _AnyFloatSub64T],
    b: onp.ToArrayStrict3D[_AnyFloatSub64T, _AnyFloatSub64T],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array2D[_AnyFloatSub64T]]: ...
@overload  # 3d +floating
def ttest_rel(
    a: onp.ToFloatStrict3D,
    b: onp.ToFloatStrict3D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array2D[np.float64 | Any]]: ...
@overload  # nd ~f64, axis=None
def ttest_rel(  # type: ignore[overload-overlap]
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    b: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    axis: None,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[np.float64]: ...
@overload  # nd ~f64, keepdims=True
def ttest_rel(  # type: ignore[overload-overlap]
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    b: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool_],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> TtestResult[onp.ArrayND[np.float64]]: ...
@overload  # nd ~T, axis=None
def ttest_rel(
    a: onp.ToArrayND[_AnyFloatSub64T, _AnyFloatSub64T],
    b: onp.ToArrayND[_AnyFloatSub64T, _AnyFloatSub64T],
    axis: None,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[_AnyFloatSub64T]: ...
@overload  # nd ~T, keepdims=True
def ttest_rel(
    a: onp.ToArrayND[_AnyFloatSub64T, _AnyFloatSub64T],
    b: onp.ToArrayND[_AnyFloatSub64T, _AnyFloatSub64T],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> TtestResult[onp.ArrayND[_AnyFloatSub64T]]: ...
@overload  # nd +floating, axis=None
def ttest_rel(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[np.float64 | Any]: ...
@overload  # nd +floating, keepdims=True
def ttest_rel(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> TtestResult[onp.ArrayND[np.float64 | Any]]: ...
@overload  # nd +floating
def ttest_rel(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[onp.ArrayND[np.float64 | Any] | np.float64 | Any]: ...

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
def fisher_exact(
    table: onp.ArrayND[_Real0D], alternative: Alternative | None = None, *, method: ResamplingMethod | None = None
) -> SignificanceResult[float]: ...

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

#
@overload  # axix: None (default)
def rankdata(
    a: onp.ToArrayND, method: _RankMethod = "average", *, axis: None = None, nan_policy: NanPolicy = "propagate"
) -> onp.Array1D[np.float64]: ...
@overload  # shape: T, axis: int
def rankdata(
    a: onp.Array[_ShapeT], method: _RankMethod = "average", *, axis: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[np.float64, _ShapeT]: ...
@overload  # shape: 1d, axis: int
def rankdata(
    a: Sequence[complex], method: _RankMethod = "average", *, axis: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array1D[np.float64]: ...
@overload  # shape: 2d, axis: int
def rankdata(
    a: Sequence[Sequence[complex]], method: _RankMethod = "average", *, axis: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.Array2D[np.float64]: ...
@overload  # shape: 3d, axis: int
def rankdata(
    a: Sequence[Sequence[Sequence[complex]]],
    method: _RankMethod = "average",
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array3D[np.float64]: ...
@overload  # shape: ?, axis: int
def rankdata(
    a: onp.ToArrayND, method: _RankMethod = "average", *, axis: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[np.float64]: ...

#
def expectile(a: onp.ToFloatND, alpha: float = 0.5, *, weights: onp.ToFloatND | None = None) -> np.float64: ...

#
@overload  # ?d, ?d
def linregress(
    x: onp.ArrayND[npc.floating | npc.integer | np.bool_, _JustAnyShape],
    y: onp.ArrayND[npc.floating | npc.integer | np.bool_, _JustAnyShape],
    alternative: Alternative = "two-sided",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> LinregressResult[np.float64 | Any]: ...
@overload  # 1d, 1d
def linregress(
    x: onp.ToFloatStrict1D,
    y: onp.ToFloatStrict1D,
    alternative: Alternative = "two-sided",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> LinregressResult[np.float64]: ...
@overload  # 2d, 2d
def linregress(
    x: onp.ToFloatStrict2D,
    y: onp.ToFloatStrict2D,
    alternative: Alternative = "two-sided",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> LinregressResult[onp.Array1D[np.float64]]: ...
@overload  # 3d, 3d
def linregress(
    x: onp.ToFloatStrict3D,
    y: onp.ToFloatStrict3D,
    alternative: Alternative = "two-sided",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> LinregressResult[onp.Array2D[np.float64]]: ...
@overload  # nd, nd
def linregress(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> LinregressResult[np.float64 | Any]: ...
@overload  # keepdims=True
def linregress(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    *,
    axis: int | None = 0,
    keepdims: L[True],
    nan_policy: NanPolicy = "propagate",
) -> LinregressResult[onp.ArrayND[np.float64]]: ...
@overload  # axis=None
def linregress(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    *,
    axis: None,
    keepdims: L[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> LinregressResult[np.float64]: ...

# NOTE: `lmoment` is currently numerically unstable after `order > 16`.
# See https://github.com/jorenham/Lmo/ for a more stable implementation that additionally supports generalized trimmed TL-moments,
# multivariate L- and TL-comoments, theoretical L- and TL-moments or `scipy.stats` distributions, and much more ;)

#
@overload  # ?d f64, order: 1d
def lmoment(
    sample: onp.ArrayND[npc.floating64 | npc.integer, _JustAnyShape],
    order: onp.ToInt1D | None = None,
    *,
    axis: int = 0,
    keepdims: bool = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d f64, order: 0d
def lmoment(
    sample: onp.ArrayND[npc.floating64 | npc.integer, _JustAnyShape],
    order: int,
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64] | np.float64: ...
@overload  # ?d f64, order: 0d, keepdims=True
def lmoment(
    sample: onp.ArrayND[npc.floating64 | npc.integer, _JustAnyShape],
    order: int,
    *,
    axis: int = 0,
    keepdims: L[True],
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d f32, order: 1d
def lmoment(
    sample: onp.ArrayND[np.float32 | np.float16, _JustAnyShape],
    order: onp.ToInt1D | None = None,
    *,
    axis: int = 0,
    keepdims: bool = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float32]: ...
@overload  # ?d f32, order: 0d
def lmoment(
    sample: onp.ArrayND[np.float32 | np.float16, _JustAnyShape],
    order: int,
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float32] | np.float32: ...
@overload  # ?d f32, order: 0d, keepdims=True
def lmoment(
    sample: onp.ArrayND[np.float32 | np.float16, _JustAnyShape],
    order: int,
    *,
    axis: int = 0,
    keepdims: L[True],
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float32]: ...
@overload  # 1d f64, order: 1d
def lmoment(
    sample: _AsFloat64_1D,
    order: onp.ToInt1D | None = None,
    *,
    axis: int | None = 0,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float64]: ...
@overload  # 1d f64, order: 1d, keepdims=True  # 8
def lmoment(
    sample: _AsFloat64_1D,
    order: onp.ToInt1D | None = None,
    *,
    axis: int | None = 0,
    keepdims: L[True],
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array2D[np.float64]: ...
@overload  # 1d f64, order: 0d
def lmoment(
    sample: _AsFloat64_1D,
    order: int,
    *,
    axis: int | None = 0,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> np.float64: ...
@overload  # 1d f64, order: 0d, keepdims=True
def lmoment(
    sample: _AsFloat64_1D,
    order: int,
    *,
    axis: int | None = 0,
    keepdims: L[True],
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float64]: ...
@overload  # 1d f32, order: 1d
def lmoment(
    sample: _AsFloat32_1D,
    order: onp.ToInt1D | None = None,
    *,
    axis: int | None = 0,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float32]: ...
@overload  # 1d f32, order: 1d, keepdims=True
def lmoment(
    sample: _AsFloat32_1D,
    order: onp.ToInt1D | None = None,
    *,
    axis: int | None = 0,
    keepdims: L[True],
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array2D[np.float32]: ...
@overload  # 1d f32, order: 0d
def lmoment(
    sample: _AsFloat32_1D,
    order: int,
    *,
    axis: int | None = 0,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> np.float32: ...
@overload  # 1d f32, order: 0d, keepdims=True
def lmoment(
    sample: _AsFloat32_1D,
    order: int,
    *,
    axis: int | None = 0,
    keepdims: L[True],
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float32]: ...
@overload  # 2d f64, order: 1d
def lmoment(
    sample: _AsFloat64_2D,
    order: onp.ToInt1D | None = None,
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array2D[np.float64]: ...
@overload  # 2d f64, order: 1d, keepdims=True
def lmoment(
    sample: _AsFloat64_2D,
    order: onp.ToInt1D | None = None,
    *,
    axis: int | None = 0,
    keepdims: L[True],
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array3D[np.float64]: ...
@overload  # 2d f64, order: 0d
def lmoment(
    sample: _AsFloat64_2D,
    order: int,
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float64]: ...
@overload  # 2d f64, order: 0d, keepdims=True
def lmoment(
    sample: _AsFloat64_2D,
    order: int,
    *,
    axis: int | None = 0,
    keepdims: L[True],
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array2D[np.float64]: ...
@overload  # 2d f32, order: 1d
def lmoment(
    sample: _AsFloat32_2D,
    order: onp.ToInt1D | None = None,
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array2D[np.float32]: ...
@overload  # 2d f32, order: 1d, keepdims=True
def lmoment(
    sample: _AsFloat32_2D,
    order: onp.ToInt1D | None = None,
    *,
    axis: int | None = 0,
    keepdims: L[True],
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array3D[np.float32]: ...
@overload  # 2d f32, order: 0d
def lmoment(
    sample: _AsFloat32_2D,
    order: int,
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float32]: ...
@overload  # 2d f32, order: 0d, keepdims=True
def lmoment(
    sample: _AsFloat32_2D,
    order: int,
    *,
    axis: int | None = 0,
    keepdims: L[True],
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array2D[np.float32]: ...
@overload  # nd f64, order: 1d
def lmoment(
    sample: _AsFloat64_ND,
    order: onp.ToInt1D | None = None,
    *,
    axis: int = 0,
    keepdims: bool = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64]: ...
@overload  # nd f64, order: 1d, axis=None
def lmoment(
    sample: _AsFloat64_ND,
    order: onp.ToInt1D | None = None,
    *,
    axis: None,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float64]: ...
@overload  # nd f64, order: 0d
def lmoment(
    sample: _AsFloat64_ND,
    order: int,
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64] | np.float64: ...
@overload  # nd f64, order: 0d, keepdims=True
def lmoment(
    sample: _AsFloat64_ND,
    order: int,
    *,
    axis: int | None = 0,
    keepdims: L[True],
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64]: ...
@overload  # nd f64, order: 0d, axis=None
def lmoment(
    sample: _AsFloat64_ND,
    order: int,
    *,
    axis: None,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> np.float64: ...
@overload  # nd f32, order: 1d
def lmoment(
    sample: _AsFloat32_ND,
    order: onp.ToInt1D | None = None,
    *,
    axis: int = 0,
    keepdims: bool = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float32]: ...
@overload  # nd f32, order: 1d, axis=None
def lmoment(
    sample: _AsFloat32_ND,
    order: onp.ToInt1D | None = None,
    *,
    axis: None,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float32]: ...
@overload  # nd f32, order: 0d
def lmoment(
    sample: _AsFloat32_ND,
    order: int,
    *,
    axis: int = 0,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float32] | np.float32: ...
@overload  # nd f32, order: 0d, keepdims=True
def lmoment(
    sample: _AsFloat32_ND,
    order: int,
    *,
    axis: int | None = 0,
    keepdims: L[True],
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float32]: ...
@overload  # nd f32, order: 0d, axis=None
def lmoment(
    sample: _AsFloat32_ND,
    order: int,
    *,
    axis: None,
    keepdims: L[False] = False,
    sorted: bool = False,
    standardize: bool = True,
    nan_policy: NanPolicy = "propagate",
) -> np.float32: ...
