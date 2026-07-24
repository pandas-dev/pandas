from collections.abc import Callable, Sequence
from dataclasses import dataclass
from types import ModuleType
from typing import Any, Generic, Literal as L, NamedTuple, Never, Protocol, Self, overload, override, type_check_only
from typing_extensions import TypeVar

import numpy as np
import numpy_typing_compat as nptc
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._resampling import BootstrapMethod, ResamplingMethod
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
    "sigmaclip",
    "skew",
    "skewtest",
    "spearmanr",
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

_FloatT = TypeVar("_FloatT", bound=npc.floating, default=np.float64 | Any)
_FloatT_co = TypeVar("_FloatT_co", bound=npc.floating, default=np.float64 | Any, covariant=True)
_RealT = TypeVar("_RealT", bound=_Real0D, default=_Real0D)
_RealT_co = TypeVar("_RealT_co", bound=_Real0D, default=np.float64 | Any, covariant=True)

_IntOrArrayT_co = TypeVar("_IntOrArrayT_co", bound=_ScalarOrND[np.intp], default=_ScalarOrND[np.intp], covariant=True)
_FloatOrArrayT_co = TypeVar(
    "_FloatOrArrayT_co",
    bound=float | npc.floating | onp.ArrayND[npc.floating, Any],
    default=float | onp.ArrayND[np.float64],
    covariant=True,
)
_IntFloatOrArrayT_co = TypeVar(
    "_IntFloatOrArrayT_co",
    bound=float | npc.integer | npc.floating | onp.ArrayND[npc.integer | npc.floating, Any],
    default=_FloatOrArrayT_co,
    covariant=True,
)
_SignOrArrayT_co = TypeVar(
    "_SignOrArrayT_co", bound=np.int8 | onp.ArrayND[np.int8], default=np.int8 | onp.ArrayND[np.int8], covariant=True
)
_FloatOrArrayT2_co = TypeVar(
    "_FloatOrArrayT2_co", bound=float | _ScalarOrND[npc.floating], default=float | onp.ArrayND[np.float64], covariant=True
)
_F64OrArrayT_co = TypeVar(
    "_F64OrArrayT_co", bound=np.float64 | onp.ArrayND[np.float64], default=np.float64 | onp.ArrayND[np.float64], covariant=True
)
_RealOrArrayT_co = TypeVar("_RealOrArrayT_co", bound=_ScalarOrND[_Real0D], default=_ScalarOrND[Any], covariant=True)

type _Real0D = npc.integer | npc.floating

type _ScalarOrND[SCT: np.generic] = SCT | onp.ArrayND[SCT]
type _FloatOrND = _ScalarOrND[npc.floating]

type _InterpolationMethod = L[
    "linear",
    "lower",
    "higher",
    "nearest",
    "midpoint",
    # from numpy.percentile:
    "inverted_cdf",
    "averaged_inverted_cdf",
    "closest_observation",
    "interpolated_inverted_cdf",
    "hazen",
    "weibull",
    "median_unbiased",
    "normal_unbiased",
]
type _QuantileInterpolation = L["fraction", "lower", "higher"]
type _PercentileInterpolation = L["rank", "weak", "strict", "mean"]
type _TrimTail = L["left", "right"]
type _KendallTauMethod = L["auto", "asymptotic", "exact"]
type _KendallTauVariant = L["b", "c"]
type _KS1TestMethod = L[_KS2TestMethod, "approx"]
type _KS2TestMethod = L["auto", "exact", "asymp"]
type _CombinePValuesMethod = L["fisher", "pearson", "tippett", "stouffer", "mudholkar_george"]
type _RankMethod = L["average", "min", "max", "dense", "ordinal"]

type _RealLimit = float | _Real0D
type _RealLimits = tuple[_RealLimit, _RealLimit]
type _ComplexLimit = complex | npc.number
type _ComplexLimits = tuple[_ComplexLimit, _ComplexLimit]

type _Weigher = Callable[[int], float | _Real0D]

type _JustAnyShape = tuple[Never, Never, Never, Never]  # workaround for https://github.com/microsoft/pyright/issues/10232
type _AsFloat64_1D = onp.ToArrayStrict1D[float, npc.floating64 | npc.integer]
type _AsFloat64_2D = onp.ToArrayStrict2D[float, npc.floating64 | npc.integer]
type _AsFloat64_ND = onp.ToArrayND[float, npc.floating64 | npc.integer]
type _AsFloat32_1D = onp.ToArrayStrict1D[np.float32, np.float32 | np.float16]
type _AsFloat32_2D = onp.ToArrayStrict2D[np.float32, np.float32 | np.float16]
type _AsFloat32_ND = onp.ToArrayND[Never, np.float32 | np.float16]

type _ToFloatStrictND = onp.ArrayND[npc.floating | npc.integer | np.bool, _JustAnyShape]

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
class _TestResultBunch(  # zuban: ignore[type-var]
    BaseBunch[_FloatOrArrayT_co, _FloatOrArrayT2_co],  # pyrefly: ignore[invalid-variance]
    Generic[_FloatOrArrayT_co, _FloatOrArrayT2_co],
):
    @property
    def statistic(self, /) -> _FloatOrArrayT_co: ...
    @property
    def pvalue(self, /) -> _FloatOrArrayT2_co: ...

    #
    @override
    def __new__(_cls, statistic: _FloatOrArrayT_co, pvalue: _FloatOrArrayT2_co) -> Self: ...  # pyrefly:ignore[bad-override]
    @override
    def __init__(self, /, statistic: _FloatOrArrayT_co, pvalue: _FloatOrArrayT2_co) -> None: ...  # pyrefly:ignore[bad-override]

###

# NOTE: On numpy<2.1, pyright reports 15 false positive incompatible overload errors here.
# pyright: reportOverlappingOverload=false

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
    nobs: np.int64
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
    binsize: np.float64
    extrapoints: np.int64

class RelfreqResult(NamedTuple):
    frequency: onp.Array1D[np.float64]
    lowerlimit: L[0] | npc.floating
    binsize: np.float64
    extrapoints: np.int64

class SigmaclipResult(NamedTuple, Generic[_RealT_co, _FloatT_co]):
    clipped: onp.Array1D[_RealT_co]
    lower: _FloatT_co
    upper: _FloatT_co

@dataclass
class AlexanderGovernResult:
    statistic: float
    pvalue: float

@dataclass
class QuantileTestResult(Generic[_FloatT]):
    statistic: _FloatT
    statistic_type: _FloatT  # 1 or 2
    pvalue: _FloatT
    _alternative: Alternative
    _x: onp.ArrayND[_FloatT]
    _p: onp.Array1D[np.float32]
    _statistic: _FloatT
    _statistic_type: onp.Array1D[_FloatT]
    _pvalue: onp.Array1D[_FloatT]
    _axis: int
    _axis_none: bool
    _keepdims: bool
    _ndim: int
    _nan_out: onp.Array1D[np.bool_]
    _xp: ModuleType

    def confidence_interval(self, /, confidence_level: float = 0.95) -> ConfidenceInterval[_FloatT_co]: ...

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

class TtestResultBase(_TestResultBunch[_FloatOrArrayT_co, _FloatOrArrayT_co], Generic[_FloatOrArrayT_co, _IntFloatOrArrayT_co]):
    @property
    def df(self, /) -> _IntFloatOrArrayT_co: ...

    #
    @override
    def __new__(_cls, statistic: _FloatOrArrayT_co, pvalue: _FloatOrArrayT_co, *, df: _FloatOrArrayT_co) -> Self: ...  # pyrefly:ignore[bad-override]
    @override
    def __init__(self, /, statistic: _FloatOrArrayT_co, pvalue: _FloatOrArrayT_co, *, df: _FloatOrArrayT_co) -> None: ...  # pyrefly:ignore[bad-override]

class TtestResult(TtestResultBase[_FloatOrArrayT_co, _IntFloatOrArrayT_co], Generic[_FloatOrArrayT_co, _IntFloatOrArrayT_co]):
    _alternative: Alternative
    _standard_error: _FloatOrArrayT_co
    _estimate: _FloatOrArrayT_co
    _statistic_np: _FloatOrArrayT_co
    _dtype: np.dtype[npc.floating]
    _xp: ModuleType

    @override
    def __init__(  # pyright: ignore[reportInconsistentConstructor]  # pyrefly:ignore[bad-override]
        self,
        /,
        statistic: _FloatOrArrayT_co,
        pvalue: _FloatOrArrayT_co,
        df: _IntFloatOrArrayT_co,
        alternative: Alternative,
        standard_error: _FloatOrArrayT_co,
        estimate: _FloatOrArrayT_co,
        statistic_np: _FloatOrArrayT_co | None = None,
        xp: ModuleType | None = None,
    ) -> None: ...
    def confidence_interval(self, /, confidence_level: float = 0.95) -> ConfidenceInterval[_FloatOrArrayT_co]: ...

class KstestResult(_TestResultBunch[_FloatOrArrayT_co, _FloatOrArrayT_co], Generic[_FloatOrArrayT_co, _SignOrArrayT_co]):
    @property
    def statistic_location(self, /) -> _FloatOrArrayT_co: ...
    @property
    def statistic_sign(self, /) -> _SignOrArrayT_co: ...

    #
    @override
    def __new__(  # pyrefly:ignore[bad-override]
        _cls,
        statistic: _FloatOrArrayT_co,
        pvalue: _FloatOrArrayT_co,
        *,
        statistic_location: _FloatOrArrayT_co,
        statistic_sign: _SignOrArrayT_co,
    ) -> Self: ...
    @override
    def __init__(  # pyrefly:ignore[bad-override]
        self,
        /,
        statistic: _FloatOrArrayT_co,
        pvalue: _FloatOrArrayT_co,
        *,
        statistic_location: _FloatOrArrayT_co,
        statistic_sign: _SignOrArrayT_co,
    ) -> None: ...

Ks_2sampResult = KstestResult

type _KstestResult0 = KstestResult[np.float64, np.int8]
# we can't use a generic shape-type here due to a variance bug in pyright
type _KstestResult1 = KstestResult[onp.Array1D[np.float64], onp.Array1D[np.int8]]
type _KstestResult2 = KstestResult[onp.Array2D[np.float64], onp.Array2D[np.int8]]
type _KstestResultN = KstestResult[onp.ArrayND[np.float64], onp.ArrayND[np.int8]]

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

# keep in sync with `hmean` and `pmean`
@overload  # ?d T@inexact
def gmean[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT, _JustAnyShape],
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # ?d i64|i32
def gmean(
    a: onp.ArrayND[npc.integer64 | npc.integer32, _JustAnyShape],
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # 1d T@inexact
def gmean[InexactT: npc.inexact](
    a: onp.Array1D[InexactT],
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloat1D | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # 1d float|i64|i32
def gmean(
    a: onp.ToArrayStrict1D[float, npc.integer64 | npc.integer32],
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloat1D | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # 1d ~complex
def gmean(
    a: list[complex],
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloat1D | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.complex128: ...
@overload  # 2d T@inexact
def gmean[InexactT: npc.inexact](
    a: onp.Array2D[InexactT],
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloat1D | onp.ToFloat2D | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[InexactT]: ...
@overload  # 2d float|i64|i32
def gmean(
    a: onp.ToArrayStrict2D[float, npc.integer64 | npc.integer32],
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloat1D | onp.ToFloat2D | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~complex
def gmean(
    a: Sequence[list[complex]],
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloat1D | onp.ToFloat2D | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[np.complex128]: ...
@overload  # Nd T@inexact
def gmean[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT],
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # Nd T@inexact, keepdims=True
def gmean[InexactT: npc.inexact, ShapeT: tuple[int, ...]](
    a: onp.ArrayND[InexactT, ShapeT],
    axis: int | None = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[InexactT, ShapeT]: ...
@overload  # Nd T@inexact, axis=None
def gmean[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT],
    axis: None,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # Nd float
def gmean(
    a: onp.SequenceND[Sequence[float]],
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64]: ...
@overload  # Nd float, keepdims=True
def gmean(
    a: onp.SequenceND[float],
    axis: int | None = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # Nd ~complex
def gmean(
    a: onp.SequenceND[list[complex]],
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.ArrayND[np.complex128]: ...
@overload  # Nd ~complex, keepdims=True
def gmean(
    a: onp.SequenceND[list[complex]] | list[complex],
    axis: int | None = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.complex128]: ...
@overload  # Nd ~complex, axis=None
def gmean(
    a: onp.SequenceND[list[complex]] | list[complex],
    axis: None,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.complex128: ...
@overload  # Nd i64|i32
def gmean(
    a: onp.ArrayND[npc.integer64 | npc.integer32],
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # Nd i64|i32, keepdims=True
def gmean[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.integer64 | npc.integer32, ShapeT],
    axis: int | None = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # Nd float|i64|i32, axis=None
def gmean(
    a: onp.ToArrayND[float, npc.integer64 | npc.integer32],
    axis: None,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # ?d, dtype=<known>
def gmean[InexactT: npc.inexact](
    a: onp.ArrayND[npc.number | np.bool, _JustAnyShape],
    axis: int = 0,
    *,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # 1d, dtype=<known>
def gmean[InexactT: npc.inexact](
    a: onp.ToComplexStrict1D,
    axis: int = 0,
    *,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # 2d, dtype=<known>
def gmean[InexactT: npc.inexact](
    a: onp.ToComplexStrict2D,
    axis: int = 0,
    *,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[InexactT]: ...
@overload  # Nd, dtype=<known>
def gmean[InexactT: npc.inexact](
    a: onp.ToComplexND,
    axis: int = 0,
    *,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # Nd, dtype=<known>, keepdims=True
def gmean[InexactT: npc.inexact](
    a: onp.ToComplexND,
    axis: int | None = 0,
    *,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[InexactT]: ...
@overload  # Nd, dtype=<known>, axis=None
def gmean[InexactT: npc.inexact](
    a: onp.ToComplexND,
    axis: None,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # dtype=? (fallback)
def gmean(
    a: onp.ToComplexND,
    axis: int | None = 0,
    *,
    dtype: str | type,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> Any: ...

# keep in sync with `gmean` and `pmean`
@overload  # ?d T@inexact
def hmean[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT, _JustAnyShape],
    axis: int = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # ?d i64|i32
def hmean(
    a: onp.ArrayND[npc.integer64 | npc.integer32, _JustAnyShape],
    axis: int = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # 1d T@inexact
def hmean[InexactT: npc.inexact](
    a: onp.Array1D[InexactT],
    axis: int = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloat1D | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # 1d float|i64|i32
def hmean(
    a: onp.ToArrayStrict1D[float, npc.integer64 | npc.integer32],
    axis: int = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloat1D | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # 1d ~complex
def hmean(
    a: list[complex],
    axis: int = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloat1D | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.complex128: ...
@overload  # 2d T@inexact
def hmean[InexactT: npc.inexact](
    a: onp.Array2D[InexactT],
    axis: int = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloat1D | onp.ToFloat2D | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[InexactT]: ...
@overload  # 2d float|i64|i32
def hmean(
    a: onp.ToArrayStrict2D[float, npc.integer64 | npc.integer32],
    axis: int = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloat1D | onp.ToFloat2D | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~complex
def hmean(
    a: Sequence[list[complex]],
    axis: int = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloat1D | onp.ToFloat2D | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[np.complex128]: ...
@overload  # Nd T@inexact
def hmean[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT],
    axis: int = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # Nd T@inexact, keepdims=True
def hmean[InexactT: npc.inexact, ShapeT: tuple[int, ...]](
    a: onp.ArrayND[InexactT, ShapeT],
    axis: int | None = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[InexactT, ShapeT]: ...
@overload  # Nd T@inexact, axis=None
def hmean[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT],
    axis: None,
    dtype: None = None,
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # Nd float
def hmean(
    a: onp.SequenceND[Sequence[float]],
    axis: int = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64]: ...
@overload  # Nd float, keepdims=True
def hmean(
    a: onp.SequenceND[float],
    axis: int | None = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # Nd ~complex
def hmean(
    a: onp.SequenceND[list[complex]],
    axis: int = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.ArrayND[np.complex128]: ...
@overload  # Nd ~complex, keepdims=True
def hmean(
    a: onp.SequenceND[list[complex]] | list[complex],
    axis: int | None = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.complex128]: ...
@overload  # Nd ~complex, axis=None
def hmean(
    a: onp.SequenceND[list[complex]] | list[complex],
    axis: None,
    dtype: None = None,
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.complex128: ...
@overload  # Nd i64|i32
def hmean(
    a: onp.ArrayND[npc.integer64 | npc.integer32],
    axis: int = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # Nd i64|i32, keepdims=True
def hmean[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.integer64 | npc.integer32, ShapeT],
    axis: int | None = 0,
    dtype: None = None,
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # Nd float|i64|i32, axis=None
def hmean(
    a: onp.ToArrayND[float, npc.integer64 | npc.integer32],
    axis: None,
    dtype: None = None,
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # ?d, dtype=<known>
def hmean[InexactT: npc.inexact](
    a: onp.ArrayND[npc.number | np.bool, _JustAnyShape],
    axis: int = 0,
    *,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # 1d, dtype=<known>
def hmean[InexactT: npc.inexact](
    a: onp.ToComplexStrict1D,
    axis: int = 0,
    *,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # 2d, dtype=<known>
def hmean[InexactT: npc.inexact](
    a: onp.ToComplexStrict2D,
    axis: int = 0,
    *,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[InexactT]: ...
@overload  # Nd, dtype=<known>
def hmean[InexactT: npc.inexact](
    a: onp.ToComplexND,
    axis: int = 0,
    *,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # Nd, dtype=<known>, keepdims=True
def hmean[InexactT: npc.inexact](
    a: onp.ToComplexND,
    axis: int | None = 0,
    *,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[InexactT]: ...
@overload  # Nd, dtype=<known>, axis=None
def hmean[InexactT: npc.inexact](
    a: onp.ToComplexND,
    axis: None,
    dtype: onp.ToDType[InexactT],
    *,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # dtype=? (fallback)
def hmean(
    a: onp.ToComplexND,
    axis: int | None = 0,
    *,
    dtype: str | type,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> Any: ...

# keep in sync with `gmean` and `hmean`
@overload  # ?d T@inexact
def pmean[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT, _JustAnyShape],
    p: float,
    *,
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # ?d i64|i32
def pmean(
    a: onp.ArrayND[npc.integer64 | npc.integer32, _JustAnyShape],
    p: float,
    *,
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # 1d T@inexact
def pmean[InexactT: npc.inexact](
    a: onp.Array1D[InexactT],
    p: float,
    *,
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloat1D | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # 1d float|i64|i32
def pmean(
    a: onp.ToArrayStrict1D[float, npc.integer64 | npc.integer32],
    p: float,
    *,
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloat1D | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # 1d ~complex
def pmean(
    a: list[complex],
    p: float,
    *,
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloat1D | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.complex128: ...
@overload  # 2d T@inexact
def pmean[InexactT: npc.inexact](
    a: onp.Array2D[InexactT],
    p: float,
    *,
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloat1D | onp.ToFloat2D | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[InexactT]: ...
@overload  # 2d float|i64|i32
def pmean(
    a: onp.ToArrayStrict2D[float, npc.integer64 | npc.integer32],
    p: float,
    *,
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloat1D | onp.ToFloat2D | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~complex
def pmean(
    a: Sequence[list[complex]],
    p: float,
    *,
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloat1D | onp.ToFloat2D | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[np.complex128]: ...
@overload  # Nd T@inexact
def pmean[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT],
    p: float,
    *,
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # Nd T@inexact, keepdims=True
def pmean[InexactT: npc.inexact, ShapeT: tuple[int, ...]](
    a: onp.ArrayND[InexactT, ShapeT],
    p: float,
    *,
    axis: int | None = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[InexactT, ShapeT]: ...
@overload  # Nd T@inexact, axis=None
def pmean[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT],
    p: float,
    *,
    axis: None,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # Nd float
def pmean(
    a: onp.SequenceND[Sequence[float]],
    p: float,
    *,
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64]: ...
@overload  # Nd float, keepdims=True
def pmean(
    a: onp.SequenceND[float],
    p: float,
    *,
    axis: int | None = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # Nd ~complex
def pmean(
    a: onp.SequenceND[list[complex]],
    p: float,
    *,
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.ArrayND[np.complex128]: ...
@overload  # Nd ~complex, keepdims=True
def pmean(
    a: onp.SequenceND[list[complex]] | list[complex],
    p: float,
    *,
    axis: int | None = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.complex128]: ...
@overload  # Nd ~complex, axis=None
def pmean(
    a: onp.SequenceND[list[complex]] | list[complex],
    p: float,
    *,
    axis: None,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.complex128: ...
@overload  # Nd i64|i32
def pmean(
    a: onp.ArrayND[npc.integer64 | npc.integer32],
    p: float,
    *,
    axis: int = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # Nd i64|i32, keepdims=True
def pmean[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.integer64 | npc.integer32, ShapeT],
    p: float,
    *,
    axis: int | None = 0,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # Nd float|i64|i32, axis=None
def pmean(
    a: onp.ToArrayND[float, npc.integer64 | npc.integer32],
    p: float,
    *,
    axis: None,
    dtype: None = None,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # ?d, dtype=<known>
def pmean[InexactT: npc.inexact](
    a: onp.ArrayND[npc.number | np.bool, _JustAnyShape],
    p: float,
    *,
    axis: int = 0,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # 1d, dtype=<known>
def pmean[InexactT: npc.inexact](
    a: onp.ToComplexStrict1D,
    p: float,
    *,
    axis: int = 0,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # 2d, dtype=<known>
def pmean[InexactT: npc.inexact](
    a: onp.ToComplexStrict2D,
    p: float,
    *,
    axis: int = 0,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[InexactT]: ...
@overload  # Nd, dtype=<known>
def pmean[InexactT: npc.inexact](
    a: onp.ToComplexND,
    p: float,
    *,
    axis: int = 0,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # Nd, dtype=<known>, keepdims=True
def pmean[InexactT: npc.inexact](
    a: onp.ToComplexND,
    p: float,
    *,
    axis: int | None = 0,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[InexactT]: ...
@overload  # Nd, dtype=<known>, axis=None
def pmean[InexactT: npc.inexact](
    a: onp.ToComplexND,
    p: float,
    *,
    axis: None,
    dtype: onp.ToDType[InexactT],
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # dtype=? (fallback)
def pmean(
    a: onp.ToComplexND,
    p: float,
    *,
    axis: int | None = 0,
    dtype: str | type,
    weights: onp.ToFloatND | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> Any: ...

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

#
@overload  # ?d T@inexact, axis=None  (default)
def tmean[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # ?d +f64, axis=None  (default)
def tmean(
    a: onp.ToArrayND[float, npc.integer | np.bool],
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # ?d ~complex, axis=None  (default)
def tmean(
    a: onp.SequenceND[list[complex]] | list[complex],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.complex128: ...
@overload  # ?d T@inexact, axis=<given>
def tmean[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT, _JustAnyShape],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    *,
    axis: int,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # ?d +integer, axis=<given>
def tmean(
    a: onp.ArrayND[npc.integer | np.bool, _JustAnyShape],
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    *,
    axis: int,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # 1d T@inexact, axis=<given>
def tmean[InexactT: npc.inexact](
    a: onp.ToArrayStrict1D[InexactT, InexactT],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    *,
    axis: int,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # 1d +float|integer, axis=<given>
def tmean(
    a: onp.ToArrayStrict1D[float, npc.integer | np.bool],
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    *,
    axis: int,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # 1d ~complex, axis=<given>
def tmean(
    a: list[complex],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    *,
    axis: int,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.complex128: ...
@overload  # 2d T@inexact, axis=<given>
def tmean[InexactT: npc.inexact](
    a: onp.ToArrayStrict2D[InexactT, InexactT],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    *,
    axis: int,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[InexactT]: ...
@overload  # 2d +float|integer, axis=<given>
def tmean(
    a: onp.ToArrayStrict2D[float, npc.integer | np.bool],
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    *,
    axis: int,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~complex, axis=<given>
def tmean(
    a: Sequence[list[complex]],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    *,
    axis: int,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[np.complex128]: ...
@overload  # S@Nd T@inexact, keepdims=True
def tmean[InexactT: npc.inexact, ShapeT: tuple[int, ...]](
    a: onp.ArrayND[InexactT, ShapeT],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[InexactT, ShapeT]: ...
@overload  # S@Nd +integer, keepdims=True
def tmean[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.integer | np.bool, ShapeT],
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # ?d +float, keepdims=True
def tmean(
    a: onp.SequenceND[float],
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d ~complex, keepdims=True
def tmean(
    a: onp.SequenceND[list[complex]] | list[complex],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.complex128]: ...

#
@overload  # ?d T@inexact
def tvar[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT, _JustAnyShape],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # ?d +integer
def tvar(
    a: onp.ArrayND[npc.integer | np.bool, _JustAnyShape],
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # 1d T@inexact
def tvar[InexactT: npc.inexact](
    a: onp.ToArrayStrict1D[InexactT, InexactT],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # 1d +float|integer
def tvar(
    a: onp.ToArrayStrict1D[float, npc.integer | np.bool],
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # 1d ~complex
def tvar(
    a: list[complex],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.complex128: ...
@overload  # 2d T@inexact
def tvar[InexactT: npc.inexact](
    a: onp.ToArrayStrict2D[InexactT, InexactT],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[InexactT]: ...
@overload  # 2d +float|integer
def tvar(
    a: onp.ToArrayStrict2D[float, npc.integer | np.bool],
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~complex
def tvar(
    a: Sequence[list[complex]],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.Array1D[np.complex128]: ...
@overload  # ?d T@inexact, axis=None
def tvar[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    *,
    axis: None,
    ddof: int = 1,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # ?d +f64, axis=None
def tvar(
    a: onp.ToArrayND[float, npc.integer | np.bool],
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    *,
    axis: None,
    ddof: int = 1,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # ?d ~complex, axis=None
def tvar(
    a: onp.SequenceND[list[complex]] | list[complex],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    *,
    axis: None,
    ddof: int = 1,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.complex128: ...
@overload  # S@Nd T@inexact, keepdims=True
def tvar[InexactT: npc.inexact, ShapeT: tuple[int, ...]](
    a: onp.ArrayND[InexactT, ShapeT],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[InexactT, ShapeT]: ...
@overload  # S@Nd +integer, keepdims=True
def tvar[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.integer | np.bool, ShapeT],
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # ?d +float, keepdims=True
def tvar(
    a: onp.SequenceND[float],
    limits: _RealLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d ~complex, keepdims=True
def tvar(
    a: onp.SequenceND[list[complex]] | list[complex],
    limits: _ComplexLimits | None = None,
    inclusive: tuple[bool, bool] = (True, True),
    axis: int = 0,
    ddof: int = 1,
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.complex128]: ...

# NOTE: These have actually different implementations, but the same signature, so we're just being lazy (and slightly incorrect)
tstd = tvar
tsem = tvar

# keep in sync with `tmax`, and structurally with `tvar` and `tmean`
@overload  # ?d T@inexact
def tmin[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT, _JustAnyShape],
    lowerlimit: _ComplexLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # ?d +integer
def tmin(
    a: onp.ArrayND[npc.integer | np.bool, _JustAnyShape],
    lowerlimit: _RealLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # 1d T@inexact
def tmin[InexactT: npc.inexact](
    a: onp.ToArrayStrict1D[InexactT, InexactT],
    lowerlimit: _ComplexLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # 1d +float|integer
def tmin(
    a: onp.ToArrayStrict1D[float, npc.integer | np.bool],
    lowerlimit: _RealLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # 1d ~complex
def tmin(
    a: list[complex],
    lowerlimit: _ComplexLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.complex128: ...
@overload  # 2d T@inexact
def tmin[InexactT: npc.inexact](
    a: onp.ToArrayStrict2D[InexactT, InexactT],
    lowerlimit: _ComplexLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array1D[InexactT]: ...
@overload  # 2d +float|integer
def tmin(
    a: onp.ToArrayStrict2D[float, npc.integer | np.bool],
    lowerlimit: _RealLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~complex
def tmin(
    a: Sequence[list[complex]],
    lowerlimit: _ComplexLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array1D[np.complex128]: ...
@overload  # ?d T@inexact, axis=None
def tmin[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT],
    lowerlimit: _ComplexLimit | None = None,
    *,
    axis: None,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # ?d +f64, axis=None
def tmin(
    a: onp.ToArrayND[float, npc.integer | np.bool],
    lowerlimit: _RealLimit | None = None,
    *,
    axis: None,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # ?d ~complex, axis=None
def tmin(
    a: onp.SequenceND[list[complex]] | list[complex],
    lowerlimit: _ComplexLimit | None = None,
    *,
    axis: None,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.complex128: ...
@overload  # S@Nd T@inexact, keepdims=True
def tmin[InexactT: npc.inexact, ShapeT: tuple[int, ...]](
    a: onp.ArrayND[InexactT, ShapeT],
    lowerlimit: _ComplexLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[InexactT, ShapeT]: ...
@overload  # S@Nd +integer, keepdims=True
def tmin[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.integer | np.bool, ShapeT],
    lowerlimit: _RealLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # ?d +float, keepdims=True
def tmin(
    a: onp.SequenceND[float],
    lowerlimit: _RealLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d ~complex, keepdims=True
def tmin(
    a: onp.SequenceND[list[complex]] | list[complex],
    lowerlimit: _ComplexLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[np.complex128]: ...

# keep in sync with `tmin`, and structurally with `tvar` and `tmean`
@overload  # ?d T@inexact
def tmax[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT, _JustAnyShape],
    upperlimit: _ComplexLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> InexactT | onp.ArrayND[InexactT]: ...
@overload  # ?d +integer
def tmax(
    a: onp.ArrayND[npc.integer | np.bool, _JustAnyShape],
    upperlimit: _RealLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # 1d T@inexact
def tmax[InexactT: npc.inexact](
    a: onp.ToArrayStrict1D[InexactT, InexactT],
    upperlimit: _ComplexLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # 1d +float|integer
def tmax(
    a: onp.ToArrayStrict1D[float, npc.integer | np.bool],
    upperlimit: _RealLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # 1d ~complex
def tmax(
    a: list[complex],
    upperlimit: _ComplexLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.complex128: ...
@overload  # 2d T@inexact
def tmax[InexactT: npc.inexact](
    a: onp.ToArrayStrict2D[InexactT, InexactT],
    upperlimit: _ComplexLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array1D[InexactT]: ...
@overload  # 2d +float|integer
def tmax(
    a: onp.ToArrayStrict2D[float, npc.integer | np.bool],
    upperlimit: _RealLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~complex
def tmax(
    a: Sequence[list[complex]],
    upperlimit: _ComplexLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array1D[np.complex128]: ...
@overload  # ?d T@inexact, axis=None
def tmax[InexactT: npc.inexact](
    a: onp.ArrayND[InexactT],
    upperlimit: _ComplexLimit | None = None,
    *,
    axis: None,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> InexactT: ...
@overload  # ?d +f64, axis=None
def tmax(
    a: onp.ToArrayND[float, npc.integer | np.bool],
    upperlimit: _RealLimit | None = None,
    *,
    axis: None,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # ?d ~complex, axis=None
def tmax(
    a: onp.SequenceND[list[complex]] | list[complex],
    upperlimit: _ComplexLimit | None = None,
    *,
    axis: None,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.complex128: ...
@overload  # S@Nd T@inexact, keepdims=True
def tmax[InexactT: npc.inexact, ShapeT: tuple[int, ...]](
    a: onp.ArrayND[InexactT, ShapeT],
    upperlimit: _ComplexLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[InexactT, ShapeT]: ...
@overload  # S@Nd +integer, keepdims=True
def tmax[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.integer | np.bool, ShapeT],
    upperlimit: _RealLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # ?d +float, keepdims=True
def tmax(
    a: onp.SequenceND[float],
    upperlimit: _RealLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d ~complex, keepdims=True
def tmax(
    a: onp.SequenceND[list[complex]] | list[complex],
    upperlimit: _ComplexLimit | None = None,
    axis: int = 0,
    inclusive: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[np.complex128]: ...

#
@overload
def gstd(
    a: _ToFloatStrictND, axis: int = 0, ddof: int = 1, *, keepdims: L[False] = False, nan_policy: NanPolicy = "propagate"
) -> np.float64 | onp.ArrayND[np.float64]: ...
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
    a: onp.ArrayND[npc.floating64 | npc.integer | np.bool, _JustAnyShape],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # ?d ~T: floating, order: 0d
def moment[FloatT: npc.floating](
    a: onp.ArrayND[FloatT, _JustAnyShape],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> FloatT | onp.ArrayND[FloatT]: ...
@overload  # 1d ~f64, order: 0d
def moment(
    a: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # 1d ~T: floating, order: 0d
def moment[FloatT: npc.floating](
    a: onp.ToArrayStrict1D[FloatT, FloatT],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> FloatT: ...
@overload  # 2d ~f64, order: 0d
def moment(
    a: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~T: floating, order: 0d
def moment[FloatT: npc.floating](
    a: onp.ToArrayStrict2D[FloatT, FloatT],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.Array1D[FloatT]: ...
@overload  # 3d ~f64, order: 0d
def moment(
    a: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.Array2D[np.float64]: ...
@overload  # 3d ~T: floating, order: 0d
def moment[FloatT: npc.floating](
    a: onp.ToArrayStrict3D[FloatT, FloatT],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.Array2D[FloatT]: ...
@overload  # nd ~f64, order: 0d
def moment(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64] | Any: ...
@overload  # nd ~f64, order: 0d, axis=None  (positional)
def moment(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    order: int,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # nd ~f64, order: 0d, axis=None  (keyword)
def moment(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    order: int = 1,
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    center: float | None = None,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # nd ~f64, order: nd
def moment(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    order: onp.ToIntND,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64]: ...
@overload  # nd ~f64, keepdims=True
def moment(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    order: int | onp.ToIntND = 1,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # nd ~T: floating, order: 0d
def moment[FloatT: npc.floating](
    a: onp.ToArrayND[FloatT, FloatT],
    order: int = 1,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.ArrayND[FloatT] | Any: ...
@overload  # nd ~T: floating, order: 0d, axis=None  (positional)
def moment[FloatT: npc.floating](
    a: onp.ToArrayND[FloatT, FloatT],
    order: int,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> FloatT: ...
@overload  # nd ~T: floating, order: 0d, axis=None  (keyword)
def moment[FloatT: npc.floating](
    a: onp.ToArrayND[FloatT, FloatT],
    order: int = 1,
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    center: float | None = None,
    keepdims: L[False] = False,
) -> FloatT: ...
@overload  # nd ~T: floating, order: nd
def moment[FloatT: npc.floating](
    a: onp.ToArrayND[FloatT, FloatT],
    order: onp.ToIntND,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[False] = False,
) -> onp.ArrayND[FloatT]: ...
@overload  # nd ~T: floating, keepdims=True
def moment[FloatT: npc.floating](
    a: onp.ToArrayND[FloatT, FloatT],
    order: int | onp.ToIntND = 1,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    center: float | None = None,
    keepdims: L[True],
) -> onp.ArrayND[FloatT]: ...
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
    a: onp.ArrayND[npc.floating64 | npc.integer | np.bool, _JustAnyShape],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # ?d ~T
def skew[FloatT: npc.floating](
    a: onp.ArrayND[FloatT, _JustAnyShape],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> FloatT | onp.ArrayND[FloatT]: ...
@overload  # 1d ~f64
def skew(
    a: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # 1d ~T
def skew[FloatT: npc.floating](
    a: onp.ToArrayStrict1D[FloatT, FloatT],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> FloatT: ...
@overload  # 2d ~f64
def skew(
    a: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~T
def skew[FloatT: npc.floating](
    a: onp.ToArrayStrict2D[FloatT, FloatT],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array1D[FloatT]: ...
@overload  # 3d ~f64
def skew(
    a: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array2D[np.float64]: ...
@overload  # 3d ~T
def skew[FloatT: npc.floating](
    a: onp.ToArrayStrict3D[FloatT, FloatT],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array2D[FloatT]: ...
@overload  # nd ~f64
def skew(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64] | Any: ...
@overload  # nd ~f64, axis=None
def skew(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    axis: None,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # nd ~f64, keepdims=True
def skew(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    axis: int | None = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # nd ~T
def skew[FloatT: npc.floating](
    a: onp.ToArrayND[FloatT, FloatT],
    axis: int = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.ArrayND[FloatT] | Any: ...
@overload  # nd ~T, axis=None
def skew[FloatT: npc.floating](
    a: onp.ToArrayND[FloatT, FloatT],
    axis: None,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> FloatT: ...
@overload  # nd ~T, keepdims=True
def skew[FloatT: npc.floating](
    a: onp.ToArrayND[FloatT, FloatT],
    axis: int | None = 0,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[FloatT]: ...
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
    a: onp.ArrayND[npc.floating64 | npc.integer | np.bool, _JustAnyShape],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # ?d ~T
def kurtosis[FloatT: npc.floating](
    a: onp.ArrayND[FloatT, _JustAnyShape],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> FloatT | onp.ArrayND[FloatT]: ...
@overload  # 1d ~f64
def kurtosis(
    a: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # 1d ~T
def kurtosis[FloatT: npc.floating](
    a: onp.ToArrayStrict1D[FloatT, FloatT],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> FloatT: ...
@overload  # 2d ~f64
def kurtosis(
    a: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~T
def kurtosis[FloatT: npc.floating](
    a: onp.ToArrayStrict2D[FloatT, FloatT],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array1D[FloatT]: ...
@overload  # 3d ~f64
def kurtosis(
    a: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array2D[np.float64]: ...
@overload  # 3d ~T
def kurtosis[FloatT: npc.floating](
    a: onp.ToArrayStrict3D[FloatT, FloatT],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.Array2D[FloatT]: ...
@overload  # nd ~f64
def kurtosis(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64] | Any: ...
@overload  # nd ~f64, axis=None
def kurtosis(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    axis: None,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # nd ~f64, keepdims=True
def kurtosis(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    axis: int | None = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # nd ~T
def kurtosis[FloatT: npc.floating](
    a: onp.ToArrayND[FloatT, FloatT],
    axis: int = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> onp.ArrayND[FloatT] | Any: ...
@overload  # nd ~T, axis=None
def kurtosis[FloatT: npc.floating](
    a: onp.ToArrayND[FloatT, FloatT],
    axis: None,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> FloatT: ...
@overload  # nd ~T, keepdims=True
def kurtosis[FloatT: npc.floating](
    a: onp.ToArrayND[FloatT, FloatT],
    axis: int | None = 0,
    fisher: bool = True,
    bias: bool = True,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> onp.ArrayND[FloatT]: ...
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
@overload  # ?d T@integer, axis=None
def describe[IntT: npc.integer](
    a: onp.ArrayND[IntT], axis: None, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[IntT, np.float64]: ...
@overload  # ?d T@floating, axis=None
def describe[FloatT: npc.floating](
    a: onp.ArrayND[FloatT], axis: None, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[FloatT, FloatT]: ...
@overload  # ?d T@integer
def describe[IntT: npc.integer](
    a: onp.ArrayND[IntT, _JustAnyShape], axis: int = 0, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[IntT | onp.ArrayND[IntT], np.float64 | onp.ArrayND[np.float64]]: ...
@overload  # ?d T@floating
def describe[FloatT: npc.floating](
    a: onp.ArrayND[FloatT, _JustAnyShape], axis: int = 0, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[FloatT | onp.ArrayND[FloatT], FloatT | onp.ArrayND[FloatT]]: ...
@overload  # 1d int
def describe(
    a: Sequence[int], axis: int | None = 0, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[np.int_, np.float64]: ...
@overload  # 1d float
def describe(
    a: list[float], axis: int | None = 0, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[np.float64, np.float64]: ...
@overload  # 1d T@integer
def describe[IntT: npc.integer](
    a: onp.Array1D[IntT], axis: int | None = 0, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[IntT, np.float64]: ...
@overload  # 1d T@floating
def describe[FloatT: npc.floating](
    a: onp.Array1D[FloatT], axis: int | None = 0, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[FloatT, FloatT]: ...
@overload  # 2d int
def describe(
    a: Sequence[Sequence[int]], axis: int = 0, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[onp.Array1D[np.int_], onp.Array1D[np.float64]]: ...
@overload  # 2d int, axis=None
def describe(
    a: Sequence[Sequence[int]], axis: None, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[np.int_, np.float64]: ...
@overload  # 2d float
def describe(
    a: Sequence[list[float]], axis: int = 0, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[onp.Array1D[np.float64], onp.Array1D[np.float64]]: ...
@overload  # 2d float, axis=None
def describe(
    a: Sequence[list[float]], axis: None, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[np.float64, np.float64]: ...
@overload  # 2d T@integer
def describe[IntT: npc.integer](
    a: onp.Array2D[IntT], axis: int = 0, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[onp.Array1D[IntT], onp.Array1D[np.float64]]: ...
@overload  # 2d T@floating
def describe[FloatT: npc.floating](
    a: onp.Array2D[FloatT], axis: int = 0, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[onp.Array1D[FloatT], onp.Array1D[FloatT]]: ...
@overload  # fallback
def describe(
    a: onp.ToFloatND, axis: int | None = 0, ddof: int = 1, bias: bool = True, nan_policy: NanPolicy = "propagate"
) -> DescribeResult[Any, Any]: ...

# keep in sync with `kurtosistest` and `normaltest`
@overload  # ?d ~f64, axis=None
def skewtest(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer],
    axis: None,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> SkewtestResult[np.float64]: ...
@overload  # ?d ~f64, axis=<given>  (default)
def skewtest(
    a: onp.ArrayND[npc.floating64 | npc.integer, _JustAnyShape],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> SkewtestResult[onp.ArrayND[np.float64] | np.float64]: ...
@overload  # 1d ~f64, axis=<given>  (default)
def skewtest(
    a: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> SkewtestResult[np.float64]: ...
@overload  # 2d ~f64, axis=<given>  (default)
def skewtest(
    a: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> SkewtestResult[onp.Array1D[np.float64]]: ...
@overload  # Nd ~f64, axis=<given>  (default)
def skewtest(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> SkewtestResult[onp.ArrayND[np.float64] | np.float64]: ...
@overload  # Nd ~f64, keepdims=True
def skewtest[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.floating64 | npc.integer, ShapeT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> SkewtestResult[onp.ArrayND[np.float64, ShapeT]]: ...
@overload  # ?d ~f64, keepdims=True
def skewtest(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> SkewtestResult[onp.ArrayND[np.float64]]: ...
@overload  # ?d ~f32, axis=None
def skewtest(
    a: onp.ToJustFloat32_ND,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> SkewtestResult[np.float32]: ...
@overload  # ?d ~f32, axis=<given>  (default)
def skewtest(
    a: onp.ArrayND[np.float32, _JustAnyShape],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> SkewtestResult[onp.ArrayND[np.float32] | np.float32]: ...
@overload  # 1d ~f32, axis=<given>  (default)
def skewtest(
    a: onp.ToJustFloat32Strict1D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> SkewtestResult[np.float32]: ...
@overload  # 2d ~f32, axis=<given>  (default)
def skewtest(
    a: onp.ToJustFloat32Strict2D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> SkewtestResult[onp.Array1D[np.float32]]: ...
@overload  # Nd ~f32, axis=<given>  (default)
def skewtest(
    a: onp.ToJustFloat32_ND,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> SkewtestResult[onp.ArrayND[np.float32] | np.float32]: ...
@overload  # Nd ~f32, keepdims=True
def skewtest[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[np.float32, ShapeT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> SkewtestResult[onp.ArrayND[np.float32, ShapeT]]: ...
@overload  # ?d ~f32, keepdims=True
def skewtest(
    a: onp.ToJustFloat32_ND,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> SkewtestResult[onp.ArrayND[np.float32]]: ...
@overload  # ?d floating, axis=None
def skewtest(
    a: onp.ToArrayND[npc.floating, npc.floating],
    axis: None,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> SkewtestResult[np.float64 | Any]: ...
@overload  # Nd floating, axis=<given>  (default)
def skewtest(
    a: onp.ToArrayND[npc.floating, npc.floating],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> SkewtestResult[onp.ArrayND[np.float64 | Any] | Any]: ...
@overload  # Nd floating, keepdims=True
def skewtest[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.floating, ShapeT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> SkewtestResult[onp.ArrayND[np.float64 | Any, ShapeT]]: ...
@overload  # ?d floating, keepdims=True
def skewtest(
    a: onp.ToArrayND[npc.floating, npc.floating],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> SkewtestResult[onp.ArrayND[np.float64 | Any]]: ...

# keep in sync with `skewtest` and `normaltest`
@overload  # ?d ~f64, axis=None
def kurtosistest(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer],
    axis: None,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> KurtosistestResult[np.float64]: ...
@overload  # ?d ~f64, axis=<given>  (default)
def kurtosistest(
    a: onp.ArrayND[npc.floating64 | npc.integer, _JustAnyShape],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> KurtosistestResult[onp.ArrayND[np.float64] | np.float64]: ...
@overload  # 1d ~f64, axis=<given>  (default)
def kurtosistest(
    a: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> KurtosistestResult[np.float64]: ...
@overload  # 2d ~f64, axis=<given>  (default)
def kurtosistest(
    a: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> KurtosistestResult[onp.Array1D[np.float64]]: ...
@overload  # Nd ~f64, axis=<given>  (default)
def kurtosistest(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> KurtosistestResult[onp.ArrayND[np.float64] | np.float64]: ...
@overload  # Nd ~f64, keepdims=True
def kurtosistest[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.floating64 | npc.integer, ShapeT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> KurtosistestResult[onp.ArrayND[np.float64, ShapeT]]: ...
@overload  # ?d ~f64, keepdims=True
def kurtosistest(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> KurtosistestResult[onp.ArrayND[np.float64]]: ...
@overload  # ?d ~f32, axis=None
def kurtosistest(
    a: onp.ToJustFloat32_ND,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> KurtosistestResult[np.float32]: ...
@overload  # ?d ~f32, axis=<given>  (default)
def kurtosistest(
    a: onp.ArrayND[np.float32, _JustAnyShape],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> KurtosistestResult[onp.ArrayND[np.float32] | np.float32]: ...
@overload  # 1d ~f32, axis=<given>  (default)
def kurtosistest(
    a: onp.ToJustFloat32Strict1D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> KurtosistestResult[np.float32]: ...
@overload  # 2d ~f32, axis=<given>  (default)
def kurtosistest(
    a: onp.ToJustFloat32Strict2D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> KurtosistestResult[onp.Array1D[np.float32]]: ...
@overload  # Nd ~f32, axis=<given>  (default)
def kurtosistest(
    a: onp.ToJustFloat32_ND,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> KurtosistestResult[onp.ArrayND[np.float32] | np.float32]: ...
@overload  # Nd ~f32, keepdims=True
def kurtosistest[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[np.float32, ShapeT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> KurtosistestResult[onp.ArrayND[np.float32, ShapeT]]: ...
@overload  # ?d ~f32, keepdims=True
def kurtosistest(
    a: onp.ToJustFloat32_ND,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> KurtosistestResult[onp.ArrayND[np.float32]]: ...
@overload  # ?d floating, axis=None
def kurtosistest(
    a: onp.ToArrayND[npc.floating, npc.floating],
    axis: None,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> KurtosistestResult[np.float64 | Any]: ...
@overload  # Nd floating, axis=<given>  (default)
def kurtosistest(
    a: onp.ToArrayND[npc.floating, npc.floating],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> KurtosistestResult[onp.ArrayND[np.float64 | Any] | Any]: ...
@overload  # Nd floating, keepdims=True
def kurtosistest[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.floating, ShapeT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> KurtosistestResult[onp.ArrayND[np.float64 | Any, ShapeT]]: ...
@overload  # ?d floating, keepdims=True
def kurtosistest(
    a: onp.ToArrayND[npc.floating, npc.floating],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> KurtosistestResult[onp.ArrayND[np.float64 | Any]]: ...

# keep in sync with `skewtest` and `kurtosistest`
@overload  # ?d ~f64, axis=None
def normaltest(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer],
    axis: None,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> NormaltestResult[np.float64]: ...
@overload  # ?d ~f64, axis=<given>  (default)
def normaltest(
    a: onp.ArrayND[npc.floating64 | npc.integer, _JustAnyShape],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> NormaltestResult[onp.ArrayND[np.float64] | np.float64]: ...
@overload  # 1d ~f64, axis=<given>  (default)
def normaltest(
    a: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> NormaltestResult[np.float64]: ...
@overload  # 2d ~f64, axis=<given>  (default)
def normaltest(
    a: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> NormaltestResult[onp.Array1D[np.float64]]: ...
@overload  # Nd ~f64, axis=<given>  (default)
def normaltest(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> NormaltestResult[onp.ArrayND[np.float64] | np.float64]: ...
@overload  # Nd ~f64, keepdims=True
def normaltest[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.floating64 | npc.integer, ShapeT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> NormaltestResult[onp.ArrayND[np.float64, ShapeT]]: ...
@overload  # ?d ~f64, keepdims=True
def normaltest(
    a: onp.ToArrayND[float, npc.floating64 | npc.integer],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[True],
) -> NormaltestResult[onp.ArrayND[np.float64]]: ...
@overload  # ?d ~f32, axis=None
def normaltest(
    a: onp.ToJustFloat32_ND, axis: None, nan_policy: NanPolicy = "propagate", *, keepdims: L[False] = False
) -> NormaltestResult[np.float32]: ...
@overload  # ?d ~f32, axis=<given>  (default)
def normaltest(
    a: onp.ArrayND[np.float32, _JustAnyShape], axis: int = 0, nan_policy: NanPolicy = "propagate", *, keepdims: L[False] = False
) -> NormaltestResult[onp.ArrayND[np.float32] | np.float32]: ...
@overload  # 1d ~f32, axis=<given>  (default)
def normaltest(
    a: onp.ToJustFloat32Strict1D, axis: int = 0, nan_policy: NanPolicy = "propagate", *, keepdims: L[False] = False
) -> NormaltestResult[np.float32]: ...
@overload  # 2d ~f32, axis=<given>  (default)
def normaltest(
    a: onp.ToJustFloat32Strict2D, axis: int = 0, nan_policy: NanPolicy = "propagate", *, keepdims: L[False] = False
) -> NormaltestResult[onp.Array1D[np.float32]]: ...
@overload  # Nd ~f32, axis=<given>  (default)
def normaltest(
    a: onp.ToJustFloat32_ND, axis: int = 0, nan_policy: NanPolicy = "propagate", *, keepdims: L[False] = False
) -> NormaltestResult[onp.ArrayND[np.float32] | np.float32]: ...
@overload  # Nd ~f32, keepdims=True
def normaltest[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[np.float32, ShapeT], axis: int | None = 0, nan_policy: NanPolicy = "propagate", *, keepdims: L[True]
) -> NormaltestResult[onp.ArrayND[np.float32, ShapeT]]: ...
@overload  # ?d ~f32, keepdims=True
def normaltest(
    a: onp.ToJustFloat32_ND, axis: int | None = 0, nan_policy: NanPolicy = "propagate", *, keepdims: L[True]
) -> NormaltestResult[onp.ArrayND[np.float32]]: ...
@overload  # ?d floating, axis=None
def normaltest(
    a: onp.ToArrayND[npc.floating, npc.floating], axis: None, nan_policy: NanPolicy = "propagate", *, keepdims: L[False] = False
) -> NormaltestResult[np.float64 | Any]: ...
@overload  # Nd floating, axis=<given>  (default)
def normaltest(
    a: onp.ToArrayND[npc.floating, npc.floating],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> NormaltestResult[onp.ArrayND[np.float64 | Any] | Any]: ...
@overload  # Nd floating, keepdims=True
def normaltest[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.floating, ShapeT], axis: int | None = 0, nan_policy: NanPolicy = "propagate", *, keepdims: L[True]
) -> NormaltestResult[onp.ArrayND[np.float64 | Any, ShapeT]]: ...
@overload  # ?d floating, keepdims=True
def normaltest(
    a: onp.ToArrayND[npc.floating, npc.floating], axis: int | None = 0, nan_policy: NanPolicy = "propagate", *, keepdims: L[True]
) -> NormaltestResult[onp.ArrayND[np.float64 | Any]]: ...

# keep in sync with `skewtest`, `kurtosistest`, and `normaltest` (but with axis=None instead of axis=0)
@overload  # ?d ~f64, axis=None (default)
def jarque_bera(
    x: onp.ToArrayND[float, npc.floating64 | npc.integer],
    *,
    axis: None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[np.float64]: ...
@overload  # ?d ~f64, axis=<given>
def jarque_bera(
    x: onp.ToArrayND[float, npc.floating64 | npc.integer],
    *,
    axis: int,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[onp.ArrayND[np.float64]]: ...
@overload  # Nd ~f64, keepdims=True
def jarque_bera[ShapeT: tuple[int, ...]](
    x: onp.ArrayND[npc.floating64 | npc.integer, ShapeT],
    *,
    axis: int | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> SignificanceResult[onp.ArrayND[np.float64, ShapeT]]: ...
@overload  # ?d ~f64, keepdims=True
def jarque_bera(
    x: onp.ToArrayND[float, npc.floating64 | npc.integer],
    *,
    axis: int | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> SignificanceResult[onp.ArrayND[np.float64]]: ...
@overload  # ?d ~f32, axis=None (default)
def jarque_bera(
    x: onp.ToJustFloat32_ND, *, axis: None = None, nan_policy: NanPolicy = "propagate", keepdims: L[False] = False
) -> SignificanceResult[np.float32]: ...
@overload  # ?d ~f32, axis=<given>
def jarque_bera(
    x: onp.ToJustFloat32_ND, *, axis: int, nan_policy: NanPolicy = "propagate", keepdims: L[False] = False
) -> SignificanceResult[onp.ArrayND[np.float32]]: ...
@overload  # Nd ~f32, keepdims=True
def jarque_bera[ShapeT: tuple[int, ...]](
    x: onp.ArrayND[np.float32, ShapeT], *, axis: int | None = None, nan_policy: NanPolicy = "propagate", keepdims: L[True]
) -> SignificanceResult[onp.ArrayND[np.float32, ShapeT]]: ...
@overload  # ?d ~f32, keepdims=True
def jarque_bera(
    x: onp.ToJustFloat32_ND, *, axis: int | None = None, nan_policy: NanPolicy = "propagate", keepdims: L[True]
) -> SignificanceResult[onp.ArrayND[np.float32]]: ...
@overload  # ?d floating, axis=None (default)
def jarque_bera(
    x: onp.ToArrayND[npc.floating, npc.floating],
    *,
    axis: None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[np.float64 | Any]: ...
@overload  # ?d floating, axis=<given>
def jarque_bera(
    x: onp.ToArrayND[npc.floating, npc.floating], *, axis: int, nan_policy: NanPolicy = "propagate", keepdims: L[False] = False
) -> SignificanceResult[onp.ArrayND[np.float64 | Any]]: ...
@overload  # Nd floating, keepdims=True
def jarque_bera[ShapeT: tuple[int, ...]](
    x: onp.ArrayND[npc.floating, ShapeT], *, axis: int | None = None, nan_policy: NanPolicy = "propagate", keepdims: L[True]
) -> SignificanceResult[onp.ArrayND[np.float64 | Any, ShapeT]]: ...
@overload  # ?d floating, keepdims=True
def jarque_bera(
    x: onp.ToArrayND[npc.floating, npc.floating],
    *,
    axis: int | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> SignificanceResult[onp.ArrayND[np.float64 | Any]]: ...

# keep in sync with `percentileofscore`
@overload
def scoreatpercentile(
    a: onp.ToFloat1D,
    per: onp.ToFloat,
    limit: _RealLimits | tuple[()] = (),
    interpolation_method: _QuantileInterpolation = "fraction",
    axis: int | None = None,
) -> np.float64: ...
@overload
def scoreatpercentile(
    a: onp.ToFloat1D,
    per: Sequence[onp.ToFloat],
    limit: _RealLimits | tuple[()] = (),
    interpolation_method: _QuantileInterpolation = "fraction",
    axis: int | None = None,
) -> onp.Array1D[np.float64]: ...
@overload
def scoreatpercentile[ShapeT: tuple[int, ...]](
    a: onp.ToFloat1D,
    per: onp.ArrayND[npc.floating | npc.integer, ShapeT],
    limit: _RealLimits | tuple[()] = (),
    interpolation_method: _QuantileInterpolation = "fraction",
    axis: int | None = None,
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload
def scoreatpercentile(
    a: onp.ToFloat1D,
    per: onp.ToFloatND,
    limit: _RealLimits | tuple[()] = (),
    interpolation_method: _QuantileInterpolation = "fraction",
    axis: int | None = None,
) -> onp.ArrayND[np.float64]: ...

# keep in sync with `scoreatpercentile`
@overload
def percentileofscore(
    a: onp.ToFloat1D, score: onp.ToFloat, kind: _PercentileInterpolation = "rank", nan_policy: NanPolicy = "propagate"
) -> np.float64: ...
@overload
def percentileofscore(
    a: onp.ToFloat1D, score: Sequence[onp.ToFloat], kind: _PercentileInterpolation = "rank", nan_policy: NanPolicy = "propagate"
) -> onp.Array1D[np.float64]: ...
@overload
def percentileofscore[ShapeT: tuple[int, ...]](
    a: onp.ToFloat1D,
    score: onp.ArrayND[npc.floating, ShapeT],
    kind: _PercentileInterpolation = "rank",
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload
def percentileofscore(
    a: onp.ToFloat1D, score: onp.ToFloatND, kind: _PercentileInterpolation = "rank", nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[np.float64]: ...

#
def cumfreq(
    a: onp.ToFloatND, numbins: int = 10, defaultreallimits: _RealLimits | None = None, weights: onp.ToFloatND | None = None
) -> CumfreqResult: ...

#
def relfreq(
    a: onp.ToFloatND, numbins: int = 10, defaultreallimits: _RealLimits | None = None, weights: onp.ToFloatND | None = None
) -> RelfreqResult: ...

#
@overload
def obrientransform(*, nan_policy: NanPolicy = "propagate") -> tuple[()]: ...
@overload
def obrientransform(x0: onp.ToFloatND, /, *, nan_policy: NanPolicy = "propagate") -> tuple[onp.ArrayND[np.float64]]: ...
@overload
def obrientransform(
    x0: onp.ToFloatND, x1: onp.ToFloatND, /, *, nan_policy: NanPolicy = "propagate"
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
@overload
def obrientransform(
    x0: onp.ToFloatND, x1: onp.ToFloatND, x2: onp.ToFloatND, /, *, nan_policy: NanPolicy = "propagate"
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
@overload
def obrientransform(
    x0: onp.ToFloatND, x1: onp.ToFloatND, x2: onp.ToFloatND, x3: onp.ToFloatND, /, *, nan_policy: NanPolicy = "propagate"
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.float64], onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
@overload
def obrientransform(*samples: onp.ToFloatND, nan_policy: NanPolicy = "propagate") -> tuple[onp.ArrayND[np.float64], ...]: ...

#
@overload  # 1d ~inexact64 | +integer, keepdims=False (default)
def sem(
    a: onp.ToArrayStrict1D[complex, npc.inexact64 | npc.integer | np.bool],
    axis: L[0, -1] | None = 0,
    ddof: int = 1,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # >1d ~inexact64 | +integer, axis: int (default)
def sem(
    a: onp.CanArray[onp.AtLeast2D, np.dtype[npc.inexact64 | npc.integer | np.bool]] | Sequence[onp.SequenceND[complex]],
    axis: int = 0,
    ddof: int = 1,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d ~inexact64 | +integer, axis=None, keepdims=False (default)
def sem(
    a: onp.ToArrayND[complex, npc.inexact64 | npc.integer | np.bool],
    axis: None,
    ddof: int = 1,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # ?d ~inexact64 | +integer, keepdims=True
def sem(
    a: onp.ToArrayND[complex, npc.inexact64 | npc.integer | np.bool],
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

# NOTE: keep in sync with `gzscore` and `zmap`
@overload  # +integer, known shape
def zscore[ShapeT: tuple[int, ...]](
    a: nptc.CanArray[ShapeT, np.dtype[npc.integer | np.bool]],
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # known inexact dtype, known shape
def zscore[ShapeT: tuple[int, ...], InexactT: npc.inexact](
    a: nptc.CanArray[ShapeT, np.dtype[InexactT]], axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[InexactT, ShapeT]: ...
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

# NOTE: keep in sync with `zscore` and `zmap`
@overload  # +integer, known shape
def gzscore[ShapeT: tuple[int, ...]](
    a: nptc.CanArray[ShapeT, np.dtype[npc.integer | np.bool]],
    *,
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # known inexact dtype, known shape
def gzscore[ShapeT: tuple[int, ...], InexactT: npc.inexact](
    a: nptc.CanArray[ShapeT, np.dtype[InexactT]], *, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[InexactT, ShapeT]: ...
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

# keep roughly in sync with `zscore` and `gzscore`
@overload  # +integer, known shape
def zmap[ShapeT: tuple[int, ...]](  # type: ignore[overload-overlap]
    scores: nptc.CanArray[ShapeT, np.dtype[npc.floating64 | npc.integer | np.bool]],
    compare: nptc.CanArray[ShapeT, np.dtype[npc.floating64 | npc.integer | np.bool]],
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # known inexact dtype, known shape
def zmap[ShapeT: tuple[int, ...], InexactT: npc.inexact](
    scores: nptc.CanArray[ShapeT, np.dtype[InexactT]],
    compare: nptc.CanArray[ShapeT, np.dtype[InexactT]],
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[InexactT, ShapeT]: ...
@overload  # float 1d
def zmap(
    scores: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool],
    compare: onp.ToFloat64Strict1D,
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float64]: ...
@overload  # float 1d
def zmap(
    scores: onp.ToFloat64Strict1D,
    compare: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool],
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float64]: ...
@overload  # float 2d
def zmap(
    scores: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool],
    compare: onp.ToFloat64Strict2D | onp.ToFloat64Strict1D,
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array2D[np.float64]: ...
@overload  # float 2d
def zmap(
    scores: onp.ToFloat64Strict2D | onp.ToFloat64Strict1D,
    compare: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool],
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array2D[np.float64]: ...
@overload  # complex 1d
def zmap(
    scores: onp.ToJustComplex128Strict1D,
    compare: onp.ToComplex128Strict1D,
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.complex128]: ...
@overload  # complex 1d
def zmap(
    scores: onp.ToComplex128Strict1D,
    compare: onp.ToJustComplex128Strict1D,
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.complex128]: ...
@overload  # complex 2d
def zmap(
    scores: onp.ToJustComplex128Strict2D,
    compare: onp.ToComplex128Strict2D | onp.ToComplex128Strict1D,
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array2D[np.complex128]: ...
@overload  # complex 2d
def zmap(
    scores: onp.ToComplex128Strict2D | onp.ToComplex128Strict1D,
    compare: onp.ToJustComplex128Strict2D,
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array2D[np.complex128]: ...
@overload  # floating fallback
def zmap(  # the weird shape-type is a workaround for a bug in pyright's overlapping overload detection on numpy<2.1
    scores: onp.ToFloatND, compare: onp.ToFloatND, axis: int | None = 0, ddof: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[npc.floating, tuple[int] | tuple[Any, ...]]: ...
@overload  # complex fallback
def zmap(
    scores: onp.ToComplexND,
    compare: onp.ToJustComplexND,
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[npc.complexfloating]: ...
@overload  # complex fallback
def zmap(
    scores: onp.ToJustComplexND,
    compare: onp.ToComplexND,
    axis: int | None = 0,
    ddof: int = 0,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[npc.complexfloating]: ...

#
@overload  # T@floating, axis=None (default)
def iqr[FloatT: npc.floating](
    x: onp.ToArrayND[FloatT, FloatT],
    axis: None = None,
    rng: tuple[float, float] = (25, 75),
    scale: L["normal"] | onp.ToFloat | onp.ToFloatND = 1.0,
    nan_policy: NanPolicy = "propagate",
    interpolation: _InterpolationMethod = "linear",
    keepdims: L[False] = False,
) -> FloatT: ...
@overload  # T@floating, keepdims=True
def iqr[FloatT: npc.floating](
    x: onp.ToArrayND[FloatT, FloatT],
    axis: int | Sequence[int] | None = None,
    rng: tuple[float, float] = (25, 75),
    scale: L["normal"] | onp.ToFloat | onp.ToFloatND = 1.0,
    nan_policy: NanPolicy = "propagate",
    interpolation: _InterpolationMethod = "linear",
    *,
    keepdims: L[True],
) -> onp.ArrayND[FloatT]: ...
@overload  # T@floating, axis=<given>
def iqr[FloatT: npc.floating](
    x: onp.ToArrayND[FloatT, FloatT],
    axis: int | Sequence[int],
    rng: tuple[float, float] = (25, 75),
    scale: L["normal"] | onp.ToFloat | onp.ToFloatND = 1.0,
    nan_policy: NanPolicy = "propagate",
    interpolation: _InterpolationMethod = "linear",
    keepdims: L[False] = False,
) -> onp.ArrayND[FloatT]: ...
@overload  # +f64, axis=None (default)
def iqr(
    x: onp.ToArrayND[float, npc.integer],
    axis: None = None,
    rng: tuple[float, float] = (25, 75),
    scale: L["normal"] | onp.ToFloat | onp.ToFloatND = 1.0,
    nan_policy: NanPolicy = "propagate",
    interpolation: _InterpolationMethod = "linear",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # +f64, keepdims=True
def iqr(
    x: onp.ToArrayND[float, npc.integer],
    axis: int | Sequence[int] | None = None,
    rng: tuple[float, float] = (25, 75),
    scale: L["normal"] | onp.ToFloat | onp.ToFloatND = 1.0,
    nan_policy: NanPolicy = "propagate",
    interpolation: _InterpolationMethod = "linear",
    *,
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # +f64, axis=<given>
def iqr(
    x: onp.ToArrayND[float, npc.integer],
    axis: int | Sequence[int],
    rng: tuple[float, float] = (25, 75),
    scale: L["normal"] | onp.ToFloat | onp.ToFloatND = 1.0,
    nan_policy: NanPolicy = "propagate",
    interpolation: _InterpolationMethod = "linear",
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64]: ...

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
@overload
def sigmaclip[IntT: npc.integer](
    a: onp.ArrayND[IntT], low: float = 4.0, high: float = 4.0, *, nan_policy: NanPolicy = "propagate"
) -> SigmaclipResult[IntT, np.float64]: ...
@overload
def sigmaclip[FloatT: npc.floating](
    a: onp.ArrayND[FloatT], low: float = 4.0, high: float = 4.0, *, nan_policy: NanPolicy = "propagate"
) -> SigmaclipResult[FloatT, FloatT]: ...
@overload
def sigmaclip(
    a: onp.SequenceND[int], low: float = 4.0, high: float = 4.0, *, nan_policy: NanPolicy = "propagate"
) -> SigmaclipResult[np.int_, np.float64]: ...
@overload
def sigmaclip(
    a: onp.SequenceND[list[float]] | list[float], low: float = 4.0, high: float = 4.0, *, nan_policy: NanPolicy = "propagate"
) -> SigmaclipResult[np.float64, np.float64]: ...
@overload
def sigmaclip(
    a: onp.ToFloatND, low: float = 4.0, high: float = 4.0, *, nan_policy: NanPolicy = "propagate"
) -> SigmaclipResult: ...

# TODO(jorenham): improve
def trimboth(a: onp.ToFloatND, proportiontocut: float, axis: int | None = 0) -> onp.ArrayND[_Real0D]: ...

# TODO(jorenham): improve
def trim1(a: onp.ToFloatND, proportiontocut: float, tail: _TrimTail = "right", axis: int | None = 0) -> onp.ArrayND[_Real0D]: ...

#
@overload
def trim_mean(
    a: _ToFloatStrictND, proportiontocut: float, axis: int = 0, *, nan_policy: NanPolicy = "propagate", keepdims: L[False] = False
) -> np.float64 | onp.ArrayND[np.float64]: ...
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
) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload
def trim_mean(
    a: onp.ToFloatND, proportiontocut: float, axis: int = 0, *, nan_policy: NanPolicy = "propagate", keepdims: L[True]
) -> onp.ArrayND[np.float64]: ...

#
@overload  # ?d, ?d|1d
def f_oneway(
    sample1: _ToFloatStrictND,
    sample2: _ToFloatStrictND | onp.ToFloatStrict1D,
    /,
    *samples: _ToFloatStrictND | onp.ToFloatStrict1D,
    equal_var: bool = True,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> F_onewayResult[np.float64 | Any]: ...
@overload  # ?d|1d, ?d
def f_oneway(
    sample1: _ToFloatStrictND | onp.ToFloatStrict1D,
    sample2: _ToFloatStrictND,
    /,
    *samples: _ToFloatStrictND | onp.ToFloatStrict1D,
    equal_var: bool = True,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> F_onewayResult[np.float64 | Any]: ...
@overload  # ?d, 2d|3d
def f_oneway(
    sample1: _ToFloatStrictND,
    sample2: onp.ToFloatStrict2D | onp.ToFloatStrict3D,
    /,
    *samples: _ToFloatStrictND | onp.ToFloatStrict1D,
    equal_var: bool = True,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> F_onewayResult[onp.ArrayND[np.float64]]: ...
@overload  # 2d|3d, ?d
def f_oneway(
    sample1: onp.ToFloatStrict2D | onp.ToFloatStrict3D,
    sample2: _ToFloatStrictND,
    /,
    *samples: _ToFloatStrictND | onp.ToFloatStrict1D,
    equal_var: bool = True,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> F_onewayResult[onp.ArrayND[np.float64]]: ...
@overload  # 1d, 1d
def f_oneway(
    sample1: onp.ToFloatStrict1D,
    sample2: onp.ToFloatStrict1D,
    /,
    *samples: onp.ToFloatStrict1D,
    equal_var: bool = True,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> F_onewayResult[np.float64]: ...
@overload  # 2d, <=2d
def f_oneway(
    sample1: onp.ToFloatStrict2D,
    sample2: onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    /,
    *samples: onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    equal_var: bool = True,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> F_onewayResult[onp.Array1D[np.float64]]: ...
@overload  # <=2d, 2d
def f_oneway(
    sample1: onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    sample2: onp.ToFloatStrict2D,
    /,
    *samples: onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    equal_var: bool = True,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> F_onewayResult[onp.Array1D[np.float64]]: ...
@overload  # 3d, <=3d
def f_oneway(
    sample1: onp.ToFloatStrict3D,
    sample2: onp.ToFloatStrict3D | onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    /,
    *samples: onp.ToFloatStrict3D | onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    equal_var: bool = True,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> F_onewayResult[onp.Array2D[np.float64]]: ...
@overload  # <=3d, 3d
def f_oneway(
    sample1: onp.ToFloatStrict3D | onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    sample2: onp.ToFloatStrict3D,
    /,
    *samples: onp.ToFloatStrict3D | onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    equal_var: bool = True,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> F_onewayResult[onp.Array2D[np.float64]]: ...
@overload  # Nd, Nd
def f_oneway(
    sample1: onp.ToFloatND,
    sample2: onp.ToFloatND,
    /,
    *samples: onp.ToFloatND,
    equal_var: bool = True,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> F_onewayResult[onp.ArrayND[np.float64] | Any]: ...
@overload  # axis=None
def f_oneway(
    sample1: onp.ToFloatND,
    sample2: onp.ToFloatND,
    /,
    *samples: onp.ToFloatND,
    equal_var: bool = True,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> F_onewayResult[np.float64]: ...
@overload  # keepdims=True
def f_oneway(
    sample1: onp.ToFloatND,
    sample2: onp.ToFloatND,
    /,
    *samples: onp.ToFloatND,
    equal_var: bool = True,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> F_onewayResult[onp.ArrayND[np.float64]]: ...

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
    a: _ToFloatStrictND,
    b: _ToFloatStrictND,
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
    x: onp.ArrayND[npc.number | np.bool, _JustAnyShape],
    y: onp.ArrayND[npc.number | np.bool, _JustAnyShape],
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
def pack_TtestResult[FloatOrArrayT: _ScalarOrND[npc.floating]](
    statistic: FloatOrArrayT,
    pvalue: FloatOrArrayT,
    df: FloatOrArrayT,
    alternative: Alternative,
    standard_error: FloatOrArrayT,
    estimate: FloatOrArrayT,
) -> TtestResult[FloatOrArrayT]: ...  # undocumented

#
def unpack_TtestResult[FloatOrArrayT: _ScalarOrND[npc.floating]](
    res: TtestResult[FloatOrArrayT], _: int
) -> tuple[
    FloatOrArrayT,  # statistic
    FloatOrArrayT,  # pvalue
    FloatOrArrayT,  # df
    Alternative,  # _alternative
    FloatOrArrayT,  # _standard_error
    FloatOrArrayT,  # _estimate
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
) -> TtestResult[np.float64, np.int_] | TtestResult[onp.ArrayND[np.float64], onp.ArrayND[np.int_]]: ...

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

# keep in sync with `ttest_rel`
@overload  # ?d ~float64
def ttest_ind(
    a: onp.ArrayND[npc.floating64 | npc.integer | np.bool, _JustAnyShape],
    b: onp.ArrayND[npc.floating64 | npc.integer | np.bool, _JustAnyShape],
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
def ttest_ind[FloatT: np.float32 | np.float16](
    a: onp.ArrayND[FloatT, _JustAnyShape],
    b: onp.ArrayND[FloatT, _JustAnyShape],
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[FloatT | Any]: ...
@overload  # 1d ~f64
def ttest_ind(
    a: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool],
    b: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool],
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
def ttest_ind[FloatT: np.float32 | np.float16](
    a: onp.ToArrayStrict1D[FloatT, FloatT],
    b: onp.ToArrayStrict1D[FloatT, FloatT],
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[FloatT]: ...
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
    a: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool],
    b: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool],
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
def ttest_ind[FloatT: np.float32 | np.float16](
    a: onp.ToArrayStrict2D[FloatT, FloatT],
    b: onp.ToArrayStrict2D[FloatT, FloatT],
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array1D[FloatT]]: ...
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
    a: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool],
    b: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool],
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
def ttest_ind[FloatT: np.float32 | np.float16](
    a: onp.ToArrayStrict3D[FloatT, FloatT],
    b: onp.ToArrayStrict3D[FloatT, FloatT],
    *,
    axis: int = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array2D[FloatT]]: ...
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
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    b: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
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
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    b: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
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
def ttest_ind[FloatT: np.float32 | np.float16](
    a: onp.ToArrayND[FloatT, FloatT],
    b: onp.ToArrayND[FloatT, FloatT],
    *,
    axis: None,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[False] = False,
) -> TtestResult[FloatT]: ...
@overload  # nd ~T, keepdims=True
def ttest_ind[FloatT: np.float32 | np.float16](
    a: onp.ToArrayND[FloatT, FloatT],
    b: onp.ToArrayND[FloatT, FloatT],
    *,
    axis: int | None = 0,
    equal_var: bool = True,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    trim: onp.ToFloat = 0,
    method: ResamplingMethod | None = None,
    keepdims: L[True],
) -> TtestResult[onp.ArrayND[FloatT]]: ...
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
    a: onp.ArrayND[npc.floating64 | npc.integer | np.bool, _JustAnyShape],
    b: onp.ArrayND[npc.floating64 | npc.integer | np.bool, _JustAnyShape],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[np.float64 | Any, np.int_ | Any]: ...
@overload  # ?d ~T
def ttest_rel[FloatT: np.float32 | np.float16](
    a: onp.ArrayND[FloatT, _JustAnyShape],
    b: onp.ArrayND[FloatT, _JustAnyShape],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[FloatT | Any, np.int_ | Any]: ...
@overload  # 1d ~f64
def ttest_rel(
    a: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool],
    b: onp.ToArrayStrict1D[float, npc.floating64 | npc.integer | np.bool],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[np.float64, np.int_]: ...
@overload  # 1d ~T
def ttest_rel[FloatT: np.float32 | np.float16](
    a: onp.ToArrayStrict1D[FloatT, FloatT],
    b: onp.ToArrayStrict1D[FloatT, FloatT],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[FloatT, np.int_]: ...
@overload  # 1d +floating
def ttest_rel(
    a: onp.ToFloatStrict1D,
    b: onp.ToFloatStrict1D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[np.float64 | Any, np.int_]: ...
@overload  # 2d ~f64
def ttest_rel(
    a: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool],
    b: onp.ToArrayStrict2D[float, npc.floating64 | npc.integer | np.bool],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array1D[np.float64], onp.Array1D[np.int_]]: ...
@overload  # 2d ~T
def ttest_rel[FloatT: np.float32 | np.float16](
    a: onp.ToArrayStrict2D[FloatT, FloatT],
    b: onp.ToArrayStrict2D[FloatT, FloatT],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array1D[FloatT], onp.Array1D[np.int_]]: ...
@overload  # 2d +floating
def ttest_rel(
    a: onp.ToFloatStrict2D,
    b: onp.ToFloatStrict2D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array1D[np.float64 | Any], onp.Array1D[np.int_]]: ...
@overload  # 3d ~f64
def ttest_rel(
    a: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool],
    b: onp.ToArrayStrict3D[float, npc.floating64 | npc.integer | np.bool],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array2D[np.float64], onp.Array2D[np.int_]]: ...
@overload  # 3d ~T
def ttest_rel[FloatT: np.float32 | np.float16](
    a: onp.ToArrayStrict3D[FloatT, FloatT],
    b: onp.ToArrayStrict3D[FloatT, FloatT],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array2D[FloatT], onp.Array1D[np.int_]]: ...
@overload  # 3d +floating
def ttest_rel(
    a: onp.ToFloatStrict3D,
    b: onp.ToFloatStrict3D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[onp.Array2D[np.float64 | Any], onp.Array2D[np.int_]]: ...
@overload  # nd ~f64, axis=None
def ttest_rel(  # type: ignore[overload-overlap]
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    b: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    axis: None,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[np.float64, np.int_]: ...
@overload  # nd ~f64, keepdims=True
def ttest_rel(  # type: ignore[overload-overlap]
    a: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    b: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> TtestResult[onp.ArrayND[np.float64], onp.ArrayND[np.int_]]: ...
@overload  # nd ~T, axis=None
def ttest_rel[FloatT: np.float32 | np.float16](
    a: onp.ToArrayND[FloatT, FloatT],
    b: onp.ToArrayND[FloatT, FloatT],
    axis: None,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[FloatT, np.int_]: ...
@overload  # nd ~T, keepdims=True
def ttest_rel[FloatT: np.float32 | np.float16](
    a: onp.ToArrayND[FloatT, FloatT],
    b: onp.ToArrayND[FloatT, FloatT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> TtestResult[onp.ArrayND[FloatT], onp.ArrayND[np.int_]]: ...
@overload  # nd +floating, axis=None
def ttest_rel(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[np.float64 | Any, np.int64]: ...
@overload  # nd +floating, keepdims=True
def ttest_rel(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[True],
) -> TtestResult[onp.ArrayND[np.float64 | Any], onp.ArrayND[np.int_]]: ...
@overload  # nd +floating
def ttest_rel(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
    *,
    keepdims: L[False] = False,
) -> TtestResult[onp.ArrayND[np.float64 | Any], onp.ArrayND[np.int_]] | TtestResult[np.float64 | Any, np.int_]: ...

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

#
@overload  # ?d, ?d|1d
def ks_2samp(
    data1: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
    data2: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape] | onp.ToFloatStrict1D,
    alternative: Alternative = "two-sided",
    method: _KS2TestMethod = "auto",
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> KstestResult[np.float64 | Any, np.int8 | Any]: ...
@overload  # ?d|1d, ?d
def ks_2samp(
    data1: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape] | onp.ToFloatStrict1D,
    data2: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
    alternative: Alternative = "two-sided",
    method: _KS2TestMethod = "auto",
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> KstestResult[np.float64 | Any, np.int8 | Any]: ...
@overload  # ?d, 2d|3d
def ks_2samp(
    data1: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
    data2: onp.ToFloatStrict2D | onp.ToFloatStrict3D,
    alternative: Alternative = "two-sided",
    method: _KS2TestMethod = "auto",
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _KstestResultN: ...
@overload  # 2d, ?d
def ks_2samp(
    data1: onp.ToFloatStrict2D | onp.ToFloatStrict3D,
    data2: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
    alternative: Alternative = "two-sided",
    method: _KS2TestMethod = "auto",
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _KstestResultN: ...
@overload  # 1d, 1d
def ks_2samp(
    data1: onp.ToFloatStrict1D,
    data2: onp.ToFloatStrict1D,
    alternative: Alternative = "two-sided",
    method: _KS2TestMethod = "auto",
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _KstestResult0: ...
@overload  # 2d, <=2d
def ks_2samp(
    data1: onp.ToFloatStrict2D,
    data2: onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    alternative: Alternative = "two-sided",
    method: _KS2TestMethod = "auto",
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _KstestResult1: ...
@overload  # <=2d, 2d
def ks_2samp(
    data1: onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    data2: onp.ToFloatStrict2D,
    alternative: Alternative = "two-sided",
    method: _KS2TestMethod = "auto",
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _KstestResult1: ...
@overload  # 3d, <=3d
def ks_2samp(
    data1: onp.ToFloatStrict3D,
    data2: onp.ToFloatStrict3D | onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    alternative: Alternative = "two-sided",
    method: _KS2TestMethod = "auto",
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _KstestResult2: ...
@overload  # <=3d, 3d
def ks_2samp(
    data1: onp.ToFloatStrict3D | onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    data2: onp.ToFloatStrict3D,
    alternative: Alternative = "two-sided",
    method: _KS2TestMethod = "auto",
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _KstestResult2: ...
@overload  # Nd
def ks_2samp(
    data1: onp.ToFloatND,
    data2: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    method: _KS2TestMethod = "auto",
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> KstestResult[np.float64 | Any, np.int8 | Any]: ...
@overload  # keepdims=True
def ks_2samp(
    data1: onp.ToFloatND,
    data2: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    method: _KS2TestMethod = "auto",
    *,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> _KstestResultN: ...
@overload  # axis=None
def ks_2samp(
    data1: onp.ToFloatND,
    data2: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    method: _KS2TestMethod = "auto",
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> _KstestResult0: ...

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

#
@overload  # ?d, ?d|1d
def kruskal(
    sample1: _ToFloatStrictND,
    sample2: _ToFloatStrictND | onp.ToFloatStrict1D,
    /,
    *samples: _ToFloatStrictND | onp.ToFloatStrict1D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> KruskalResult[np.float64 | Any]: ...
@overload  # ?d|1d, ?d
def kruskal(
    sample1: _ToFloatStrictND | onp.ToFloatStrict1D,
    sample2: _ToFloatStrictND,
    /,
    *samples: _ToFloatStrictND | onp.ToFloatStrict1D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> KruskalResult[np.float64 | Any]: ...
@overload  # ?d, 2d|3d
def kruskal(
    sample1: _ToFloatStrictND,
    sample2: onp.ToFloatStrict2D | onp.ToFloatStrict3D,
    /,
    *samples: _ToFloatStrictND | onp.ToFloatStrict1D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> KruskalResult[onp.ArrayND[np.float64]]: ...
@overload  # 2d|3d, ?d
def kruskal(
    sample1: onp.ToFloatStrict2D | onp.ToFloatStrict3D,
    sample2: _ToFloatStrictND,
    /,
    *samples: _ToFloatStrictND | onp.ToFloatStrict1D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> KruskalResult[onp.ArrayND[np.float64]]: ...
@overload  # 1d, 1d
def kruskal(
    sample1: onp.ToFloatStrict1D,
    sample2: onp.ToFloatStrict1D,
    /,
    *samples: onp.ToFloatStrict1D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> KruskalResult[np.float64]: ...
@overload  # 2d, <=2d
def kruskal(
    sample1: onp.ToFloatStrict2D,
    sample2: onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    /,
    *samples: onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> KruskalResult[onp.Array1D[np.float64]]: ...
@overload  # <=2d, 2d
def kruskal(
    sample1: onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    sample2: onp.ToFloatStrict2D,
    /,
    *samples: onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> KruskalResult[onp.Array1D[np.float64]]: ...
@overload  # 3d, <=3d
def kruskal(
    sample1: onp.ToFloatStrict3D,
    sample2: onp.ToFloatStrict3D | onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    /,
    *samples: onp.ToFloatStrict3D | onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> KruskalResult[onp.Array2D[np.float64]]: ...
@overload  # <=3d, 3d
def kruskal(
    sample1: onp.ToFloatStrict3D | onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    sample2: onp.ToFloatStrict3D,
    /,
    *samples: onp.ToFloatStrict3D | onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> KruskalResult[onp.Array2D[np.float64]]: ...
@overload  # Nd, Nd
def kruskal(
    sample1: onp.ToFloatND,
    sample2: onp.ToFloatND,
    /,
    *samples: onp.ToFloatND,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> KruskalResult[onp.ArrayND[np.float64] | Any]: ...
@overload  # axis=None
def kruskal(
    sample1: onp.ToFloatND,
    sample2: onp.ToFloatND,
    /,
    *samples: onp.ToFloatND,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> KruskalResult[np.float64]: ...
@overload  # keepdims=True
def kruskal(
    sample1: onp.ToFloatND,
    sample2: onp.ToFloatND,
    /,
    *samples: onp.ToFloatND,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> KruskalResult[onp.ArrayND[np.float64]]: ...

# TODO(jorenham): improve
def friedmanchisquare(
    *samples: onp.ToFloatND, axis: int | None = 0, nan_policy: NanPolicy = "propagate", keepdims: bool = False
) -> FriedmanchisquareResult: ...

#
@overload  # ?d, ?d|1d
def brunnermunzel(
    x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
    y: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape] | onp.ToFloatStrict1D,
    alternative: Alternative = "two-sided",
    distribution: L["t", "normal"] = "t",
    nan_policy: NanPolicy = "propagate",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
) -> BrunnerMunzelResult[np.float64 | Any]: ...
@overload  # ?d|1d, ?d
def brunnermunzel(
    x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape] | onp.ToFloatStrict1D,
    y: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
    alternative: Alternative = "two-sided",
    distribution: L["t", "normal"] = "t",
    nan_policy: NanPolicy = "propagate",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
) -> BrunnerMunzelResult[np.float64 | Any]: ...
@overload  # ?d, 2d|3d
def brunnermunzel(
    x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
    y: onp.ToFloatStrict2D | onp.ToFloatStrict3D,
    alternative: Alternative = "two-sided",
    distribution: L["t", "normal"] = "t",
    nan_policy: NanPolicy = "propagate",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
) -> BrunnerMunzelResult[onp.ArrayND[np.float64]]: ...
@overload  # 2d|3d, ?d
def brunnermunzel(
    x: onp.ToFloatStrict2D | onp.ToFloatStrict3D,
    y: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
    alternative: Alternative = "two-sided",
    distribution: L["t", "normal"] = "t",
    nan_policy: NanPolicy = "propagate",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
) -> BrunnerMunzelResult[onp.ArrayND[np.float64]]: ...
@overload  # 1d, 1d
def brunnermunzel(
    x: onp.ToFloatStrict1D,
    y: onp.ToFloatStrict1D,
    alternative: Alternative = "two-sided",
    distribution: L["t", "normal"] = "t",
    nan_policy: NanPolicy = "propagate",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
) -> BrunnerMunzelResult[np.float64]: ...
@overload  # 2d, <=2d
def brunnermunzel(
    x: onp.ToFloatStrict2D,
    y: onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    alternative: Alternative = "two-sided",
    distribution: L["t", "normal"] = "t",
    nan_policy: NanPolicy = "propagate",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
) -> BrunnerMunzelResult[onp.Array1D[np.float64]]: ...
@overload  # <=2d, 2d
def brunnermunzel(
    x: onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    y: onp.ToFloatStrict2D,
    alternative: Alternative = "two-sided",
    distribution: L["t", "normal"] = "t",
    nan_policy: NanPolicy = "propagate",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
) -> BrunnerMunzelResult[onp.Array1D[np.float64]]: ...
@overload  # 3d, <=3d
def brunnermunzel(
    x: onp.ToFloatStrict3D,
    y: onp.ToFloatStrict3D | onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    alternative: Alternative = "two-sided",
    distribution: L["t", "normal"] = "t",
    nan_policy: NanPolicy = "propagate",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
) -> BrunnerMunzelResult[onp.Array2D[np.float64]]: ...
@overload  # <=3d, 3d
def brunnermunzel(
    x: onp.ToFloatStrict3D | onp.ToFloatStrict2D | onp.ToFloatStrict1D,
    y: onp.ToFloatStrict3D,
    alternative: Alternative = "two-sided",
    distribution: L["t", "normal"] = "t",
    nan_policy: NanPolicy = "propagate",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
) -> BrunnerMunzelResult[onp.Array2D[np.float64]]: ...
@overload  # Nd, Nd
def brunnermunzel(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    distribution: L["t", "normal"] = "t",
    nan_policy: NanPolicy = "propagate",
    *,
    axis: int = 0,
    keepdims: L[False] = False,
) -> BrunnerMunzelResult[onp.ArrayND[np.float64] | Any]: ...
@overload  # axis=None
def brunnermunzel(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    distribution: L["t", "normal"] = "t",
    nan_policy: NanPolicy = "propagate",
    *,
    axis: None,
    keepdims: L[False] = False,
) -> BrunnerMunzelResult[np.float64]: ...
@overload  # keepdims=True
def brunnermunzel(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    alternative: Alternative = "two-sided",
    distribution: L["t", "normal"] = "t",
    nan_policy: NanPolicy = "propagate",
    *,
    axis: int | None = 0,
    keepdims: L[True],
) -> BrunnerMunzelResult[onp.ArrayND[np.float64]]: ...

#
@overload  # ?d T@floating
def combine_pvalues[FloatT: npc.floating](
    pvalues: onp.ArrayND[FloatT, _JustAnyShape],
    method: _CombinePValuesMethod = "fisher",
    weights: onp.ToFloatND | None = None,
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[FloatT | onp.ArrayND[FloatT]]: ...
@overload  # 1d float
def combine_pvalues(
    pvalues: Sequence[float],
    method: _CombinePValuesMethod = "fisher",
    weights: onp.ToFloat1D | None = None,
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[np.float64]: ...
@overload  # 1d T@floating
def combine_pvalues[FloatT: npc.floating](
    pvalues: onp.Array1D[FloatT],
    method: _CombinePValuesMethod = "fisher",
    weights: onp.ToFloat1D | None = None,
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[FloatT]: ...
@overload  # 2d float
def combine_pvalues(
    pvalues: Sequence[Sequence[float]],
    method: _CombinePValuesMethod = "fisher",
    weights: onp.ToFloat1D | onp.ToFloat2D | None = None,
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[onp.Array1D[np.float64]]: ...
@overload  # 2d T@floating
def combine_pvalues[FloatT: npc.floating](
    pvalues: onp.Array2D[FloatT],
    method: _CombinePValuesMethod = "fisher",
    weights: onp.ToFloat1D | onp.ToFloat2D | None = None,
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[onp.Array1D[FloatT]]: ...
@overload  # Nd T@floating
def combine_pvalues[FloatT: npc.floating](
    pvalues: onp.ArrayND[FloatT],
    method: _CombinePValuesMethod = "fisher",
    weights: onp.ToFloatND | None = None,
    *,
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[FloatT | onp.ArrayND[FloatT]]: ...
@overload  # Nd T@floating, axis=None
def combine_pvalues[FloatT: npc.floating](
    pvalues: onp.ArrayND[FloatT],
    method: _CombinePValuesMethod = "fisher",
    weights: onp.ToFloatND | None = None,
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[FloatT]: ...
@overload  # Nd float, axis=None
def combine_pvalues(
    pvalues: onp.SequenceND[float],
    method: _CombinePValuesMethod = "fisher",
    weights: onp.ToFloatND | None = None,
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[np.float64]: ...
@overload  # Nd floating, axis=None
def combine_pvalues(
    pvalues: onp.ToFloatND,
    method: _CombinePValuesMethod = "fisher",
    weights: onp.ToFloatND | None = None,
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> SignificanceResult[np.float64 | Any]: ...
@overload  # Nd T@floating, keepdims=True
def combine_pvalues[FloatT: npc.floating, ShapeT: tuple[int, ...]](
    pvalues: onp.ArrayND[FloatT, ShapeT],
    method: _CombinePValuesMethod = "fisher",
    weights: onp.ToFloatND | None = None,
    *,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> SignificanceResult[onp.ArrayND[FloatT, ShapeT]]: ...
@overload  # Nd float, keepdims=True
def combine_pvalues(
    pvalues: onp.SequenceND[float],
    method: _CombinePValuesMethod = "fisher",
    weights: onp.ToFloatND | None = None,
    *,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> SignificanceResult[onp.ArrayND[np.float64]]: ...
@overload  # ?d floating, keepdims=True
def combine_pvalues(
    pvalues: onp.ToFloatND,
    method: _CombinePValuesMethod = "fisher",
    weights: onp.ToFloatND | None = None,
    *,
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> SignificanceResult[onp.ArrayND[np.float64 | Any]]: ...

#
def fisher_exact(
    table: onp.ArrayND[_Real0D], alternative: Alternative | None = None, *, method: ResamplingMethod | None = None
) -> SignificanceResult[float]: ...

# undocumented
def quantile_test_iv(
    x: onp.ToFloatND,
    q: float | _Real0D,
    p: float | npc.floating,
    alternative: Alternative,
    axis: int | None,
    keepdims: bool | None,
) -> tuple[onp.ArrayND[_Real0D], _Real0D, npc.floating, Alternative]: ...

#
@overload
def quantile_test[FloatT: npc.floating](
    x: onp.ArrayND[FloatT],
    *,
    q: float | _Real0D = 0.0,
    p: float | npc.floating = 0.5,
    alternative: Alternative = "two-sided",
    axis: int | None = 0,
    keepdims: bool | None = None,
) -> QuantileTestResult[FloatT]: ...
@overload
def quantile_test(
    x: onp.ToIntND | onp.ToJustFloat64_ND,
    *,
    q: float | _Real0D = 0.0,
    p: float | npc.floating = 0.5,
    alternative: Alternative = "two-sided",
    axis: int | None = 0,
    keepdims: bool | None = None,
) -> QuantileTestResult[np.float64]: ...
@overload
def quantile_test(
    x: onp.ToFloatND,
    *,
    q: float | _Real0D = 0.0,
    p: float | npc.floating = 0.5,
    alternative: Alternative = "two-sided",
    axis: int | None = 0,
    keepdims: bool | None = None,
) -> QuantileTestResult[np.float64 | Any]: ...

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
def rankdata[ShapeT: tuple[int, ...]](
    a: onp.Array[ShapeT], method: _RankMethod = "average", *, axis: int = 0, nan_policy: NanPolicy = "propagate"
) -> onp.ArrayND[np.float64, ShapeT]: ...
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
@overload  # axis=None (default)
def expectile(
    a: onp.ToFloatND,
    alpha: float = 0.5,
    *,
    weights: onp.ToFloatND | None = None,
    axis: None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> np.float64: ...
@overload  # axis=<given>
def expectile(
    a: onp.ToFloatND,
    alpha: float = 0.5,
    *,
    weights: onp.ToFloatND | None = None,
    axis: int,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[False] = False,
) -> onp.ArrayND[np.float64]: ...
@overload  # known shape, keepdims=True
def expectile[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[npc.floating | npc.integer, ShapeT],
    alpha: float = 0.5,
    *,
    weights: onp.ToFloatND | None = None,
    axis: int | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # known shape, keepdims=True
def expectile(
    a: onp.ToFloatND,
    alpha: float = 0.5,
    *,
    weights: onp.ToFloatND | None = None,
    axis: int | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: L[True],
) -> onp.ArrayND[np.float64]: ...

#
@overload  # ?d, ?d
def linregress(
    x: _ToFloatStrictND,
    y: _ToFloatStrictND,
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
