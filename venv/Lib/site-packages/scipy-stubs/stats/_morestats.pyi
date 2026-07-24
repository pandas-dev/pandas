from collections.abc import Callable
from types import ModuleType
from typing import (
    Any,
    Generic,
    Literal,
    NamedTuple,
    Never,
    Protocol,
    Self,
    SupportsIndex,
    final,
    overload,
    override,
    type_check_only,
)
from typing_extensions import TypeVar, deprecated

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._continuous_distns import gengamma_gen, invgamma_gen, norm_gen, t_gen
from ._distn_infrastructure import rv_continuous_frozen
from ._fit import FitResult
from ._resampling import MonteCarloMethod, PermutationMethod
from ._stats_py import SignificanceResult
from ._typing import Alternative, BaseBunch, NanPolicy
from scipy.optimize import OptimizeResult

__all__ = [
    "anderson",
    "anderson_ksamp",
    "ansari",
    "bartlett",
    "bayes_mvs",
    "boxcox",
    "boxcox_llf",
    "boxcox_normmax",
    "boxcox_normplot",
    "circmean",
    "circstd",
    "circvar",
    "directional_stats",
    "false_discovery_control",
    "fligner",
    "kstat",
    "kstatvar",
    "levene",
    "median_test",
    "mood",
    "mvsdist",
    "ppcc_max",
    "ppcc_plot",
    "probplot",
    "shapiro",
    "wilcoxon",
    "yeojohnson",
    "yeojohnson_llf",
    "yeojohnson_normmax",
    "yeojohnson_normplot",
]

###

_NDT_co = TypeVar(
    "_NDT_co",
    covariant=True,
    bound=np.float64 | onp.ArrayND[np.float64],
    default=np.float64 | onp.ArrayND[np.float64],
)  # fmt: skip

_DirectionT_co = TypeVar("_DirectionT_co", bound=onp.ArrayND[npc.inexact], default=onp.ArrayND[np.float64 | Any], covariant=True)
_LengthT_co = TypeVar(
    "_LengthT_co", bound=npc.floating | onp.ArrayND[npc.floating], default=onp.ArrayND[np.float64 | Any], covariant=True
)

type _JustAnyShape = tuple[Never, Never, Never, Never]  # workaround for https://github.com/microsoft/pyright/issues/10232
type _Tuple2[T] = tuple[T, T]
type _Tuple3[T] = tuple[T, T, T]
type _Float1D = onp.Array1D[np.float64]

type _KStatOrder = Literal[1, 2, 3, 4]
type _CenterMethod = Literal["mean", "median", "trimmed"]
type _RVCAnderson = Literal["norm", "expon", "logistic", "extreme1", "gumbel", "gumbel_l", "gumbel_r", "weibull_min"]
type _RVC0 = Literal[
    "anglit",
    "arcsine",
    "cauchy",
    "cosine",
    "expon",
    "gibrat",
    "gumbel_l",
    "gumbel_r",
    "halfcauchy",
    "halflogistic",
    "halfnorm",
    "hypsecant",
    "kstwobign",
    "laplace",
    "levy",
    "levy_l",
    "logistic",
    "maxwell",
    "moyal",
    "norm",
    "rayleigh",
    "semicircular",
    "uniform",
    "wald",
]
type _RVC1 = Literal[
    "alpha",
    "argus",
    "bradford",
    "chi",
    "chi2",
    "erlang",
    "exponnorm",
    "exponpow",
    "fatiguelife",
    "fisk",
    "gamma",
    "genextreme",
    "genlogistic",
    "gennorm",
    "genpareto",
    "gompertz",
    "invgamma",
    "invgauss",
    "invweibull",
    "loggamma",
    "loglaplace",
    "lognorm",
    "lomax",
    "nakagami",
    "pareto",
    "pearson3",
    "powerlaw",
    "powernorm",
    "rdist",
    "recipinvgauss",
    "rel_breitwigner",
    "rice",
    "skewcauchy",
    "skewnorm",
    "t",
    "triang",
    "tukeylambda",
    "vonmises",
    "vonmises_line",
    "weibull_max",
    "weibull_min",
    "wrapcauchy",
]
type _AnsariMethod = Literal["auto", "asymptotic", "exact"] | PermutationMethod

type _ObjFun1D = Callable[[float], float | npc.floating]
type _MinFun1D = Callable[[_ObjFun1D], _HasX] | Callable[[_ObjFun1D], OptimizeResult]

type _AndersonResult = FitResult[Callable[[onp.ToFloat, onp.ToFloat], np.float64]]

@type_check_only
class _TestResult(NamedTuple, Generic[_NDT_co]):
    statistic: _NDT_co
    pvalue: _NDT_co

@type_check_only
class _ConfidenceInterval(NamedTuple):
    statistic: float
    minmax: tuple[float, float]

# represents the e.g. `matplotlib.pyplot` module and a `matplotlib.axes.Axes` object with a `plot` and `text` method
@type_check_only
class _CanPlotText(Protocol):
    # NOTE: `Any` is required as return type because it's covariant, and not shouldn't be `Never`.
    def plot(self, /, *args: float | onp.ToFloatND | str, **kwargs: object) -> Any: ...
    def text(self, /, x: float, y: float, s: str, fontdict: dict[str, Any] | None = None, **kwargs: object) -> Any: ...

@type_check_only
class _CanPPF(Protocol):
    def ppf(self, q: onp.ArrayND[np.float64], /) -> onp.ArrayND[np.float64]: ...

@type_check_only
class _HasX(Protocol):
    x: float | npc.floating

###

@final
class _BigFloat: ...

@final
class DirectionalStats(Generic[_DirectionT_co, _LengthT_co]):
    mean_direction: _DirectionT_co
    mean_resultant_length: _LengthT_co

    def __init__(self, /, mean_direction: _DirectionT_co, mean_resultant_length: _LengthT_co) -> None: ...

class ShapiroResult(_TestResult[_NDT_co], Generic[_NDT_co]): ...
class AnsariResult(_TestResult[_NDT_co], Generic[_NDT_co]): ...
class BartlettResult(_TestResult[_NDT_co], Generic[_NDT_co]): ...
class LeveneResult(_TestResult[_NDT_co], Generic[_NDT_co]): ...
class FlignerResult(_TestResult[_NDT_co], Generic[_NDT_co]): ...

#
class Mean(_ConfidenceInterval): ...
class Variance(_ConfidenceInterval): ...
class Std_dev(_ConfidenceInterval): ...

class AndersonResult(BaseBunch[np.float64, _Float1D, _Float1D]):
    @property
    def statistic(self, /) -> np.float64: ...
    @property
    def critical_values(self, /) -> _Float1D: ...
    @property
    def significance_level(self, /) -> _Float1D: ...
    @property
    def fit_result(self, /) -> _AndersonResult: ...

    #
    @override
    def __new__(  # pyrefly:ignore[bad-override]
        _cls, statistic: np.float64, critical_values: _Float1D, significance_level: _Float1D, *, fit_result: _AndersonResult
    ) -> Self: ...
    @override
    def __init__(  # pyrefly:ignore[bad-override]
        self, /, statistic: np.float64, critical_values: _Float1D, significance_level: _Float1D, *, fit_result: _AndersonResult
    ) -> None: ...

class Anderson_ksampResult(BaseBunch[np.float64, _Float1D, np.float64]):
    @property
    def statistic(self, /) -> np.float64: ...
    @property
    @deprecated("Present only when `variant` is unspecified.")
    def critical_values(self, /) -> _Float1D: ...
    @property
    def pvalue(self, /) -> np.float64: ...

    #
    @override
    def __new__(_cls, statistic: np.float64, critical_values: _Float1D, pvalue: np.float64) -> Self: ...  # pyrefly:ignore[bad-override]
    @override
    def __init__(self, /, statistic: np.float64, critical_values: _Float1D, pvalue: np.float64) -> None: ...  # pyrefly:ignore[bad-override]

class WilcoxonResult(BaseBunch[_NDT_co, _NDT_co], Generic[_NDT_co]):  # pyright: ignore[reportInvalidTypeArguments]  # pyrefly: ignore[invalid-variance]  # zuban: ignore[type-var]
    zstatistic: _NDT_co  # might not be set (depends on `method`)

    @property
    def statistic(self, /) -> _NDT_co: ...
    @property
    def pvalue(self, /) -> _NDT_co: ...

    #
    @override
    def __new__(_cls, statistic: _NDT_co, pvalue: _NDT_co) -> Self: ...  # pyrefly:ignore[bad-override]
    @override
    def __init__(self, /, statistic: _NDT_co, pvalue: _NDT_co) -> None: ...  # pyrefly:ignore[bad-override]

class MedianTestResult(BaseBunch[np.float64, np.float64, np.float64, onp.Array2D[np.float64]]):
    @property
    def statistic(self, /) -> np.float64: ...
    @property
    def pvalue(self, /) -> np.float64: ...
    @property
    def median(self, /) -> np.float64: ...
    @property
    def table(self, /) -> onp.Array2D[np.int_]: ...

    #
    @override
    def __new__(_cls, statistic: np.float64, pvalue: np.float64, median: np.float64, table: onp.Array2D[np.float64]) -> Self: ...  # pyrefly:ignore[bad-override]
    @override
    def __init__(  # pyrefly:ignore[bad-override]
        self, /, statistic: np.float64, pvalue: np.float64, median: np.float64, table: onp.Array2D[np.float64]
    ) -> None: ...

def bayes_mvs(data: onp.ToFloatND, alpha: onp.ToFloat = 0.9) -> tuple[Mean, Variance, Std_dev]: ...

#
def mvsdist(
    data: onp.ToFloatND,
) -> (
    tuple[
        rv_continuous_frozen[norm_gen, np.float64],
        rv_continuous_frozen[norm_gen, np.float64],
        rv_continuous_frozen[norm_gen, np.float64],
    ]
    | tuple[
        rv_continuous_frozen[t_gen, np.float64],
        rv_continuous_frozen[invgamma_gen, np.float64],
        rv_continuous_frozen[gengamma_gen, np.float64],
    ]
): ...

#
@overload
def kstat(
    data: onp.ToFloatND,
    n: _KStatOrder = 2,
    *,
    axis: None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload
def kstat(
    data: onp.ToFloatND,
    n: _KStatOrder = 2,
    *,
    axis: SupportsIndex | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.ArrayND[np.float64]: ...
@overload
def kstat(
    data: onp.ToFloatND,
    n: _KStatOrder = 2,
    *,
    axis: SupportsIndex | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...

#
@overload
def kstatvar(
    data: onp.ToFloatND,
    n: _KStatOrder = 2,
    *,
    axis: None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload
def kstatvar(
    data: onp.ToFloatND,
    n: _KStatOrder = 2,
    *,
    axis: SupportsIndex | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> onp.ArrayND[np.float64]: ...
@overload
def kstatvar(
    data: onp.ToFloatND,
    n: _KStatOrder = 2,
    *,
    axis: SupportsIndex | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...

#
@overload
def probplot(
    x: onp.ToFloat | onp.ToFloatND,
    sparams: tuple[()] = (),
    dist: _RVC0 | _CanPPF = "norm",
    fit: Literal[True] = True,
    plot: _CanPlotText | ModuleType | None = None,
    rvalue: bool = False,
) -> tuple[_Tuple2[onp.ArrayND[np.float64]], _Tuple3[np.float64]]: ...
@overload
def probplot(
    x: onp.ToFloat | onp.ToFloatND,
    sparams: tuple[()] = (),
    dist: _RVC0 | _CanPPF = "norm",
    *,
    fit: Literal[False],
    plot: _CanPlotText | ModuleType | None = None,
    rvalue: bool = False,
) -> _Tuple2[onp.ArrayND[np.float64]]: ...
@overload
def probplot(
    x: onp.ToFloat | onp.ToFloatND,
    sparams: tuple[onp.ToFloat, ...],
    dist: str | _CanPPF = "norm",
    fit: Literal[True] = True,
    plot: _CanPlotText | ModuleType | None = None,
    rvalue: bool = False,
) -> tuple[_Tuple2[onp.ArrayND[np.float64]], _Tuple3[np.float64]]: ...
@overload
def probplot(
    x: onp.ToFloat | onp.ToFloatND,
    sparams: tuple[onp.ToFloat],
    dist: str | _CanPPF = "norm",
    *,
    fit: Literal[False],
    plot: _CanPlotText | ModuleType | None = None,
    rvalue: bool = False,
) -> _Tuple2[onp.ArrayND[np.float64]]: ...

#
def ppcc_max(
    x: onp.ToFloat | onp.ToFloatND,
    brack: _Tuple2[onp.ToFloat] | _Tuple3[onp.ToFloat] = (0.0, 1.0),
    dist: _RVC1 | _CanPPF = "tukeylambda",
) -> np.float64: ...

#
def ppcc_plot(
    x: onp.ToFloat | onp.ToFloatND,
    a: onp.ToFloat,
    b: onp.ToFloat,
    dist: _RVC1 | _CanPPF = "tukeylambda",
    plot: _CanPlotText | ModuleType | None = None,
    N: int = 80,
) -> _Tuple2[onp.ArrayND[np.float64]]: ...

# technically this also supports conplex data, but since boxcox is only for real data, we limit to real here as well
@overload
def boxcox_llf(
    lmb: float | np.float64,
    data: onp.ArrayND[npc.integer, _JustAnyShape],
    *,
    axis: int = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64] | np.float64: ...
@overload
def boxcox_llf(
    lmb: float | np.float64,
    data: onp.ToArrayStrict1D[float, npc.integer],
    *,
    axis: int | None = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> np.float64: ...
@overload
def boxcox_llf(
    lmb: float | np.float64,
    data: onp.ToArrayStrict2D[float, npc.integer],
    *,
    axis: int = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float64]: ...
@overload
def boxcox_llf(
    lmb: float | np.float64,
    data: onp.ToArrayND[float, npc.integer],
    *,
    axis: None,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> np.float64: ...
@overload
def boxcox_llf(
    lmb: float | np.float64,
    data: onp.ToArrayND[float, npc.integer],
    *,
    axis: int | None = 0,
    keepdims: Literal[True],
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64]: ...
@overload
def boxcox_llf[FloatingT: npc.floating](
    lmb: float | np.float64,
    data: onp.ArrayND[FloatingT, _JustAnyShape],
    *,
    axis: int = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[FloatingT] | FloatingT: ...
@overload
def boxcox_llf[FloatingT: npc.floating](
    lmb: float | np.float64,
    data: onp.ToArrayStrict1D[FloatingT, FloatingT],
    *,
    axis: int | None = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> FloatingT: ...
@overload
def boxcox_llf[FloatingT: npc.floating](
    lmb: float | np.float64,
    data: onp.ToArrayStrict2D[FloatingT, FloatingT],
    *,
    axis: int = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[FloatingT]: ...
@overload
def boxcox_llf[FloatingT: npc.floating](
    lmb: float | np.float64,
    data: onp.ToArrayND[FloatingT, FloatingT],
    *,
    axis: None,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> FloatingT: ...
@overload
def boxcox_llf[FloatingT: npc.floating](
    lmb: float | np.float64,
    data: onp.ToArrayND[FloatingT, FloatingT],
    *,
    axis: int | None = 0,
    keepdims: Literal[True],
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[FloatingT]: ...
@overload
def boxcox_llf(
    lmb: float | np.float64,
    data: onp.ToFloatStrict1D,
    *,
    axis: int | None = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> np.float64 | Any: ...
@overload
def boxcox_llf(
    lmb: float | np.float64,
    data: onp.ToFloatND,
    *,
    axis: None,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> np.float64 | Any: ...
@overload
def boxcox_llf(
    lmb: float | np.float64,
    data: onp.ToFloatND,
    *,
    axis: int | None = 0,
    keepdims: Literal[True],
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64 | Any]: ...
@overload
def boxcox_llf(
    lmb: float | np.float64,
    data: onp.ToFloatND,
    *,
    axis: int | None = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64 | Any] | Any: ...

# keep in sync with `boxcox_llf`
@overload
def yeojohnson_llf(
    lmb: float | np.float64,
    data: onp.ArrayND[npc.integer, _JustAnyShape],
    *,
    axis: int = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64] | np.float64: ...
@overload
def yeojohnson_llf(
    lmb: float | np.float64,
    data: onp.ToArrayStrict1D[float, npc.integer],
    *,
    axis: int | None = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> np.float64: ...
@overload
def yeojohnson_llf(
    lmb: float | np.float64,
    data: onp.ToArrayStrict2D[float, npc.integer],
    *,
    axis: int = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float64]: ...
@overload
def yeojohnson_llf(
    lmb: float | np.float64,
    data: onp.ToArrayND[float, npc.integer],
    *,
    axis: None,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> np.float64: ...
@overload
def yeojohnson_llf(
    lmb: float | np.float64,
    data: onp.ToArrayND[float, npc.integer],
    *,
    axis: int | None = 0,
    keepdims: Literal[True],
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64]: ...
@overload
def yeojohnson_llf[FloatingT: npc.floating](
    lmb: float | np.float64,
    data: onp.ArrayND[FloatingT, _JustAnyShape],
    *,
    axis: int = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[FloatingT] | FloatingT: ...
@overload
def yeojohnson_llf[FloatingT: npc.floating](
    lmb: float | np.float64,
    data: onp.ToArrayStrict1D[FloatingT, FloatingT],
    *,
    axis: int | None = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> FloatingT: ...
@overload
def yeojohnson_llf[FloatingT: npc.floating](
    lmb: float | np.float64,
    data: onp.ToArrayStrict2D[FloatingT, FloatingT],
    *,
    axis: int = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[FloatingT]: ...
@overload
def yeojohnson_llf[FloatingT: npc.floating](
    lmb: float | np.float64,
    data: onp.ToArrayND[FloatingT, FloatingT],
    *,
    axis: None,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> FloatingT: ...
@overload
def yeojohnson_llf[FloatingT: npc.floating](
    lmb: float | np.float64,
    data: onp.ToArrayND[FloatingT, FloatingT],
    *,
    axis: int | None = 0,
    keepdims: Literal[True],
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[FloatingT]: ...
@overload
def yeojohnson_llf(
    lmb: float | np.float64,
    data: onp.ToFloatStrict1D,
    *,
    axis: int | None = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> np.float64 | Any: ...
@overload
def yeojohnson_llf(
    lmb: float | np.float64,
    data: onp.ToFloatND,
    *,
    axis: None,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> np.float64 | Any: ...
@overload
def yeojohnson_llf(
    lmb: float | np.float64,
    data: onp.ToFloatND,
    *,
    axis: int | None = 0,
    keepdims: Literal[True],
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64 | Any]: ...
@overload
def yeojohnson_llf(
    lmb: float | np.float64,
    data: onp.ToFloatND,
    *,
    axis: int | None = 0,
    keepdims: Literal[False] = False,
    nan_policy: NanPolicy = "propagate",
) -> onp.ArrayND[np.float64 | Any] | Any: ...

#
@overload
def boxcox(
    x: onp.ToFloat1D,
    lmbda: None = None,
    alpha: None = None,
    optimizer: _MinFun1D | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
) -> tuple[_Float1D, np.float64]: ...
@overload
def boxcox(
    x: onp.ToFloat1D,
    lmbda: onp.ToFloat,
    alpha: float | None = None,
    optimizer: _MinFun1D | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
) -> _Float1D: ...
@overload
def boxcox(
    x: onp.ToFloat1D, lmbda: None, alpha: float, optimizer: _MinFun1D | None = None, *, nan_policy: NanPolicy = "propagate"
) -> tuple[_Float1D, np.float64, _Tuple2[float]]: ...
@overload
def boxcox(
    x: onp.ToFloat1D, lmbda: None = None, *, alpha: float, optimizer: _MinFun1D | None = None, nan_policy: NanPolicy = "propagate"
) -> tuple[_Float1D, np.float64, _Tuple2[float]]: ...

#
@overload
def yeojohnson(x: onp.ToFloat1D, lmbda: None = None, *, nan_policy: NanPolicy = "propagate") -> tuple[_Float1D, np.float64]: ...
@overload
def yeojohnson(x: onp.ToFloat1D, lmbda: onp.ToFloat, *, nan_policy: NanPolicy = "propagate") -> _Float1D: ...

#
@overload
def boxcox_normmax(
    x: onp.ToFloat1D,
    brack: _Tuple2[float] | None = None,
    method: Literal["pearsonr", "mle"] = "pearsonr",
    optimizer: _MinFun1D | None = None,
    *,
    ymax: onp.ToFloat | _BigFloat = ...,
    nan_policy: NanPolicy = "propagate",
) -> np.float64: ...
@overload
def boxcox_normmax(
    x: onp.ToFloat1D,
    brack: _Tuple2[float] | None,
    method: Literal["all"],
    optimizer: _MinFun1D | None = None,
    *,
    ymax: onp.ToFloat | _BigFloat = ...,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float64]: ...
@overload
def boxcox_normmax(
    x: onp.ToFloat1D,
    brack: _Tuple2[float] | None = None,
    *,
    method: Literal["all"],
    optimizer: _MinFun1D | None = None,
    ymax: onp.ToFloat | _BigFloat = ...,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float64]: ...

#
@overload
def yeojohnson_normmax(
    x: onp.ArrayND[npc.floating | npc.integer, _JustAnyShape],
    brack: _Tuple2[onp.ToFloat] | None = None,
    *,
    nan_policy: NanPolicy = "propagate",
) -> onp.Array1D[np.float64] | np.float64: ...
@overload
def yeojohnson_normmax(
    x: onp.ToFloatStrict1D, brack: _Tuple2[onp.ToFloat] | None = None, *, nan_policy: NanPolicy = "propagate"
) -> np.float64: ...
@overload
def yeojohnson_normmax(
    x: onp.ToFloatStrict2D, brack: _Tuple2[onp.ToFloat] | None = None, *, nan_policy: NanPolicy = "propagate"
) -> onp.Array1D[np.float64]: ...
@overload
def yeojohnson_normmax(
    x: onp.ToFloatND, brack: _Tuple2[onp.ToFloat] | None = None, *, nan_policy: NanPolicy = "propagate"
) -> onp.Array1D[np.float64] | np.float64: ...

#
def boxcox_normplot(
    x: onp.ToFloat1D, la: float, lb: float, plot: _CanPlotText | ModuleType | None = None, N: int = 80
) -> _Tuple2[onp.Array1D[np.float64]]: ...

#
def yeojohnson_normplot(
    x: onp.ToFloat1D, la: float, lb: float, plot: _CanPlotText | ModuleType | None = None, N: int = 80
) -> _Tuple2[onp.Array1D[np.float64]]: ...

#
@overload
@deprecated(
    "As of SciPy 1.17, users must choose a p-value calculation method by providing the `method` parameter. "
    "`method='interpolate'` interpolates the p-value from pre-calculated tables; `method` may also be an instance of "
    "`MonteCarloMethod` to approximate the p-value via Monte Carlo simulation. "
    "When `method` is specified, the result object will include a `pvalue` attribute and not attributes `critical_value`, "
    "`significance_level`, or `fit_result`. "
    "Beginning in 1.19.0, these other attributes will no longer be available, "
    "and a p-value will always be computed according to one of the available `method` options.",
    category=FutureWarning,
)
def anderson(x: onp.ToFloatND, dist: _RVCAnderson = "norm", *, method: None = None) -> AndersonResult: ...
@overload
def anderson(
    x: onp.ToFloatND, dist: _RVCAnderson = "norm", *, method: MonteCarloMethod | Literal["interpolated"]
) -> AndersonResult: ...

#
@overload
@deprecated(
    "Parameter `variant` has been introduced to replace `midrank`; "
    "`midrank` will be removed in SciPy 1.19.0. Specify `variant` to silence this warning. "
    "Note that the returned object will no longer be unpackable as a tuple, and `critical_values` will be omitted."
)
def anderson_ksamp(
    samples: onp.ToFloatND, midrank: bool, *, variant: op.JustObject = ..., method: PermutationMethod | None = None
) -> Anderson_ksampResult: ...
@overload
def anderson_ksamp(
    samples: onp.ToFloatND,
    midrank: op.JustObject = ...,
    *,
    variant: Literal["midrank", "right", "continuous"] | op.JustObject = ...,
    method: PermutationMethod | None = None,
) -> Anderson_ksampResult: ...

#
@overload
def shapiro(
    x: onp.ToFloat | onp.ToFloatND, *, axis: None = None, nan_policy: NanPolicy = "propagate", keepdims: Literal[False] = False
) -> ShapiroResult[np.float64]: ...
@overload
def shapiro(
    x: onp.ToFloat | onp.ToFloatND,
    *,
    axis: SupportsIndex | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> ShapiroResult[onp.ArrayND[np.float64]]: ...
@overload
def shapiro(
    x: onp.ToFloat | onp.ToFloatND,
    *,
    axis: SupportsIndex | None = None,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> ShapiroResult: ...

#
@overload
def ansari(
    x: onp.ToFloat | onp.ToFloatND,
    y: onp.ToFloat | onp.ToFloatND,
    alternative: Alternative = "two-sided",
    *,
    axis: None,
    method: _AnsariMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> AnsariResult[np.float64]: ...
@overload
def ansari(
    x: onp.ToFloat | onp.ToFloatND,
    y: onp.ToFloat | onp.ToFloatND,
    alternative: Alternative = "two-sided",
    *,
    axis: SupportsIndex | None = 0,
    method: _AnsariMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> AnsariResult[onp.ArrayND[np.float64]]: ...
@overload
def ansari(
    x: onp.ToFloat | onp.ToFloatND,
    y: onp.ToFloat | onp.ToFloatND,
    alternative: Alternative = "two-sided",
    *,
    axis: SupportsIndex | None = 0,
    method: _AnsariMethod = "auto",
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> AnsariResult: ...

#
@overload
def bartlett(
    *samples: onp.ToFloatND, axis: None, nan_policy: NanPolicy = "propagate", keepdims: Literal[False] = False
) -> BartlettResult[np.float64]: ...
@overload
def bartlett(
    *samples: onp.ToFloatND, axis: SupportsIndex | None = 0, nan_policy: NanPolicy = "propagate", keepdims: Literal[True]
) -> BartlettResult[onp.ArrayND[np.float64]]: ...
@overload
def bartlett(
    *samples: onp.ToFloatND, axis: SupportsIndex | None = 0, nan_policy: NanPolicy = "propagate", keepdims: bool = False
) -> BartlettResult: ...

#
@overload
def levene(
    *samples: onp.ToFloatND,
    center: _CenterMethod = "median",
    proportiontocut: onp.ToFloat = 0.05,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> LeveneResult[np.float64]: ...
@overload
def levene(
    *samples: onp.ToFloatND,
    center: _CenterMethod = "median",
    proportiontocut: onp.ToFloat = 0.05,
    axis: SupportsIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> LeveneResult[onp.ArrayND[np.float64]]: ...
@overload
def levene(
    *samples: onp.ToFloatND,
    center: _CenterMethod = "median",
    proportiontocut: onp.ToFloat = 0.05,
    axis: SupportsIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> LeveneResult: ...

#
@overload
def fligner(
    *samples: onp.ToFloatND,
    center: _CenterMethod = "median",
    proportiontocut: onp.ToFloat = 0.05,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> FlignerResult[np.float64]: ...
@overload
def fligner(
    *samples: onp.ToFloatND,
    center: _CenterMethod = "median",
    proportiontocut: onp.ToFloat = 0.05,
    axis: SupportsIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> FlignerResult[onp.ArrayND[np.float64]]: ...
@overload
def fligner(
    *samples: onp.ToFloatND,
    center: _CenterMethod = "median",
    proportiontocut: onp.ToFloat = 0.05,
    axis: SupportsIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> FlignerResult: ...

#
@overload
def mood(
    x: onp.ToFloat | onp.ToFloatND,
    y: onp.ToFloat | onp.ToFloatND,
    axis: None,
    alternative: Alternative = "two-sided",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> SignificanceResult[np.float64]: ...
@overload
def mood(
    x: onp.ToFloat | onp.ToFloatND,
    y: onp.ToFloat | onp.ToFloatND,
    axis: SupportsIndex | None = 0,
    alternative: Alternative = "two-sided",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> SignificanceResult[onp.ArrayND[np.float64]]: ...
@overload
def mood(
    x: onp.ToFloat | onp.ToFloatND,
    y: onp.ToFloat | onp.ToFloatND,
    axis: SupportsIndex | None = 0,
    alternative: Alternative = "two-sided",
    *,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> SignificanceResult[np.float64 | onp.ArrayND[np.float64]]: ...

#
@overload
def wilcoxon(
    x: onp.ToFloat | onp.ToFloatND,
    y: onp.ToFloat | onp.ToFloatND | None = None,
    zero_method: Literal["wilcox", "pratt", "zsplit"] = "wilcox",
    correction: bool = False,
    alternative: Alternative = "two-sided",
    method: Literal["auto", "exact", "approx"] | PermutationMethod = "auto",
    *,
    axis: None,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[False] = False,
) -> WilcoxonResult[np.float64]: ...
@overload
def wilcoxon(
    x: onp.ToFloat | onp.ToFloatND,
    y: onp.ToFloat | onp.ToFloatND | None = None,
    zero_method: Literal["wilcox", "pratt", "zsplit"] = "wilcox",
    correction: bool = False,
    alternative: Alternative = "two-sided",
    method: Literal["auto", "exact", "approx", "asymptotic"] | PermutationMethod = "auto",
    *,
    axis: SupportsIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: Literal[True],
) -> WilcoxonResult[onp.ArrayND[np.float64]]: ...
@overload
def wilcoxon(
    x: onp.ToFloat | onp.ToFloatND,
    y: onp.ToFloat | onp.ToFloatND | None = None,
    zero_method: Literal["wilcox", "pratt", "zsplit"] = "wilcox",
    correction: bool = False,
    alternative: Alternative = "two-sided",
    method: Literal["auto", "exact", "approx"] | PermutationMethod = "auto",
    *,
    axis: SupportsIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    keepdims: bool = False,
) -> WilcoxonResult: ...

#
def wilcoxon_result_object(
    statistic: np.float64, pvalue: np.float64, zstatistic: np.float64 | None = None
) -> WilcoxonResult: ...  # undocumented
def wilcoxon_result_unpacker(res: WilcoxonResult, _: int) -> _Tuple2[np.float64] | _Tuple3[np.float64]: ...  # undocumented
def wilcoxon_outputs(kwds: dict[str, str]) -> Literal[2, 3]: ...  # undocumented

#
def median_test(
    *samples: onp.ToFloatND,
    ties: Literal["below", "above", "ignore"] = "below",
    correction: bool = True,
    lambda_: onp.ToFloat | str = 1,
    nan_policy: NanPolicy = "propagate",
) -> MedianTestResult: ...

#
@overload
def circmean(
    samples: onp.ToFloatND,
    high: onp.ToFloat = 6.283_185_307_179_586,  # 2 * pi
    low: onp.ToFloat = 0,
    axis: None = None,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload
def circmean(
    samples: onp.ToFloatND,
    high: onp.ToFloat = 6.283_185_307_179_586,  # 2 * pi
    low: onp.ToFloat = 0,
    axis: SupportsIndex | None = None,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: Literal[True],
) -> onp.ArrayND[np.float64]: ...
@overload
def circmean(
    samples: onp.ToFloatND,
    high: onp.ToFloat = 6.283_185_307_179_586,  # 2 * pi
    low: onp.ToFloat = 0,
    axis: SupportsIndex | None = None,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: bool = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...

#
@overload
def circvar(
    samples: onp.ToFloatND,
    high: onp.ToFloat = 6.283_185_307_179_586,  # 2 * pi
    low: onp.ToFloat = 0,
    axis: None = None,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload
def circvar(
    samples: onp.ToFloatND,
    high: onp.ToFloat = 6.283_185_307_179_586,  # 2 * pi
    low: onp.ToFloat = 0,
    axis: SupportsIndex | None = None,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: Literal[True],
) -> onp.ArrayND[np.float64]: ...
@overload
def circvar(
    samples: onp.ToFloatND,
    high: onp.ToFloat = 6.283_185_307_179_586,  # 2 * pi
    low: onp.ToFloat = 0,
    axis: SupportsIndex | None = None,
    nan_policy: NanPolicy = "propagate",
    *,
    keepdims: bool = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...

#
@overload
def circstd(
    samples: onp.ToFloatND,
    high: onp.ToFloat = 6.283_185_307_179_586,  # 2 * pi
    low: onp.ToFloat = 0,
    axis: None = None,
    nan_policy: NanPolicy = "propagate",
    *,
    normalize: bool = False,
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload
def circstd(
    samples: onp.ToFloatND,
    high: onp.ToFloat = 6.283_185_307_179_586,  # 2 * pi
    low: onp.ToFloat = 0,
    axis: SupportsIndex | None = None,
    nan_policy: NanPolicy = "propagate",
    *,
    normalize: bool = False,
    keepdims: Literal[True],
) -> onp.ArrayND[np.float64]: ...
@overload
def circstd(
    samples: onp.ToFloatND,
    high: onp.ToFloat = 6.283_185_307_179_586,  # 2 * pi
    low: onp.ToFloat = 0,
    axis: SupportsIndex | None = None,
    nan_policy: NanPolicy = "propagate",
    *,
    normalize: bool = False,
    keepdims: bool = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...

#
@overload  # ?d +T@floating
def directional_stats[ScalarT: npc.floating](
    samples: onp.ArrayND[ScalarT, _JustAnyShape], *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.ArrayND[ScalarT], ScalarT | onp.ArrayND[ScalarT]]: ...
@overload  # ?d +integer
def directional_stats(
    samples: onp.ArrayND[npc.integer | np.bool, _JustAnyShape], *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.ArrayND[np.float64], np.float64 | onp.ArrayND[np.float64]]: ...
@overload  # ?d ~c128
def directional_stats(
    samples: onp.ArrayND[npc.complexfloating128, _JustAnyShape], *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.ArrayND[np.complex128], np.float64 | onp.ArrayND[np.float64]]: ...
@overload  # ?d ~c64
def directional_stats(
    samples: onp.ArrayND[npc.complexfloating64, _JustAnyShape], *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.ArrayND[np.complex64], np.float32 | onp.ArrayND[np.float32]]: ...
@overload  # ?d ~c160
def directional_stats(
    samples: onp.ArrayND[npc.complexfloating160, _JustAnyShape], *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.ArrayND[np.clongdouble], np.longdouble | onp.ArrayND[np.longdouble]]: ...
@overload  # 2d +T@floating
def directional_stats[ScalarT: npc.floating](
    samples: onp.Array2D[ScalarT], *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.Array1D[ScalarT], ScalarT]: ...
@overload  # 2d +integer | +float
def directional_stats(
    samples: onp.ToArrayStrict2D[float, npc.integer | np.bool], *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.Array1D[np.float64], np.float64]: ...
@overload  # 2d ~c128
def directional_stats(
    samples: onp.ToJustComplex128Strict2D, *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.Array1D[np.complex128], np.float64]: ...
@overload  # 2d ~c64
def directional_stats(
    samples: onp.ToJustComplex64Strict2D, *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.Array1D[np.complex64], np.float32]: ...
@overload  # 2d ~c160
def directional_stats(
    samples: onp.ToJustCLongDoubleStrict2D, *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.Array1D[np.clongdouble], np.longdouble]: ...
@overload  # 3d +T@floating
def directional_stats[ScalarT: npc.floating](
    samples: onp.Array3D[ScalarT], *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.Array2D[ScalarT], onp.Array1D[ScalarT]]: ...
@overload  # 3d +integer | +float
def directional_stats(
    samples: onp.ToArrayStrict3D[float, npc.integer | np.bool], *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.Array2D[np.float64], onp.Array1D[np.float64]]: ...
@overload  # 3d ~c128
def directional_stats(
    samples: onp.ToJustComplex128Strict3D, *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.Array2D[np.complex128], onp.Array1D[np.float64]]: ...
@overload  # 3d ~c64
def directional_stats(
    samples: onp.ToJustComplex64Strict3D, *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.Array2D[np.complex64], onp.Array1D[np.float32]]: ...
@overload  # 3d ~c160
def directional_stats(
    samples: onp.ToJustCLongDoubleStrict3D, *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.Array2D[np.clongdouble], onp.Array1D[np.longdouble]]: ...
@overload  # Nd +T@floating
def directional_stats[ScalarT: npc.floating](
    samples: onp.ArrayND[ScalarT], *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.ArrayND[ScalarT], ScalarT | onp.ArrayND[ScalarT]]: ...
@overload  # Nd +integer | +float
def directional_stats(
    samples: onp.ToArrayND[float, npc.integer | np.bool], *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.ArrayND[np.float64], np.float64 | onp.ArrayND[np.float64]]: ...
@overload  # Nd ~c128
def directional_stats(
    samples: onp.ToJustComplex128_ND, *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.ArrayND[np.complex128], np.float64 | onp.ArrayND[np.float64]]: ...
@overload  # Nd ~c64
def directional_stats(
    samples: onp.ToJustComplex64_ND, *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.ArrayND[np.complex64], np.float32 | onp.ArrayND[np.float32]]: ...
@overload  # Nd ~c160
def directional_stats(
    samples: onp.ToJustCLongDoubleND, *, axis: SupportsIndex = 0, normalize: bool = True
) -> DirectionalStats[onp.ArrayND[np.clongdouble], np.longdouble | onp.ArrayND[np.longdouble]]: ...

#
@overload
def false_discovery_control(
    ps: onp.ToFloat, *, axis: SupportsIndex | None = 0, method: Literal["bh", "by"] = "bh"
) -> np.float64: ...
@overload
def false_discovery_control(
    ps: onp.ToFloatND, *, axis: SupportsIndex | None = 0, method: Literal["bh", "by"] = "bh"
) -> onp.ArrayND[np.float64]: ...
