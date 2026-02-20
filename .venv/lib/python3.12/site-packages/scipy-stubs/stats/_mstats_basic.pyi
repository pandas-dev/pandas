from collections.abc import Callable
from typing import Concatenate, Final, Generic, Literal, NamedTuple, Self, TypeAlias, TypedDict, overload, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc
from numpy._typing import _ArrayLike

from ._stats_mstats_common import SiegelslopesResult, TheilslopesResult
from ._stats_py import KstestResult, LinregressResult, SignificanceResult
from ._typing import Alternative, BaseBunch, NanPolicy

__all__ = [
    "argstoarray",
    "brunnermunzel",
    "count_tied_groups",
    "describe",
    "f_oneway",
    "find_repeats",
    "friedmanchisquare",
    "kendalltau",
    "kendalltau_seasonal",
    "kruskal",
    "kruskalwallis",
    "ks_1samp",
    "ks_2samp",
    "ks_twosamp",
    "kstest",
    "kurtosis",
    "kurtosistest",
    "linregress",
    "mannwhitneyu",
    "meppf",
    "mode",
    "moment",
    "mquantiles",
    "msign",
    "normaltest",
    "obrientransform",
    "pearsonr",
    "plotting_positions",
    "pointbiserialr",
    "rankdata",
    "scoreatpercentile",
    "sem",
    "sen_seasonal_slopes",
    "siegelslopes",
    "skew",
    "skewtest",
    "spearmanr",
    "theilslopes",
    "tmax",
    "tmean",
    "tmin",
    "trim",
    "trima",
    "trimboth",
    "trimmed_mean",
    "trimmed_std",
    "trimmed_stde",
    "trimmed_var",
    "trimr",
    "trimtail",
    "tsem",
    "ttest_1samp",
    "ttest_ind",
    "ttest_onesamp",
    "ttest_rel",
    "tvar",
    "variation",
    "winsorize",
]

_SCT = TypeVar("_SCT", bound=np.generic, default=np.float64)
_SCT_f = TypeVar("_SCT_f", bound=npc.floating, default=np.float64)
_SCT_bifc = TypeVar("_SCT_bifc", bound=npc.number | np.bool_, default=np.float64)
_SCT_bifcmO = TypeVar("_SCT_bifcmO", bound=npc.number | np.timedelta64 | np.bool_ | np.object_)

_MArrayOrND: TypeAlias = _SCT | onp.MArray[_SCT]

_NDT_f_co = TypeVar("_NDT_f_co", covariant=True, bound=float | _MArrayOrND[npc.floating], default=onp.MArray[np.float64])
_NDT_fc_co = TypeVar(
    "_NDT_fc_co",
    covariant=True,
    bound=complex | _MArrayOrND[npc.inexact],
    default=_MArrayOrND[np.float64 | np.complex128],
)  # fmt: skip

_KendallTauMethod: TypeAlias = Literal["auto", "asymptotic", "exact"]
_TheilSlopesMethod: TypeAlias = Literal["joint", "separate"]
_SiegelSlopesMethod: TypeAlias = Literal["hierarchical", "separate"]

_KSMethod: TypeAlias = Literal["auto", "exact", "asymp"]
_KTestMethod: TypeAlias = Literal[_KSMethod, "approx"]

@type_check_only
class _TestResult(NamedTuple, Generic[_NDT_f_co, _NDT_fc_co]):
    statistic: _NDT_fc_co
    pvalue: _NDT_f_co

_KendallTauSeasonalResult = TypedDict(
    "_KendallTauSeasonalResult",
    {
        "seasonal tau": _MArrayOrND[np.float64],
        "global tau": np.float64,
        "global tau (alt)": np.float64,
        "seasonal p-value": onp.ArrayND[np.float64],
        "global p-value (indep)": np.float64,
        "global p-value (dep)": np.float64,
        "chi2 total": onp.MArray[np.float64],
        "chi2 trend": onp.MArray[np.float64],
    },
)

###

trimdoc: Final[str] = ...

class ModeResult(NamedTuple):
    mode: onp.MArray[np.float64]
    count: onp.MArray[np.float64]  # type: ignore[assignment]  # pyright: ignore[reportIncompatibleMethodOverride]

class DescribeResult(NamedTuple):
    nobs: np.int_ | onp.ArrayND[np.int_]
    minmax: tuple[onp.MArray[npc.floating | npc.integer], onp.MArray[npc.floating | npc.integer]]
    mean: npc.floating
    variance: npc.floating
    skewness: npc.floating
    kurtosis: npc.floating

class PointbiserialrResult(NamedTuple):
    correlation: np.float64
    pvalue: np.float64

class Ttest_relResult(_TestResult[_NDT_f_co, _NDT_fc_co], Generic[_NDT_f_co, _NDT_fc_co]): ...
class Ttest_indResult(_TestResult[_NDT_f_co, _NDT_fc_co], Generic[_NDT_f_co, _NDT_fc_co]): ...
class Ttest_1sampResult(_TestResult[_NDT_f_co, _NDT_fc_co], Generic[_NDT_f_co, _NDT_fc_co]): ...
class SkewtestResult(_TestResult[_NDT_f_co, _NDT_fc_co], Generic[_NDT_f_co, _NDT_fc_co]): ...
class KurtosistestResult(_TestResult[_NDT_f_co, _NDT_fc_co], Generic[_NDT_f_co, _NDT_fc_co]): ...
class NormaltestResult(_TestResult[_NDT_f_co, _NDT_f_co], Generic[_NDT_f_co]): ...
class MannwhitneyuResult(_TestResult[np.float64, np.float64]): ...
class F_onewayResult(_TestResult[np.float64, np.float64]): ...
class KruskalResult(_TestResult[np.float64, np.float64]): ...
class FriedmanchisquareResult(_TestResult[np.float64, np.float64]): ...
class BrunnerMunzelResult(_TestResult[np.float64, np.float64]): ...

class SenSeasonalSlopesResult(BaseBunch[onp.MArray[np.float64], np.float64]):
    def __new__(_cls, intra_slope: float, inter_slope: float) -> Self: ...
    def __init__(self, /, intra_slope: float, inter_slope: float) -> None: ...
    @property
    def intra_slope(self, /) -> onp.MArray[np.float64]: ...
    @property
    def inter_slope(self, /) -> float: ...

# TODO(jorenham): Overloads for scalar vs. array
# TODO(jorenham): Overloads for specific dtypes

def argstoarray(*args: onp.ToFloatND) -> onp.MArray[np.float64]: ...
def find_repeats(arr: onp.ToFloatND) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.intp]]: ...
def count_tied_groups(x: onp.ToFloatND, use_missing: bool = False) -> dict[np.intp, np.intp]: ...
def rankdata(data: onp.ToFloatND, axis: op.CanIndex | None = None, use_missing: bool = False) -> onp.ArrayND[np.float64]: ...
def mode(a: onp.ToFloatND, axis: op.CanIndex | None = 0) -> ModeResult: ...

#
@overload
def msign(x: _ArrayLike[_SCT_bifcmO]) -> onp.ArrayND[_SCT_bifcmO]: ...
@overload
def msign(x: onp.ToComplexND) -> onp.ArrayND[npc.number | np.timedelta64 | np.bool_ | np.object_]: ...

#
def pearsonr(x: onp.ToFloatND, y: onp.ToFloatND) -> tuple[np.float64, np.float64]: ...
def spearmanr(
    x: onp.ToFloatND,
    y: onp.ToFloatND | None = None,
    use_ties: bool = True,
    axis: op.CanIndex | None = None,
    nan_policy: NanPolicy = "propagate",
    alternative: Alternative = "two-sided",
) -> SignificanceResult: ...
def kendalltau(
    x: onp.ToFloatND,
    y: onp.ToFloatND,
    use_ties: bool = True,
    use_missing: bool = False,
    method: _KendallTauMethod = "auto",
    alternative: Alternative = "two-sided",
) -> SignificanceResult: ...
def kendalltau_seasonal(x: onp.ToFloatND) -> _KendallTauSeasonalResult: ...
def pointbiserialr(x: onp.ToFloatND, y: onp.ToFloatND) -> PointbiserialrResult: ...
def linregress(x: onp.ToFloatND, y: onp.ToFloatND | None = None) -> LinregressResult: ...
def theilslopes(
    y: onp.ToFloatND, x: onp.ToFloatND | None = None, alpha: onp.ToJustFloat = 0.95, method: _TheilSlopesMethod = "separate"
) -> TheilslopesResult: ...
def siegelslopes(
    y: onp.ToFloatND, x: onp.ToFloatND | None = None, method: _SiegelSlopesMethod = "hierarchical"
) -> SiegelslopesResult: ...
def sen_seasonal_slopes(x: onp.ToFloatND) -> SenSeasonalSlopesResult: ...

#
def ttest_1samp(
    a: onp.ToFloatND, popmean: onp.ToFloat | onp.ToFloatND, axis: op.CanIndex | None = 0, alternative: Alternative = "two-sided"
) -> Ttest_1sampResult: ...
def ttest_ind(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    axis: op.CanIndex | None = 0,
    equal_var: op.CanBool = True,
    alternative: Alternative = "two-sided",
) -> Ttest_indResult: ...
def ttest_rel(
    a: onp.ToFloatND, b: onp.ToFloatND, axis: op.CanIndex | None = 0, alternative: Alternative = "two-sided"
) -> Ttest_relResult: ...
def mannwhitneyu(x: onp.ToFloatND, y: onp.ToFloatND, use_continuity: op.CanBool = True) -> MannwhitneyuResult: ...
def kruskal(arg0: onp.ToFloatND, arg1: onp.ToFloatND, /, *args: onp.ToFloatND) -> KruskalResult: ...

#
@overload
def ks_1samp(
    x: onp.ToFloatND,
    cdf: str | Callable[[float], onp.ToFloat],
    args: tuple[()] = (),
    alternative: Alternative = "two-sided",
    method: _KSMethod = "auto",
) -> KstestResult: ...
@overload
def ks_1samp(
    x: onp.ToFloatND,
    cdf: str | Callable[Concatenate[float, ...], onp.ToFloat],
    args: tuple[object, ...],
    alternative: Alternative = "two-sided",
    method: _KSMethod = "auto",
) -> KstestResult: ...

#
def ks_2samp(
    data1: onp.ToFloatND, data2: onp.ToFloatND, alternative: Alternative = "two-sided", method: _KSMethod = "auto"
) -> KstestResult: ...

#
@overload
def kstest(
    data1: onp.ToFloatND,
    data2: onp.ToFloatND | str | Callable[[float], onp.ToFloat],
    args: tuple[()] = (),
    alternative: Alternative = "two-sided",
    method: _KTestMethod = "auto",
) -> KstestResult: ...
@overload
def kstest(
    data1: onp.ToFloatND,
    data2: Callable[Concatenate[float, ...], onp.ToFloat],
    args: tuple[object, ...],
    alternative: Alternative = "two-sided",
    method: _KTestMethod = "auto",
) -> KstestResult: ...

#
@overload
def trima(
    a: onp.SequenceND[bool],
    limits: tuple[onp.ToInt, onp.ToInt] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
) -> onp.MArray[np.bool_]: ...
@overload
def trima(
    a: onp.SequenceND[op.JustInt],
    limits: tuple[onp.ToInt, onp.ToInt] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
) -> onp.MArray[np.int_]: ...
@overload
def trima(
    a: onp.SequenceND[float],
    limits: tuple[onp.ToFloat, onp.ToFloat] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
) -> onp.MArray[np.float64 | np.int_ | np.bool_]: ...
@overload
def trima(
    a: onp.SequenceND[complex],
    limits: tuple[onp.ToComplex, onp.ToComplex] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
) -> onp.MArray[np.complex128 | np.float64 | np.int_ | np.bool_]: ...
@overload
def trima(
    a: _ArrayLike[_SCT_bifc],
    limits: tuple[onp.ToComplex, onp.ToComplex] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
) -> onp.MArray[_SCT_bifc]: ...

#
@overload
def trimr(
    a: onp.SequenceND[op.JustInt | np.int_],
    limits: tuple[onp.ToFloat, onp.ToFloat] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = None,
) -> onp.MArray[np.int_]: ...
@overload
def trimr(
    a: onp.SequenceND[float],
    limits: tuple[onp.ToFloat, onp.ToFloat] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = None,
) -> onp.MArray[np.float64 | np.int_]: ...
@overload
def trimr(
    a: onp.SequenceND[complex],
    limits: tuple[onp.ToComplex, onp.ToComplex] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = None,
) -> onp.MArray[np.complex128 | np.float64 | np.int_]: ...
@overload
def trimr(
    a: _ArrayLike[_SCT_bifc],
    limits: tuple[onp.ToComplex, onp.ToComplex] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = None,
) -> onp.MArray[_SCT_bifc]: ...

#
@overload
def trim(
    a: onp.SequenceND[op.JustInt | np.int_],
    limits: tuple[onp.ToFloat, onp.ToFloat] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    relative: op.CanBool = False,
    axis: op.CanIndex | None = None,
) -> onp.MArray[np.int_]: ...
@overload
def trim(
    a: onp.SequenceND[float],
    limits: tuple[onp.ToFloat, onp.ToFloat] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    relative: op.CanBool = False,
    axis: op.CanIndex | None = None,
) -> onp.MArray[np.float64 | np.int_]: ...
@overload
def trim(
    a: onp.SequenceND[complex],
    limits: tuple[onp.ToComplex, onp.ToComplex] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    relative: op.CanBool = False,
    axis: op.CanIndex | None = None,
) -> onp.MArray[np.complex128 | np.float64 | np.int_]: ...
@overload
def trim(
    a: _ArrayLike[_SCT_bifc],
    limits: tuple[onp.ToComplex, onp.ToComplex] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    relative: op.CanBool = False,
    axis: op.CanIndex | None = None,
) -> onp.MArray[_SCT_bifc]: ...

#
@overload
def trimboth(
    data: onp.SequenceND[op.JustInt | np.int_],
    proportiontocut: float | npc.floating = 0.2,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = None,
) -> onp.MArray[np.int_]: ...
@overload
def trimboth(
    data: onp.SequenceND[float],
    proportiontocut: float | npc.floating = 0.2,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = None,
) -> onp.MArray[np.float64 | np.int_]: ...
@overload
def trimboth(
    data: onp.SequenceND[complex],
    proportiontocut: float | npc.floating = 0.2,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = None,
) -> onp.MArray[np.complex128 | np.float64 | np.int_]: ...
@overload
def trimboth(
    data: _ArrayLike[_SCT_bifc],
    proportiontocut: float | npc.floating = 0.2,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = None,
) -> onp.MArray[_SCT_bifc]: ...

#
@overload
def trimtail(
    data: onp.SequenceND[op.JustInt | np.int_],
    proportiontocut: float | npc.floating = 0.2,
    tail: Literal["left", "right"] = "left",
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = None,
) -> onp.MArray[np.int_]: ...
@overload
def trimtail(
    data: onp.SequenceND[float],
    proportiontocut: float | npc.floating = 0.2,
    tail: Literal["left", "right"] = "left",
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = None,
) -> onp.MArray[np.float64 | np.int_]: ...
@overload
def trimtail(
    data: onp.SequenceND[complex],
    proportiontocut: float | npc.floating = 0.2,
    tail: Literal["left", "right"] = "left",
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = None,
) -> onp.MArray[np.complex128 | np.float64 | np.int_]: ...
@overload
def trimtail(
    data: _ArrayLike[_SCT_bifc],
    proportiontocut: float | npc.floating = 0.2,
    tail: Literal["left", "right"] = "left",
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = None,
) -> onp.MArray[_SCT_bifc]: ...

#
@overload
def trimmed_mean(
    a: onp.ToFloatND,
    limits: tuple[onp.ToFloat, onp.ToFloat] = (0.1, 0.1),
    inclusive: tuple[op.CanBool, op.CanBool] = (1, 1),
    relative: op.CanBool = True,
    axis: op.CanIndex | None = None,
) -> _MArrayOrND[npc.floating]: ...
@overload
def trimmed_mean(
    a: onp.ToComplexND,
    limits: tuple[onp.ToComplex, onp.ToComplex] = (0.1, 0.1),
    inclusive: tuple[op.CanBool, op.CanBool] = (1, 1),
    relative: op.CanBool = True,
    axis: op.CanIndex | None = None,
) -> _MArrayOrND[npc.floating | np.complex128]: ...

#
def trimmed_var(
    a: onp.ToComplexND,
    limits: tuple[onp.ToFloat, onp.ToFloat] = (0.1, 0.1),
    inclusive: tuple[op.CanBool, op.CanBool] = (1, 1),
    relative: op.CanBool = True,
    axis: op.CanIndex | None = None,
    ddof: onp.ToInt = 0,
) -> _MArrayOrND[np.float64]: ...

#
def trimmed_std(
    a: onp.ToComplexND,
    limits: tuple[onp.ToFloat, onp.ToFloat] = (0.1, 0.1),
    inclusive: tuple[op.CanBool, op.CanBool] = (1, 1),
    relative: op.CanBool = True,
    axis: op.CanIndex | None = None,
    ddof: onp.ToInt = 0,
) -> _MArrayOrND[np.float64]: ...

#
def trimmed_stde(
    a: onp.ToComplexND,
    limits: tuple[onp.ToFloat, onp.ToFloat] = (0.1, 0.1),
    inclusive: tuple[op.CanBool, op.CanBool] = (1, 1),
    axis: op.CanIndex | None = None,
) -> _MArrayOrND[np.float64]: ...

#
@overload
def tmean(
    a: onp.ToFloatND,
    limits: tuple[onp.ToFloat, onp.ToFloat] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = None,
) -> _MArrayOrND[npc.floating]: ...
@overload
def tmean(
    a: onp.ToComplexND,
    limits: tuple[onp.ToComplex, onp.ToComplex] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = None,
) -> _MArrayOrND[npc.inexact]: ...

#
def tvar(
    a: onp.MArray[npc.floating | npc.integer],
    limits: tuple[onp.ToFloat, onp.ToFloat] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = 0,
    ddof: onp.ToInt = 1,
) -> _MArrayOrND[npc.floating]: ...

#
@overload
def tmin(
    a: onp.SequenceND[op.JustInt | np.int_],
    lowerlimit: onp.ToFloat | None = None,
    axis: op.CanIndex | None = 0,
    inclusive: onp.ToBool = True,
) -> _MArrayOrND[np.int_]: ...
@overload
def tmin(
    a: onp.SequenceND[float], lowerlimit: onp.ToFloat | None = None, axis: op.CanIndex | None = 0, inclusive: onp.ToBool = True
) -> _MArrayOrND[np.float64 | np.int_]: ...
@overload
def tmin(
    a: onp.SequenceND[complex],
    lowerlimit: onp.ToComplex | None = None,
    axis: op.CanIndex | None = 0,
    inclusive: onp.ToBool = True,
) -> _MArrayOrND[np.complex128 | np.float64 | np.int_]: ...
@overload
def tmin(
    a: _ArrayLike[_SCT_bifc], lowerlimit: onp.ToComplex | None = None, axis: op.CanIndex | None = 0, inclusive: onp.ToBool = True
) -> _MArrayOrND[_SCT_bifc]: ...

#
@overload
def tmax(
    a: onp.SequenceND[op.JustInt | np.int_],
    upperlimit: onp.ToFloat | None = None,
    axis: op.CanIndex | None = 0,
    inclusive: onp.ToBool = True,
) -> _MArrayOrND[np.int_]: ...
@overload
def tmax(
    a: onp.SequenceND[float], upperlimit: onp.ToFloat | None = None, axis: op.CanIndex | None = 0, inclusive: onp.ToBool = True
) -> _MArrayOrND[np.float64 | np.int_]: ...
@overload
def tmax(
    a: onp.SequenceND[complex],
    upperlimit: onp.ToComplex | None = None,
    axis: op.CanIndex | None = 0,
    inclusive: onp.ToBool = True,
) -> _MArrayOrND[np.complex128 | np.float64 | np.int_]: ...
@overload
def tmax(
    a: _ArrayLike[_SCT_bifc], upperlimit: onp.ToComplex | None = None, axis: op.CanIndex | None = 0, inclusive: onp.ToBool = True
) -> _MArrayOrND[_SCT_bifc]: ...

#
def tsem(
    a: onp.ToComplexND,
    limits: tuple[onp.ToFloat, onp.ToFloat] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    axis: op.CanIndex | None = 0,
    ddof: onp.ToInt = 1,
) -> _MArrayOrND: ...

#
@overload
def winsorize(
    a: onp.ToIntND,
    limits: tuple[onp.ToFloat, onp.ToFloat] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    inplace: onp.ToBool = False,
    axis: op.CanIndex | None = None,
    nan_policy: NanPolicy = "propagate",
) -> onp.MArray[np.int_]: ...
@overload
def winsorize(
    a: _ArrayLike[_SCT_f],
    limits: tuple[onp.ToFloat, onp.ToFloat] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    inplace: onp.ToBool = False,
    axis: op.CanIndex | None = None,
    nan_policy: NanPolicy = "propagate",
) -> onp.MArray[_SCT_f]: ...
@overload
def winsorize(
    a: onp.ToFloatND,
    limits: tuple[onp.ToFloat, onp.ToFloat] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    inplace: onp.ToBool = False,
    axis: op.CanIndex | None = None,
    nan_policy: NanPolicy = "propagate",
) -> onp.MArray[npc.floating | np.int_]: ...
@overload
def winsorize(
    a: onp.ToComplexND,
    limits: tuple[onp.ToComplex, onp.ToComplex] | None = None,
    inclusive: tuple[op.CanBool, op.CanBool] = (True, True),
    inplace: onp.ToBool = False,
    axis: op.CanIndex | None = None,
    nan_policy: NanPolicy = "propagate",
) -> onp.MArray[np.complex128 | npc.floating | np.int_]: ...

# TODO(jorenham): Overloads for complex array-likes
def moment(a: onp.ToFloatND, moment: onp.ToInt | onp.ToIntND = 1, axis: op.CanIndex | None = 0) -> _MArrayOrND[npc.floating]: ...
def variation(a: onp.ToFloatND, axis: op.CanIndex | None = 0, ddof: onp.ToInt = 0) -> _MArrayOrND[npc.floating]: ...
def skew(a: onp.ToFloatND, axis: op.CanIndex | None = 0, bias: op.CanBool = True) -> _MArrayOrND[npc.floating]: ...
def kurtosis(
    a: onp.ToFloatND, axis: op.CanIndex | None = 0, fisher: op.CanBool = True, bias: op.CanBool = True
) -> _MArrayOrND[npc.floating]: ...
def describe(a: onp.ToFloatND, axis: op.CanIndex | None = 0, ddof: onp.ToInt = 0, bias: op.CanBool = True) -> DescribeResult: ...

#
@overload
def stde_median(data: onp.ToFloatND, axis: op.CanIndex | None = None) -> _MArrayOrND[npc.floating]: ...
@overload
def stde_median(data: onp.ToComplexND, axis: op.CanIndex | None = None) -> _MArrayOrND[npc.inexact]: ...

#
@overload
def skewtest(
    a: onp.ToFloatND, axis: op.CanIndex | None = 0, alternative: Alternative = "two-sided"
) -> SkewtestResult[_MArrayOrND[np.float64], _MArrayOrND[np.float64]]: ...
@overload
def skewtest(
    a: onp.ToComplexND, axis: op.CanIndex | None = 0, alternative: Alternative = "two-sided"
) -> SkewtestResult[_MArrayOrND[np.float64], _MArrayOrND[np.float64 | np.complex128]]: ...

#
@overload
def kurtosistest(
    a: onp.ToFloatND, axis: op.CanIndex | None = 0, alternative: Alternative = "two-sided"
) -> KurtosistestResult[_MArrayOrND[np.float64], _MArrayOrND[np.float64]]: ...
@overload
def kurtosistest(
    a: onp.ToComplexND, axis: op.CanIndex | None = 0, alternative: Alternative = "two-sided"
) -> KurtosistestResult[_MArrayOrND[np.float64], _MArrayOrND[np.float64 | np.complex128]]: ...

#
def normaltest(a: onp.ToFloatND, axis: op.CanIndex | None = 0) -> NormaltestResult[_MArrayOrND[np.float64]]: ...

#
def mquantiles(
    a: onp.ToFloatND,
    prob: onp.ToFloatND = (0.25, 0.5, 0.75),
    alphap: onp.ToFloat = 0.4,
    betap: onp.ToFloat = 0.4,
    axis: op.CanIndex | None = None,
    limit: tuple[onp.ToFloat, onp.ToFloat] | tuple[()] = (),
) -> onp.MArray[np.float64]: ...

#
def scoreatpercentile(
    data: onp.ToFloatND,
    per: onp.ToFloat,
    limit: tuple[onp.ToFloat, onp.ToFloat] | tuple[()] = (),
    alphap: onp.ToFloat = 0.4,
    betap: onp.ToFloat = 0.4,
) -> onp.MArray[np.float64]: ...

#
def plotting_positions(data: onp.ToFloatND, alpha: onp.ToFloat = 0.4, beta: onp.ToFloat = 0.4) -> onp.MArray[np.float64]: ...

#
def obrientransform(arg0: onp.ToFloatND, /, *args: onp.ToFloatND) -> onp.MArray[np.float64]: ...

#
def sem(a: onp.ToFloatND, axis: op.CanIndex | None = 0, ddof: onp.ToInt = 1) -> np.float64 | onp.MArray[np.float64]: ...

#
def f_oneway(arg0: onp.ToFloatND, arg1: onp.ToFloatND, /, *args: onp.ToFloatND) -> F_onewayResult: ...

#
def friedmanchisquare(arg0: onp.ToFloatND, *args: onp.ToFloatND) -> FriedmanchisquareResult: ...

#
def brunnermunzel(
    x: onp.ToFloatND, y: onp.ToFloatND, alternative: Alternative = "two-sided", distribution: Literal["t", "normal"] = "t"
) -> BrunnerMunzelResult: ...

#
ttest_onesamp = ttest_1samp
kruskalwallis = kruskal
ks_twosamp = ks_2samp
trim1 = trimtail
meppf = plotting_positions
