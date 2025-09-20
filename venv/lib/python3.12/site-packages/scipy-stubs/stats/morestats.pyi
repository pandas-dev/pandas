# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from . import distributions

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
    "chi2_contingency",
    "circmean",
    "circstd",
    "circvar",
    "distributions",
    "find_repeats",
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

# mstats_basic
@deprecated("will be removed in SciPy v2.0.0")
def find_repeats(arr: object) -> object: ...

# contingency
@deprecated("will be removed in SciPy v2.0.0")
def chi2_contingency(observed: object, correction: object = ..., lambda_: object = ..., *, method: object = ...) -> object: ...

# morestats
@deprecated("will be removed in SciPy v2.0.0")
def bayes_mvs(data: object, alpha: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def mvsdist(data: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kstat(data: object, n: object = ..., *, axis: object = ..., nan_policy: object = ..., keepdims: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kstatvar(
    data: object, n: object = ..., *, axis: object = ..., nan_policy: object = ..., keepdims: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def probplot(
    x: object, sparams: object = ..., dist: object = ..., fit: object = ..., plot: object = ..., rvalue: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ppcc_max(x: object, brack: object = ..., dist: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ppcc_plot(x: object, a: object, b: object, dist: object = ..., plot: object = ..., N: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def boxcox_llf(lmb: object, data: object, *, axis: object = ..., keepdims: object = ..., nan_policy: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def boxcox(x: object, lmbda: object = ..., alpha: object = ..., optimizer: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def boxcox_normmax(
    x: object, brack: object = ..., method: object = ..., optimizer: object = ..., *, ymax: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def boxcox_normplot(x: object, la: object, lb: object, plot: object = ..., N: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def yeojohnson(x: object, lmbda: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def yeojohnson_llf(lmb: object, data: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def yeojohnson_normmax(x: object, brack: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def yeojohnson_normplot(x: object, la: object, lb: object, plot: object = ..., N: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def shapiro(x: object, *, axis: object = ..., nan_policy: object = ..., keepdims: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def anderson(x: object, dist: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def anderson_ksamp(samples: object, midrank: object = ..., *, method: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ansari(
    x: object, y: object, alternative: object = ..., *, axis: object = ..., nan_policy: object = ..., keepdims: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bartlett(*samples: object, axis: object = ..., nan_policy: object = ..., keepdims: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def levene(
    *samples: object,
    center: object = ...,
    proportiontocut: object = ...,
    axis: object = ...,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fligner(
    *samples: object,
    center: object = ...,
    proportiontocut: object = ...,
    axis: object = ...,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def mood(
    x: object, y: object, axis: object = ..., alternative: object = ..., *, nan_policy: object = ..., keepdims: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def wilcoxon(
    x: object,
    y: object = ...,
    zero_method: object = ...,
    correction: object = ...,
    alternative: object = ...,
    method: object = ...,
    *,
    axis: object = ...,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def median_test(
    *samples: object, ties: object = ..., correction: object = ..., lambda_: object = ..., nan_policy: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def circmean(
    samples: object,
    high: object = ...,
    low: object = ...,
    axis: object = ...,
    nan_policy: object = ...,
    *,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def circvar(
    samples: object,
    high: object = ...,
    low: object = ...,
    axis: object = ...,
    nan_policy: object = ...,
    *,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def circstd(
    samples: object,
    high: object = ...,
    low: object = ...,
    axis: object = ...,
    nan_policy: object = ...,
    *,
    normalize: object = ...,
    keepdims: object = ...,
) -> object: ...
