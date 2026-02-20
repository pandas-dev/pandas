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

# contingency
@deprecated("will be removed in SciPy v2.0.0")
def chi2_contingency(observed: object, correction: object = True, lambda_: object = None, *, method: object = None) -> object: ...

# morestats
@deprecated("will be removed in SciPy v2.0.0")
def bayes_mvs(data: object, alpha: object = 0.9) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def mvsdist(data: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kstat(
    data: object, n: object = 2, *, axis: object = None, nan_policy: object = "propagate", keepdims: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kstatvar(
    data: object, n: object = 2, *, axis: object = None, nan_policy: object = "propagate", keepdims: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def probplot(
    x: object, sparams: object = (), dist: object = "norm", fit: object = True, plot: object = None, rvalue: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ppcc_max(x: object, brack: object = (0.0, 1.0), dist: object = "tukeylambda") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ppcc_plot(x: object, a: object, b: object, dist: object = "tukeylambda", plot: object = None, N: object = 80) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def boxcox_llf(
    lmb: object, data: object, *, axis: object = 0, keepdims: object = False, nan_policy: object = "propagate"
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def boxcox(x: object, lmbda: object = None, alpha: object = None, optimizer: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def boxcox_normmax(
    x: object, brack: object = None, method: object = "pearsonr", optimizer: object = None, *, ymax: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def boxcox_normplot(x: object, la: object, lb: object, plot: object = None, N: object = 80) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def yeojohnson(x: object, lmbda: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def yeojohnson_llf(
    lmb: object, data: object, *, axis: object = 0, nan_policy: object = "propagate", keepdims: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def yeojohnson_normmax(x: object, brack: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def yeojohnson_normplot(x: object, la: object, lb: object, plot: object = None, N: object = 80) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def shapiro(x: object, *, axis: object = None, nan_policy: object = "propagate", keepdims: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def anderson(x: object, dist: object = "norm", *, method: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def anderson_ksamp(samples: object, midrank: object = ..., *, variant: object = ..., method: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ansari(
    x: object,
    y: object,
    alternative: object = "two-sided",
    *,
    axis: object = 0,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bartlett(*samples: object, axis: object = 0, nan_policy: object = "propagate", keepdims: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def levene(
    *samples: object,
    center: object = "median",
    proportiontocut: object = 0.05,
    axis: object = 0,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fligner(
    *samples: object,
    center: object = "median",
    proportiontocut: object = 0.05,
    axis: object = 0,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def mood(
    x: object,
    y: object,
    axis: object = 0,
    alternative: object = "two-sided",
    *,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def wilcoxon(
    x: object,
    y: object = None,
    zero_method: object = "wilcox",
    correction: object = False,
    alternative: object = "two-sided",
    method: object = "auto",
    *,
    axis: object = 0,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def median_test(
    *samples: object, ties: object = "below", correction: object = True, lambda_: object = 1, nan_policy: object = "propagate"
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def circmean(
    samples: object,
    high: object = ...,
    low: object = 0,
    axis: object = None,
    nan_policy: object = "propagate",
    *,
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def circvar(
    samples: object,
    high: object = ...,
    low: object = 0,
    axis: object = None,
    nan_policy: object = "propagate",
    *,
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def circstd(
    samples: object,
    high: object = ...,
    low: object = 0,
    axis: object = None,
    nan_policy: object = "propagate",
    *,
    normalize: object = False,
    keepdims: object = False,
) -> object: ...
