# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from . import distributions, mstats_basic

__all__ = [
    "alexandergovern",
    "brunnermunzel",
    "chisquare",
    "combine_pvalues",
    "cumfreq",
    "describe",
    "distributions",
    "energy_distance",
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
    "median_abs_deviation",
    "mode",
    "moment",
    "mstats_basic",
    "multiscale_graphcorr",
    "normaltest",
    "obrientransform",
    "pearsonr",
    "percentileofscore",
    "pmean",
    "pointbiserialr",
    "power_divergence",
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
    "weightedtau",
    "zmap",
    "zscore",
]

# mgc
@deprecated("will be removed in SciPy v2.0.0")
def multiscale_graphcorr(
    x: object,
    y: object,
    compute_distance: object = ...,
    reps: object = 1_000,
    workers: object = 1,
    is_twosamp: object = False,
    random_state: object = None,
) -> object: ...

# _stats
@deprecated("will be removed in SciPy v2.0.0")
def gmean(
    a: object,
    axis: object = 0,
    dtype: object = None,
    weights: object = None,
    *,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hmean(
    a: object,
    axis: object = 0,
    dtype: object = None,
    *,
    weights: object = None,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pmean(
    a: object,
    p: object,
    *,
    dtype: object = None,
    weights: object = None,
    axis: object = 0,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def mode(a: object, axis: object = 0, nan_policy: object = "propagate", keepdims: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tmean(
    a: object,
    limits: object = None,
    inclusive: object = (True, True),
    axis: object = None,
    *,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tvar(
    a: object,
    limits: object = None,
    inclusive: object = (True, True),
    axis: object = 0,
    ddof: object = 1,
    *,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tmin(
    a: object,
    lowerlimit: object = None,
    axis: object = 0,
    inclusive: object = True,
    nan_policy: object = "propagate",
    *,
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tmax(
    a: object,
    upperlimit: object = None,
    axis: object = 0,
    inclusive: object = True,
    nan_policy: object = "propagate",
    *,
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tstd(
    a: object,
    limits: object = None,
    inclusive: object = (True, True),
    axis: object = 0,
    ddof: object = 1,
    *,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tsem(
    a: object,
    limits: object = None,
    inclusive: object = (True, True),
    axis: object = 0,
    ddof: object = 1,
    *,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def moment(
    a: object,
    order: object = 1,
    axis: object = 0,
    nan_policy: object = "propagate",
    *,
    center: object = None,
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def skew(
    a: object, axis: object = 0, bias: object = True, nan_policy: object = "propagate", *, keepdims: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kurtosis(
    a: object,
    axis: object = 0,
    fisher: object = True,
    bias: object = True,
    nan_policy: object = "propagate",
    *,
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def describe(a: object, axis: object = 0, ddof: object = 1, bias: object = True, nan_policy: object = "propagate") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def skewtest(
    a: object, axis: object = 0, nan_policy: object = "propagate", alternative: object = "two-sided", *, keepdims: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kurtosistest(
    a: object, axis: object = 0, nan_policy: object = "propagate", alternative: object = "two-sided", *, keepdims: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def normaltest(a: object, axis: object = 0, nan_policy: object = "propagate", *, keepdims: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def jarque_bera(x: object, *, axis: object = None, nan_policy: object = "propagate", keepdims: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def scoreatpercentile(
    a: object, per: object, limit: object = (), interpolation_method: object = "fraction", axis: object = None
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def percentileofscore(a: object, score: object, kind: object = "rank", nan_policy: object = "propagate") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def cumfreq(a: object, numbins: object = 10, defaultreallimits: object = None, weights: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def relfreq(a: object, numbins: object = 10, defaultreallimits: object = None, weights: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def obrientransform(*samples: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sem(
    a: object, axis: object = 0, ddof: object = 1, nan_policy: object = "propagate", *, keepdims: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def zscore(a: object, axis: object = 0, ddof: object = 0, nan_policy: object = "propagate") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gzscore(a: object, *, axis: object = 0, ddof: object = 0, nan_policy: object = "propagate") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def zmap(scores: object, compare: object, axis: object = 0, ddof: object = 0, nan_policy: object = "propagate") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gstd(
    a: object, axis: object = 0, ddof: object = 1, *, keepdims: object = False, nan_policy: object = "propagate"
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def iqr(
    x: object,
    axis: object = None,
    rng: object = (25, 75),
    scale: object = 1.0,
    nan_policy: object = "propagate",
    interpolation: object = "linear",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def median_abs_deviation(
    x: object,
    axis: object = 0,
    center: object = ...,
    scale: object = 1.0,
    nan_policy: object = "propagate",
    *,
    keepdims: bool = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sigmaclip(a: object, low: object = 4.0, high: object = 4.0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trimboth(a: object, proportiontocut: object, axis: object = 0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trim1(a: object, proportiontocut: object, tail: object = "right", axis: object = 0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trim_mean(
    a: object, proportiontocut: object, axis: object = 0, *, nan_policy: object = "propagate", keepdims: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def f_oneway(
    *samples: object, axis: object = 0, equal_var: object = True, nan_policy: object = "propagate", keepdims: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def alexandergovern(*samples: object, nan_policy: object = "propagate", axis: object = 0, keepdims: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pearsonr(x: object, y: object, *, alternative: object = "two-sided", method: object = None, axis: object = 0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fisher_exact(table: object, alternative: object = None, *, method: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def spearmanr(
    a: object, b: object = None, axis: object = 0, nan_policy: object = "propagate", alternative: object = "two-sided"
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pointbiserialr(
    x: object, y: object, *, axis: object = 0, nan_policy: object = "propagate", keepdims: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kendalltau(
    x: object,
    y: object,
    *,
    method: object = "auto",
    variant: object = "b",
    alternative: object = "two-sided",
    axis: object = None,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def weightedtau(
    x: object,
    y: object,
    rank: object = True,
    weigher: object = None,
    additive: object = True,
    *,
    axis: object = None,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ttest_1samp(
    a: object,
    popmean: object,
    axis: object = 0,
    nan_policy: object = "propagate",
    alternative: object = "two-sided",
    *,
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ttest_ind_from_stats(
    mean1: object,
    std1: object,
    nobs1: object,
    mean2: object,
    std2: object,
    nobs2: object,
    equal_var: object = True,
    alternative: object = "two-sided",
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ttest_ind(
    a: object,
    b: object,
    *,
    axis: object = 0,
    equal_var: object = True,
    nan_policy: object = "propagate",
    alternative: object = "two-sided",
    trim: object = 0,
    method: object = None,
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ttest_rel(
    a: object,
    b: object,
    axis: object = 0,
    nan_policy: object = "propagate",
    alternative: object = "two-sided",
    *,
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def power_divergence(
    f_obs: object,
    f_exp: object = None,
    ddof: object = 0,
    axis: object = 0,
    lambda_: object = None,
    *,
    keepdims: object = False,
    nan_policy: object = "propagate",
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def chisquare(
    f_obs: object,
    f_exp: object = None,
    ddof: object = 0,
    axis: object = 0,
    *,
    sum_check: bool = True,
    keepdims: object = False,
    nan_policy: object = "propagate",
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ks_1samp(
    x: object,
    cdf: object,
    args: object = (),
    alternative: object = "two-sided",
    method: object = "auto",
    *,
    axis: object = 0,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ks_2samp(
    data1: object,
    data2: object,
    alternative: object = "two-sided",
    method: object = "auto",
    *,
    axis: object = 0,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kstest(
    rvs: object,
    cdf: object,
    args: object = (),
    N: object = 20,
    alternative: object = "two-sided",
    method: object = "auto",
    *,
    axis: object = 0,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tiecorrect(rankvals: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ranksums(
    x: object,
    y: object,
    alternative: object = "two-sided",
    *,
    axis: object = 0,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kruskal(*samples: object, nan_policy: object = "propagate", axis: object = 0, keepdims: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def friedmanchisquare(
    *samples: object, axis: object = 0, nan_policy: object = "propagate", keepdims: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def brunnermunzel(
    x: object,
    y: object,
    alternative: object = "two-sided",
    distribution: object = "t",
    nan_policy: object = "propagate",
    *,
    axis: object = 0,
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def combine_pvalues(
    pvalues: object,
    method: object = "fisher",
    weights: object = None,
    *,
    axis: object = 0,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def wasserstein_distance(u_values: object, v_values: object, u_weights: object = None, v_weights: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def energy_distance(u_values: object, v_values: object, u_weights: object = None, v_weights: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rankdata(a: object, method: object = "average", *, axis: object = None, nan_policy: object = "propagate") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def linregress(
    x: object,
    y: object,
    alternative: object = "two-sided",
    *,
    axis: object = 0,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...

# mstats_basic
@deprecated("will be removed in SciPy v2.0.0")
def theilslopes(
    y: object,
    x: object = None,
    alpha: object = 0.95,
    method: object = "separate",
    *,
    axis: object = None,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def siegelslopes(
    y: object,
    x: object = None,
    method: object = "hierarchical",
    *,
    axis: object = None,
    nan_policy: object = "propagate",
    keepdims: object = False,
) -> object: ...
