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
    reps: object = ...,
    workers: object = ...,
    is_twosamp: object = ...,
    random_state: object = ...,
) -> object: ...

# _stats
@deprecated("will be removed in SciPy v2.0.0")
def gmean(
    a: object, axis: object = ..., dtype: object = ..., weights: object = ..., *, nan_policy: object = ..., keepdims: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hmean(
    a: object, axis: object = ..., dtype: object = ..., *, weights: object = ..., nan_policy: object = ..., keepdims: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pmean(
    a: object,
    p: object,
    *,
    dtype: object = ...,
    weights: object = ...,
    axis: object = ...,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def mode(a: object, axis: object = ..., nan_policy: object = ..., keepdims: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tmean(
    a: object,
    limits: object = ...,
    inclusive: object = ...,
    axis: object = ...,
    *,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tvar(
    a: object,
    limits: object = ...,
    inclusive: object = ...,
    axis: object = ...,
    ddof: object = ...,
    *,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tmin(
    a: object,
    lowerlimit: object = ...,
    axis: object = ...,
    inclusive: object = ...,
    nan_policy: object = ...,
    *,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tmax(
    a: object,
    upperlimit: object = ...,
    axis: object = ...,
    inclusive: object = ...,
    nan_policy: object = ...,
    *,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tstd(
    a: object,
    limits: object = ...,
    inclusive: object = ...,
    axis: object = ...,
    ddof: object = ...,
    *,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tsem(
    a: object,
    limits: object = ...,
    inclusive: object = ...,
    axis: object = ...,
    ddof: object = ...,
    *,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def moment(
    a: object, order: object = ..., axis: object = ..., nan_policy: object = ..., *, center: object = ..., keepdims: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def skew(a: object, axis: object = ..., bias: object = ..., nan_policy: object = ..., *, keepdims: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kurtosis(
    a: object, axis: object = ..., fisher: object = ..., bias: object = ..., nan_policy: object = ..., *, keepdims: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def describe(a: object, axis: object = ..., ddof: object = ..., bias: object = ..., nan_policy: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def skewtest(
    a: object, axis: object = ..., nan_policy: object = ..., alternative: object = ..., *, keepdims: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kurtosistest(
    a: object, axis: object = ..., nan_policy: object = ..., alternative: object = ..., *, keepdims: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def normaltest(a: object, axis: object = ..., nan_policy: object = ..., *, keepdims: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def jarque_bera(x: object, *, axis: object = ..., nan_policy: object = ..., keepdims: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def scoreatpercentile(
    a: object, per: object, limit: object = ..., interpolation_method: object = ..., axis: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def percentileofscore(a: object, score: object, kind: object = ..., nan_policy: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def cumfreq(a: object, numbins: object = ..., defaultreallimits: object = ..., weights: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def relfreq(a: object, numbins: object = ..., defaultreallimits: object = ..., weights: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def obrientransform(*samples: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sem(a: object, axis: object = ..., ddof: object = ..., nan_policy: object = ..., *, keepdims: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def zscore(a: object, axis: object = ..., ddof: object = ..., nan_policy: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gzscore(a: object, *, axis: object = ..., ddof: object = ..., nan_policy: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def zmap(scores: object, compare: object, axis: object = ..., ddof: object = ..., nan_policy: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gstd(a: object, axis: object = ..., ddof: object = ..., *, keepdims: object = ..., nan_policy: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def iqr(
    x: object,
    axis: object = ...,
    rng: object = ...,
    scale: object = ...,
    nan_policy: object = ...,
    interpolation: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def median_abs_deviation(
    x: object, axis: object = ..., center: object = ..., scale: object = ..., nan_policy: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sigmaclip(a: object, low: object = ..., high: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trimboth(a: object, proportiontocut: object, axis: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trim1(a: object, proportiontocut: object, tail: object = ..., axis: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trim_mean(a: object, proportiontocut: object, axis: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def f_oneway(
    *samples: object, axis: object = ..., equal_var: object = ..., nan_policy: object = ..., keepdims: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def alexandergovern(*samples: object, nan_policy: object = ..., axis: object = ..., keepdims: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pearsonr(x: object, y: object, *, alternative: object = ..., method: object = ..., axis: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fisher_exact(table: object, alternative: object = ..., *, method: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def spearmanr(a: object, b: object = ..., axis: object = ..., nan_policy: object = ..., alternative: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pointbiserialr(x: object, y: object, *, axis: object = ..., nan_policy: object = ..., keepdims: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kendalltau(
    x: object,
    y: object,
    *,
    method: object = ...,
    variant: object = ...,
    alternative: object = ...,
    axis: object = ...,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def weightedtau(
    x: object,
    y: object,
    rank: object = ...,
    weigher: object = ...,
    additive: object = ...,
    *,
    axis: object = ...,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ttest_1samp(
    a: object, popmean: object, axis: object = ..., nan_policy: object = ..., alternative: object = ..., *, keepdims: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ttest_ind_from_stats(
    mean1: object,
    std1: object,
    nobs1: object,
    mean2: object,
    std2: object,
    nobs2: object,
    equal_var: object = ...,
    alternative: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ttest_ind(
    a: object,
    b: object,
    *,
    axis: object = ...,
    equal_var: object = ...,
    nan_policy: object = ...,
    permutations: object = ...,
    random_state: object = ...,
    alternative: object = ...,
    trim: object = ...,
    method: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ttest_rel(
    a: object, b: object, axis: object = ..., nan_policy: object = ..., alternative: object = ..., *, keepdims: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def power_divergence(
    f_obs: object,
    f_exp: object = ...,
    ddof: object = ...,
    axis: object = ...,
    lambda_: object = ...,
    *,
    keepdims: object = ...,
    nan_policy: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def chisquare(
    f_obs: object,
    f_exp: object = ...,
    ddof: object = ...,
    axis: object = ...,
    *,
    sum_check: bool = ...,
    keepdims: object = ...,
    nan_policy: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ks_1samp(
    x: object,
    cdf: object,
    args: object = ...,
    alternative: object = ...,
    method: object = ...,
    *,
    axis: object = ...,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ks_2samp(
    data1: object,
    data2: object,
    alternative: object = ...,
    method: object = ...,
    *,
    axis: object = ...,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kstest(
    rvs: object,
    cdf: object,
    args: object = ...,
    N: object = ...,
    alternative: object = ...,
    method: object = ...,
    *,
    axis: object = ...,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tiecorrect(rankvals: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ranksums(
    x: object, y: object, alternative: object = ..., *, axis: object = ..., nan_policy: object = ..., keepdims: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kruskal(*samples: object, nan_policy: object = ..., axis: object = ..., keepdims: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def friedmanchisquare(*samples: object, axis: object = ..., nan_policy: object = ..., keepdims: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def brunnermunzel(
    x: object,
    y: object,
    alternative: object = ...,
    distribution: object = ...,
    nan_policy: object = ...,
    *,
    axis: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def combine_pvalues(
    pvalues: object,
    method: object = ...,
    weights: object = ...,
    *,
    axis: object = ...,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def wasserstein_distance(u_values: object, v_values: object, u_weights: object = ..., v_weights: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def energy_distance(u_values: object, v_values: object, u_weights: object = ..., v_weights: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def find_repeats(arr: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rankdata(a: object, method: object = ..., *, axis: object = ..., nan_policy: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def linregress(
    x: object, y: object, alternative: object = ..., *, axis: object = ..., nan_policy: object = ..., keepdims: object = ...
) -> object: ...

# mstats_basic
@deprecated("will be removed in SciPy v2.0.0")
def theilslopes(
    y: object,
    x: object = ...,
    alpha: object = ...,
    method: object = ...,
    *,
    axis: object = ...,
    nan_policy: object = ...,
    keepdims: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def siegelslopes(
    y: object, x: object = ..., method: object = ..., *, axis: object = ..., nan_policy: object = ..., keepdims: object = ...
) -> object: ...
