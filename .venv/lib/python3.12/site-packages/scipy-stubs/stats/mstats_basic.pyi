# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

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

@deprecated("will be removed in SciPy v2.0.0")
def argstoarray(*args: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def find_repeats(arr: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def count_tied_groups(x: object, use_missing: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def rankdata(data: object, axis: object = None, use_missing: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def mode(a: object, axis: object = 0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def msign(x: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pearsonr(x: object, y: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def spearmanr(
    x: object,
    y: object = None,
    use_ties: object = True,
    axis: object = None,
    nan_policy: object = "propagate",
    alternative: object = "two-sided",
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kendalltau(
    x: object,
    y: object,
    use_ties: object = True,
    use_missing: object = False,
    method: object = "auto",
    alternative: object = "two-sided",
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kendalltau_seasonal(x: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def pointbiserialr(x: object, y: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def linregress(x: object, y: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def theilslopes(y: object, x: object = None, alpha: object = 0.95, method: object = "separate") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def siegelslopes(y: object, x: object = None, method: object = "hierarchical") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sen_seasonal_slopes(x: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ttest_onesamp(a: object, popmean: object, axis: object = 0, alternative: object = "two-sided") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ttest_1samp(a: object, popmean: object, axis: object = 0, alternative: object = "two-sided") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ttest_ind(a: object, b: object, axis: object = 0, equal_var: object = True, alternative: object = "two-sided") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ttest_rel(a: object, b: object, axis: object = 0, alternative: object = "two-sided") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def mannwhitneyu(x: object, y: object, use_continuity: object = True) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kruskal(*args: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kruskalwallis(*args: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ks_1samp(x: object, cdf: object, args: object = (), alternative: object = "two-sided", method: object = "auto") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ks_2samp(data1: object, data2: object, alternative: object = "two-sided", method: object = "auto") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ks_twosamp(data1: object, data2: object, alternative: object = "two-sided", method: object = "auto") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kstest(
    data1: object, data2: object, args: object = (), alternative: object = "two-sided", method: object = "auto"
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trima(a: object, limits: object = None, inclusive: object = (True, True)) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trimr(a: object, limits: object = None, inclusive: object = (True, True), axis: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trim(
    a: object, limits: object = None, inclusive: object = (True, True), relative: object = False, axis: object = None
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trimboth(data: object, proportiontocut: object = 0.2, inclusive: object = (True, True), axis: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trimtail(
    data: object, proportiontocut: object = 0.2, tail: object = "left", inclusive: object = (True, True), axis: object = None
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trimmed_mean(
    a: object, limits: object = (0.1, 0.1), inclusive: object = (1, 1), relative: object = True, axis: object = None
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trimmed_var(
    a: object,
    limits: object = (0.1, 0.1),
    inclusive: object = (1, 1),
    relative: object = True,
    axis: object = None,
    ddof: object = 0,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trimmed_std(
    a: object,
    limits: object = (0.1, 0.1),
    inclusive: object = (1, 1),
    relative: object = True,
    axis: object = None,
    ddof: object = 0,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def trimmed_stde(a: object, limits: object = (0.1, 0.1), inclusive: object = (1, 1), axis: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tmean(a: object, limits: object = None, inclusive: object = (True, True), axis: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tvar(a: object, limits: object = None, inclusive: object = (True, True), axis: object = 0, ddof: object = 1) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tmin(a: object, lowerlimit: object = None, axis: object = 0, inclusive: object = True) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tmax(a: object, upperlimit: object = None, axis: object = 0, inclusive: object = True) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tsem(a: object, limits: object = None, inclusive: object = (True, True), axis: object = 0, ddof: object = 1) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def winsorize(
    a: object,
    limits: object = None,
    inclusive: object = (True, True),
    inplace: object = False,
    axis: object = None,
    nan_policy: object = "propagate",
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def moment(a: object, moment: object = 1, axis: object = 0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def variation(a: object, axis: object = 0, ddof: object = 0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def skew(a: object, axis: object = 0, bias: object = True) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kurtosis(a: object, axis: object = 0, fisher: object = True, bias: object = True) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def describe(a: object, axis: object = 0, ddof: object = 0, bias: object = True) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def skewtest(a: object, axis: object = 0, alternative: object = "two-sided") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kurtosistest(a: object, axis: object = 0, alternative: object = "two-sided") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def normaltest(a: object, axis: object = 0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def mquantiles(
    a: object,
    prob: object = (0.25, 0.5, 0.75),
    alphap: object = 0.4,
    betap: object = 0.4,
    axis: object = None,
    limit: object = (),
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def scoreatpercentile(data: object, per: object, limit: object = (), alphap: object = 0.4, betap: object = 0.4) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def meppf(data: object, alpha: object = 0.4, beta: object = 0.4) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def plotting_positions(data: object, alpha: object = 0.4, beta: object = 0.4) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def obrientransform(*args: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sem(a: object, axis: object = 0, ddof: object = 1) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def f_oneway(*args: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def friedmanchisquare(*args: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def brunnermunzel(x: object, y: object, alternative: object = "two-sided", distribution: object = "t") -> object: ...
