import numpy as np
import scipy._external.array_api_extra as xpx
from scipy._lib._array_api import xp_capabilities, array_namespace, xp_promote, is_jax
from scipy.optimize.elementwise import find_root
from scipy.special import ndtri
from scipy.special import _ufuncs as scu
from ._common import ConfidenceInterval


class BinomTestResult:
    """
    Result of `scipy.stats.binomtest`.

    Attributes
    ----------
    k : int
        The number of successes (copied from `binomtest` input).
    n : int
        The number of trials (copied from `binomtest` input).
    alternative : str
        Indicates the alternative hypothesis specified in the input
        to `binomtest`.  It will be one of ``'two-sided'``, ``'greater'``,
        or ``'less'``.
    statistic: float
        The estimate of the proportion of successes.
    pvalue : float
        The p-value of the hypothesis test.

    """
    def __init__(self, k, n, alternative, statistic, pvalue, xp):
        self.k = k
        self.n = n
        self.alternative = alternative
        self.statistic = statistic
        self.pvalue = pvalue
        self._xp = xp

        # add alias for backward compatibility
        self.proportion_estimate = statistic

    def __repr__(self):
        s = ("BinomTestResult("
             f"k={self.k}, "
             f"n={self.n}, "
             f"alternative={self.alternative!r}, "
             f"statistic={self.statistic}, "
             f"pvalue={self.pvalue})")
        return s

    def proportion_ci(self, confidence_level=0.95, method='exact'):
        """
        Compute the confidence interval for ``statistic``.

        Parameters
        ----------
        confidence_level : float, optional
            Confidence level for the computed confidence interval
            of the estimated proportion. Default is 0.95.
        method : {'exact', 'wilson', 'wilsoncc'}, optional
            Selects the method used to compute the confidence interval
            for the estimate of the proportion:

            'exact' :
                Use the Clopper-Pearson exact method [1]_.
            'wilson' :
                Wilson's method, without continuity correction ([2]_, [3]_).
            'wilsoncc' :
                Wilson's method, with continuity correction ([2]_, [3]_).

            Default is ``'exact'``.

        Returns
        -------
        ci : ``ConfidenceInterval`` object
            The object has attributes ``low`` and ``high`` that hold the
            lower and upper bounds of the confidence interval.

        References
        ----------
        .. [1] C. J. Clopper and E. S. Pearson, The use of confidence or
               fiducial limits illustrated in the case of the binomial,
               Biometrika, Vol. 26, No. 4, pp 404-413 (Dec. 1934).
        .. [2] E. B. Wilson, Probable inference, the law of succession, and
               statistical inference, J. Amer. Stat. Assoc., 22, pp 209-212
               (1927).
        .. [3] Robert G. Newcombe, Two-sided confidence intervals for the
               single proportion: comparison of seven methods, Statistics
               in Medicine, 17, pp 857-872 (1998).

        Examples
        --------
        >>> from scipy.stats import binomtest
        >>> result = binomtest(k=7, n=50, p=0.1)
        >>> result.statistic
        0.14
        >>> result.proportion_ci()
        ConfidenceInterval(low=0.05819170033997342, high=0.26739600249700846)
        """
        if method not in ('exact', 'wilson', 'wilsoncc'):
            raise ValueError(f"method ('{method}') must be one of 'exact', "
                             "'wilson' or 'wilsoncc'.")
        if not (0 <= confidence_level <= 1):
            raise ValueError(f'confidence_level ({confidence_level}) must be in '
                             'the interval [0, 1].')

        xp = self._xp
        k, n, confidence_level = xp_promote(self.k, self.n, confidence_level, xp=xp)
        if method == 'exact':
            low, high = _binom_exact_conf_int(
                k, n, confidence_level, self.alternative, xp=xp)
        else:
            # method is 'wilson' or 'wilsoncc'
            correction = method == 'wilsoncc'
            low, high = _binom_wilson_conf_int(
                k, n, confidence_level, self.alternative, correction, xp=xp)
        low = low[()] if low.ndim == 0 else low
        high = high[()] if high.ndim == 0 else high
        return ConfidenceInterval(low=low, high=high)


def _binom_exact_conf_int(k, n, confidence_level, alternative, *, xp):
    """
    Compute the estimate and confidence interval for the binomial test.

    Returns proportion, prop_low, prop_high
    """
    init = (xp.zeros_like(k), xp.ones_like(k))
    args = (k, n)
    alpha = ((1 - confidence_level) / 2 if alternative == 'two-sided'
             else 1 - confidence_level)

    # I think using the private methods here is fine, since we will only evaluate with
    # valid `p`, `k`, and `n` (or all NaNs). One exception is when `k=0` and
    # `binom._sf(k-1, n, p)`: evaluates to NaN, but that's not a problem because
    # `plow` has a special case for `k=0` below.
    plow = (xp.zeros_like(k) if alternative == 'less' else
            find_root(lambda p, k, n: _SimpleBinomial(n, p).sf(k-1) - alpha,
                      init, args=args).x)
    phigh = (xp.ones_like(k) if alternative == 'greater' else
             find_root(lambda p, k, n: _SimpleBinomial(n, p).cdf(k) - alpha,
                       init, args=args).x)

    plow = xp.where(k == 0, 0.0, plow)
    phigh = xp.where(k == n, 1.0, phigh)
    return plow, phigh


def _binom_wilson_conf_int(k, n, confidence_level, alternative, correction, *, xp):
    # This function assumes that the arguments have already been validated.
    # In particular, `alternative` must be one of 'two-sided', 'less' or
    # 'greater'.
    p = k / n
    if alternative == 'two-sided':
        z = ndtri(0.5 + 0.5*confidence_level)
    else:
        z = ndtri(confidence_level)

    # For reference, the formulas implemented here are from
    # Newcombe (1998) (ref. [3] in the proportion_ci docstring).
    denom = 2*(n + z**2)
    center = (2*n*p + z**2)/denom
    q = 1 - p

    if correction:
        with np.errstate(divide='ignore', invalid='ignore'):
            dlo = (1 + z*xp.sqrt(z**2 - 2 - 1/n + 4*p*(n*q + 1)))/denom
            dhi = (1 + z*xp.sqrt(z**2 + 2 - 1/n + 4*p*(n*q - 1)))/denom
    else:
        delta = z / denom * xp.sqrt(4*n*p*q + z**2)
        dlo, dhi = delta, delta

    lo = xp.where((k == 0) | (alternative == 'less'), 0.0, center - dlo)
    hi = xp.where((k == n) | (alternative == 'greater'), 1.0, center + dhi)
    return lo, hi


@xp_capabilities(skip_backends=[('dask.array', "")], cpu_only=True,
                 reason="binomial distribution ufuncs only available for NumPy",
                 extra_note="`alternative='two-sided'` is incompatible with JAX arrays.")
def binomtest(k, n, p=0.5, alternative='two-sided'):
    """
    Perform a test that the probability of success is p.

    The binomial test [1]_ is a test of the null hypothesis that the
    probability of success in a Bernoulli experiment is `p`.

    Details of the test can be found in many texts on statistics, such
    as section 24.5 of [2]_.

    The documentation is written as though the function accepts and returns Python
    scalars, but the function is vectorized to work elementwise with NumPy arrays.

    Parameters
    ----------
    k : int
        The number of successes.
    n : int
        The number of trials.
    p : float, optional
        The hypothesized probability of success, i.e. the expected
        proportion of successes.  The value must be in the interval
        ``0 <= p <= 1``. The default value is ``p = 0.5``.
    alternative : {'two-sided', 'greater', 'less'}, optional
        Indicates the alternative hypothesis. The default value is
        'two-sided'.

    Returns
    -------
    result : `~scipy.stats._result_classes.BinomTestResult` instance
        The return value is an object with the following attributes:

        k : int
            The number of successes (copied from `binomtest` input).
        n : int
            The number of trials (copied from `binomtest` input).
        alternative : str
            Indicates the alternative hypothesis specified in the input
            to `binomtest`.  It will be one of ``'two-sided'``, ``'greater'``,
            or ``'less'``.
        statistic : float
            The estimate of the proportion of successes.
        pvalue : float
            The p-value of the hypothesis test.

        The object has the following methods:

        proportion_ci(confidence_level=0.95, method='exact') :
            Compute the confidence interval for ``statistic``.

    Notes
    -----
    .. versionadded:: 1.7.0

    References
    ----------
    .. [1] Binomial test, https://en.wikipedia.org/wiki/Binomial_test
    .. [2] Jerrold H. Zar, Biostatistical Analysis (fifth edition),
           Prentice Hall, Upper Saddle River, New Jersey USA (2010)

    Examples
    --------
    >>> from scipy.stats import binomtest

    A car manufacturer claims that no more than 10% of their cars are unsafe.
    15 cars are inspected for safety, 3 were found to be unsafe. Test the
    manufacturer's claim:

    >>> result = binomtest(3, n=15, p=0.1, alternative='greater')
    >>> result.pvalue
    0.18406106910639114

    The null hypothesis cannot be rejected at the 5% level of significance
    because the returned p-value is greater than the critical value of 5%.

    The test statistic is equal to the estimated proportion, which is simply
    ``3/15``:

    >>> result.statistic
    0.2

    We can use the `proportion_ci()` method of the result to compute the
    confidence interval of the estimate:

    >>> result.proportion_ci(confidence_level=0.95)
    ConfidenceInterval(low=0.05684686759024681, high=1.0)

    """
    xp = array_namespace(k, n, p)
    k, n, p = xp_promote(k, n, p, force_floating=True, broadcast=True, xp=xp)
    k_valid = (k >= 0) & (k <= n) & (k == xp.floor(k))
    n_valid = (n >= 1) & (n == xp.floor(n))
    p_valid = (p >= 0) & (p <= 1)
    valid = k_valid & n_valid & p_valid
    k = xp.where(valid, k, xp.nan)
    n = xp.where(valid, n, xp.nan)
    p = xp.where(valid, p, xp.nan)

    if alternative not in ('two-sided', 'less', 'greater'):
        raise ValueError(f"alternative ('{alternative}') not recognized; \n"
                         "must be 'two-sided', 'less' or 'greater'")

    B = _SimpleBinomial(n, p)
    if alternative == 'less':
        pval = B.cdf(k)
    elif alternative == 'greater':
        pval = B.sf(k - 1)
    else:
        if is_jax(xp):
            message = "`alternative='two-sided'` is incompatible with JAX arrays."
            raise ValueError(message)

        # alternative is 'two-sided'
        d = B.pmf(k)
        rerr = 1 + 1e-7

        def k_lt_pn(d, k, p, n):
            B = _SimpleBinomial(n, p)
            ix = _binary_search_for_binom_tst(lambda x1: -B.pmf(x1), -d*rerr,
                                              xp.ceil(p * n), n, xp=xp)
            # y is the number of terms between mode and n that are <= d*rerr.
            # ix gave us the first term where a(ix) <= d*rerr < a(ix-1)
            # if the first equality doesn't hold, y=n-ix. Otherwise, we
            # need to include ix as well as the equality holds. Note that
            # the equality will hold in very very rare situations due to rerr.
            y = n - ix + xp.asarray(d*rerr == B.pmf(ix), dtype=ix.dtype)
            pval = B.cdf(k) + B.sf(n - y)
            return pval

        def k_gte_pn(d, k, p, n):
            B = _SimpleBinomial(n, p)
            ix = _binary_search_for_binom_tst(B.pmf, d*rerr,
                                              xp.zeros_like(n), xp.floor(p * n), xp=xp)
            # y is the number of terms between 0 and mode that are <= d*rerr.
            # we need to add a 1 to account for the 0 index.
            # For comparing this with old behavior, see
            # tst_binary_srch_for_binom_tst method in test_morestats.
            y = ix + 1
            pval = B.cdf(y-1) + B.sf(k-1)
            return pval

        pval = xpx.apply_where(k < p*n, (d, k, p, n), k_lt_pn,  k_gte_pn)
        # xp.minimum(1.0, pval) but for data-apis/array-api-compat#271
        pval = xp.minimum(xp.asarray(1.0, dtype=pval.dtype), pval)

    statistic = xp.where(valid, k/n, xp.nan)
    pval = xp.where(valid, pval, xp.nan)
    if statistic.ndim == 0:
        k, n, statistic, pval = k[()], n[()], statistic[()], pval[()]

    result = BinomTestResult(k=k, n=n, alternative=alternative,
                             statistic=statistic, pvalue=pval, xp=xp)
    return result


def _binary_search_for_binom_tst(a, d, lo, hi, *, xp):
    """
    Conducts an implicit binary search on a function specified by `a`.

    Meant to be used on the binomial PMF for the case of two-sided tests
    to obtain the value on the other side of the mode where the tail
    probability should be computed. The values on either side of
    the mode are always in order, meaning binary search is applicable.

    Parameters
    ----------
    a : callable
      The function over which to perform binary search. Its values
      for inputs lo and hi should be in ascending order.
    d : float
      The value to search.
    lo : int
      The lower end of range to search.
    hi : int
      The higher end of the range to search.

    Returns
    -------
    int
      The index, i between lo and hi
      such that a(i)<=d<a(i+1)
    """
    d = xp.asarray(d, copy=True)
    lo = xp.asarray(lo, copy=True)
    hi = xp.asarray(hi, copy=True)
    while xp.any(lo < hi):
        mid = lo + (hi-lo)//2
        midval = a(mid)

        i_lt = midval < d
        lo = xpx.at(lo)[i_lt].set(mid[i_lt] + 1)

        i_gt = midval > d
        hi = xpx.at(hi)[i_gt].set(mid[i_gt] - 1)

        i_eq = (midval == d)
        mid_i_eq = mid[i_eq]
        lo = xpx.at(lo)[i_eq].set(mid_i_eq)
        hi = xpx.at(hi)[i_eq].set(mid_i_eq)

    return xp.where(a(lo) <= d, lo, lo-1)


class _SimpleBinomial:
    # A very simple, array-API compatible binomial distribution for use in
    # hypothesis tests. May be replaced by new infrastructure Binomial
    # distribution in due time.
    def __init__(self, n, p, *, xp=None):
        xp = array_namespace(n, p) if xp is None else xp
        self.n = n
        self.p = p
        self.xp = xp

    def f(self, x, fun):
        xp = self.xp
        return xpx.lazy_apply(fun, x, self.n, self.p, as_numpy=True, xp=xp)

    def cdf(self, x):
        return self.xp.where(x >= 0, self.f(self.xp.floor(x), scu._binom_cdf), 0.0)

    def sf(self, x):
        return self.xp.where(x >= 0, self.f(self.xp.floor(x), scu._binom_sf), 1.0)

    def ppf(self, x):
        return self.f(x, scu._binom_ppf)

    def isf(self, x):
        return self.f(x, scu._binom_isf)

    def pmf(self, x):
        return self.f(self.xp.floor(x), scu._binom_pmf)
