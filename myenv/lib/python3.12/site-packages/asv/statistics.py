# Licensed under a 3-clause BSD style license - see LICENSE.rst

# Author: Pauli Virtanen, 2016

import math
from operator import index


def get_weight(stats):
    """
    Return a data point weight for the result.
    """
    if stats is None or 'ci_99_a' not in stats or 'ci_99_b' not in stats:
        return None

    try:
        a = stats['ci_99_a']
        b = stats['ci_99_b']

        if math.isinf(a) or math.isinf(b):
            # Infinite interval is due to too few samples --- consider
            # weight as missing
            return None

        return 2 / abs(b - a)
    except ZeroDivisionError:
        return None


def is_different(samples_a, samples_b, stats_a, stats_b, p_threshold=0.002):
    """Check whether the samples are statistically different.

    If sample data is not provided, or the sample is too small, falls
    back to a pessimistic CI-based check. If it returns True, then the
    difference is statistically significant. If it returns False, it
    might or might not be statistically significant.

    Parameters
    ----------
    samples_a, samples_b
        Input samples
    stats_a, stats_b
        Input stats data

    """

    if samples_a is not None and samples_b is not None:
        # Raw data present: Mann-Whitney U test, but only if there's
        # enough data so that the test can return True
        a = [x for x in samples_a if not math.isnan(x)]
        b = [x for x in samples_b if not math.isnan(x)]

        p_min = 1 / binom(len(a) + len(b), min(len(a), len(b)))
        if p_min < p_threshold:
            _, p = mann_whitney_u(a, b)
            return p < p_threshold

    # If confidence intervals overlap, reject.
    # Corresponds to a test with ill-specified threshold p-value,
    # which generally can be significantly smaller than p <= 0.01
    # depending on the actual data. For normal test (known variance),
    # 0.00027 <= p <= 0.01.
    ci_a = (stats_a['ci_99_a'], stats_a['ci_99_b'])
    ci_b = (stats_b['ci_99_a'], stats_b['ci_99_b'])

    if ci_a[1] >= ci_b[0] and ci_a[0] <= ci_b[1]:
        return False

    return True


_mann_whitney_u_memo = {}


def mann_whitney_u(x, y, method='auto'):
    """
    Mann-Whitney U test

    Ties are handled conservatively, returning the least significant
    tie breaking.

    Parameters
    ----------
    x, y : list of float
        Samples to test
    method : {'auto', 'exact', 'normal'}
        Whether to compute p-value exactly of via normal approximation.
        The option 'auto' switches to approximation for sample size > 20.

    Returns
    -------
    u : int
        U-statistic
    p : float
        p-value for two-sided alternative

    References
    ----------
    .. [1] Mann & Whitney, Ann. Math. Statist. 18, 50 (1947).
    .. [2] Gibbons & Chakraborti, "Nonparametric statistical inference". (2003)

    """
    memo = _mann_whitney_u_memo
    if len(memo) > 100000:
        memo.clear()

    m = len(x)
    n = len(y)

    if method == 'auto':
        if max(m, n) > 20:
            method = 'normal'
        else:
            method = 'exact'

    u, ties = mann_whitney_u_u(x, y)

    # Conservative tie breaking
    if u <= m * n // 2 and u + ties >= m * n // 2:
        ties = m * n // 2 - u

    ux1 = min(u, m * n - u)
    ux2 = min(u + ties, m * n - (u + ties))

    if ux1 >= ux2:
        ux = ux1
    else:
        u = u + ties
        ux = ux2

    # Get p-value
    if method == 'exact':
        p1 = mann_whitney_u_cdf(m, n, ux, memo)
        p2 = 1.0 - mann_whitney_u_cdf(m, n, max(m * n // 2, m * n - ux - 1), memo)
        p = p1 + p2
    elif method == 'normal':
        N = m + n
        var = m * n * (N + 1) / 12
        z = (ux - m * n / 2) / math.sqrt(var)
        cdf = 0.5 * math.erfc(-z / math.sqrt(2))
        p = 2 * cdf
    else:
        raise ValueError(f"Unknown method {repr(method)}")

    return u, p


def mann_whitney_u_u(x, y):
    u = 0
    ties = 0
    for xx in x:
        for yy in y:
            if xx > yy:
                u += 1
            elif xx == yy:
                ties += 1
    return u, ties


def mann_whitney_u_cdf(m, n, u, memo=None):
    if memo is None:
        memo = {}
    cdf = 0
    for uu in range(u + 1):
        cdf += mann_whitney_u_pmf(m, n, uu, memo)
    return cdf


def mann_whitney_u_pmf(m, n, u, memo=None):
    if memo is None:
        memo = {}
    return mann_whitney_u_r(m, n, u, memo) / binom(m + n, m)


def mann_whitney_u_r(m, n, u, memo=None):
    """
    Number of orderings in Mann-Whitney U test.

    The PMF of U for samples of sizes (m, n) is given by
    p(u) = r(m, n, u) / binom(m + n, m).

    References
    ----------
    .. [1] Mann & Whitney, Ann. Math. Statist. 18, 50 (1947).
    """
    if u < 0:
        value = 0
    elif m == 0 or n == 0:
        value = 1 if u == 0 else 0
    else:
        # Don't bother figuring out table construction, memoization
        # sorts it out
        if memo is None:
            memo = {}
        key = (m, n, u)
        value = memo.get(key)
        if value is not None:
            return value

        value = (mann_whitney_u_r(m, n - 1, u, memo) +
                 mann_whitney_u_r(m - 1, n, u - n, memo))

        memo[key] = value
    return value


def binom(n, k):
    """
    Binomial coefficient (n over k)
    """
    n = index(n)
    k = index(k)
    if not 0 <= k <= n:
        return 0
    m = n + 1
    num = 1
    den = 1
    for j in range(1, min(k, n - k) + 1):
        num *= m - j
        den *= j
    return num // den
