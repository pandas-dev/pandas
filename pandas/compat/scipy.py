"""
Shipping functions from SciPy to reduce dependency on having SciPy installed
"""

import numpy as np


def scoreatpercentile(a, per, limit=(), interpolation_method='fraction'):
    """
    Calculate the score at the given `per` percentile of the sequence `a`.

    For example, the score at `per=50` is the median. If the desired quantile
    lies between two data points, we interpolate between them, according to
    the value of `interpolation`. If the parameter `limit` is provided, it
    should be a tuple (lower, upper) of two values. Values of `a` outside
    this (closed) interval will be ignored.

    The `interpolation_method` parameter supports three values, namely
    `fraction` (default), `lower` and `higher`. Interpolation is done only,
    if the desired quantile lies between two data points `i` and `j`. For
    `fraction`, the result is an interpolated value between `i` and `j`;
    for `lower`, the result is `i`, for `higher` the result is `j`.

    Parameters
    ----------
    a : ndarray
        Values from which to extract score.
    per : scalar
        Percentile at which to extract score.
    limit : tuple, optional
        Tuple of two scalars, the lower and upper limits within which to
        compute the percentile.
    interpolation_method : {'fraction', 'lower', 'higher'}, optional
        This optional parameter specifies the interpolation method to use,
        when the desired quantile lies between two data points `i` and `j`:

        - fraction: `i + (j - i)*fraction`, where `fraction` is the
                    fractional part of the index surrounded by `i` and `j`.
        - lower: `i`.
        - higher: `j`.

    Returns
    -------
    score : float
        Score at percentile.

    See Also
    --------
    percentileofscore

    Examples
    --------
    >>> from scipy import stats
    >>> a = np.arange(100)
    >>> stats.scoreatpercentile(a, 50)
    49.5

    """
    # TODO: this should be a simple wrapper around a well-written quantile
    # function.  GNU R provides 9 quantile algorithms (!), with differing
    # behaviour at, for example, discontinuities.
    values = np.sort(a, axis=0)
    if limit:
        values = values[(limit[0] <= values) & (values <= limit[1])]

    idx = per / 100. * (values.shape[0] - 1)
    if (idx % 1 == 0):
        score = values[idx]
    else:
        if interpolation_method == 'fraction':
            score = _interpolate(values[int(idx)], values[int(idx) + 1],
                                 idx % 1)
        elif interpolation_method == 'lower':
            score = values[np.floor(idx)]
        elif interpolation_method == 'higher':
            score = values[np.ceil(idx)]
        else:
            raise ValueError("interpolation_method can only be 'fraction', "
                             "'lower' or 'higher'")

    return score


def _interpolate(a, b, fraction):
    """Returns the point at the given fraction between a and b, where
    'fraction' must be between 0 and 1.
    """
    return a + (b - a) * fraction


def rankdata(a):
    """
    Ranks the data, dealing with ties appropriately.

    Equal values are assigned a rank that is the average of the ranks that
    would have been otherwise assigned to all of the values within that set.
    Ranks begin at 1, not 0.

    Parameters
    ----------
    a : array_like
        This array is first flattened.

    Returns
    -------
    rankdata : ndarray
         An array of length equal to the size of `a`, containing rank scores.

    Examples
    --------
    >>> stats.rankdata([0, 2, 2, 3])
    array([ 1. ,  2.5,  2.5,  4. ])

    """
    a = np.ravel(a)
    n = len(a)
    svec, ivec = fastsort(a)
    sumranks = 0
    dupcount = 0
    newarray = np.zeros(n, float)
    for i in xrange(n):
        sumranks += i
        dupcount += 1
        if i == n - 1 or svec[i] != svec[i + 1]:
            averank = sumranks / float(dupcount) + 1
            for j in xrange(i - dupcount + 1, i + 1):
                newarray[ivec[j]] = averank
            sumranks = 0
            dupcount = 0
    return newarray


def fastsort(a):
    """
    Sort an array and provide the argsort.

    Parameters
    ----------
    a : array_like
        Input array.

    Returns
    -------
    fastsort : ndarray of type int
        sorted indices into the original array

    """
    # TODO: the wording in the docstring is nonsense.
    it = np.argsort(a)
    as_ = a[it]
    return as_, it


def percentileofscore(a, score, kind='rank'):
    '''
    The percentile rank of a score relative to a list of scores.

    A `percentileofscore` of, for example, 80% means that 80% of the
    scores in `a` are below the given score. In the case of gaps or
    ties, the exact definition depends on the optional keyword, `kind`.

    Parameters
    ----------
    a: array like
        Array of scores to which `score` is compared.
    score: int or float
        Score that is compared to the elements in `a`.
    kind: {'rank', 'weak', 'strict', 'mean'}, optional
        This optional parameter specifies the interpretation of the
        resulting score:

        - "rank": Average percentage ranking of score.  In case of
                  multiple matches, average the percentage rankings of
                  all matching scores.
        - "weak": This kind corresponds to the definition of a cumulative
                  distribution function.  A percentileofscore of 80%
                  means that 80% of values are less than or equal
                  to the provided score.
        - "strict": Similar to "weak", except that only values that are
                    strictly less than the given score are counted.
        - "mean": The average of the "weak" and "strict" scores, often used in
                  testing.  See

                  http://en.wikipedia.org/wiki/Percentile_rank

    Returns
    -------
    pcos : float
        Percentile-position of score (0-100) relative to `a`.

    Examples
    --------
    Three-quarters of the given values lie below a given score:

    >>> percentileofscore([1, 2, 3, 4], 3)
    75.0

    With multiple matches, note how the scores of the two matches, 0.6
    and 0.8 respectively, are averaged:

    >>> percentileofscore([1, 2, 3, 3, 4], 3)
    70.0

    Only 2/5 values are strictly less than 3:

    >>> percentileofscore([1, 2, 3, 3, 4], 3, kind='strict')
    40.0

    But 4/5 values are less than or equal to 3:

    >>> percentileofscore([1, 2, 3, 3, 4], 3, kind='weak')
    80.0

    The average between the weak and the strict scores is

    >>> percentileofscore([1, 2, 3, 3, 4], 3, kind='mean')
    60.0

    '''
    a = np.array(a)
    n = len(a)

    if kind == 'rank':
        if not(np.any(a == score)):
            a = np.append(a, score)
            a_len = np.array(range(len(a)))
        else:
            a_len = np.array(range(len(a))) + 1.0

        a = np.sort(a)
        idx = [a == score]
        pct = (np.mean(a_len[idx]) / n) * 100.0
        return pct

    elif kind == 'strict':
        return sum(a < score) / float(n) * 100
    elif kind == 'weak':
        return sum(a <= score) / float(n) * 100
    elif kind == 'mean':
        return (sum(a < score) + sum(a <= score)) * 50 / float(n)
    else:
        raise ValueError("kind can only be 'rank', 'strict', 'weak' or 'mean'")
