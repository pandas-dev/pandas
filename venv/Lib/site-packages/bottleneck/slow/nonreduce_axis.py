import numpy as np
from numpy import partition, argpartition

__all__ = ["rankdata", "nanrankdata", "partition", "argpartition", "push"]


def rankdata(a, axis=None):
    "Slow rankdata function used for unaccelerated dtypes."
    return _rank(scipy_rankdata, a, axis)


def nanrankdata(a, axis=None):
    "Slow nanrankdata function used for unaccelerated dtypes."
    return _rank(_nanrankdata_1d, a, axis)


def _rank(func1d, a, axis):
    a = np.asarray(a)
    if axis is None:
        a = a.ravel()
        axis = 0
    if a.size == 0:
        y = a.astype(np.float64, copy=True)
    else:
        y = np.apply_along_axis(func1d, axis, a)
        if a.dtype != np.float64:
            y = y.astype(np.float64)
    return y


def _nanrankdata_1d(a):
    y = np.empty(a.shape, dtype=np.float64)
    y.fill(np.nan)
    idx = ~np.isnan(a)
    y[idx] = scipy_rankdata(a[idx])
    return y


def push(a, n=None, axis=-1):
    "Slow push used for unaccelerated dtypes."
    if n is None:
        n = np.inf
    y = np.array(a)
    ndim = y.ndim
    if axis != -1 or axis != ndim - 1:
        y = np.rollaxis(y, axis, ndim)
    if ndim == 1:
        y = y[None, :]
    elif ndim == 0:
        return y
    fidx = ~np.isnan(y)
    recent = np.empty(y.shape[:-1])
    count = np.empty(y.shape[:-1])
    recent.fill(np.nan)
    count.fill(np.nan)
    with np.errstate(invalid="ignore"):
        for i in range(y.shape[-1]):
            idx = (i - count) > n
            recent[idx] = np.nan
            idx = ~fidx[..., i]
            y[idx, i] = recent[idx]
            idx = fidx[..., i]
            count[idx] = i
            recent[idx] = y[idx, i]
    if axis != -1 or axis != ndim - 1:
        y = np.rollaxis(y, ndim - 1, axis)
    if ndim == 1:
        return y[0]
    return y


# ---------------------------------------------------------------------------
#
# SciPy
#
# Local copy of SciPy's rankdata to avoid a SciPy dependency. The SciPy
# license is included in the Bottleneck license file, which is distributed
# with Bottleneck.
#
# Code taken from scipy master branch on Aug 31, 2016.


def scipy_rankdata(a, method="average"):
    """
    rankdata(a, method='average')
    Assign ranks to data, dealing with ties appropriately.
    Ranks begin at 1.  The `method` argument controls how ranks are assigned
    to equal values.  See [1]_ for further discussion of ranking methods.
    Parameters
    ----------
    a : array_like
        The array of values to be ranked.  The array is first flattened.
    method : str, optional
        The method used to assign ranks to tied elements.
        The options are 'average', 'min', 'max', 'dense' and 'ordinal'.
        'average':
            The average of the ranks that would have been assigned to
            all the tied values is assigned to each value.
        'min':
            The minimum of the ranks that would have been assigned to all
            the tied values is assigned to each value.  (This is also
            referred to as "competition" ranking.)
        'max':
            The maximum of the ranks that would have been assigned to all
            the tied values is assigned to each value.
        'dense':
            Like 'min', but the rank of the next highest element is assigned
            the rank immediately after those assigned to the tied elements.
        'ordinal':
            All values are given a distinct rank, corresponding to the order
            that the values occur in `a`.
        The default is 'average'.
    Returns
    -------
    ranks : ndarray
         An array of length equal to the size of `a`, containing rank
         scores.
    References
    ----------
    .. [1] "Ranking", http://en.wikipedia.org/wiki/Ranking
    Examples
    --------
    >>> from scipy.stats import rankdata
    >>> rankdata([0, 2, 3, 2])
    array([ 1. ,  2.5,  4. ,  2.5])
    >>> rankdata([0, 2, 3, 2], method='min')
    array([ 1,  2,  4,  2])
    >>> rankdata([0, 2, 3, 2], method='max')
    array([ 1,  3,  4,  3])
    >>> rankdata([0, 2, 3, 2], method='dense')
    array([ 1,  2,  3,  2])
    >>> rankdata([0, 2, 3, 2], method='ordinal')
    array([ 1,  2,  4,  3])
    """
    if method not in ("average", "min", "max", "dense", "ordinal"):
        raise ValueError('unknown method "{0}"'.format(method))

    a = np.ravel(np.asarray(a))
    algo = "mergesort" if method == "ordinal" else "quicksort"
    sorter = np.argsort(a, kind=algo)

    inv = np.empty(sorter.size, dtype=np.intp)
    inv[sorter] = np.arange(sorter.size, dtype=np.intp)

    if method == "ordinal":
        return inv + 1

    a = a[sorter]
    obs = np.r_[True, a[1:] != a[:-1]]
    dense = obs.cumsum()[inv]

    if method == "dense":
        return dense

    # cumulative counts of each unique value
    count = np.r_[np.nonzero(obs)[0], len(obs)]

    if method == "max":
        return count[dense]

    if method == "min":
        return count[dense - 1] + 1

    # average method
    return 0.5 * (count[dense] + count[dense - 1] + 1)
