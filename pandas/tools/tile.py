"""
Quantilization functions and related stuff
"""

from pandas.core.api import DataFrame, Series
from pandas.core.categorical import Categorical
from pandas.core.index import _ensure_index
from pandas.core.interval import IntervalIndex, Interval
import pandas.core.algorithms as algos
import pandas.core.common as com
import pandas.core.nanops as nanops
from pandas.compat import zip

import numpy as np

import warnings


def cut(x, bins, right=True, labels=None, retbins=False, precision=3,
        include_lowest=False):
    """
    Return indices of half-open bins to which each value of `x` belongs.

    Parameters
    ----------
    x : array-like
        Input array to be binned. It has to be 1-dimensional.
    bins : int or sequence of scalars
        If `bins` is an int, it defines the number of equal-width bins in the
        range of `x`. However, in this case, the range of `x` is extended
        by .1% on each side to include the min or max values of `x`. If
        `bins` is a sequence it defines the bin edges allowing for
        non-uniform bin width. No extension of the range of `x` is done in
        this case.
    right : bool, optional
        Indicates whether the bins include the rightmost edge or not. If
        right == True (the default), then the bins [1,2,3,4] indicate
        (1,2], (2,3], (3,4].
    labels : array or boolean, default None
        Used as labels for the resulting bins. Must be of the same length as the resulting
        bins. If False, return only integer indicators of the bins.
    retbins : bool, optional
        Whether to return the bins or not. Can be useful if bins is given
        as a scalar.
    precision : int, optional
        The precision at which to store and display the bins labels
    include_lowest : bool, optional
        Whether the first interval should be left-inclusive or not.

    Returns
    -------
    out : Categorical or Series or array of integers if labels is False
        The return type (Categorical or Series) depends on the input: a Series of type category if
        input is a Series else Categorical. Bins are represented as categories when categorical
        data is returned.
    bins : ndarray of floats
        Returned only if `retbins` is True.

    Notes
    -----
    The `cut` function can be useful for going from a continuous variable to
    a categorical variable. For example, `cut` could convert ages to groups
    of age ranges.

    Any NA values will be NA in the result.  Out of bounds values will be NA in
    the resulting Categorical object


    Examples
    --------
    >>> pd.cut(np.array([.2, 1.4, 2.5, 6.2, 9.7, 2.1]), 3, retbins=True)
    ([(0.191, 3.367], (0.191, 3.367], (0.191, 3.367], (3.367, 6.533], (6.533, 9.7], (0.191, 3.367]]
    Categories (3, object): [(0.191, 3.367] < (3.367, 6.533] < (6.533, 9.7]],
    array([ 0.1905    ,  3.36666667,  6.53333333,  9.7       ]))
    >>> pd.cut(np.array([.2, 1.4, 2.5, 6.2, 9.7, 2.1]), 3, labels=["good","medium","bad"])
    [good, good, good, medium, bad, good]
    Categories (3, object): [good < medium < bad]
    >>> pd.cut(np.ones(5), 4, labels=False)
    array([1, 1, 1, 1, 1], dtype=int64)
    """
    # NOTE: this binning code is changed a bit from histogram for var(x) == 0
    if not np.iterable(bins):
        if np.isscalar(bins) and bins < 1:
            raise ValueError("`bins` should be a positive integer.")
        try:  # for array-like
            sz = x.size
        except AttributeError:
            x = np.asarray(x)
            sz = x.size

        if sz == 0:
            raise ValueError('Cannot cut empty array')

        rng = (nanops.nanmin(x), nanops.nanmax(x))
        mn, mx = [mi + 0.0 for mi in rng]

        if mn == mx:  # adjust end points before binning
            mn -= .001 * mn
            mx += .001 * mx
            bins = np.linspace(mn, mx, bins + 1, endpoint=True)
        else:  # adjust end points after binning
            bins = np.linspace(mn, mx, bins + 1, endpoint=True)
            adj = (mx - mn) * 0.001  # 0.1% of the range
            if right:
                bins[0] -= adj
            else:
                bins[-1] += adj

    else:
        bins = np.asarray(bins)
        if (np.diff(bins) < 0).any():
            raise ValueError('bins must increase monotonically.')

    return _bins_to_cuts(x, bins, right=right, labels=labels, retbins=retbins,
                         precision=precision, include_lowest=include_lowest)


def qcut(x, q, labels=None, retbins=False, precision=3):
    """
    Quantile-based discretization function. Discretize variable into
    equal-sized buckets based on rank or based on sample quantiles. For example
    1000 values for 10 quantiles would produce a Categorical object indicating
    quantile membership for each data point.

    Parameters
    ----------
    x : ndarray or Series
    q : integer or array of quantiles
        Number of quantiles. 10 for deciles, 4 for quartiles, etc. Alternately
        array of quantiles, e.g. [0, .25, .5, .75, 1.] for quartiles
    labels : array or boolean, default None
        Used as labels for the resulting bins. Must be of the same length as the resulting
        bins. If False, return only integer indicators of the bins.
    retbins : bool, optional
        Whether to return the bins or not. Can be useful if bins is given
        as a scalar.
    precision : int, optional
        The precision at which to store and display the bins labels

    Returns
    -------
    out : Categorical or Series or array of integers if labels is False
        The return type (Categorical or Series) depends on the input: a Series of type category if
        input is a Series else Categorical. Bins are represented as categories when categorical
        data is returned.
    bins : ndarray of floats
        Returned only if `retbins` is True.

    Notes
    -----
    Out of bounds values will be NA in the resulting Categorical object

    Examples
    --------
    >>> pd.qcut(range(5), 4)
    [[0, 1], [0, 1], (1, 2], (2, 3], (3, 4]]
    Categories (4, object): [[0, 1] < (1, 2] < (2, 3] < (3, 4]]
    >>> pd.qcut(range(5), 3, labels=["good","medium","bad"])
    [good, good, medium, bad, bad]
    Categories (3, object): [good < medium < bad]
    >>> pd.qcut(range(5), 4, labels=False)
    array([0, 0, 1, 2, 3], dtype=int64)
    """
    if com.is_integer(q):
        quantiles = np.linspace(0, 1, q + 1)
    else:
        quantiles = q
    bins = algos.quantile(x, quantiles)
    return _bins_to_cuts(x, bins, labels=labels, retbins=retbins,
                         precision=precision, include_lowest=True)


def _bins_to_cuts(x, bins, right=True, labels=None, retbins=False,
                  precision=None, name=None, include_lowest=False):
    x_is_series = isinstance(x, Series)
    series_index = None

    if x_is_series:
        series_index = x.index
        if name is None:
            name = x.name

    x = np.asarray(x)

    side = 'left' if right else 'right'
    ids = bins.searchsorted(x, side=side)

    if len(algos.unique(bins)) < len(bins):
        raise ValueError('Bin edges must be unique: %s' % repr(bins))

    if include_lowest:
        ids[x == bins[0]] = 1

    na_mask = com.isnull(x) | (ids == len(bins)) | (ids == 0)
    has_nas = na_mask.any()

    if labels is not False:
        if labels is None:
            closed = 'right' if right else 'left'
            precision = _infer_precision(precision, bins)
            breaks = [_round_frac(b, precision) for b in bins]
            labels = IntervalIndex.from_breaks(breaks, closed=closed).values

            if right and include_lowest:
                labels[0] = Interval(labels[0].left, labels[0].right,
                                     closed='both')

        else:
            if len(labels) != len(bins) - 1:
                raise ValueError('Bin labels must be one fewer than '
                                 'the number of bin edges')

        if not com.is_categorical(labels):
            labels = np.asarray(labels)

        np.putmask(ids, na_mask, 0)
        result = com.take_nd(labels, ids - 1)

    else:
        result = ids - 1
        if has_nas:
            result = result.astype(np.float64)
            np.putmask(result, na_mask, np.nan)

    if x_is_series:
        result = Series(result, index=series_index, name=name)

    if not retbins:
        return result

    return result, bins


def _round_frac(x, precision):
    """Round the fractional part of the given number
    """
    if not np.isfinite(x) or x == 0:
        return x
    else:
        frac, whole = np.modf(x)
        if whole == 0:
            digits = -int(np.floor(np.log10(abs(frac)))) - 1 + precision
        else:
            digits = precision
        return np.around(x, digits)


def _infer_precision(base_precision, bins):
    """Infer an appropriate precision for _round_frac
    """
    for precision in range(base_precision, 20):
        levels = [_round_frac(b, precision) for b in bins]
        if algos.unique(levels).size == bins.size:
            return precision
    return base_precision  # default
