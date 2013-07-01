"""
Quantilization functions and related stuff
"""

from pandas.core.api import DataFrame, Series
from pandas.core.categorical import Categorical
from pandas.core.index import _ensure_index
import pandas.core.algorithms as algos
import pandas.core.common as com
import pandas.core.nanops as nanops

import numpy as np


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
        Labels to use for bin edges, or False to return integer bin labels
    retbins : bool, optional
        Whether to return the bins or not. Can be useful if bins is given
        as a scalar.

    Returns
    -------
    out : Categorical or array of integers if labels is False
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
    >>> cut(np.array([.2, 1.4, 2.5, 6.2, 9.7, 2.1]), 3, retbins=True)
    (array([(0.191, 3.367], (0.191, 3.367], (0.191, 3.367], (3.367, 6.533],
           (6.533, 9.7], (0.191, 3.367]], dtype=object),
     array([ 0.1905    ,  3.36666667,  6.53333333,  9.7       ]))
    >>> cut(np.ones(5), 4, labels=False)
    array([2, 2, 2, 2, 2])
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
            # handle empty arrays. Can't determine range, so use 0-1.
            # rng = (0, 1)
        else:
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

    return _bins_to_cuts(x, bins, right=right, labels=labels,
                         retbins=retbins, precision=precision,
                         include_lowest=include_lowest)


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
        Labels to use for bin edges, or False to return integer bin labels
    retbins : bool, optional
        Whether to return the bins or not. Can be useful if bins is given
        as a scalar.

    Returns
    -------
    cat : Categorical

    Notes
    -----
    Out of bounds values will be NA in the resulting Categorical object

    Examples
    --------
    """
    if com.is_integer(q):
        quantiles = np.linspace(0, 1, q + 1)
    else:
        quantiles = q
    bins = algos.quantile(x, quantiles)
    return _bins_to_cuts(x, bins, labels=labels, retbins=retbins,
                         precision=precision, include_lowest=True)


def _bins_to_cuts(x, bins, right=True, labels=None, retbins=False,
                  precision=3, name=None, include_lowest=False):
    if name is None and isinstance(x, Series):
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
            increases = 0
            while True:
                try:
                    levels = _format_levels(bins, precision, right=right,
                                            include_lowest=include_lowest)
                except ValueError:
                    increases += 1
                    precision += 1
                    if increases >= 20:
                        raise
                else:
                    break

        else:
            if len(labels) != len(bins) - 1:
                raise ValueError('Bin labels must be one fewer than '
                                 'the number of bin edges')
            levels = labels

        levels = np.asarray(levels, dtype=object)
        np.putmask(ids, na_mask, 0)
        fac = Categorical(ids - 1, levels, name=name)
    else:
        fac = ids - 1
        if has_nas:
            fac = fac.astype(np.float64)
            np.putmask(fac, na_mask, np.nan)

    if not retbins:
        return fac

    return fac, bins


def _format_levels(bins, prec, right=True,
                   include_lowest=False):
    fmt = lambda v: _format_label(v, precision=prec)
    if right:
        levels = []
        for a, b in zip(bins, bins[1:]):
            fa, fb = fmt(a), fmt(b)

            if a != b and fa == fb:
                raise ValueError('precision too low')

            formatted = '(%s, %s]' % (fa, fb)

            levels.append(formatted)

        if include_lowest:
            levels[0] = '[' + levels[0][1:]
    else:
        levels = ['[%s, %s)' % (fmt(a), fmt(b))
                  for a, b in zip(bins, bins[1:])]

    return levels


def _format_label(x, precision=3):
    fmt_str = '%%.%dg' % precision
    if com.is_float(x):
        frac, whole = np.modf(x)
        sgn = '-' if x < 0 else ''
        whole = abs(whole)
        if frac != 0.0:
            val = fmt_str % frac

            # rounded up or down
            if '.' not in val:
                if x < 0:
                    return '%d' % (-whole - 1)
                else:
                    return '%d' % (whole + 1)

            if 'e' in val:
                return _trim_zeros(fmt_str % x)
            else:
                val = _trim_zeros(val)
                if '.' in val:
                    return sgn + '.'.join(('%d' % whole, val.split('.')[1]))
                else:  # pragma: no cover
                    return sgn + '.'.join(('%d' % whole, val))
        else:
            return sgn + '%d' % whole
    else:
        return str(x)


def _trim_zeros(x):
    while len(x) > 1 and x[-1] == '0':
        x = x[:-1]
    if len(x) > 1 and x[-1] == '.':
        x = x[:-1]
    return x
