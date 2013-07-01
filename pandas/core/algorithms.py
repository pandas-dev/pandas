"""
Generic data algorithms. This module is experimental at the moment and not
intended for public consumption
"""

import numpy as np

import pandas.core.common as com
import pandas.algos as algos
import pandas.hashtable as htable


def match(to_match, values, na_sentinel=-1):
    """
    Compute locations of to_match into values

    Parameters
    ----------
    to_match : array-like
        values to find positions of
    values : array-like
        Unique set of values
    na_sentinel : int, default -1
        Value to mark "not found"

    Examples
    --------

    Returns
    -------
    match : ndarray of integers
    """
    values = com._asarray_tuplesafe(values)
    if issubclass(values.dtype.type, basestring):
        values = np.array(values, dtype='O')

    f = lambda htype, caster: _match_generic(to_match, values, htype, caster)
    return _hashtable_algo(f, values.dtype)


def unique(values):
    """
    Compute unique values (not necessarily sorted) efficiently from input array
    of values

    Parameters
    ----------
    values : array-like

    Returns
    -------
    uniques
    """
    values = com._asarray_tuplesafe(values)
    f = lambda htype, caster: _unique_generic(values, htype, caster)
    return _hashtable_algo(f, values.dtype)


# def count(values, uniques=None):
#     f = lambda htype, caster: _count_generic(values, htype, caster)

#     if uniques is not None:
#         raise NotImplementedError
#     else:
#         return _hashtable_algo(f, values.dtype)


def _hashtable_algo(f, dtype):
    """
    f(HashTable, type_caster) -> result
    """
    if com.is_float_dtype(dtype):
        return f(htable.Float64HashTable, com._ensure_float64)
    elif com.is_integer_dtype(dtype):
        return f(htable.Int64HashTable, com._ensure_int64)
    else:
        return f(htable.PyObjectHashTable, com._ensure_object)


def _count_generic(values, table_type, type_caster):
    from pandas.core.series import Series

    values = type_caster(values)
    table = table_type(min(len(values), 1000000))
    uniques, labels = table.factorize(values)

    return Series(counts, index=uniques)


def _match_generic(values, index, table_type, type_caster):
    values = type_caster(values)
    index = type_caster(index)
    table = table_type(min(len(index), 1000000))
    table.map_locations(index)
    return table.lookup(values)


def _unique_generic(values, table_type, type_caster):
    values = type_caster(values)
    table = table_type(min(len(values), 1000000))
    uniques = table.unique(values)
    return type_caster(uniques)


def factorize(values, sort=False, order=None, na_sentinel=-1):
    """
    Encode input values as an enumerated type or categorical variable

    Parameters
    ----------
    values : sequence
    sort :
    order :

    Returns
    -------
    """
    from pandas.tseries.period import PeriodIndex
    vals = np.asarray(values)
    is_datetime = com.is_datetime64_dtype(vals)
    (hash_klass, vec_klass), vals = _get_data_algo(vals, _hashtables)

    table = hash_klass(len(vals))
    uniques = vec_klass()
    labels = table.get_labels(vals, uniques, 0, na_sentinel)

    labels = com._ensure_platform_int(labels)

    uniques = uniques.to_array()

    if sort and len(uniques) > 0:
        sorter = uniques.argsort()
        reverse_indexer = np.empty(len(sorter), dtype=np.int_)
        reverse_indexer.put(sorter, np.arange(len(sorter)))

        mask = labels < 0
        labels = reverse_indexer.take(labels)
        np.putmask(labels, mask, -1)

        uniques = uniques.take(sorter)

    if is_datetime:
        uniques = uniques.view('M8[ns]')
    if isinstance(values, PeriodIndex):
        uniques = PeriodIndex(ordinal=uniques, freq=values.freq)

    return labels, uniques


def value_counts(values, sort=True, ascending=False, normalize=False):
    """
    Compute a histogram of the counts of non-null values

    Parameters
    ----------
    values : ndarray (1-d)
    sort : boolean, default True
        Sort by values
    ascending : boolean, default False
        Sort in ascending order
    normalize: boolean, default False
        If True then compute a relative histogram

    Returns
    -------
    value_counts : Series
    """
    from pandas.core.series import Series

    values = np.asarray(values)

    if com.is_integer_dtype(values.dtype):
        values = com._ensure_int64(values)
        keys, counts = htable.value_count_int64(values)
    elif issubclass(values.dtype.type, (np.datetime64,np.timedelta64)):

        dtype = values.dtype
        values = values.view(np.int64)
        keys, counts = htable.value_count_int64(values)

        # convert the keys back to the dtype we came in
        keys = Series(keys,dtype=dtype)
    else:
        mask = com.isnull(values)
        values = com._ensure_object(values)
        keys, counts = htable.value_count_object(values, mask)

    result = Series(counts, index=keys)

    if sort:
        result.sort()
        if not ascending:
            result = result[::-1]

    if normalize:
        result = result / float(values.size)

    return result


def rank(values, axis=0, method='average', na_option='keep',
         ascending=True):
    """

    """
    if values.ndim == 1:
        f, values = _get_data_algo(values, _rank1d_functions)
        ranks = f(values, ties_method=method, ascending=ascending,
                  na_option=na_option)
    elif values.ndim == 2:
        f, values = _get_data_algo(values, _rank2d_functions)
        ranks = f(values, axis=axis, ties_method=method,
                  ascending=ascending, na_option=na_option)
    return ranks


def quantile(x, q, interpolation_method='fraction'):
    """
    Compute sample quantile or quantiles of the input array. For example, q=0.5
    computes the median.

    The `interpolation_method` parameter supports three values, namely
    `fraction` (default), `lower` and `higher`. Interpolation is done only,
    if the desired quantile lies between two data points `i` and `j`. For
    `fraction`, the result is an interpolated value between `i` and `j`;
    for `lower`, the result is `i`, for `higher` the result is `j`.

    Parameters
    ----------
    x : ndarray
        Values from which to extract score.
    q : scalar or array
        Percentile at which to extract score.
    interpolation_method : {'fraction', 'lower', 'higher'}, optional
        This optional parameter specifies the interpolation method to use,
        when the desired quantile lies between two data points `i` and `j`:

        - fraction: `i + (j - i)*fraction`, where `fraction` is the
                    fractional part of the index surrounded by `i` and `j`.
        -lower: `i`.
        - higher: `j`.

    Returns
    -------
    score : float
        Score at percentile.

    Examples
    --------
    >>> from scipy import stats
    >>> a = np.arange(100)
    >>> stats.scoreatpercentile(a, 50)
    49.5

    """
    x = np.asarray(x)
    mask = com.isnull(x)

    x = x[-mask]

    values = np.sort(x)

    def _get_score(at):
        if len(values) == 0:
            return np.nan

        idx = at * (len(values) - 1)
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
                raise ValueError("interpolation_method can only be 'fraction' "
                                 ", 'lower' or 'higher'")

        return score

    if np.isscalar(q):
        return _get_score(q)
    else:
        q = np.asarray(q, np.float64)
        return algos.arrmap_float64(q, _get_score)


def _interpolate(a, b, fraction):
    """Returns the point at the given fraction between a and b, where
    'fraction' must be between 0 and 1.
    """
    return a + (b - a) * fraction


def _get_data_algo(values, func_map):
    if com.is_float_dtype(values):
        f = func_map['float64']
        values = com._ensure_float64(values)
    elif com.is_datetime64_dtype(values):
        f = func_map['int64']
        values = values.view('i8')
    elif com.is_integer_dtype(values):
        f = func_map['int64']
        values = com._ensure_int64(values)
    else:
        f = func_map['generic']
        values = com._ensure_object(values)
    return f, values


def group_position(*args):
    """
    Get group position
    """
    from collections import defaultdict
    table = defaultdict(int)

    result = []
    for tup in zip(*args):
        result.append(table[tup])
        table[tup] += 1

    return result


_rank1d_functions = {
    'float64': algos.rank_1d_float64,
    'int64': algos.rank_1d_int64,
    'generic': algos.rank_1d_generic
}

_rank2d_functions = {
    'float64': algos.rank_2d_float64,
    'int64': algos.rank_2d_int64,
    'generic': algos.rank_2d_generic
}

_hashtables = {
    'float64': (htable.Float64HashTable, htable.Float64Vector),
    'int64': (htable.Int64HashTable, htable.Int64Vector),
    'generic': (htable.PyObjectHashTable, htable.ObjectVector)
}
