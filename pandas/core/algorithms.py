"""
Generic data algorithms. This module is experimental at the moment and not
intended for public consumption
"""

import numpy as np

import pandas.core.common as com
import pandas.lib as lib
import pandas._algos as _algos

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


def count(values, uniques=None):
    f = lambda htype, caster: _count_generic(values, htype, caster)

    if uniques is not None:
        raise NotImplementedError
    else:
        return _hashtable_algo(f, values.dtype)

def _hashtable_algo(f, dtype):
    """
    f(HashTable, type_caster) -> result
    """
    if com.is_float_dtype(dtype):
        return f(lib.Float64HashTable, com._ensure_float64)
    elif com.is_integer_dtype(dtype):
        return f(lib.Int64HashTable, com._ensure_int64)
    else:
        return f(lib.PyObjectHashTable, com._ensure_object)


def _count_generic(values, table_type, type_caster):
    from pandas.core.series import Series

    values = type_caster(values)
    table = table_type(min(len(values), 1000000))
    uniques, labels, counts = table.factorize(values)

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
    values = np.asarray(values)
    is_datetime = com.is_datetime64_dtype(values)
    hash_klass, values = _get_data_algo(values, _hashtables)

    uniques = []
    table = hash_klass(len(values))
    labels, counts = table.get_labels(values, uniques, 0, na_sentinel)

    labels = com._ensure_platform_int(labels)

    uniques = com._asarray_tuplesafe(uniques)
    if sort and len(counts) > 0:
        sorter = uniques.argsort()
        reverse_indexer = np.empty(len(sorter), dtype=np.int_)
        reverse_indexer.put(sorter, np.arange(len(sorter)))

        mask = labels < 0
        labels = reverse_indexer.take(labels)
        np.putmask(labels, mask, -1)

        uniques = uniques.take(sorter)
        counts = counts.take(sorter)

    if is_datetime:
        uniques = np.array(uniques, dtype='M8[ns]')

    return labels, uniques, counts

def value_counts(values, sort=True, ascending=False):
    """
    Compute a histogram of the counts of non-null values

    Parameters
    ----------
    values : ndarray (1-d)
    sort : boolean, default True
        Sort by values
    ascending : boolean, default False
        Sort in ascending order

    Returns
    -------
    value_counts : Series
    """
    from pandas.core.series import Series
    from collections import defaultdict

    values = np.asarray(values)

    if com.is_integer_dtype(values.dtype):
        values = com._ensure_int64(values)
        keys, counts = lib.value_count_int64(values)
        result = Series(counts, index=keys)
    else:
        counter = defaultdict(lambda: 0)
        values = values[com.notnull(values)]
        for value in values:
            counter[value] += 1
        result = Series(counter)

    if sort:
        result.sort()
        if not ascending:
            result = result[::-1]

    return result


def rank(values, axis=0, method='average', na_option='keep',
         ascending=True):
    """

    """
    if values.ndim == 1:
        f, values = _get_data_algo(values, _rank1d_functions)
        ranks = f(values, ties_method=method, ascending=ascending)
    elif values.ndim == 2:
        f, values = _get_data_algo(values, _rank2d_functions)
        ranks = f(values, axis=axis, ties_method=method,
                  ascending=ascending)
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
    a : ndarray
        Values from which to extract score.
    q : scalar or array
        Percentile at which to extract score.
    interpolation : {'fraction', 'lower', 'higher'}, optional
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
                raise ValueError("interpolation_method can only be 'fraction', " \
                                 "'lower' or 'higher'")

        return score

    if np.isscalar(q):
        return _get_score(q)
    else:
        q = np.asarray(q, np.float64)
        return _algos.arrmap_float64(q, _get_score)

def _interpolate(a, b, fraction):
    """Returns the point at the given fraction between a and b, where
    'fraction' must be between 0 and 1.
    """
    return a + (b - a)*fraction


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
    'float64' : lib.rank_1d_float64,
    'int64' : lib.rank_1d_int64,
    'generic' : lib.rank_1d_generic
}

_rank2d_functions = {
    'float64' : lib.rank_2d_float64,
    'int64' : lib.rank_2d_int64,
    'generic' : lib.rank_2d_generic
}

_hashtables = {
    'float64' : lib.Float64HashTable,
    'int64' : lib.Int64HashTable,
    'generic' : lib.PyObjectHashTable
}
