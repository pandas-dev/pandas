"""
Generic data algorithms. This module is experimental at the moment and not
intended for public consumption
"""

import numpy as np

from pandas.core.series import Series
import pandas.core.common as com
import pandas._tseries as lib

def match(values, index):
    """


    Parameters
    ----------

    Returns
    -------
    match : ndarray
    """
    if com.is_float_dtype(index):
        return _match_generic(values, index, lib.Float64HashTable,
                              com._ensure_float64)
    elif com.is_integer_dtype(index):
        return _match_generic(values, index, lib.Int64HashTable,
                              com._ensure_int64)
    else:
        return _match_generic(values, index, lib.PyObjectHashTable,
                              com._ensure_object)

def count(values, uniques=None):
    if uniques is not None:
        raise NotImplementedError
    else:
        if com.is_float_dtype(values):
            return _count_generic(values, lib.Float64HashTable,
                                  com._ensure_float64)
        elif com.is_integer_dtype(values):
            return _count_generic(values, lib.Int64HashTable,
                                  com._ensure_int64)
        else:
            return _count_generic(values, lib.PyObjectHashTable,
                                  com._ensure_object)

def _count_generic(values, table_type, type_caster):
    values = type_caster(values)
    table = table_type(len(values))
    uniques, labels, counts = table.factorize(values)

    return Series(counts, index=uniques)

def _match_generic(values, index, table_type, type_caster):
    values = type_caster(values)
    index = type_caster(index)
    table = table_type(len(index))
    table.map_locations(index)
    return table.lookup(values)

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
    hash_klass, values = _get_data_algo(values, _hashtables)

    uniques = []
    table = hash_klass(len(values))
    labels, counts = table.get_labels(values, uniques, 0, na_sentinel)

    uniques = com._asarray_tuplesafe(uniques)
    if sort and len(counts) > 0:
        sorter = uniques.argsort()
        reverse_indexer = np.empty(len(sorter), dtype=np.int32)
        reverse_indexer.put(sorter, np.arange(len(sorter)))

        mask = labels < 0
        labels = reverse_indexer.take(labels)
        np.putmask(labels, mask, -1)

        uniques = uniques.take(sorter)
        counts = counts.take(sorter)

    return labels, uniques, counts

def value_counts(values, sort=True, ascending=False):
    """
    Compute a histogram of the counts of non-null values

    Returns
    -------
    value_counts : Series
    """
    from collections import defaultdict
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


def _get_data_algo(values, func_map):
    if com.is_float_dtype(values):
        f = func_map['float64']
        values = com._ensure_float64(values)
    elif com.is_integer_dtype(values):
        f = func_map['int64']
        values = com._ensure_int64(values)
    else:
        f = func_map['generic']
        values = com._ensure_object(values)
    return f, values

_rank1d_functions = {
    'float64' : lib.rank_1d_float64,
    'int64' : lib.rank_1d_int64,
    'generic' : lib.rank_1d_generic
}

_rank2d_functions = {
    'float64' : lib.rank_2d_float64,
    'generic' : lib.rank_2d_generic
}

_hashtables = {
    'float64' : lib.Float64HashTable,
    'int64' : lib.Int64HashTable,
    'generic' : lib.PyObjectHashTable
}

def unique(values):
    """

    """
    pass
