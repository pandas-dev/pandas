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

def _get_hash_table_and_cast(values):
    if com.is_float_dtype(values):
        klass = lib.Float64HashTable
        values = com._ensure_float64(values)
    elif com.is_integer_dtype(values):
        klass = lib.Int64HashTable
        values = com._ensure_int64(values)
    else:
        klass = lib.PyObjectHashTable
        values = com._ensure_object(values)
    return klass, values

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
    hash_klass, values = _get_hash_table_and_cast(values)

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

def unique(values):
    """

    """
    pass
