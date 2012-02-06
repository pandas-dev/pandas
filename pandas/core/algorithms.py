"""
Generic data algorithms
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
                              _ensure_float64)
    elif com.is_integer_dtype(index):
        return _match_generic(values, index, lib.Int64HashTable,
                              _ensure_int64)
    else:
        return _match_generic(values, index, lib.PyObjectHashTable,
                              _ensure_object)


def count(values, uniques=None):
    if uniques is not None:
        raise NotImplementedError
    else:
        if com.is_float_dtype(values):
            return _count_generic(values, lib.Float64HashTable,
                                  _ensure_float64)
        elif com.is_integer_dtype(values):
            return _count_generic(values, lib.Int64HashTable,
                                  _ensure_int64)
        else:
            return _count_generic(values, lib.PyObjectHashTable,
                                  _ensure_object)

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

def factorize(values):
    pass

def unique(values):
    pass

def _ensure_float64(arr):
    if arr.dtype != np.float64:
        arr = arr.astype(np.float64)
    return arr

def _ensure_int64(arr):
    if arr.dtype != np.int64:
        arr = arr.astype(np.int64)
    return arr

def _ensure_object(arr):
    if arr.dtype != np.object_:
        arr = arr.astype('O')
    return arr
