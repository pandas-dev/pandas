"""
Utilities for making working with rpy2 more user- and
developer-friendly.
"""
from __future__ import print_function

from pandas.compat import zip, range
import numpy as np

import pandas as pd
import pandas.core.common as com
import pandas.util.testing as _test

from rpy2.robjects.packages import importr
from rpy2.robjects import r
import rpy2.robjects as robj

import itertools as IT


__all__ = ['convert_robj', 'load_data', 'convert_to_r_dataframe',
           'convert_to_r_matrix']


def load_data(name, package=None, convert=True):
    if package:
        importr(package)

    r.data(name)

    robj = r[name]

    if convert:
        return convert_robj(robj)
    else:
        return robj


def _rclass(obj):
    """
    Return R class name for input object
    """
    return r['class'](obj)[0]


def _is_null(obj):
    return _rclass(obj) == 'NULL'


def _convert_list(obj):
    """
    Convert named Vector to dict, factors to list
    """
    try:
        values = [convert_robj(x) for x in obj]
        keys = r['names'](obj)
        return dict(zip(keys, values))
    except TypeError:
        # For state.division and state.region
        factors = list(r['factor'](obj))
        level = list(r['levels'](obj))
        result = [level[index-1] for index in factors]
        return result


def _convert_array(obj):
    """
    Convert Array to DataFrame
    """
    def _list(item):
        try:
            return list(item)
        except TypeError:
            return []
        
    # For iris3, HairEyeColor, UCBAdmissions, Titanic
    dim = list(obj.dim)
    values = np.array(list(obj))
    names = r['dimnames'](obj)
    try:
        columns = list(r['names'](names))[::-1]
    except TypeError:
        columns = ['X{:d}'.format(i) for i in range(len(names))][::-1]
    columns.append('value')
    name_list = [(_list(x) or range(d)) for x, d in zip(names, dim)][::-1]
    arr = np.array(list(IT.product(*name_list)))
    arr = np.column_stack([arr,values])
    df = pd.DataFrame(arr, columns=columns)
    return df


def _convert_vector(obj):
    if isinstance(obj, robj.IntVector):
        return _convert_int_vector(obj)
    elif isinstance(obj, robj.StrVector):
        return _convert_str_vector(obj)
    # Check if the vector has extra information attached to it that can be used
    # as an index
    try:
        attributes = set(r['attributes'](obj).names)
    except AttributeError:
        return list(obj)
    if 'names' in attributes:
        return pd.Series(list(obj), index=r['names'](obj)) 
    elif 'tsp' in attributes:
        return pd.Series(list(obj), index=r['time'](obj)) 
    elif 'labels' in attributes:
        return pd.Series(list(obj), index=r['labels'](obj))
    if _rclass(obj) == 'dist':
        # For 'eurodist'. WARNING: This results in a DataFrame, not a Series or list.
        matrix = r['as.matrix'](obj)
        return convert_robj(matrix)
    else:
        return list(obj)

NA_INTEGER = -2147483648


def _convert_int_vector(obj):
    arr = np.asarray(obj)
    mask = arr == NA_INTEGER
    if mask.any():
        arr = arr.astype(float)
        arr[mask] = np.nan
    return arr


def _convert_str_vector(obj):
    arr = np.asarray(obj, dtype=object)
    mask = arr == robj.NA_Character
    if mask.any():
        arr[mask] = np.nan
    return arr


def _convert_DataFrame(rdf):
    columns = list(rdf.colnames)
    rows = np.array(rdf.rownames)

    data = {}
    for i, col in enumerate(columns):
        vec = rdf.rx2(i + 1)
        values = _convert_vector(vec)

        if isinstance(vec, robj.FactorVector):
            levels = np.asarray(vec.levels)
            if com.is_float_dtype(values):
                mask = np.isnan(values)
                notmask = -mask
                result = np.empty(len(values), dtype=object)
                result[mask] = np.nan

                locs = (values[notmask] - 1).astype(np.int_)
                result[notmask] = levels.take(locs)
                values = result
            else:
                values = np.asarray(vec.levels).take(values - 1)

        data[col] = values

    return pd.DataFrame(data, index=_check_int(rows), columns=columns)


def _convert_Matrix(mat):
    columns = mat.colnames
    rows = mat.rownames

    columns = None if _is_null(columns) else list(columns)
    index = r['time'](mat) if _is_null(rows) else list(rows)
    return pd.DataFrame(np.array(mat), index=_check_int(index),
                        columns=columns)


def _check_int(vec):
    try:
        # R observation numbers come through as strings
        vec = vec.astype(int)
    except Exception:
        pass

    return vec

_pandas_converters = [
    (robj.DataFrame, _convert_DataFrame),
    (robj.Matrix, _convert_Matrix),
    (robj.StrVector, _convert_vector),
    (robj.FloatVector, _convert_vector),
    (robj.Array, _convert_array),
    (robj.Vector, _convert_list),
]

_converters = [
    (robj.DataFrame, lambda x: _convert_DataFrame(x).toRecords(index=False)),
    (robj.Matrix, lambda x: _convert_Matrix(x).toRecords(index=False)),
    (robj.IntVector, _convert_vector),
    (robj.StrVector, _convert_vector),
    (robj.FloatVector, _convert_vector),
    (robj.Array, _convert_array),
    (robj.Vector, _convert_list),
]


def convert_robj(obj, use_pandas=True):
    """
    Convert rpy2 object to a pandas-friendly form

    Parameters
    ----------
    obj : rpy2 object

    Returns
    -------
    Non-rpy data structure, mix of NumPy and pandas objects
    """
    if not isinstance(obj, robj.RObjectMixin):
        return obj

    converters = _pandas_converters if use_pandas else _converters

    for rpy_type, converter in converters:
        if isinstance(obj, rpy_type):
            return converter(obj)

    raise TypeError('Do not know what to do with %s object' % type(obj))


def convert_to_r_posixct(obj):
    """
    Convert DatetimeIndex or np.datetime array to R POSIXct using
    m8[s] format.

    Parameters
    ----------
    obj : source pandas object (one of [DatetimeIndex, np.datetime])

    Returns
    -------
    An R POSIXct vector (rpy2.robjects.vectors.POSIXct)

    """
    import time
    from rpy2.rinterface import StrSexpVector

    # convert m8[ns] to m8[s]
    vals = robj.vectors.FloatSexpVector(obj.values.view('i8') / 1E9)
    as_posixct = robj.baseenv.get('as.POSIXct')
    origin = StrSexpVector([time.strftime("%Y-%m-%d",
                                          time.gmtime(0)), ])

    # We will be sending ints as UTC
    tz = obj.tz.zone if hasattr(
        obj, 'tz') and hasattr(obj.tz, 'zone') else 'UTC'
    tz = StrSexpVector([tz])
    utc_tz = StrSexpVector(['UTC'])

    posixct = as_posixct(vals, origin=origin, tz=utc_tz)
    posixct.do_slot_assign('tzone', tz)
    return posixct


VECTOR_TYPES = {np.float64: robj.FloatVector,
                np.float32: robj.FloatVector,
                np.float: robj.FloatVector,
                np.int: robj.IntVector,
                np.int32: robj.IntVector,
                np.int64: robj.IntVector,
                np.object_: robj.StrVector,
                np.str: robj.StrVector,
                np.bool: robj.BoolVector}

NA_TYPES = {np.float64: robj.NA_Real,
            np.float32: robj.NA_Real,
            np.float: robj.NA_Real,
            np.int: robj.NA_Integer,
            np.int32: robj.NA_Integer,
            np.int64: robj.NA_Integer,
            np.object_: robj.NA_Character,
            np.str: robj.NA_Character,
            np.bool: robj.NA_Logical}


def convert_to_r_dataframe(df, strings_as_factors=False):
    """
    Convert a pandas DataFrame to a R data.frame.

    Parameters
    ----------
    df: The DataFrame being converted
    strings_as_factors: Whether to turn strings into R factors (default: False)

    Returns
    -------
    A R data.frame

    """

    import rpy2.rlike.container as rlc

    columns = rlc.OrdDict()

    # FIXME: This doesn't handle MultiIndex

    for column in df:
        value = df[column]
        value_type = value.dtype.type

        if value_type == np.datetime64:
            value = convert_to_r_posixct(value)
        else:
            value = [item if pd.notnull(item) else NA_TYPES[value_type]
                     for item in value]

            value = VECTOR_TYPES[value_type](value)

            if not strings_as_factors:
                I = robj.baseenv.get("I")
                value = I(value)

        columns[column] = value

    r_dataframe = robj.DataFrame(columns)

    del columns

    r_dataframe.rownames = robj.StrVector(df.index)

    return r_dataframe


def convert_to_r_matrix(df, strings_as_factors=False):

    """
    Convert a pandas DataFrame to a R matrix.

    Parameters
    ----------
    df: The DataFrame being converted
    strings_as_factors: Whether to turn strings into R factors (default: False)

    Returns
    -------
    A R matrix

    """

    if df._is_mixed_type:
        raise TypeError("Conversion to matrix only possible with non-mixed "
                        "type DataFrames")

    r_dataframe = convert_to_r_dataframe(df, strings_as_factors)
    as_matrix = robj.baseenv.get("as.matrix")
    r_matrix = as_matrix(r_dataframe)

    return r_matrix

if __name__ == '__main__':
    pass
