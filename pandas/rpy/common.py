"""
Utilities for making working with rpy2 more user- and
developer-friendly.
"""

import numpy as np

import pandas as pd
import pandas.core.common as com
import pandas.util.testing as _test

from rpy2.robjects.packages import importr
from rpy2.robjects import r
import rpy2.robjects as robj

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
    Convert named Vector to dict
    """
    values = [convert_robj(x) for x in obj]
    return dict(zip(obj.names, values))


def _convert_array(obj):
    """
    Convert Array to ndarray
    """
    # this royally sucks. "Matrices" (arrays) with dimension > 3 in R aren't
    # really matrices-- things come out Fortran order in the first two
    # dimensions. Maybe I'm wrong?

    dim = list(obj.dim)
    values = np.array(list(obj))

    if len(dim) == 3:
        arr = values.reshape(dim[-1:] + dim[:-1]).swapaxes(1, 2)

    if obj.names is not None:
        name_list = [list(x) for x in obj.names]
        if len(dim) == 2:
            return pd.DataFrame(arr, index=name_list[0], columns=name_list[1])
        elif len(dim) == 3:
            return pd.Panel(arr, items=name_list[2],
                            major_axis=name_list[0],
                            minor_axis=name_list[1])
        else:
            print ('Cannot handle dim=%d' % len(dim))
    else:
        return arr


def _convert_vector(obj):
    if isinstance(obj, robj.IntVector):
        return _convert_int_vector(obj)
    elif isinstance(obj, robj.StrVector):
        return _convert_str_vector(obj)

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
    index = None if _is_null(rows) else list(rows)

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

    raise Exception('Do not know what to do with %s object' % type(obj))


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


def test_convert_list():
    obj = r('list(a=1, b=2, c=3)')

    converted = convert_robj(obj)
    expected = {'a': [1], 'b': [2], 'c': [3]}

    _test.assert_dict_equal(converted, expected)


def test_convert_nested_list():
    obj = r('list(a=list(foo=1, bar=2))')

    converted = convert_robj(obj)
    expected = {'a': {'foo': [1], 'bar': [2]}}

    _test.assert_dict_equal(converted, expected)


def test_convert_frame():
    # built-in dataset
    df = r['faithful']

    converted = convert_robj(df)

    assert np.array_equal(converted.columns, ['eruptions', 'waiting'])
    assert np.array_equal(converted.index, np.arange(1, 273))


def _test_matrix():
    r('mat <- matrix(rnorm(9), ncol=3)')
    r('colnames(mat) <- c("one", "two", "three")')
    r('rownames(mat) <- c("a", "b", "c")')

    return r['mat']


def test_convert_matrix():
    mat = _test_matrix()

    converted = convert_robj(mat)

    assert np.array_equal(converted.index, ['a', 'b', 'c'])
    assert np.array_equal(converted.columns, ['one', 'two', 'three'])


def test_convert_r_dataframe():

    is_na = robj.baseenv.get("is.na")

    seriesd = _test.getSeriesData()
    frame = pd.DataFrame(seriesd, columns=['D', 'C', 'B', 'A'])

    # Null data
    frame["E"] = [np.nan for item in frame["A"]]
    # Some mixed type data
    frame["F"] = ["text" if item % 2 == 0 else np.nan for item in range(30)]

    r_dataframe = convert_to_r_dataframe(frame)

    assert np.array_equal(convert_robj(r_dataframe.rownames), frame.index)
    assert np.array_equal(convert_robj(r_dataframe.colnames), frame.columns)
    assert all(is_na(item) for item in r_dataframe.rx2("E"))

    for column in frame[["A", "B", "C", "D"]]:
        coldata = r_dataframe.rx2(column)
        original_data = frame[column]
        assert np.array_equal(convert_robj(coldata), original_data)

    for column in frame[["D", "E"]]:
        for original, converted in zip(frame[column],
                                       r_dataframe.rx2(column)):

            if pd.isnull(original):
                assert is_na(converted)
            else:
                assert original == converted


def test_convert_r_matrix():

    is_na = robj.baseenv.get("is.na")

    seriesd = _test.getSeriesData()
    frame = pd.DataFrame(seriesd, columns=['D', 'C', 'B', 'A'])
    # Null data
    frame["E"] = [np.nan for item in frame["A"]]

    r_dataframe = convert_to_r_matrix(frame)

    assert np.array_equal(convert_robj(r_dataframe.rownames), frame.index)
    assert np.array_equal(convert_robj(r_dataframe.colnames), frame.columns)
    assert all(is_na(item) for item in r_dataframe.rx(True, "E"))

    for column in frame[["A", "B", "C", "D"]]:
        coldata = r_dataframe.rx(True, column)
        original_data = frame[column]
        assert np.array_equal(convert_robj(coldata),
                              original_data)

    # Pandas bug 1282
    frame["F"] = ["text" if item % 2 == 0 else np.nan for item in range(30)]

    # FIXME: Ugly, this whole module needs to be ported to nose/unittest
    try:
        wrong_matrix = convert_to_r_matrix(frame)
    except TypeError:
        pass
    except Exception:
        raise


if __name__ == '__main__':
    pass
