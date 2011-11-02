"""
Utilities for making working with rpy2 more user- and
developer-friendly.
"""

import numpy as np

import pandas as pn
import pandas.util.testing as _test

from rpy2.robjects.packages import importr
from rpy2.robjects import r
import rpy2.robjects as robj

__all__ = ['convert_robj', 'load_data']

def load_data(name, package=None, convert=True):
    if package:
        pack = importr(package)

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
            return pn.DataFrame(arr, index=name_list[0], columns=name_list[1])
        elif len(dim) == 3:
            return pn.Panel(arr, items=name_list[2],
                            major_axis=name_list[0],
                            minor_axis=name_list[1])
        else:
            print 'Cannot handle dim=%d' % len(dim)
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
            values = np.asarray(vec.levels).take(values - 1)

        data[col] = values

    return pn.DataFrame(data, index=_check_int(rows), columns=columns)

def _convert_Matrix(mat):
    columns = mat.colnames
    rows = mat.rownames

    columns = None if _is_null(columns) else list(columns)
    index = None if _is_null(rows) else list(rows)

    return pn.DataFrame(np.array(mat), index=_check_int(index),
                        columns=columns)

def _check_int(vec):
    try:
        # R observation numbers come through as strings
        vec = vec.astype(int)
    except Exception:
        pass

    return vec

_pandas_converters = [
    (robj.DataFrame , _convert_DataFrame),
    (robj.Matrix , _convert_Matrix),
    (robj.StrVector, _convert_vector),
    (robj.FloatVector, _convert_vector),
    (robj.Array, _convert_array),
    (robj.Vector, _convert_list),
]

_converters = [
    (robj.DataFrame , lambda x: _convert_DataFrame(x).toRecords(index=False)),
    (robj.Matrix , lambda x: _convert_Matrix(x).toRecords(index=False)),
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

def test_convert_list():
    obj = r('list(a=1, b=2, c=3)')

    converted = convert_robj(obj)
    expected = {'a' : [1], 'b' : [2], 'c' : [3]}

    _test.assert_dict_equal(converted, expected)

def test_convert_nested_list():
    obj = r('list(a=list(foo=1, bar=2))')

    converted = convert_robj(obj)
    expected = {'a' : {'foo' : [1], 'bar' : [2]}}

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


if __name__ == '__main__':
    pass
