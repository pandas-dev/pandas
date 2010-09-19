"""
Utilities for making working with rpy2 more user- and
developer-friendly.
"""

import numpy as np

from pandas import DataFrame, DataMatrix

from rpy2.robjects.packages import importr
from rpy2.robjects import r
import rpy2.robjects as robj

__all__ = ['convert_robj', 'load_data']

def load_data(name, package=None):
    if package:
        pack = importr(package)

    r.data(name)
    return convert_robj(r[name])

def _rclass(obj):
    """
    Return R class name for input object
    """
    return r['class'](obj)[0]

def _is_null(obj):
    return _rclass(obj) == 'NULL'

def _convert_list(obj):
    pass

def _convert_named_list(obj):
    pass

def _convert_DataFrame(rdf):
    columns = list(rdf.colnames)
    rows = np.array(rdf.rownames)

    data = {}
    for i, col in enumerate(columns):
        vec = rdf.rx2(i + 1)
        data[col] = list(vec)

    return DataFrame(data, index=rows)

def _convert_Matrix(mat):
    columns = mat.colnames
    rows = mat.rownames

    columns = None if _is_null(columns) else list(columns)
    index = None if _is_null(index) else list(index)

    return DataMatrix(np.array(mat), index=index, columns=columns)

def _check_int(vec):
    try:
        # R observation numbers come through as strings
        vec = vec.astype(int)
    except Exception:
        pass

    return vec

_converters = [
    (robj.DataFrame , _convert_DataFrame),
    (robj.Matrix , _convert_Matrix),
]

def convert_robj(obj):
    """
    Convert rpy2 object to a pandas-friendly form

    Parameters
    ----------
    obj : rpy2 object

    Returns
    -------
    Non-rpy data structure, mix of NumPy and pandas objects
    """
    if not isinstance(obj, orbj.RObjectMixin):
        return obj

    for rpy_type, converter in _converters:
        if isinstance(obj, rpy_type):
            return converter(obj)

    raise Exception('Do not know what to do with %s object' % klass)


import pandas.util.testing as _test

def test_convert_list():
    obj = r('list(a=1, b=2, c=3)')
    converted = convert_robj(obj)

    _test.assert_dict_equal

def test_convert_frame():
    # built-in dataset
    df = r['faithful']

    converted = convert_robj(obj)

def _named_matrix():
    r('mat <- matrix(rnorm(9), ncol=3)')
    r('colnames(mat) <- c("one", "two", "three")')
    r('rownames(mat) <- c("a", "b", "c")')

    return r['mat']

def test_convert_matrix():
    mat = _named_matrix()

    converted = convert_robj(mat)

    assert np.array_equal(converted.index, ['a', 'b', 'c'])
    assert np.array_equal(converted.columns, ['one', 'two', 'three'])

def test_convert_nested():
    pass


