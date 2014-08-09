"""
Misc tools for implementing data structures
"""

import re
import collections
import numbers
import codecs
import csv
import types
from datetime import datetime, timedelta
from functools import partial

from numpy.lib.format import read_array, write_array
import numpy as np

import pandas as pd
import pandas.algos as algos
import pandas.lib as lib
import pandas.tslib as tslib
from pandas import compat
from pandas.compat import StringIO, BytesIO, range, long, u, zip, map, string_types, iteritems
from pandas.core.dtypes import CategoricalDtype, CategoricalDtypeType, DatetimeTZDtype, DatetimeTZDtypeType
from pandas.core.config import get_option

class PandasError(Exception):
    pass


class SettingWithCopyError(ValueError):
    pass


class SettingWithCopyWarning(Warning):
    pass


class AmbiguousIndexError(PandasError, KeyError):
    pass


class AbstractMethodError(NotImplementedError):
    """Raise this error instead of NotImplementedError for abstract methods
    while keeping compatibility with Python 2 and Python 3.
    """
    def __init__(self, class_instance):
        self.class_instance = class_instance

    def __str__(self):
        return "This method must be defined in the concrete class of " \
               + self.class_instance.__class__.__name__

_POSSIBLY_CAST_DTYPES = set([np.dtype(t).name
                             for t in ['O', 'int8',
                                       'uint8', 'int16', 'uint16', 'int32',
                                       'uint32', 'int64', 'uint64']])

_NS_DTYPE = np.dtype('M8[ns]')
_TD_DTYPE = np.dtype('m8[ns]')
_INT64_DTYPE = np.dtype(np.int64)
_DATELIKE_DTYPES = set([np.dtype(t) for t in ['M8[ns]', '<M8[ns]', '>M8[ns]',
                                              'm8[ns]', '<m8[ns]', '>m8[ns]']])
_int8_max = np.iinfo(np.int8).max
_int16_max = np.iinfo(np.int16).max
_int32_max = np.iinfo(np.int32).max
_int64_max = np.iinfo(np.int64).max

# define abstract base classes to enable isinstance type checking on our
# objects
def create_pandas_abc_type(name, attr, comp):
    @classmethod
    def _check(cls, inst):
        return getattr(inst, attr, '_typ') in comp
    dct = dict(__instancecheck__=_check,
               __subclasscheck__=_check)
    meta = type("ABCBase", (type,), dct)
    return meta(name, tuple(), dct)


ABCIndex = create_pandas_abc_type("ABCIndex", "_typ", ("index",))
ABCInt64Index = create_pandas_abc_type("ABCInt64Index", "_typ", ("int64index",))
ABCFloat64Index = create_pandas_abc_type("ABCFloat64Index", "_typ", ("float64index",))
ABCMultiIndex = create_pandas_abc_type("ABCMultiIndex", "_typ", ("multiindex",))
ABCDatetimeIndex = create_pandas_abc_type("ABCDatetimeIndex", "_typ", ("datetimeindex",))
ABCTimedeltaIndex = create_pandas_abc_type("ABCTimedeltaIndex", "_typ", ("timedeltaindex",))
ABCPeriodIndex = create_pandas_abc_type("ABCPeriodIndex", "_typ", ("periodindex",))
ABCCategoricalIndex = create_pandas_abc_type("ABCCategoricalIndex", "_typ", ("categoricalindex",))
ABCIndexClass = create_pandas_abc_type("ABCIndexClass", "_typ", ("index",
                                                                 "int64index",
                                                                 "float64index",
                                                                 "multiindex",
                                                                 "datetimeindex",
                                                                 "timedeltaindex",
                                                                 "periodindex",
                                                                 "categoricalindex"))

ABCSeries = create_pandas_abc_type("ABCSeries", "_typ", ("series",))
ABCDataFrame = create_pandas_abc_type("ABCDataFrame", "_typ", ("dataframe",))
ABCPanel = create_pandas_abc_type("ABCPanel", "_typ", ("panel",))
ABCSparseSeries = create_pandas_abc_type("ABCSparseSeries", "_subtyp",
                                         ('sparse_series',
                                          'sparse_time_series'))
ABCSparseArray = create_pandas_abc_type("ABCSparseArray", "_subtyp",
                                        ('sparse_array', 'sparse_series'))
ABCCategorical = create_pandas_abc_type("ABCCategorical","_typ",("categorical"))
ABCPeriod = create_pandas_abc_type("ABCPeriod", "_typ", ("period",))

class _ABCGeneric(type):

    def __instancecheck__(cls, inst):
        return hasattr(inst, "_data")


ABCGeneric = _ABCGeneric("ABCGeneric", tuple(), {})


def bind_method(cls, name, func):
    """Bind a method to class, python 2 and python 3 compatible.

    Parameters
    ----------

    cls : type
        class to receive bound method
    name : basestring
        name of method on class instance
    func : function
        function to be bound as method


    Returns
    -------
    None
    """
    # only python 2 has bound/unbound method issue
    if not compat.PY3:
        setattr(cls, name, types.MethodType(func, None, cls))
    else:
        setattr(cls, name, func)

def isnull(obj):
    """Detect missing values (NaN in numeric arrays, None/NaN in object arrays)

    Parameters
    ----------
    arr : ndarray or object value
        Object to check for null-ness

    Returns
    -------
    isnulled : array-like of bool or bool
        Array or bool indicating whether an object is null or if an array is
        given which of the element is null.

    See also
    --------
    pandas.notnull: boolean inverse of pandas.isnull
    """
    return _isnull(obj)


def _isnull_new(obj):
    if lib.isscalar(obj):
        return lib.checknull(obj)
    # hack (for now) because MI registers as ndarray
    elif isinstance(obj, pd.MultiIndex):
        raise NotImplementedError("isnull is not defined for MultiIndex")
    elif isinstance(obj, (ABCSeries, np.ndarray, pd.Index)):
        return _isnull_ndarraylike(obj)
    elif isinstance(obj, ABCGeneric):
        return obj._constructor(obj._data.isnull(func=isnull))
    elif isinstance(obj, list) or hasattr(obj, '__array__'):
        return _isnull_ndarraylike(np.asarray(obj))
    else:
        return obj is None


def _isnull_old(obj):
    """Detect missing values. Treat None, NaN, INF, -INF as null.

    Parameters
    ----------
    arr: ndarray or object value

    Returns
    -------
    boolean ndarray or boolean
    """
    if lib.isscalar(obj):
        return lib.checknull_old(obj)
    # hack (for now) because MI registers as ndarray
    elif isinstance(obj, pd.MultiIndex):
        raise NotImplementedError("isnull is not defined for MultiIndex")
    elif isinstance(obj, (ABCSeries, np.ndarray, pd.Index)):
        return _isnull_ndarraylike_old(obj)
    elif isinstance(obj, ABCGeneric):
        return obj._constructor(obj._data.isnull(func=_isnull_old))
    elif isinstance(obj, list) or hasattr(obj, '__array__'):
        return _isnull_ndarraylike_old(np.asarray(obj))
    else:
        return obj is None

_isnull = _isnull_new


def _use_inf_as_null(key):
    """Option change callback for null/inf behaviour
    Choose which replacement for numpy.isnan / -numpy.isfinite is used.

    Parameters
    ----------
    flag: bool
        True means treat None, NaN, INF, -INF as null (old way),
        False means None and NaN are null, but INF, -INF are not null
        (new way).

    Notes
    -----
    This approach to setting global module values is discussed and
    approved here:

    * http://stackoverflow.com/questions/4859217/
      programmatically-creating-variables-in-python/4859312#4859312
    """
    flag = get_option(key)
    if flag:
        globals()['_isnull'] = _isnull_old
    else:
        globals()['_isnull'] = _isnull_new


def _isnull_ndarraylike(obj):

    values = getattr(obj, 'values', obj)
    dtype = values.dtype

    if dtype.kind in ('O', 'S', 'U'):
        if is_categorical_dtype(values):
            from pandas import Categorical
            if not isinstance(values, Categorical):
                values = values.values
            result = values.isnull()
        else:

            # Working around NumPy ticket 1542
            shape = values.shape

            if dtype.kind in ('S', 'U'):
                result = np.zeros(values.shape, dtype=bool)
            else:
                result = np.empty(shape, dtype=bool)
                vec = lib.isnullobj(values.ravel())
                result[...] = vec.reshape(shape)

    elif is_datetimelike(obj):
        # this is the NaT pattern
        result = values.view('i8') == tslib.iNaT
    else:
        result = np.isnan(values)

    # box
    if isinstance(obj, ABCSeries):
        from pandas import Series
        result = Series(result, index=obj.index, name=obj.name, copy=False)

    return result

def _isnull_ndarraylike_old(obj):
    values = getattr(obj, 'values', obj)
    dtype = values.dtype

    if dtype.kind in ('O', 'S', 'U'):
        # Working around NumPy ticket 1542
        shape = values.shape

        if values.dtype.kind in ('S', 'U'):
            result = np.zeros(values.shape, dtype=bool)
        else:
            result = np.empty(shape, dtype=bool)
            vec = lib.isnullobj_old(values.ravel())
            result[:] = vec.reshape(shape)

    elif dtype in _DATELIKE_DTYPES:
        # this is the NaT pattern
        result = values.view('i8') == tslib.iNaT
    else:
        result = ~np.isfinite(values)

    # box
    if isinstance(obj, ABCSeries):
        from pandas import Series
        result = Series(result, index=obj.index, name=obj.name, copy=False)

    return result


def notnull(obj):
    """Replacement for numpy.isfinite / -numpy.isnan which is suitable for use
    on object arrays.

    Parameters
    ----------
    arr : ndarray or object value
        Object to check for *not*-null-ness

    Returns
    -------
    isnulled : array-like of bool or bool
        Array or bool indicating whether an object is *not* null or if an array
        is given which of the element is *not* null.

    See also
    --------
    pandas.isnull : boolean inverse of pandas.notnull
    """
    res = isnull(obj)
    if np.isscalar(res):
        return not res
    return ~res

def is_null_datelike_scalar(other):
    """ test whether the object is a null datelike, e.g. Nat
    but guard against passing a non-scalar """
    if other is pd.NaT or other is None:
        return True
    elif np.isscalar(other):

        # a timedelta
        if hasattr(other,'dtype'):
            return other.view('i8') == tslib.iNaT
        elif is_integer(other) and other == tslib.iNaT:
            return True
        return isnull(other)
    return False

def array_equivalent(left, right, strict_nan=False):
    """
    True if two arrays, left and right, have equal non-NaN elements, and NaNs in
    corresponding locations.  False otherwise. It is assumed that left and right
    are NumPy arrays of the same dtype. The behavior of this function
    (particularly with respect to NaNs) is not defined if the dtypes are
    different.

    Parameters
    ----------
    left, right : ndarrays
    strict_nan : bool, default False
        If True, consider NaN and None to be different.

    Returns
    -------
    b : bool
        Returns True if the arrays are equivalent.

    Examples
    --------
    >>> array_equivalent(
    ...     np.array([1, 2, np.nan]),
    ...     np.array([1, 2, np.nan]))
    True
    >>> array_equivalent(
    ...     np.array([1, np.nan, 2]),
    ...     np.array([1, 2, np.nan]))
    False
    """

    left, right = np.asarray(left), np.asarray(right)
    if left.shape != right.shape: return False

    # Object arrays can contain None, NaN and NaT.
    if issubclass(left.dtype.type, np.object_) or issubclass(right.dtype.type, np.object_):

        if not strict_nan:
            # pd.isnull considers NaN and None to be equivalent.
            return lib.array_equivalent_object(_ensure_object(left.ravel()),
                                               _ensure_object(right.ravel()))

        for left_value, right_value in zip(left, right):
            if left_value is tslib.NaT and right_value is not tslib.NaT:
                return False

            elif isinstance(left_value, float) and np.isnan(left_value):
                if not isinstance(right_value, float) or not np.isnan(right_value):
                    return False
            else:
                if left_value != right_value:
                    return False
        return True

    # NaNs can occur in float and complex arrays.
    if issubclass(left.dtype.type, (np.floating, np.complexfloating)):
        return ((left == right) | (np.isnan(left) & np.isnan(right))).all()

    # numpy will will not allow this type of datetimelike vs integer comparison
    elif is_datetimelike_v_numeric(left, right):
        return False

    # NaNs cannot occur otherwise.
    return np.array_equal(left, right)

def _iterable_not_string(x):
    return (isinstance(x, collections.Iterable) and
            not isinstance(x, compat.string_types))


def flatten(l):
    """Flatten an arbitrarily nested sequence.

    Parameters
    ----------
    l : sequence
        The non string sequence to flatten

    Notes
    -----
    This doesn't consider strings sequences.

    Returns
    -------
    flattened : generator
    """
    for el in l:
        if _iterable_not_string(el):
            for s in flatten(el):
                yield s
        else:
            yield el


def mask_missing(arr, values_to_mask):
    """
    Return a masking array of same size/shape as arr
    with entries equaling any member of values_to_mask set to True
    """
    if not isinstance(values_to_mask, (list, np.ndarray)):
        values_to_mask = [values_to_mask]

    try:
        values_to_mask = np.array(values_to_mask, dtype=arr.dtype)
    except Exception:
        values_to_mask = np.array(values_to_mask, dtype=object)

    na_mask = isnull(values_to_mask)
    nonna = values_to_mask[~na_mask]

    mask = None
    for x in nonna:
        if mask is None:
            mask = arr == x

            # if x is a string and arr is not, then we get False and we must
            # expand the mask to size arr.shape
            if np.isscalar(mask):
                mask = np.zeros(arr.shape, dtype=bool)
        else:
            mask |= arr == x

    if na_mask.any():
        if mask is None:
            mask = isnull(arr)
        else:
            mask |= isnull(arr)

    return mask


def _pickle_array(arr):
    arr = arr.view(np.ndarray)

    buf = BytesIO()
    write_array(buf, arr)

    return buf.getvalue()


def _unpickle_array(bytes):
    arr = read_array(BytesIO(bytes))

    # All datetimes should be stored as M8[ns].  When unpickling with
    # numpy1.6, it will read these as M8[us].  So this ensures all
    # datetime64 types are read as MS[ns]
    if is_datetime64_dtype(arr):
        arr = arr.view(_NS_DTYPE)

    return arr


def _view_wrapper(f, arr_dtype=None, out_dtype=None, fill_wrap=None):
    def wrapper(arr, indexer, out, fill_value=np.nan):
        if arr_dtype is not None:
            arr = arr.view(arr_dtype)
        if out_dtype is not None:
            out = out.view(out_dtype)
        if fill_wrap is not None:
            fill_value = fill_wrap(fill_value)
        f(arr, indexer, out, fill_value=fill_value)
    return wrapper


def _convert_wrapper(f, conv_dtype):
    def wrapper(arr, indexer, out, fill_value=np.nan):
        arr = arr.astype(conv_dtype)
        f(arr, indexer, out, fill_value=fill_value)
    return wrapper


def _take_2d_multi_generic(arr, indexer, out, fill_value, mask_info):
    # this is not ideal, performance-wise, but it's better than raising
    # an exception (best to optimize in Cython to avoid getting here)
    row_idx, col_idx = indexer
    if mask_info is not None:
        (row_mask, col_mask), (row_needs, col_needs) = mask_info
    else:
        row_mask = row_idx == -1
        col_mask = col_idx == -1
        row_needs = row_mask.any()
        col_needs = col_mask.any()
    if fill_value is not None:
        if row_needs:
            out[row_mask, :] = fill_value
        if col_needs:
            out[:, col_mask] = fill_value
    for i in range(len(row_idx)):
        u_ = row_idx[i]
        for j in range(len(col_idx)):
            v = col_idx[j]
            out[i, j] = arr[u_, v]


def _take_nd_generic(arr, indexer, out, axis, fill_value, mask_info):
    if mask_info is not None:
        mask, needs_masking = mask_info
    else:
        mask = indexer == -1
        needs_masking = mask.any()
    if arr.dtype != out.dtype:
        arr = arr.astype(out.dtype)
    if arr.shape[axis] > 0:
        arr.take(_ensure_platform_int(indexer), axis=axis, out=out)
    if needs_masking:
        outindexer = [slice(None)] * arr.ndim
        outindexer[axis] = mask
        out[tuple(outindexer)] = fill_value


_take_1d_dict = {
    ('int8', 'int8'): algos.take_1d_int8_int8,
    ('int8', 'int32'): algos.take_1d_int8_int32,
    ('int8', 'int64'): algos.take_1d_int8_int64,
    ('int8', 'float64'): algos.take_1d_int8_float64,
    ('int16', 'int16'): algos.take_1d_int16_int16,
    ('int16', 'int32'): algos.take_1d_int16_int32,
    ('int16', 'int64'): algos.take_1d_int16_int64,
    ('int16', 'float64'): algos.take_1d_int16_float64,
    ('int32', 'int32'): algos.take_1d_int32_int32,
    ('int32', 'int64'): algos.take_1d_int32_int64,
    ('int32', 'float64'): algos.take_1d_int32_float64,
    ('int64', 'int64'): algos.take_1d_int64_int64,
    ('int64', 'float64'): algos.take_1d_int64_float64,
    ('float32', 'float32'): algos.take_1d_float32_float32,
    ('float32', 'float64'): algos.take_1d_float32_float64,
    ('float64', 'float64'): algos.take_1d_float64_float64,
    ('object', 'object'): algos.take_1d_object_object,
    ('bool', 'bool'):
    _view_wrapper(algos.take_1d_bool_bool, np.uint8, np.uint8),
    ('bool', 'object'):
    _view_wrapper(algos.take_1d_bool_object, np.uint8, None),
    ('datetime64[ns]', 'datetime64[ns]'):
    _view_wrapper(algos.take_1d_int64_int64, np.int64, np.int64, np.int64)
}


_take_2d_axis0_dict = {
    ('int8', 'int8'): algos.take_2d_axis0_int8_int8,
    ('int8', 'int32'): algos.take_2d_axis0_int8_int32,
    ('int8', 'int64'): algos.take_2d_axis0_int8_int64,
    ('int8', 'float64'): algos.take_2d_axis0_int8_float64,
    ('int16', 'int16'): algos.take_2d_axis0_int16_int16,
    ('int16', 'int32'): algos.take_2d_axis0_int16_int32,
    ('int16', 'int64'): algos.take_2d_axis0_int16_int64,
    ('int16', 'float64'): algos.take_2d_axis0_int16_float64,
    ('int32', 'int32'): algos.take_2d_axis0_int32_int32,
    ('int32', 'int64'): algos.take_2d_axis0_int32_int64,
    ('int32', 'float64'): algos.take_2d_axis0_int32_float64,
    ('int64', 'int64'): algos.take_2d_axis0_int64_int64,
    ('int64', 'float64'): algos.take_2d_axis0_int64_float64,
    ('float32', 'float32'): algos.take_2d_axis0_float32_float32,
    ('float32', 'float64'): algos.take_2d_axis0_float32_float64,
    ('float64', 'float64'): algos.take_2d_axis0_float64_float64,
    ('object', 'object'): algos.take_2d_axis0_object_object,
    ('bool', 'bool'):
    _view_wrapper(algos.take_2d_axis0_bool_bool, np.uint8, np.uint8),
    ('bool', 'object'):
    _view_wrapper(algos.take_2d_axis0_bool_object, np.uint8, None),
    ('datetime64[ns]', 'datetime64[ns]'):
    _view_wrapper(algos.take_2d_axis0_int64_int64, np.int64, np.int64,
                  fill_wrap=np.int64)
}


_take_2d_axis1_dict = {
    ('int8', 'int8'): algos.take_2d_axis1_int8_int8,
    ('int8', 'int32'): algos.take_2d_axis1_int8_int32,
    ('int8', 'int64'): algos.take_2d_axis1_int8_int64,
    ('int8', 'float64'): algos.take_2d_axis1_int8_float64,
    ('int16', 'int16'): algos.take_2d_axis1_int16_int16,
    ('int16', 'int32'): algos.take_2d_axis1_int16_int32,
    ('int16', 'int64'): algos.take_2d_axis1_int16_int64,
    ('int16', 'float64'): algos.take_2d_axis1_int16_float64,
    ('int32', 'int32'): algos.take_2d_axis1_int32_int32,
    ('int32', 'int64'): algos.take_2d_axis1_int32_int64,
    ('int32', 'float64'): algos.take_2d_axis1_int32_float64,
    ('int64', 'int64'): algos.take_2d_axis1_int64_int64,
    ('int64', 'float64'): algos.take_2d_axis1_int64_float64,
    ('float32', 'float32'): algos.take_2d_axis1_float32_float32,
    ('float32', 'float64'): algos.take_2d_axis1_float32_float64,
    ('float64', 'float64'): algos.take_2d_axis1_float64_float64,
    ('object', 'object'): algos.take_2d_axis1_object_object,
    ('bool', 'bool'):
    _view_wrapper(algos.take_2d_axis1_bool_bool, np.uint8, np.uint8),
    ('bool', 'object'):
    _view_wrapper(algos.take_2d_axis1_bool_object, np.uint8, None),
    ('datetime64[ns]', 'datetime64[ns]'):
    _view_wrapper(algos.take_2d_axis1_int64_int64, np.int64, np.int64,
                  fill_wrap=np.int64)
}


_take_2d_multi_dict = {
    ('int8', 'int8'): algos.take_2d_multi_int8_int8,
    ('int8', 'int32'): algos.take_2d_multi_int8_int32,
    ('int8', 'int64'): algos.take_2d_multi_int8_int64,
    ('int8', 'float64'): algos.take_2d_multi_int8_float64,
    ('int16', 'int16'): algos.take_2d_multi_int16_int16,
    ('int16', 'int32'): algos.take_2d_multi_int16_int32,
    ('int16', 'int64'): algos.take_2d_multi_int16_int64,
    ('int16', 'float64'): algos.take_2d_multi_int16_float64,
    ('int32', 'int32'): algos.take_2d_multi_int32_int32,
    ('int32', 'int64'): algos.take_2d_multi_int32_int64,
    ('int32', 'float64'): algos.take_2d_multi_int32_float64,
    ('int64', 'int64'): algos.take_2d_multi_int64_int64,
    ('int64', 'float64'): algos.take_2d_multi_int64_float64,
    ('float32', 'float32'): algos.take_2d_multi_float32_float32,
    ('float32', 'float64'): algos.take_2d_multi_float32_float64,
    ('float64', 'float64'): algos.take_2d_multi_float64_float64,
    ('object', 'object'): algos.take_2d_multi_object_object,
    ('bool', 'bool'):
    _view_wrapper(algos.take_2d_multi_bool_bool, np.uint8, np.uint8),
    ('bool', 'object'):
    _view_wrapper(algos.take_2d_multi_bool_object, np.uint8, None),
    ('datetime64[ns]', 'datetime64[ns]'):
    _view_wrapper(algos.take_2d_multi_int64_int64, np.int64, np.int64,
                  fill_wrap=np.int64)
}


def _get_take_nd_function(ndim, arr_dtype, out_dtype, axis=0, mask_info=None):
    if ndim <= 2:
        tup = (arr_dtype.name, out_dtype.name)
        if ndim == 1:
            func = _take_1d_dict.get(tup, None)
        elif ndim == 2:
            if axis == 0:
                func = _take_2d_axis0_dict.get(tup, None)
            else:
                func = _take_2d_axis1_dict.get(tup, None)
        if func is not None:
            return func

        tup = (out_dtype.name, out_dtype.name)
        if ndim == 1:
            func = _take_1d_dict.get(tup, None)
        elif ndim == 2:
            if axis == 0:
                func = _take_2d_axis0_dict.get(tup, None)
            else:
                func = _take_2d_axis1_dict.get(tup, None)
        if func is not None:
            func = _convert_wrapper(func, out_dtype)
            return func

    def func(arr, indexer, out, fill_value=np.nan):
        indexer = _ensure_int64(indexer)
        _take_nd_generic(arr, indexer, out, axis=axis,
                         fill_value=fill_value, mask_info=mask_info)
    return func


def take_nd(arr, indexer, axis=0, out=None, fill_value=np.nan,
            mask_info=None, allow_fill=True):
    """
    Specialized Cython take which sets NaN values in one pass

    Parameters
    ----------
    arr : ndarray
        Input array
    indexer : ndarray
        1-D array of indices to take, subarrays corresponding to -1 value
        indicies are filed with fill_value
    axis : int, default 0
        Axis to take from
    out : ndarray or None, default None
        Optional output array, must be appropriate type to hold input and
        fill_value together, if indexer has any -1 value entries; call
        common._maybe_promote to determine this type for any fill_value
    fill_value : any, default np.nan
        Fill value to replace -1 values with
    mask_info : tuple of (ndarray, boolean)
        If provided, value should correspond to:
            (indexer != -1, (indexer != -1).any())
        If not provided, it will be computed internally if necessary
    allow_fill : boolean, default True
        If False, indexer is assumed to contain no -1 values so no filling
        will be done.  This short-circuits computation of a mask.  Result is
        undefined if allow_fill == False and -1 is present in indexer.
    """

    # dispatch to internal type takes
    if is_categorical(arr):
        return arr.take_nd(indexer, fill_value=fill_value,
                           allow_fill=allow_fill)
    elif is_datetimetz(arr):
        return arr.take(indexer, fill_value=fill_value, allow_fill=allow_fill)

    if indexer is None:
        indexer = np.arange(arr.shape[axis], dtype=np.int64)
        dtype, fill_value = arr.dtype, arr.dtype.type()
    else:
        indexer = _ensure_int64(indexer)
        if not allow_fill:
            dtype, fill_value = arr.dtype, arr.dtype.type()
            mask_info = None, False
        else:
            # check for promotion based on types only (do this first because
            # it's faster than computing a mask)
            dtype, fill_value = _maybe_promote(arr.dtype, fill_value)
            if dtype != arr.dtype and (out is None or out.dtype != dtype):
                # check if promotion is actually required based on indexer
                if mask_info is not None:
                    mask, needs_masking = mask_info
                else:
                    mask = indexer == -1
                    needs_masking = mask.any()
                    mask_info = mask, needs_masking
                if needs_masking:
                    if out is not None and out.dtype != dtype:
                        raise TypeError('Incompatible type for fill_value')
                else:
                    # if not, then depromote, set fill_value to dummy
                    # (it won't be used but we don't want the cython code
                    # to crash when trying to cast it to dtype)
                    dtype, fill_value = arr.dtype, arr.dtype.type()

    flip_order = False
    if arr.ndim == 2:
        if arr.flags.f_contiguous:
            flip_order = True

    if flip_order:
        arr = arr.T
        axis = arr.ndim - axis - 1
        if out is not None:
            out = out.T

    # at this point, it's guaranteed that dtype can hold both the arr values
    # and the fill_value
    if out is None:
        out_shape = list(arr.shape)
        out_shape[axis] = len(indexer)
        out_shape = tuple(out_shape)
        if arr.flags.f_contiguous and axis == arr.ndim - 1:
            # minor tweak that can make an order-of-magnitude difference
            # for dataframes initialized directly from 2-d ndarrays
            # (s.t. df.values is c-contiguous and df._data.blocks[0] is its
            # f-contiguous transpose)
            out = np.empty(out_shape, dtype=dtype, order='F')
        else:
            out = np.empty(out_shape, dtype=dtype)

    func = _get_take_nd_function(arr.ndim, arr.dtype, out.dtype,
                                 axis=axis, mask_info=mask_info)
    indexer = _ensure_int64(indexer)
    func(arr, indexer, out, fill_value)

    if flip_order:
        out = out.T
    return out


take_1d = take_nd


def take_2d_multi(arr, indexer, out=None, fill_value=np.nan,
                  mask_info=None, allow_fill=True):
    """
    Specialized Cython take which sets NaN values in one pass
    """
    if indexer is None or (indexer[0] is None and indexer[1] is None):
        row_idx = np.arange(arr.shape[0], dtype=np.int64)
        col_idx = np.arange(arr.shape[1], dtype=np.int64)
        indexer = row_idx, col_idx
        dtype, fill_value = arr.dtype, arr.dtype.type()
    else:
        row_idx, col_idx = indexer
        if row_idx is None:
            row_idx = np.arange(arr.shape[0], dtype=np.int64)
        else:
            row_idx = _ensure_int64(row_idx)
        if col_idx is None:
            col_idx = np.arange(arr.shape[1], dtype=np.int64)
        else:
            col_idx = _ensure_int64(col_idx)
        indexer = row_idx, col_idx
        if not allow_fill:
            dtype, fill_value = arr.dtype, arr.dtype.type()
            mask_info = None, False
        else:
            # check for promotion based on types only (do this first because
            # it's faster than computing a mask)
            dtype, fill_value = _maybe_promote(arr.dtype, fill_value)
            if dtype != arr.dtype and (out is None or out.dtype != dtype):
                # check if promotion is actually required based on indexer
                if mask_info is not None:
                    (row_mask, col_mask), (row_needs, col_needs) = mask_info
                else:
                    row_mask = row_idx == -1
                    col_mask = col_idx == -1
                    row_needs = row_mask.any()
                    col_needs = col_mask.any()
                    mask_info = (row_mask, col_mask), (row_needs, col_needs)
                if row_needs or col_needs:
                    if out is not None and out.dtype != dtype:
                        raise TypeError('Incompatible type for fill_value')
                else:
                    # if not, then depromote, set fill_value to dummy
                    # (it won't be used but we don't want the cython code
                    # to crash when trying to cast it to dtype)
                    dtype, fill_value = arr.dtype, arr.dtype.type()

    # at this point, it's guaranteed that dtype can hold both the arr values
    # and the fill_value
    if out is None:
        out_shape = len(row_idx), len(col_idx)
        out = np.empty(out_shape, dtype=dtype)

    func = _take_2d_multi_dict.get((arr.dtype.name, out.dtype.name), None)
    if func is None and arr.dtype != out.dtype:
        func = _take_2d_multi_dict.get((out.dtype.name, out.dtype.name), None)
        if func is not None:
            func = _convert_wrapper(func, out.dtype)
    if func is None:
        def func(arr, indexer, out, fill_value=np.nan):
            _take_2d_multi_generic(arr, indexer, out,
                                   fill_value=fill_value, mask_info=mask_info)
    func(arr, indexer, out=out, fill_value=fill_value)
    return out

_diff_special = {
    'float64': algos.diff_2d_float64,
    'float32': algos.diff_2d_float32,
    'int64': algos.diff_2d_int64,
    'int32': algos.diff_2d_int32,
    'int16': algos.diff_2d_int16,
    'int8': algos.diff_2d_int8,
}

def diff(arr, n, axis=0):
    """ difference of n between self,
        analagoust to s-s.shift(n) """

    n = int(n)
    na = np.nan
    dtype = arr.dtype
    is_timedelta = False
    if needs_i8_conversion(arr):
        dtype = np.float64
        arr = arr.view('i8')
        na = tslib.iNaT
        is_timedelta = True
    elif issubclass(dtype.type, np.integer):
        dtype = np.float64
    elif issubclass(dtype.type, np.bool_):
        dtype = np.object_

    dtype = np.dtype(dtype)
    out_arr = np.empty(arr.shape, dtype=dtype)

    na_indexer = [slice(None)] * arr.ndim
    na_indexer[axis] = slice(None, n) if n >= 0 else slice(n, None)
    out_arr[tuple(na_indexer)] = na

    if arr.ndim == 2 and arr.dtype.name in _diff_special:
        f = _diff_special[arr.dtype.name]
        f(arr, out_arr, n, axis)
    else:
        res_indexer = [slice(None)] * arr.ndim
        res_indexer[axis] = slice(n, None) if n >= 0 else slice(None, n)
        res_indexer = tuple(res_indexer)

        lag_indexer = [slice(None)] * arr.ndim
        lag_indexer[axis] = slice(None, -n) if n > 0 else slice(-n, None)
        lag_indexer = tuple(lag_indexer)

        # need to make sure that we account for na for datelike/timedelta
        # we don't actually want to subtract these i8 numbers
        if is_timedelta:
            res = arr[res_indexer]
            lag = arr[lag_indexer]

            mask = (arr[res_indexer] == na) | (arr[lag_indexer] == na)
            if mask.any():
                res = res.copy()
                res[mask] = 0
                lag = lag.copy()
                lag[mask] = 0

            result = res - lag
            result[mask] = na
            out_arr[res_indexer] = result
        else:
            out_arr[res_indexer] = arr[res_indexer] - arr[lag_indexer]

    if is_timedelta:
        from pandas import TimedeltaIndex
        out_arr = TimedeltaIndex(out_arr.ravel().astype('int64')).asi8.reshape(out_arr.shape).astype('timedelta64[ns]')

    return out_arr

def _coerce_indexer_dtype(indexer, categories):
    """ coerce the indexer input array to the smallest dtype possible """
    l = len(categories)
    if l < _int8_max:
        return _ensure_int8(indexer)
    elif l < _int16_max:
        return _ensure_int16(indexer)
    elif l < _int32_max:
        return _ensure_int32(indexer)
    return _ensure_int64(indexer)

def _coerce_to_dtypes(result, dtypes):
    """ given a dtypes and a result set, coerce the result elements to the
    dtypes
    """
    if len(result) != len(dtypes):
        raise AssertionError("_coerce_to_dtypes requires equal len arrays")

    from pandas.tseries.timedeltas import _coerce_scalar_to_timedelta_type

    def conv(r, dtype):
        try:
            if isnull(r):
                pass
            elif dtype == _NS_DTYPE:
                r = lib.Timestamp(r)
            elif dtype == _TD_DTYPE:
                r = _coerce_scalar_to_timedelta_type(r)
            elif dtype == np.bool_:
                # messy. non 0/1 integers do not get converted.
                if is_integer(r) and r not in [0,1]:
                    return int(r)
                r = bool(r)
            elif dtype.kind == 'f':
                r = float(r)
            elif dtype.kind == 'i':
                r = int(r)
        except:
            pass

        return r

    return [conv(r, dtype) for r, dtype in zip(result, dtypes)]


def _infer_fill_value(val):
    """
    infer the fill value for the nan/NaT from the provided scalar/ndarray/list-like
    if we are a NaT, return the correct dtyped element to provide proper block construction

    """

    if not is_list_like(val):
        val = [val]
    val = np.array(val,copy=False)
    if is_datetimelike(val):
        return np.array('NaT',dtype=val.dtype)
    elif is_object_dtype(val.dtype):
        dtype = lib.infer_dtype(_ensure_object(val))
        if dtype in ['datetime','datetime64']:
            return np.array('NaT',dtype=_NS_DTYPE)
        elif dtype in ['timedelta','timedelta64']:
            return np.array('NaT',dtype=_TD_DTYPE)
    return np.nan


def _infer_dtype_from_scalar(val):
    """ interpret the dtype from a scalar, upcast floats and ints
        return the new value and the dtype """

    dtype = np.object_

    # a 1-element ndarray
    if isinstance(val, np.ndarray):
        if val.ndim != 0:
            raise ValueError(
                "invalid ndarray passed to _infer_dtype_from_scalar")

        dtype = val.dtype
        val = val.item()

    elif isinstance(val, compat.string_types):

        # If we create an empty array using a string to infer
        # the dtype, NumPy will only allocate one character per entry
        # so this is kind of bad. Alternately we could use np.repeat
        # instead of np.empty (but then you still don't want things
        # coming out as np.str_!

        dtype = np.object_

    elif isinstance(val, (np.datetime64, datetime)) and getattr(val,'tz',None) is None:
        val = lib.Timestamp(val).value
        dtype = np.dtype('M8[ns]')

    elif isinstance(val, (np.timedelta64, timedelta)):
        val = tslib.convert_to_timedelta(val,'ns')
        dtype = np.dtype('m8[ns]')

    elif is_bool(val):
        dtype = np.bool_

    # provide implicity upcast on scalars
    elif is_integer(val):
        dtype = np.int64

    elif is_float(val):
        dtype = np.float64

    elif is_complex(val):
        dtype = np.complex_

    return dtype, val


def _maybe_promote(dtype, fill_value=np.nan):

    # if we passed an array here, determine the fill value by dtype
    if isinstance(fill_value, np.ndarray):
        if issubclass(fill_value.dtype.type, (np.datetime64, np.timedelta64)):
            fill_value = tslib.iNaT
        else:

            # we need to change to object type as our
            # fill_value is of object type
            if fill_value.dtype == np.object_:
                dtype = np.dtype(np.object_)
            fill_value = np.nan

    # returns tuple of (dtype, fill_value)
    if issubclass(dtype.type, (np.datetime64, np.timedelta64)):
        # for now: refuse to upcast datetime64
        # (this is because datetime64 will not implicitly upconvert
        #  to object correctly as of numpy 1.6.1)
        if isnull(fill_value):
            fill_value = tslib.iNaT
        else:
            if issubclass(dtype.type, np.datetime64):
                try:
                    fill_value = lib.Timestamp(fill_value).value
                except:
                    # the proper thing to do here would probably be to upcast
                    # to object (but numpy 1.6.1 doesn't do this properly)
                    fill_value = tslib.iNaT
            else:
                fill_value = tslib.iNaT
    elif is_float(fill_value):
        if issubclass(dtype.type, np.bool_):
            dtype = np.object_
        elif issubclass(dtype.type, np.integer):
            dtype = np.float64
    elif is_bool(fill_value):
        if not issubclass(dtype.type, np.bool_):
            dtype = np.object_
    elif is_integer(fill_value):
        if issubclass(dtype.type, np.bool_):
            dtype = np.object_
        elif issubclass(dtype.type, np.integer):
            # upcast to prevent overflow
            arr = np.asarray(fill_value)
            if arr != arr.astype(dtype):
                dtype = arr.dtype
    elif is_complex(fill_value):
        if issubclass(dtype.type, np.bool_):
            dtype = np.object_
        elif issubclass(dtype.type, (np.integer, np.floating)):
            dtype = np.complex128
    else:
        dtype = np.object_

    # in case we have a string that looked like a number
    if is_categorical_dtype(dtype):
        dtype = dtype
    elif issubclass(np.dtype(dtype).type, compat.string_types):
        dtype = np.object_

    return dtype, fill_value


def _maybe_upcast_putmask(result, mask, other):
    """
    A safe version of putmask that potentially upcasts the result

    Parameters
    ----------
    result : ndarray
        The destination array. This will be mutated in-place if no upcasting is
        necessary.
    mask : boolean ndarray
    other : ndarray or scalar
        The source array or value

    Returns
    -------
    result : ndarray
    changed : boolean
        Set to true if the result array was upcasted
    """

    if mask.any():
        # Two conversions for date-like dtypes that can't be done automatically
        # in np.place:
        #   NaN -> NaT
        #   integer or integer array -> date-like array
        if result.dtype in _DATELIKE_DTYPES:
            if lib.isscalar(other):
                if isnull(other):
                    other = result.dtype.type('nat')
                elif is_integer(other):
                    other = np.array(other, dtype=result.dtype)
            elif is_integer_dtype(other):
                other = np.array(other, dtype=result.dtype)

        def changeit():

            # try to directly set by expanding our array to full
            # length of the boolean
            try:
                om = other[mask]
                om_at = om.astype(result.dtype)
                if (om == om_at).all():
                    new_result = result.values.copy()
                    new_result[mask] = om_at
                    result[:] = new_result
                    return result, False
            except:
                pass

            # we are forced to change the dtype of the result as the input
            # isn't compatible
            r, _ = _maybe_upcast(result, fill_value=other, copy=True)
            np.place(r, mask, other)

            return r, True

        # we want to decide whether place will work
        # if we have nans in the False portion of our mask then we need to
        # upcast (possibly), otherwise we DON't want to upcast (e.g. if we
        # have values, say integers, in the success portion then it's ok to not
        # upcast)
        new_dtype, _ = _maybe_promote(result.dtype, other)
        if new_dtype != result.dtype:

            # we have a scalar or len 0 ndarray
            # and its nan and we are changing some values
            if (np.isscalar(other) or
                    (isinstance(other, np.ndarray) and other.ndim < 1)):
                if isnull(other):
                    return changeit()

            # we have an ndarray and the masking has nans in it
            else:

                if isnull(other[mask]).any():
                    return changeit()

        try:
            np.place(result, mask, other)
        except:
            return changeit()

    return result, False


def _maybe_upcast(values, fill_value=np.nan, dtype=None, copy=False):
    """ provide explict type promotion and coercion

    Parameters
    ----------
    values : the ndarray that we want to maybe upcast
    fill_value : what we want to fill with
    dtype : if None, then use the dtype of the values, else coerce to this type
    copy : if True always make a copy even if no upcast is required
    """

    if is_internal_type(values):
        if copy:
            values = values.copy()
    else:
        if dtype is None:
            dtype = values.dtype
        new_dtype, fill_value = _maybe_promote(dtype, fill_value)
        if new_dtype != values.dtype:
            values = values.astype(new_dtype)
        elif copy:
            values = values.copy()

    return values, fill_value


def _possibly_cast_item(obj, item, dtype):
    chunk = obj[item]

    if chunk.values.dtype != dtype:
        if dtype in (np.object_, np.bool_):
            obj[item] = chunk.astype(np.object_)
        elif not issubclass(dtype, (np.integer, np.bool_)):  # pragma: no cover
            raise ValueError("Unexpected dtype encountered: %s" % dtype)


def _possibly_downcast_to_dtype(result, dtype):
    """ try to cast to the specified dtype (e.g. convert back to bool/int
    or could be an astype of float64->float32
    """

    if np.isscalar(result):
        return result

    trans = lambda x: x
    if isinstance(dtype, compat.string_types):
        if dtype == 'infer':
            inferred_type = lib.infer_dtype(_ensure_object(result.ravel()))
            if inferred_type == 'boolean':
                dtype = 'bool'
            elif inferred_type == 'integer':
                dtype = 'int64'
            elif inferred_type == 'datetime64':
                dtype = 'datetime64[ns]'
            elif inferred_type == 'timedelta64':
                dtype = 'timedelta64[ns]'

            # try to upcast here
            elif inferred_type == 'floating':
                dtype = 'int64'
                if issubclass(result.dtype.type, np.number):
                    trans = lambda x: x.round()

            else:
                dtype = 'object'

    if isinstance(dtype, compat.string_types):
        dtype = np.dtype(dtype)

    try:

        # don't allow upcasts here (except if empty)
        if dtype.kind == result.dtype.kind:
            if result.dtype.itemsize <= dtype.itemsize and np.prod(result.shape):
                return result

        if issubclass(dtype.type, np.floating):
            return result.astype(dtype)
        elif dtype == np.bool_ or issubclass(dtype.type, np.integer):

            # if we don't have any elements, just astype it
            if not np.prod(result.shape):
                return trans(result).astype(dtype)

            # do a test on the first element, if it fails then we are done
            r = result.ravel()
            arr = np.array([r[0]])
            if not np.allclose(arr, trans(arr).astype(dtype)):
                return result

            # a comparable, e.g. a Decimal may slip in here
            elif not isinstance(r[0], (np.integer, np.floating, np.bool, int,
                                       float, bool)):
                return result

            if (issubclass(result.dtype.type, (np.object_, np.number)) and
                    notnull(result).all()):
                new_result = trans(result).astype(dtype)
                try:
                    if np.allclose(new_result, result):
                        return new_result
                except:

                    # comparison of an object dtype with a number type could
                    # hit here
                    if (new_result == result).all():
                        return new_result

        # a datetimelike
        elif dtype.kind in ['M','m'] and result.dtype.kind in ['i']:
            try:
                result = result.astype(dtype)
            except:
                pass

    except:
        pass

    return result


def _maybe_convert_string_to_object(values):
    """
    Convert string-like and string-like array to convert object dtype.
    This is to avoid numpy to handle the array as str dtype.
    """
    if isinstance(values, string_types):
        values = np.array([values], dtype=object)
    elif (isinstance(values, np.ndarray) and
        issubclass(values.dtype.type, (np.string_, np.unicode_))):
        values = values.astype(object)
    return values


def _lcd_dtypes(a_dtype, b_dtype):
    """ return the lcd dtype to hold these types """

    if is_datetime64_dtype(a_dtype) or is_datetime64_dtype(b_dtype):
        return _NS_DTYPE
    elif is_timedelta64_dtype(a_dtype) or is_timedelta64_dtype(b_dtype):
        return _TD_DTYPE
    elif is_complex_dtype(a_dtype):
        if is_complex_dtype(b_dtype):
            return a_dtype
        return np.float64
    elif is_integer_dtype(a_dtype):
        if is_integer_dtype(b_dtype):
            if a_dtype.itemsize == b_dtype.itemsize:
                return a_dtype
            return np.int64
        return np.float64
    elif is_float_dtype(a_dtype):
        if is_float_dtype(b_dtype):
            if a_dtype.itemsize == b_dtype.itemsize:
                return a_dtype
            else:
                return np.float64
        elif is_integer(b_dtype):
            return np.float64
    return np.object


def _fill_zeros(result, x, y, name, fill):
    """
    if this is a reversed op, then flip x,y

    if we have an integer value (or array in y)
    and we have 0's, fill them with the fill,
    return the result

    mask the nan's from x
    """
    if fill is None or is_float_dtype(result):
        return result

    if name.startswith(('r', '__r')):
        x,y = y,x

    is_typed_variable = (hasattr(y, 'dtype') or hasattr(y,'type'))
    is_scalar = lib.isscalar(y)

    if not is_typed_variable and not is_scalar:
        return result

    if is_scalar:
        y = np.array(y)

    if is_integer_dtype(y):

        if (y == 0).any():

            # GH 7325, mask and nans must be broadcastable (also: PR 9308)
            # Raveling and then reshaping makes np.putmask faster
            mask = ((y == 0) & ~np.isnan(result)).ravel()

            shape = result.shape
            result = result.astype('float64', copy=False).ravel()

            np.putmask(result, mask, fill)

            # if we have a fill of inf, then sign it correctly
            # (GH 6178 and PR 9308)
            if np.isinf(fill):
                signs = np.sign(y if name.startswith(('r', '__r')) else x)
                negative_inf_mask = (signs.ravel() < 0) & mask
                np.putmask(result, negative_inf_mask, -fill)

            if "floordiv" in name:  # (PR 9308)
                nan_mask = ((y == 0) & (x == 0)).ravel()
                np.putmask(result, nan_mask, np.nan)

            result = result.reshape(shape)

    return result


def _interp_wrapper(f, wrap_dtype, na_override=None):
    def wrapper(arr, mask, limit=None):
        view = arr.view(wrap_dtype)
        f(view, mask, limit=limit)
    return wrapper


_pad_1d_datetime = _interp_wrapper(algos.pad_inplace_int64, np.int64)
_pad_2d_datetime = _interp_wrapper(algos.pad_2d_inplace_int64, np.int64)
_backfill_1d_datetime = _interp_wrapper(algos.backfill_inplace_int64,
                                        np.int64)
_backfill_2d_datetime = _interp_wrapper(algos.backfill_2d_inplace_int64,
                                        np.int64)


def pad_1d(values, limit=None, mask=None, dtype=None):

    if dtype is None:
        dtype = values.dtype
    _method = None
    if is_float_dtype(values):
        _method = getattr(algos, 'pad_inplace_%s' % dtype.name, None)
    elif dtype in _DATELIKE_DTYPES or is_datetime64_dtype(values):
        _method = _pad_1d_datetime
    elif is_integer_dtype(values):
        values = _ensure_float64(values)
        _method = algos.pad_inplace_float64
    elif values.dtype == np.object_:
        _method = algos.pad_inplace_object

    if _method is None:
        raise ValueError('Invalid dtype for pad_1d [%s]' % dtype.name)

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)
    _method(values, mask, limit=limit)
    return values


def backfill_1d(values, limit=None, mask=None, dtype=None):

    if dtype is None:
        dtype = values.dtype
    _method = None
    if is_float_dtype(values):
        _method = getattr(algos, 'backfill_inplace_%s' % dtype.name, None)
    elif dtype in _DATELIKE_DTYPES or is_datetime64_dtype(values):
        _method = _backfill_1d_datetime
    elif is_integer_dtype(values):
        values = _ensure_float64(values)
        _method = algos.backfill_inplace_float64
    elif values.dtype == np.object_:
        _method = algos.backfill_inplace_object

    if _method is None:
        raise ValueError('Invalid dtype for backfill_1d [%s]' % dtype.name)

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)

    _method(values, mask, limit=limit)
    return values


def pad_2d(values, limit=None, mask=None, dtype=None):

    if dtype is None:
        dtype = values.dtype
    _method = None
    if is_float_dtype(values):
        _method = getattr(algos, 'pad_2d_inplace_%s' % dtype.name, None)
    elif dtype in _DATELIKE_DTYPES or is_datetime64_dtype(values):
        _method = _pad_2d_datetime
    elif is_integer_dtype(values):
        values = _ensure_float64(values)
        _method = algos.pad_2d_inplace_float64
    elif values.dtype == np.object_:
        _method = algos.pad_2d_inplace_object

    if _method is None:
        raise ValueError('Invalid dtype for pad_2d [%s]' % dtype.name)

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)

    if np.all(values.shape):
        _method(values, mask, limit=limit)
    else:
        # for test coverage
        pass
    return values


def backfill_2d(values, limit=None, mask=None, dtype=None):

    if dtype is None:
        dtype = values.dtype
    _method = None
    if is_float_dtype(values):
        _method = getattr(algos, 'backfill_2d_inplace_%s' % dtype.name, None)
    elif dtype in _DATELIKE_DTYPES or is_datetime64_dtype(values):
        _method = _backfill_2d_datetime
    elif is_integer_dtype(values):
        values = _ensure_float64(values)
        _method = algos.backfill_2d_inplace_float64
    elif values.dtype == np.object_:
        _method = algos.backfill_2d_inplace_object

    if _method is None:
        raise ValueError('Invalid dtype for backfill_2d [%s]' % dtype.name)

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)

    if np.all(values.shape):
        _method(values, mask, limit=limit)
    else:
        # for test coverage
        pass
    return values


def _clean_interp_method(method, **kwargs):
    order = kwargs.get('order')
    valid = ['linear', 'time', 'index', 'values', 'nearest', 'zero', 'slinear',
             'quadratic', 'cubic', 'barycentric', 'polynomial',
             'krogh', 'piecewise_polynomial',
             'pchip', 'spline']
    if method in ('spline', 'polynomial') and order is None:
        raise ValueError("You must specify the order of the spline or "
                         "polynomial.")
    if method not in valid:
        raise ValueError("method must be one of {0}."
                         "Got '{1}' instead.".format(valid, method))
    return method


def interpolate_1d(xvalues, yvalues, method='linear', limit=None,
                   limit_direction='forward',
                   fill_value=None, bounds_error=False, order=None, **kwargs):
    """
    Logic for the 1-d interpolation.  The result should be 1-d, inputs
    xvalues and yvalues will each be 1-d arrays of the same length.

    Bounds_error is currently hardcoded to False since non-scipy ones don't
    take it as an argumnet.
    """
    # Treat the original, non-scipy methods first.

    invalid = isnull(yvalues)
    valid = ~invalid

    if not valid.any():
        # have to call np.asarray(xvalues) since xvalues could be an Index
        # which cant be mutated
        result = np.empty_like(np.asarray(xvalues), dtype=np.float64)
        result.fill(np.nan)
        return result

    if valid.all():
        return yvalues

    if method == 'time':
        if not getattr(xvalues, 'is_all_dates', None):
        # if not issubclass(xvalues.dtype.type, np.datetime64):
            raise ValueError('time-weighted interpolation only works '
                             'on Series or DataFrames with a '
                             'DatetimeIndex')
        method = 'values'

    def _interp_limit(invalid, fw_limit, bw_limit):
        "Get idx of values that won't be filled b/c they exceed the limits."
        for x in np.where(invalid)[0]:
            if invalid[max(0, x - fw_limit):x + bw_limit + 1].all():
                yield x

    valid_limit_directions = ['forward', 'backward', 'both']
    limit_direction = limit_direction.lower()
    if limit_direction not in valid_limit_directions:
        msg = 'Invalid limit_direction: expecting one of %r, got %r.' % (
            valid_limit_directions, limit_direction)
        raise ValueError(msg)

    from pandas import Series
    ys = Series(yvalues)
    start_nans = set(range(ys.first_valid_index()))
    end_nans = set(range(1 + ys.last_valid_index(), len(valid)))

    # This is a list of the indexes in the series whose yvalue is currently NaN,
    # but whose interpolated yvalue will be overwritten with NaN after computing
    # the interpolation. For each index in this list, one of these conditions is
    # true of the corresponding NaN in the yvalues:
    #
    # a) It is one of a chain of NaNs at the beginning of the series, and either
    #    limit is not specified or limit_direction is 'forward'.
    # b) It is one of a chain of NaNs at the end of the series, and limit is
    #    specified and limit_direction is 'backward' or 'both'.
    # c) Limit is nonzero and it is further than limit from the nearest non-NaN
    #    value (with respect to the limit_direction setting).
    #
    # The default behavior is to fill forward with no limit, ignoring NaNs at
    # the beginning (see issues #9218 and #10420)
    violate_limit = sorted(start_nans)

    if limit:
        if limit_direction == 'forward':
            violate_limit = sorted(start_nans | set(_interp_limit(invalid, limit, 0)))
        if limit_direction == 'backward':
            violate_limit = sorted(end_nans | set(_interp_limit(invalid, 0, limit)))
        if limit_direction == 'both':
            violate_limit = sorted(_interp_limit(invalid, limit, limit))

    xvalues = getattr(xvalues, 'values', xvalues)
    yvalues = getattr(yvalues, 'values', yvalues)
    result = yvalues.copy()

    if method in ['linear', 'time', 'index', 'values']:
        if method in ('values', 'index'):
            inds = np.asarray(xvalues)
            # hack for DatetimeIndex, #1646
            if issubclass(inds.dtype.type, np.datetime64):
                inds = inds.view(np.int64)
            if inds.dtype == np.object_:
                inds = lib.maybe_convert_objects(inds)
        else:
            inds = xvalues
        result[invalid] = np.interp(inds[invalid], inds[valid], yvalues[valid])
        result[violate_limit] = np.nan
        return result

    sp_methods = ['nearest', 'zero', 'slinear', 'quadratic', 'cubic',
                  'barycentric', 'krogh', 'spline', 'polynomial',
                  'piecewise_polynomial', 'pchip']
    if method in sp_methods:
        inds = np.asarray(xvalues)
        # hack for DatetimeIndex, #1646
        if issubclass(inds.dtype.type, np.datetime64):
            inds = inds.view(np.int64)
        result[invalid] = _interpolate_scipy_wrapper(
            inds[valid], yvalues[valid], inds[invalid], method=method,
            fill_value=fill_value,
            bounds_error=bounds_error, order=order, **kwargs)
        result[violate_limit] = np.nan
        return result


def _interpolate_scipy_wrapper(x, y, new_x, method, fill_value=None,
                               bounds_error=False, order=None, **kwargs):
    """
    passed off to scipy.interpolate.interp1d. method is scipy's kind.
    Returns an array interpolated at new_x.  Add any new methods to
    the list in _clean_interp_method
    """
    try:
        from scipy import interpolate
        from pandas import DatetimeIndex
    except ImportError:
        raise ImportError('{0} interpolation requires Scipy'.format(method))

    new_x = np.asarray(new_x)

    # ignores some kwargs that could be passed along.
    alt_methods = {
        'barycentric': interpolate.barycentric_interpolate,
        'krogh': interpolate.krogh_interpolate,
        'piecewise_polynomial': interpolate.piecewise_polynomial_interpolate,
    }

    if getattr(x, 'is_all_dates', False):
        # GH 5975, scipy.interp1d can't hande datetime64s
        x, new_x = x._values.astype('i8'), new_x.astype('i8')

    try:
        alt_methods['pchip'] = interpolate.pchip_interpolate
    except AttributeError:
        if method == 'pchip':
            raise ImportError("Your version of scipy does not support "
                              "PCHIP interpolation.")

    interp1d_methods = ['nearest', 'zero', 'slinear', 'quadratic', 'cubic',
                        'polynomial']
    if method in interp1d_methods:
        if method == 'polynomial':
            method = order
        terp = interpolate.interp1d(x, y, kind=method, fill_value=fill_value,
                                    bounds_error=bounds_error)
        new_y = terp(new_x)
    elif method == 'spline':
        # GH #10633
        if not order:
            raise ValueError("order needs to be specified and greater than 0")
        terp = interpolate.UnivariateSpline(x, y, k=order, **kwargs)
        new_y = terp(new_x)
    else:
        # GH 7295: need to be able to write for some reason
        # in some circumstances: check all three
        if not x.flags.writeable:
            x = x.copy()
        if not y.flags.writeable:
            y = y.copy()
        if not new_x.flags.writeable:
            new_x = new_x.copy()
        method = alt_methods[method]
        new_y = method(x, y, new_x, **kwargs)
    return new_y


def interpolate_2d(values, method='pad', axis=0, limit=None, fill_value=None, dtype=None):
    """ perform an actual interpolation of values, values will be make 2-d if
    needed fills inplace, returns the result
    """

    transf = (lambda x: x) if axis == 0 else (lambda x: x.T)

    # reshape a 1 dim if needed
    ndim = values.ndim
    if values.ndim == 1:
        if axis != 0:  # pragma: no cover
            raise AssertionError("cannot interpolate on a ndim == 1 with "
                                 "axis != 0")
        values = values.reshape(tuple((1,) + values.shape))

    if fill_value is None:
        mask = None
    else:  # todo create faster fill func without masking
        mask = mask_missing(transf(values), fill_value)

    method = _clean_fill_method(method)
    if method == 'pad':
        values = transf(pad_2d(transf(values), limit=limit, mask=mask, dtype=dtype))
    else:
        values = transf(backfill_2d(transf(values), limit=limit, mask=mask, dtype=dtype))

    # reshape back
    if ndim == 1:
        values = values[0]

    return values


def _consensus_name_attr(objs):
    name = objs[0].name
    for obj in objs[1:]:
        if obj.name != name:
            return None
    return name


_fill_methods = {'pad': pad_1d, 'backfill': backfill_1d}


def _get_fill_func(method):
    method = _clean_fill_method(method)
    return _fill_methods[method]


#----------------------------------------------------------------------
# Lots of little utilities

def _validate_date_like_dtype(dtype):
    try:
        typ = np.datetime_data(dtype)[0]
    except ValueError as e:
        raise TypeError('%s' % e)
    if typ != 'generic' and typ != 'ns':
        raise ValueError('%r is too specific of a frequency, try passing %r'
                         % (dtype.name, dtype.type.__name__))


def _invalidate_string_dtypes(dtype_set):
    """Change string like dtypes to object for ``DataFrame.select_dtypes()``."""
    non_string_dtypes = dtype_set - _string_dtypes
    if non_string_dtypes != dtype_set:
        raise TypeError("string dtypes are not allowed, use 'object' instead")


def _get_dtype_from_object(dtype):
    """Get a numpy dtype.type-style object. This handles the
       datetime64[ns] and datetime64[ns, TZ] compat

    Notes
    -----
    If nothing can be found, returns ``object``.
    """
    # type object from a dtype
    if isinstance(dtype, type) and issubclass(dtype, np.generic):
        return dtype
    elif is_categorical(dtype):
        return CategoricalDtype().type
    elif is_datetimetz(dtype):
        return DatetimeTZDtype(dtype).type
    elif isinstance(dtype, np.dtype):  # dtype object
        try:
            _validate_date_like_dtype(dtype)
        except TypeError:
            # should still pass if we don't have a datelike
            pass
        return dtype.type
    elif isinstance(dtype, compat.string_types):
        if dtype == 'datetime' or dtype == 'timedelta':
            dtype += '64'

        try:
            return _get_dtype_from_object(getattr(np,dtype))
        except AttributeError:
            # handles cases like _get_dtype(int)
            # i.e., python objects that are valid dtypes (unlike user-defined
            # types, in general)
            # further handle internal types
            pass

    return _get_dtype_from_object(np.dtype(dtype))


def _get_info_slice(obj, indexer):
    """Slice the info axis of `obj` with `indexer`."""
    if not hasattr(obj, '_info_axis_number'):
        raise TypeError('object of type %r has no info axis' %
                        type(obj).__name__)
    slices = [slice(None)] * obj.ndim
    slices[obj._info_axis_number] = indexer
    return tuple(slices)


def _maybe_box(indexer, values, obj, key):

    # if we have multiples coming back, box em
    if isinstance(values, np.ndarray):
        return obj[indexer.get_loc(key)]

    # return the value
    return values


def _maybe_box_datetimelike(value):
    # turn a datetime like into a Timestamp/timedelta as needed

    if isinstance(value, (np.datetime64, datetime)):
        value = tslib.Timestamp(value)
    elif isinstance(value, (np.timedelta64, timedelta)):
        value = tslib.Timedelta(value)

    return value

_values_from_object = lib.values_from_object


def _possibly_castable(arr):
    # return False to force a non-fastpath

    # check datetime64[ns]/timedelta64[ns] are valid
    # otherwise try to coerce
    kind = arr.dtype.kind
    if kind == 'M' or kind == 'm':
        return arr.dtype in _DATELIKE_DTYPES

    return arr.dtype.name not in _POSSIBLY_CAST_DTYPES


def _possibly_convert_platform(values):
    """ try to do platform conversion, allow ndarray or list here """

    if isinstance(values, (list, tuple)):
        values = lib.list_to_object_array(values)
    if getattr(values, 'dtype', None) == np.object_:
        if hasattr(values, '_values'):
            values = values._values
        values = lib.maybe_convert_objects(values)

    return values


def _possibly_cast_to_datetime(value, dtype, errors='raise'):
    """ try to cast the array/value to a datetimelike dtype, converting float
    nan to iNaT
    """
    from pandas.tseries.timedeltas import to_timedelta
    from pandas.tseries.tools import to_datetime

    if dtype is not None:
        if isinstance(dtype, compat.string_types):
            dtype = np.dtype(dtype)

        is_datetime64 = is_datetime64_dtype(dtype)
        is_datetime64tz = is_datetime64tz_dtype(dtype)
        is_timedelta64 = is_timedelta64_dtype(dtype)

        if is_datetime64 or is_datetime64tz or is_timedelta64:

            # force the dtype if needed
            if is_datetime64 and not is_dtype_equal(dtype,_NS_DTYPE):
                if dtype.name == 'datetime64[ns]':
                    dtype = _NS_DTYPE
                else:
                    raise TypeError(
                        "cannot convert datetimelike to dtype [%s]" % dtype)
            elif is_datetime64tz:
                pass
            elif is_timedelta64 and not is_dtype_equal(dtype,_TD_DTYPE):
                if dtype.name == 'timedelta64[ns]':
                    dtype = _TD_DTYPE
                else:
                    raise TypeError(
                        "cannot convert timedeltalike to dtype [%s]" % dtype)

            if np.isscalar(value):
                if value == tslib.iNaT or isnull(value):
                    value = tslib.iNaT
            else:
                value = np.array(value,copy=False)

                # have a scalar array-like (e.g. NaT)
                if value.ndim == 0:
                    value = tslib.iNaT

                # we have an array of datetime or timedeltas & nulls
                elif np.prod(value.shape) and not is_dtype_equal(value.dtype, dtype):
                    try:
                        if is_datetime64:
                            value = to_datetime(value, errors=errors)._values
                        elif is_datetime64tz:

                            # input has to be UTC at this point, so just localize
                            value = to_datetime(value, errors=errors).tz_localize(dtype.tz)
                        elif is_timedelta64:
                            value = to_timedelta(value, errors=errors)._values
                    except (AttributeError, ValueError):
                        pass

        # coerce datetimelike to object
        elif is_datetime64_dtype(value) and not is_datetime64_dtype(dtype):
            if is_object_dtype(dtype):
                ints = np.asarray(value).view('i8')
                return tslib.ints_to_pydatetime(ints)

            # we have a non-castable dtype that was passed
            raise TypeError('Cannot cast datetime64 to %s' % dtype)

    else:

        is_array = isinstance(value, np.ndarray)

        # catch a datetime/timedelta that is not of ns variety
        # and no coercion specified
        if is_array and value.dtype.kind in ['M', 'm']:
            dtype = value.dtype

            if dtype.kind == 'M' and dtype != _NS_DTYPE:
                value = value.astype(_NS_DTYPE)

            elif dtype.kind == 'm' and dtype != _TD_DTYPE:
                value = to_timedelta(value)

        # only do this if we have an array and the dtype of the array is not
        # setup already we are not an integer/object, so don't bother with this
        # conversion
        elif not (is_array and not (issubclass(value.dtype.type, np.integer) or
                                    value.dtype == np.object_)):
            value = _possibly_infer_to_datetimelike(value)

    return value


def _possibly_infer_to_datetimelike(value, convert_dates=False):
    """
    we might have a array (or single object) that is datetime like,
    and no dtype is passed don't change the value unless we find a
    datetime/timedelta set

    this is pretty strict in that a datetime/timedelta is REQUIRED
    in addition to possible nulls/string likes

    ONLY strings are NOT datetimelike

    Parameters
    ----------
    value : np.array / Series / Index / list-like
    convert_dates : boolean, default False
       if True try really hard to convert dates (such as datetime.date), other
       leave inferred dtype 'date' alone

    """

    if isinstance(value, (ABCDatetimeIndex, ABCPeriodIndex)):
        return value
    elif isinstance(value, ABCSeries):
        if isinstance(value._values, ABCDatetimeIndex):
            return value._values

    v = value
    if not is_list_like(v):
        v = [v]
    v = np.array(v,copy=False)
    shape = v.shape
    if not v.ndim == 1:
        v = v.ravel()

    if len(v):

        def _try_datetime(v):
            # safe coerce to datetime64
            try:
                v = tslib.array_to_datetime(v, errors='raise')
            except ValueError:

                # we might have a sequence of the same-datetimes with tz's
                # if so coerce to a DatetimeIndex; if they are not the same, then
                # these stay as object dtype
                try:
                    from pandas import to_datetime
                    return to_datetime(v)
                except:
                    pass

            except:
                pass

            return v.reshape(shape)

        def _try_timedelta(v):
            # safe coerce to timedelta64

            # will try first with a string & object conversion
            from pandas.tseries.timedeltas import to_timedelta
            try:
                return to_timedelta(v)._values.reshape(shape)
            except:
                return v

        # do a quick inference for perf
        sample = v[:min(3,len(v))]
        inferred_type = lib.infer_dtype(sample)

        if inferred_type in ['datetime', 'datetime64'] or (convert_dates and inferred_type in ['date']):
            value = _try_datetime(v)
        elif inferred_type in ['timedelta', 'timedelta64']:
            value = _try_timedelta(v)

        # its possible to have nulls intermixed within the datetime or timedelta
        # these will in general have an inferred_type of 'mixed', so have to try
        # both datetime and timedelta

        # try timedelta first to avoid spurious datetime conversions
        # e.g. '00:00:01' is a timedelta but technically is also a datetime
        elif inferred_type in ['mixed']:

            if lib.is_possible_datetimelike_array(_ensure_object(v)):
                value = _try_timedelta(v)
                if lib.infer_dtype(value) in ['mixed']:
                    value = _try_datetime(v)

    return value


def is_bool_indexer(key):
    if isinstance(key, (ABCSeries, np.ndarray)):
        if key.dtype == np.object_:
            key = np.asarray(_values_from_object(key))

            if not lib.is_bool_array(key):
                if isnull(key).any():
                    raise ValueError('cannot index with vector containing '
                                     'NA / NaN values')
                return False
            return True
        elif key.dtype == np.bool_:
            return True
    elif isinstance(key, list):
        try:
            arr = np.asarray(key)
            return arr.dtype == np.bool_ and len(arr) == len(key)
        except TypeError:  # pragma: no cover
            return False

    return False


def _default_index(n):
    from pandas.core.index import Int64Index
    values = np.arange(n, dtype=np.int64)
    result = Int64Index(values,name=None)
    result.is_unique = True
    return result


def ensure_float(arr):
    if issubclass(arr.dtype.type, (np.integer, np.bool_)):
        arr = arr.astype(float)
    return arr


def _mut_exclusive(**kwargs):
    item1, item2 = kwargs.items()
    label1, val1 = item1
    label2, val2 = item2
    if val1 is not None and val2 is not None:
        raise TypeError('mutually exclusive arguments: %r and %r' %
                        (label1, label2))
    elif val1 is not None:
        return val1
    else:
        return val2


def _not_none(*args):
    return (arg for arg in args if arg is not None)

def _any_none(*args):
    for arg in args:
        if arg is None:
            return True
    return False


def _all_not_none(*args):
    for arg in args:
        if arg is None:
            return False
    return True


def _try_sort(iterable):
    listed = list(iterable)
    try:
        return sorted(listed)
    except Exception:
        return listed


def _count_not_none(*args):
    return sum(x is not None for x in args)

#------------------------------------------------------------------------------
# miscellaneous python tools




def adjoin(space, *lists, **kwargs):
    """
    Glues together two sets of strings using the amount of space requested.
    The idea is to prettify.

    ----------
    space : int
        number of spaces for padding
    lists : str
        list of str which being joined
    strlen : callable
        function used to calculate the length of each str. Needed for unicode
        handling.
    justfunc : callable
        function used to justify str. Needed for unicode handling.
    """
    strlen = kwargs.pop('strlen', len)
    justfunc = kwargs.pop('justfunc', _justify)

    out_lines = []
    newLists = []
    lengths = [max(map(strlen, x)) + space for x in lists[:-1]]
    # not the last one
    lengths.append(max(map(len, lists[-1])))
    maxLen = max(map(len, lists))
    for i, lst in enumerate(lists):
        nl = justfunc(lst, lengths[i], mode='left')
        nl.extend([' ' * lengths[i]] * (maxLen - len(lst)))
        newLists.append(nl)
    toJoin = zip(*newLists)
    for lines in toJoin:
        out_lines.append(_join_unicode(lines))
    return _join_unicode(out_lines, sep='\n')

def _justify(texts, max_len, mode='right'):
    """
    Perform ljust, center, rjust against string or list-like
    """
    if mode == 'left':
        return [x.ljust(max_len) for x in texts]
    elif mode == 'center':
        return [x.center(max_len) for x in texts]
    else:
        return [x.rjust(max_len) for x in texts]

def _join_unicode(lines, sep=''):
    try:
        return sep.join(lines)
    except UnicodeDecodeError:
        sep = compat.text_type(sep)
        return sep.join([x.decode('utf-8') if isinstance(x, str) else x
                         for x in lines])


def iterpairs(seq):
    """
    Parameters
    ----------
    seq: sequence

    Returns
    -------
    iterator returning overlapping pairs of elements

    Examples
    --------
    >>> list(iterpairs([1, 2, 3, 4]))
    [(1, 2), (2, 3), (3, 4)]
    """
    # input may not be sliceable
    seq_it = iter(seq)
    seq_it_next = iter(seq)
    next(seq_it_next)

    return zip(seq_it, seq_it_next)


def split_ranges(mask):
    """ Generates tuples of ranges which cover all True value in mask

    >>> list(split_ranges([1,0,0,1,0]))
    [(0, 1), (3, 4)]
    """
    ranges = [(0, len(mask))]

    for pos, val in enumerate(mask):
        if not val:  # this pos should be ommited, split off the prefix range
            r = ranges.pop()
            if pos > r[0]:  # yield non-zero range
                yield (r[0], pos)
            if pos + 1 < len(mask):  # save the rest for processing
                ranges.append((pos + 1, len(mask)))
    if ranges:
        yield ranges[-1]


def indent(string, spaces=4):
    dent = ' ' * spaces
    return '\n'.join([dent + x for x in string.split('\n')])


def banner(message):
    """
    Return 80-char width message declaration with = bars on top and bottom.
    """
    bar = '=' * 80
    return '%s\n%s\n%s' % (bar, message, bar)


def _long_prod(vals):
    result = long(1)
    for x in vals:
        result *= x
    return result


class groupby(dict):

    """
    A simple groupby different from the one in itertools.

    Does not require the sequence elements to be sorted by keys,
    however it is slower.
    """

    def __init__(self, seq, key=lambda x: x):
        for value in seq:
            k = key(value)
            self.setdefault(k, []).append(value)
    try:
        __iter__ = dict.iteritems
    except AttributeError:  # pragma: no cover
        # Python 3
        def __iter__(self):
            return iter(dict.items(self))


def map_indices_py(arr):
    """
    Returns a dictionary with (element, index) pairs for each element in the
    given array/list
    """
    return dict([(x, i) for i, x in enumerate(arr)])


def union(*seqs):
    result = set([])
    for seq in seqs:
        if not isinstance(seq, set):
            seq = set(seq)
        result |= seq
    return type(seqs[0])(list(result))


def difference(a, b):
    return type(a)(list(set(a) - set(b)))


def intersection(*seqs):
    result = set(seqs[0])
    for seq in seqs:
        if not isinstance(seq, set):
            seq = set(seq)
        result &= seq
    return type(seqs[0])(list(result))


def _asarray_tuplesafe(values, dtype=None):
    from pandas.core.index import Index

    if not (isinstance(values, (list, tuple))
            or hasattr(values, '__array__')):
        values = list(values)
    elif isinstance(values, Index):
        return values.values

    if isinstance(values, list) and dtype in [np.object_, object]:
        return lib.list_to_object_array(values)

    result = np.asarray(values, dtype=dtype)

    if issubclass(result.dtype.type, compat.string_types):
        result = np.asarray(values, dtype=object)

    if result.ndim == 2:
        if isinstance(values, list):
            return lib.list_to_object_array(values)
        else:
            # Making a 1D array that safely contains tuples is a bit tricky
            # in numpy, leading to the following
            try:
                result = np.empty(len(values), dtype=object)
                result[:] = values
            except ValueError:
                # we have a list-of-list
                result[:] = [tuple(x) for x in values]

    return result


def _index_labels_to_array(labels):
    if isinstance(labels, (compat.string_types, tuple)):
        labels = [labels]

    if not isinstance(labels, (list, np.ndarray)):
        try:
            labels = list(labels)
        except TypeError:  # non-iterable
            labels = [labels]

    labels = _asarray_tuplesafe(labels)

    return labels


def _maybe_make_list(obj):
    if obj is not None and not isinstance(obj, (tuple, list)):
        return [obj]
    return obj

########################
##### TYPE TESTING #####
########################

is_bool = lib.is_bool


is_integer = lib.is_integer


is_float = lib.is_float


is_complex = lib.is_complex


def is_iterator(obj):
    # python 3 generators have __next__ instead of next
    return hasattr(obj, 'next') or hasattr(obj, '__next__')


def is_number(obj):
    return isinstance(obj, (numbers.Number, np.number))

def is_period_arraylike(arr):
    """ return if we are period arraylike / PeriodIndex """
    if isinstance(arr, pd.PeriodIndex):
        return True
    elif isinstance(arr, (np.ndarray, ABCSeries)):
        return arr.dtype == object and lib.infer_dtype(arr) == 'period'
    return getattr(arr, 'inferred_type', None) == 'period'

def is_datetime_arraylike(arr):
    """ return if we are datetime arraylike / DatetimeIndex """
    if isinstance(arr, ABCDatetimeIndex):
        return True
    elif isinstance(arr, (np.ndarray, ABCSeries)):
        return arr.dtype == object and lib.infer_dtype(arr) == 'datetime'
    return getattr(arr, 'inferred_type', None) == 'datetime'

def is_datetimelike(arr):
    return arr.dtype in _DATELIKE_DTYPES or isinstance(arr, ABCPeriodIndex) or is_datetimetz(arr)

def _coerce_to_dtype(dtype):
    """ coerce a string / np.dtype to a dtype """
    if is_categorical_dtype(dtype):
        dtype = CategoricalDtype()
    elif is_datetime64tz_dtype(dtype):
        dtype = DatetimeTZDtype(dtype)
    else:
        dtype = np.dtype(dtype)
    return dtype

def _get_dtype(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        return arr_or_dtype
    elif isinstance(arr_or_dtype, type):
        return np.dtype(arr_or_dtype)
    elif isinstance(arr_or_dtype, CategoricalDtype):
        return arr_or_dtype
    elif isinstance(arr_or_dtype, DatetimeTZDtype):
        return arr_or_dtype
    elif isinstance(arr_or_dtype, compat.string_types):
        if is_categorical_dtype(arr_or_dtype):
            return CategoricalDtype.construct_from_string(arr_or_dtype)
        elif is_datetime64tz_dtype(arr_or_dtype):
            return DatetimeTZDtype.construct_from_string(arr_or_dtype)

    if hasattr(arr_or_dtype, 'dtype'):
        arr_or_dtype = arr_or_dtype.dtype
    return np.dtype(arr_or_dtype)

def _get_dtype_type(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        return arr_or_dtype.type
    elif isinstance(arr_or_dtype, type):
        return np.dtype(arr_or_dtype).type
    elif isinstance(arr_or_dtype, CategoricalDtype):
        return CategoricalDtypeType
    elif isinstance(arr_or_dtype, DatetimeTZDtype):
        return DatetimeTZDtypeType
    elif isinstance(arr_or_dtype, compat.string_types):
        if is_categorical_dtype(arr_or_dtype):
            return CategoricalDtypeType
        elif is_datetime64tz_dtype(arr_or_dtype):
            return DatetimeTZDtypeType
        return _get_dtype_type(np.dtype(arr_or_dtype))
    try:
        return arr_or_dtype.dtype.type
    except AttributeError:
        return type(None)

def is_dtype_equal(source, target):
    """ return a boolean if the dtypes are equal """
    try:
        source = _get_dtype(source)
        target = _get_dtype(target)
        return source == target
    except (TypeError, AttributeError):

        # invalid comparison
        # object == category will hit this
        return False

def is_any_int_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.integer)


def is_integer_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return (issubclass(tipo, np.integer) and
            not issubclass(tipo, (np.datetime64, np.timedelta64)))

def is_int64_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.int64)


def is_int_or_datetime_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return (issubclass(tipo, np.integer) or
            issubclass(tipo, (np.datetime64, np.timedelta64)))

def is_datetime64_dtype(arr_or_dtype):
    try:
        tipo = _get_dtype_type(arr_or_dtype)
    except TypeError:
        return False
    return issubclass(tipo, np.datetime64)

def is_datetime64tz_dtype(arr_or_dtype):
    return DatetimeTZDtype.is_dtype(arr_or_dtype)

def is_datetime64_any_dtype(arr_or_dtype):
    return is_datetime64_dtype(arr_or_dtype) or is_datetime64tz_dtype(arr_or_dtype)

def is_datetime64_ns_dtype(arr_or_dtype):
    try:
        tipo = _get_dtype(arr_or_dtype)
    except TypeError:
        return False
    return tipo == _NS_DTYPE

def is_timedelta64_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.timedelta64)


def is_timedelta64_ns_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return tipo == _TD_DTYPE


def is_datetime_or_timedelta_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, (np.datetime64, np.timedelta64))


def is_datetimelike_v_numeric(a, b):
    # return if we have an i8 convertible and numeric comparision
    if not hasattr(a,'dtype'):
        a = np.asarray(a)
    if not hasattr(b, 'dtype'):
        b = np.asarray(b)
    is_numeric = lambda x: is_integer_dtype(x) or is_float_dtype(x)
    is_datetimelike = needs_i8_conversion
    return (is_datetimelike(a) and is_numeric(b)) or (
        is_datetimelike(b) and is_numeric(a))

def is_datetimelike_v_object(a, b):
    # return if we have an i8 convertible and object comparision
    if not hasattr(a,'dtype'):
        a = np.asarray(a)
    if not hasattr(b, 'dtype'):
        b = np.asarray(b)
    f = lambda x: is_object_dtype(x)
    is_object = lambda x: is_integer_dtype(x) or is_float_dtype(x)
    is_datetimelike = needs_i8_conversion
    return (is_datetimelike(a) and is_object(b)) or (
        is_datetimelike(b) and is_object(a))

needs_i8_conversion = lambda arr_or_dtype: is_datetime_or_timedelta_dtype(arr_or_dtype) or \
                      is_datetime64tz_dtype(arr_or_dtype)

def i8_boxer(arr_or_dtype):
    """ return the scalar boxer for the dtype """
    if is_datetime64_dtype(arr_or_dtype) or is_datetime64tz_dtype(arr_or_dtype):
        return lib.Timestamp
    elif is_timedelta64_dtype(arr_or_dtype):
        return lambda x: lib.Timedelta(x,unit='ns')
    raise ValueError("cannot find a scalar boxer for {0}".format(arr_or_dtype))

def is_numeric_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return (issubclass(tipo, (np.number, np.bool_))
            and not issubclass(tipo, (np.datetime64, np.timedelta64)))


def is_float_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.floating)


def is_floating_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return isinstance(tipo, np.floating)


def is_bool_dtype(arr_or_dtype):
    try:
        tipo = _get_dtype_type(arr_or_dtype)
    except ValueError:
        # this isn't even a dtype
        return False
    return issubclass(tipo, np.bool_)

def is_sparse(array):
    """ return if we are a sparse array """
    return isinstance(array, (ABCSparseArray, ABCSparseSeries))

def is_datetimetz(array):
    """ return if we are a datetime with tz array """
    return (isinstance(array, ABCDatetimeIndex) and getattr(array,'tz',None) is not None) or is_datetime64tz_dtype(array)

def is_internal_type(value):
    """
    if we are a klass that is preserved by the internals
    these are internal klasses that we represent (and don't use a np.array)
    """
    if is_categorical(value):
        return True
    elif is_sparse(value):
        return True
    elif is_datetimetz(value):
        return True
    return False

def is_categorical(array):
    """ return if we are a categorical possibility """
    return isinstance(array, ABCCategorical) or is_categorical_dtype(array)

def is_categorical_dtype(arr_or_dtype):
    return CategoricalDtype.is_dtype(arr_or_dtype)

def is_complex_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.complexfloating)


def is_object_dtype(arr_or_dtype):
    tipo = _get_dtype_type(arr_or_dtype)
    return issubclass(tipo, np.object_)


def is_re(obj):
    return isinstance(obj, re._pattern_type)


def is_re_compilable(obj):
    try:
        re.compile(obj)
    except TypeError:
        return False
    else:
        return True


def is_list_like(arg):
     return (hasattr(arg, '__iter__') and
            not isinstance(arg, compat.string_and_binary_types))

def is_null_slice(obj):
    """ we have a null slice """
    return (isinstance(obj, slice) and obj.start is None and
            obj.stop is None and obj.step is None)

def is_full_slice(obj, l):
    """ we have a full length slice """
    return (isinstance(obj, slice) and obj.start == 0 and
            obj.stop == l and obj.step is None)

def is_hashable(arg):
    """Return True if hash(arg) will succeed, False otherwise.

    Some types will pass a test against collections.Hashable but fail when they
    are actually hashed with hash().

    Distinguish between these and other types by trying the call to hash() and
    seeing if they raise TypeError.

    Examples
    --------
    >>> a = ([],)
    >>> isinstance(a, collections.Hashable)
    True
    >>> is_hashable(a)
    False
    """
    # unfortunately, we can't use isinstance(arg, collections.Hashable), which
    # can be faster than calling hash, because numpy scalars on Python 3 fail
    # this test

    # reconsider this decision once this numpy bug is fixed:
    # https://github.com/numpy/numpy/issues/5562

    try:
        hash(arg)
    except TypeError:
        return False
    else:
        return True


def is_sequence(x):
    try:
        iter(x)
        len(x)  # it has a length
        return not isinstance(x, compat.string_and_binary_types)
    except (TypeError, AttributeError):
        return False


def _get_callable_name(obj):
    # typical case has name
    if hasattr(obj, '__name__'):
        return getattr(obj, '__name__')
    # some objects don't; could recurse
    if isinstance(obj, partial):
        return _get_callable_name(obj.func)
    # fall back to class name
    if hasattr(obj, '__call__'):
        return obj.__class__.__name__
    # everything failed (probably because the argument
    # wasn't actually callable); we return None
    # instead of the empty string in this case to allow
    # distinguishing between no name and a name of ''
    return None

_string_dtypes = frozenset(map(_get_dtype_from_object, (compat.binary_type,
                                                        compat.text_type)))


_ensure_float64 = algos.ensure_float64
_ensure_float32 = algos.ensure_float32
_ensure_int64 = algos.ensure_int64
_ensure_int32 = algos.ensure_int32
_ensure_int16 = algos.ensure_int16
_ensure_int8 = algos.ensure_int8
_ensure_platform_int = algos.ensure_platform_int
_ensure_object = algos.ensure_object


def _astype_nansafe(arr, dtype, copy=True):
    """ return a view if copy is False, but
        need to be very careful as the result shape could change! """
    if not isinstance(dtype, np.dtype):
        dtype = _coerce_to_dtype(dtype)

    if issubclass(dtype.type, compat.text_type):
        # in Py3 that's str, in Py2 that's unicode
        return lib.astype_unicode(arr.ravel()).reshape(arr.shape)
    elif issubclass(dtype.type, compat.string_types):
        return lib.astype_str(arr.ravel()).reshape(arr.shape)
    elif is_datetime64_dtype(arr):
        if dtype == object:
            return tslib.ints_to_pydatetime(arr.view(np.int64))
        elif dtype == np.int64:
            return arr.view(dtype)
        elif dtype != _NS_DTYPE:
            raise TypeError("cannot astype a datetimelike from [%s] to [%s]" %
                            (arr.dtype, dtype))
        return arr.astype(_NS_DTYPE)
    elif is_timedelta64_dtype(arr):
        if dtype == np.int64:
            return arr.view(dtype)
        elif dtype == object:
            return tslib.ints_to_pytimedelta(arr.view(np.int64))

        # in py3, timedelta64[ns] are int64
        elif ((compat.PY3 and dtype not in [_INT64_DTYPE, _TD_DTYPE]) or
                (not compat.PY3 and dtype != _TD_DTYPE)):

            # allow frequency conversions
            if dtype.kind == 'm':
                mask = isnull(arr)
                result = arr.astype(dtype).astype(np.float64)
                result[mask] = np.nan
                return result

            raise TypeError("cannot astype a timedelta from [%s] to [%s]" %
                            (arr.dtype, dtype))

        return arr.astype(_TD_DTYPE)
    elif (np.issubdtype(arr.dtype, np.floating) and
          np.issubdtype(dtype, np.integer)):

        if np.isnan(arr).any():
            raise ValueError('Cannot convert NA to integer')
    elif arr.dtype == np.object_ and np.issubdtype(dtype.type, np.integer):
        # work around NumPy brokenness, #1987
        return lib.astype_intsafe(arr.ravel(), dtype).reshape(arr.shape)

    if copy:
        return arr.astype(dtype)
    return arr.view(dtype)


def _clean_fill_method(method, allow_nearest=False):
    if method is None:
        return None
    method = method.lower()
    if method == 'ffill':
        method = 'pad'
    if method == 'bfill':
        method = 'backfill'

    valid_methods = ['pad', 'backfill']
    expecting = 'pad (ffill) or backfill (bfill)'
    if allow_nearest:
        valid_methods.append('nearest')
        expecting = 'pad (ffill), backfill (bfill) or nearest'
    if method not in valid_methods:
        msg = ('Invalid fill method. Expecting %s. Got %s'
               % (expecting, method))
        raise ValueError(msg)
    return method


def _clean_reindex_fill_method(method):
    return _clean_fill_method(method, allow_nearest=True)


def _all_none(*args):
    for arg in args:
        if arg is not None:
            return False
    return True


class UTF8Recoder:

    """
    Iterator that reads an encoded stream and reencodes the input to UTF-8
    """

    def __init__(self, f, encoding):
        self.reader = codecs.getreader(encoding)(f)

    def __iter__(self):
        return self

    def read(self, bytes=-1):
        return self.reader.read(bytes).encode('utf-8')

    def readline(self):
        return self.reader.readline().encode('utf-8')

    def next(self):
        return next(self.reader).encode("utf-8")

    # Python 3 iterator
    __next__ = next


def _get_handle(path, mode, encoding=None, compression=None):
    """Gets file handle for given path and mode.
    NOTE: Under Python 3.2, getting a compressed file handle means reading in
    the entire file, decompressing it and decoding it to ``str`` all at once
    and then wrapping it in a StringIO.
    """
    if compression is not None:
        if encoding is not None and not compat.PY3:
            msg = 'encoding + compression not yet supported in Python 2'
            raise ValueError(msg)

        if compression == 'gzip':
            import gzip
            f = gzip.GzipFile(path, 'rb')
        elif compression == 'bz2':
            import bz2

            f = bz2.BZ2File(path, 'rb')
        else:
            raise ValueError('Unrecognized compression type: %s' %
                             compression)
        if compat.PY3:
            from io import TextIOWrapper
            f = TextIOWrapper(f, encoding=encoding)
        return f
    else:
        if compat.PY3:
            if encoding:
                f = open(path, mode, encoding=encoding)
            else:
                f = open(path, mode, errors='replace')
        else:
            f = open(path, mode)

    return f


if compat.PY3:  # pragma: no cover
    def UnicodeReader(f, dialect=csv.excel, encoding="utf-8", **kwds):
        # ignore encoding
        return csv.reader(f, dialect=dialect, **kwds)

    def UnicodeWriter(f, dialect=csv.excel, encoding="utf-8", **kwds):
        return csv.writer(f, dialect=dialect, **kwds)
else:
    class UnicodeReader:

        """
        A CSV reader which will iterate over lines in the CSV file "f",
        which is encoded in the given encoding.

        On Python 3, this is replaced (below) by csv.reader, which handles
        unicode.
        """

        def __init__(self, f, dialect=csv.excel, encoding="utf-8", **kwds):
            f = UTF8Recoder(f, encoding)
            self.reader = csv.reader(f, dialect=dialect, **kwds)

        def next(self):
            row = next(self.reader)
            return [compat.text_type(s, "utf-8") for s in row]

        # python 3 iterator
        __next__ = next

        def __iter__(self):  # pragma: no cover
            return self

    class UnicodeWriter:

        """
        A CSV writer which will write rows to CSV file "f",
        which is encoded in the given encoding.
        """

        def __init__(self, f, dialect=csv.excel, encoding="utf-8", **kwds):
            # Redirect output to a queue
            self.queue = StringIO()
            self.writer = csv.writer(self.queue, dialect=dialect, **kwds)
            self.stream = f
            self.encoder = codecs.getincrementalencoder(encoding)()
            self.quoting = kwds.get("quoting", None)

        def writerow(self, row):
            def _check_as_is(x):
                return (self.quoting == csv.QUOTE_NONNUMERIC and
                        is_number(x)) or isinstance(x, str)

            row = [x if _check_as_is(x)
                   else pprint_thing(x).encode('utf-8') for x in row]

            self.writer.writerow([s for s in row])
            # Fetch UTF-8 output from the queue ...
            data = self.queue.getvalue()
            data = data.decode("utf-8")
            # ... and reencode it into the target encoding
            data = self.encoder.encode(data)
            # write to the target stream
            self.stream.write(data)
            # empty queue
            self.queue.truncate(0)

        def writerows(self, rows):
            def _check_as_is(x):
                return (self.quoting == csv.QUOTE_NONNUMERIC and
                        is_number(x)) or isinstance(x, str)

            for i, row in enumerate(rows):
                rows[i] = [x if _check_as_is(x)
                           else pprint_thing(x).encode('utf-8') for x in row]

            self.writer.writerows([[s for s in row] for row in rows])
            # Fetch UTF-8 output from the queue ...
            data = self.queue.getvalue()
            data = data.decode("utf-8")
            # ... and reencode it into the target encoding
            data = self.encoder.encode(data)
            # write to the target stream
            self.stream.write(data)
            # empty queue
            self.queue.truncate(0)


def get_dtype_kinds(l):
    """
    Parameters
    ----------
    l : list of arrays

    Returns
    -------
    a set of kinds that exist in this list of arrays
    """

    typs = set()
    for arr in l:

        dtype = arr.dtype
        if is_categorical_dtype(dtype):
            typ = 'category'
        elif is_sparse(arr):
            typ = 'sparse'
        elif is_datetimetz(arr):
            typ = 'datetimetz'
        elif is_datetime64_dtype(dtype):
            typ = 'datetime'
        elif is_timedelta64_dtype(dtype):
            typ = 'timedelta'
        elif is_object_dtype(dtype):
            typ = 'object'
        elif is_bool_dtype(dtype):
            typ = 'bool'
        else:
            typ = dtype.kind
        typs.add(typ)
    return typs

def _concat_compat(to_concat, axis=0):
    """
    provide concatenation of an array of arrays each of which is a single
    'normalized' dtypes (in that for example, if its object, then it is a non-datetimelike
    provde a combined dtype for the resulting array the preserves the overall dtype if possible)

    Parameters
    ----------
    to_concat : array of arrays
    axis : axis to provide concatenation

    Returns
    -------
    a single array, preserving the combined dtypes
    """

    # filter empty arrays
    # 1-d dtypes always are included here
    def is_nonempty(x):
        try:
            return x.shape[axis] > 0
        except Exception:
            return True
    nonempty = [x for x in to_concat if is_nonempty(x)]

    # If all arrays are empty, there's nothing to convert, just short-cut to
    # the concatenation, #3121.
    #
    # Creating an empty array directly is tempting, but the winnings would be
    # marginal given that it would still require shape & dtype calculation and
    # np.concatenate which has them both implemented is compiled.

    typs = get_dtype_kinds(to_concat)

    # these are mandated to handle empties as well
    if 'datetime' in typs or 'datetimetz' in typs or 'timedelta' in typs:
        from pandas.tseries.common import _concat_compat
        return _concat_compat(to_concat, axis=axis)

    elif 'sparse' in typs:
        from pandas.sparse.array import _concat_compat
        return _concat_compat(to_concat, axis=axis)

    elif 'category' in typs:
        from pandas.core.categorical import _concat_compat
        return _concat_compat(to_concat, axis=axis)

    if not nonempty:

        # we have all empties, but may need to coerce the result dtype to object if we
        # have non-numeric type operands (numpy would otherwise cast this to float)
        typs = get_dtype_kinds(to_concat)
        if len(typs) != 1:

            if not len(typs-set(['i','u','f'])) or not len(typs-set(['bool','i','u'])):
                # let numpy coerce
                pass
            else:
                # coerce to object
                to_concat = [ x.astype('object') for x in to_concat ]

    return np.concatenate(to_concat,axis=axis)

def _where_compat(mask, arr1, arr2):
    if arr1.dtype == _NS_DTYPE and arr2.dtype == _NS_DTYPE:
        new_vals = np.where(mask, arr1.view('i8'), arr2.view('i8'))
        return new_vals.view(_NS_DTYPE)

    import pandas.tslib as tslib
    if arr1.dtype == _NS_DTYPE:
        arr1 = tslib.ints_to_pydatetime(arr1.view('i8'))
    if arr2.dtype == _NS_DTYPE:
        arr2 = tslib.ints_to_pydatetime(arr2.view('i8'))

    return np.where(mask, arr1, arr2)

def _dict_compat(d):
    """
    Helper function to convert datetimelike-keyed dicts to Timestamp-keyed dict

    Parameters
    ----------
    d: dict like object

    Returns
    -------
    dict

    """
    return dict((_maybe_box_datetimelike(key), value) for key, value in iteritems(d))

def sentinel_factory():

    class Sentinel(object):
        pass

    return Sentinel()

def in_interactive_session():
    """ check if we're running in an interactive shell

    returns True if running under python/ipython interactive shell
    """
    def check_main():
        import __main__ as main
        return (not hasattr(main, '__file__') or
                get_option('mode.sim_interactive'))

    try:
        return __IPYTHON__ or check_main()
    except:
        return check_main()


def in_qtconsole():
    """
    check if we're inside an IPython qtconsole

    DEPRECATED: This is no longer needed, or working, in IPython 3 and above.
    """
    try:
        ip = get_ipython()
        front_end = (
            ip.config.get('KernelApp', {}).get('parent_appname', "") or
            ip.config.get('IPKernelApp', {}).get('parent_appname', "")
        )
        if 'qtconsole' in front_end.lower():
            return True
    except:
        return False
    return False


def in_ipnb():
    """
    check if we're inside an IPython Notebook

    DEPRECATED: This is no longer used in pandas, and won't work in IPython 3
    and above.
    """
    try:
        ip = get_ipython()
        front_end = (
            ip.config.get('KernelApp', {}).get('parent_appname', "") or
            ip.config.get('IPKernelApp', {}).get('parent_appname', "")
        )
        if 'notebook' in front_end.lower():
            return True
    except:
        return False
    return False


def in_ipython_frontend():
    """
    check if we're inside an an IPython zmq frontend
    """
    try:
        ip = get_ipython()
        return 'zmq' in str(type(ip)).lower()
    except:
        pass

    return False

# Unicode consolidation
# ---------------------
#
# pprinting utility functions for generating Unicode text or
# bytes(3.x)/str(2.x) representations of objects.
# Try to use these as much as possible rather then rolling your own.
#
# When to use
# -----------
#
# 1) If you're writing code internal to pandas (no I/O directly involved),
#    use pprint_thing().
#
#    It will always return unicode text which can handled by other
#    parts of the package without breakage.
#
# 2) If you need to send something to the console, use console_encode().
#
#    console_encode() should (hopefully) choose the right encoding for you
#    based on the encoding set in option "display.encoding"
#
# 3) if you need to write something out to file, use
#    pprint_thing_encoded(encoding).
#
#    If no encoding is specified, it defaults to utf-8. Since encoding pure
#    ascii with utf-8 is a no-op you can safely use the default utf-8 if you're
#    working with straight ascii.


def _pprint_seq(seq, _nest_lvl=0, max_seq_items=None, **kwds):
    """
    internal. pprinter for iterables. you should probably use pprint_thing()
    rather then calling this directly.

    bounds length of printed sequence, depending on options
    """
    if isinstance(seq, set):
        fmt = u("set([%s])")
    else:
        fmt = u("[%s]") if hasattr(seq, '__setitem__') else u("(%s)")

    if max_seq_items is False:
        nitems = len(seq)
    else:
        nitems = max_seq_items or get_option("max_seq_items") or len(seq)

    s = iter(seq)
    r = []
    for i in range(min(nitems, len(seq))):  # handle sets, no slicing
        r.append(pprint_thing(next(s), _nest_lvl + 1, max_seq_items=max_seq_items, **kwds))
    body = ", ".join(r)

    if nitems < len(seq):
        body += ", ..."
    elif isinstance(seq, tuple) and len(seq) == 1:
        body += ','

    return fmt % body


def _pprint_dict(seq, _nest_lvl=0, max_seq_items=None, **kwds):
    """
    internal. pprinter for iterables. you should probably use pprint_thing()
    rather then calling this directly.
    """
    fmt = u("{%s}")
    pairs = []

    pfmt = u("%s: %s")

    if max_seq_items is False:
        nitems = len(seq)
    else:
        nitems = max_seq_items or get_option("max_seq_items") or len(seq)

    for k, v in list(seq.items())[:nitems]:
        pairs.append(pfmt % (pprint_thing(k, _nest_lvl + 1, max_seq_items=max_seq_items, **kwds),
                             pprint_thing(v, _nest_lvl + 1, max_seq_items=max_seq_items, **kwds)))

    if nitems < len(seq):
        return fmt % (", ".join(pairs) + ", ...")
    else:
        return fmt % ", ".join(pairs)


def pprint_thing(thing, _nest_lvl=0, escape_chars=None, default_escapes=False,
                 quote_strings=False, max_seq_items=None):
    """
    This function is the sanctioned way of converting objects
    to a unicode representation.

    properly handles nested sequences containing unicode strings
    (unicode(object) does not)

    Parameters
    ----------
    thing : anything to be formatted
    _nest_lvl : internal use only. pprint_thing() is mutually-recursive
        with pprint_sequence, this argument is used to keep track of the
        current nesting level, and limit it.
    escape_chars : list or dict, optional
        Characters to escape. If a dict is passed the values are the
        replacements
    default_escapes : bool, default False
        Whether the input escape characters replaces or adds to the defaults
    max_seq_items : False, int, default None
        Pass thru to other pretty printers to limit sequence printing

    Returns
    -------
    result - unicode object on py2, str on py3. Always Unicode.

    """
    def as_escaped_unicode(thing, escape_chars=escape_chars):
        # Unicode is fine, else we try to decode using utf-8 and 'replace'
        # if that's not it either, we have no way of knowing and the user
        # should deal with it himself.

        try:
            result = compat.text_type(thing)  # we should try this first
        except UnicodeDecodeError:
            # either utf-8 or we replace errors
            result = str(thing).decode('utf-8', "replace")

        translate = {'\t': r'\t',
                     '\n': r'\n',
                     '\r': r'\r',
                     }
        if isinstance(escape_chars, dict):
            if default_escapes:
                translate.update(escape_chars)
            else:
                translate = escape_chars
            escape_chars = list(escape_chars.keys())
        else:
            escape_chars = escape_chars or tuple()
        for c in escape_chars:
            result = result.replace(c, translate[c])

        return compat.text_type(result)

    if (compat.PY3 and hasattr(thing, '__next__')) or hasattr(thing, 'next'):
        return compat.text_type(thing)
    elif (isinstance(thing, dict) and
          _nest_lvl < get_option("display.pprint_nest_depth")):
        result = _pprint_dict(thing, _nest_lvl, quote_strings=True, max_seq_items=max_seq_items)
    elif is_sequence(thing) and _nest_lvl < \
            get_option("display.pprint_nest_depth"):
        result = _pprint_seq(thing, _nest_lvl, escape_chars=escape_chars,
                             quote_strings=quote_strings, max_seq_items=max_seq_items)
    elif isinstance(thing, compat.string_types) and quote_strings:
        if compat.PY3:
            fmt = "'%s'"
        else:
            fmt = "u'%s'"
        result = fmt % as_escaped_unicode(thing)
    else:
        result = as_escaped_unicode(thing)

    return compat.text_type(result)  # always unicode


def pprint_thing_encoded(object, encoding='utf-8', errors='replace', **kwds):
    value = pprint_thing(object)  # get unicode representation of object
    return value.encode(encoding, errors, **kwds)


def console_encode(object, **kwds):
    """
    this is the sanctioned way to prepare something for
    sending *to the console*, it delegates to pprint_thing() to get
    a unicode representation of the object relies on the global encoding
    set in display.encoding. Use this everywhere
    where you output to the console.
    """
    return pprint_thing_encoded(object,
                                get_option("display.encoding"))

def _maybe_match_name(a, b):
    a_has = hasattr(a, 'name')
    b_has = hasattr(b, 'name')
    if a_has and b_has:
        if a.name == b.name:
            return a.name
        else:
            return None
    elif a_has:
        return a.name
    elif b_has:
        return b.name
    return None

def _random_state(state=None):
    """
    Helper function for processing random_state arguments.

    Parameters
    ----------
    state : int, np.random.RandomState, None.
        If receives an int, passes to np.random.RandomState() as seed.
        If receives an np.random.RandomState object, just returns object.
        If receives `None`, returns an np.random.RandomState object.
        If receives anything else, raises an informative ValueError.
        Default None.

    Returns
    -------
    np.random.RandomState
    """

    if is_integer(state):
        return np.random.RandomState(state)
    elif isinstance(state, np.random.RandomState):
        return state
    elif state is None:
        return np.random.RandomState()
    else:
        raise ValueError("random_state must be an integer, a numpy RandomState, or None")
