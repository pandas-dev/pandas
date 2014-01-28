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

from numpy.lib.format import read_array, write_array
import numpy as np

import pandas as pd
import pandas.algos as algos
import pandas.lib as lib
import pandas.tslib as tslib
from pandas import compat
from pandas.compat import StringIO, BytesIO, range, long, u, zip, map

from pandas.core.config import get_option
from pandas.core import array as pa


class PandasError(Exception):
    pass


class SettingWithCopyError(ValueError):
    pass


class SettingWithCopyWarning(Warning):
    pass


class AmbiguousIndexError(PandasError, KeyError):
    pass


_POSSIBLY_CAST_DTYPES = set([np.dtype(t).name
                             for t in ['O', 'int8',
                                       'uint8', 'int16', 'uint16', 'int32',
                                       'uint32', 'int64', 'uint64']])

_NS_DTYPE = np.dtype('M8[ns]')
_TD_DTYPE = np.dtype('m8[ns]')
_INT64_DTYPE = np.dtype(np.int64)
_DATELIKE_DTYPES = set([np.dtype(t) for t in ['M8[ns]', '<M8[ns]', '>M8[ns]',
                                              'm8[ns]', '<m8[ns]', '>m8[ns]']])


# define abstract base classes to enable isinstance type checking on our
# objects
def create_pandas_abc_type(name, attr, comp):
    @classmethod
    def _check(cls, inst):
        return getattr(inst, attr, None) in comp
    dct = dict(__instancecheck__=_check,
               __subclasscheck__=_check)
    meta = type("ABCBase", (type,), dct)
    return meta(name, tuple(), dct)


ABCSeries = create_pandas_abc_type("ABCSeries", "_typ", ("series",))
ABCDataFrame = create_pandas_abc_type("ABCDataFrame", "_typ", ("dataframe",))
ABCPanel = create_pandas_abc_type("ABCPanel", "_typ", ("panel",))
ABCSparseSeries = create_pandas_abc_type("ABCSparseSeries", "_subtyp",
                                         ('sparse_series',
                                          'sparse_time_series'))
ABCSparseArray = create_pandas_abc_type("ABCSparseArray", "_subtyp",
                                        ('sparse_array', 'sparse_series'))


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
    """
    return _isnull(obj)


def _isnull_new(obj):
    if lib.isscalar(obj):
        return lib.checknull(obj)
    # hack (for now) because MI registers as ndarray
    elif isinstance(obj, pd.MultiIndex):
        raise NotImplementedError("isnull is not defined for MultiIndex")
    elif isinstance(obj, (ABCSeries, np.ndarray)):
        return _isnull_ndarraylike(obj)
    elif isinstance(obj, ABCGeneric):
        return obj._constructor(obj._data.apply(lambda x: isnull(x.values)))
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
    elif isinstance(obj, (ABCSeries, np.ndarray)):
        return _isnull_ndarraylike_old(obj)
    elif isinstance(obj, ABCGeneric):
        return obj._constructor(obj._data.apply(
            lambda x: _isnull_old(x.values)))
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
        # Working around NumPy ticket 1542
        shape = values.shape

        if dtype.kind in ('S', 'U'):
            result = np.zeros(values.shape, dtype=bool)
        else:
            result = np.empty(shape, dtype=bool)
            vec = lib.isnullobj(values.ravel())
            result[:] = vec.reshape(shape)

    elif dtype in _DATELIKE_DTYPES:
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
        result = -np.isfinite(values)

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
    """
    res = isnull(obj)
    if np.isscalar(res):
        return not res
    return -res

def _is_null_datelike_scalar(other):
    """ test whether the object is a null datelike, e.g. Nat
    but guard against passing a non-scalar """
    return (np.isscalar(other) and (isnull(other) or other == tslib.iNaT)) or other is pd.NaT or other is None

def array_equivalent(left, right):
    """
    True if two arrays, left and right, have equal non-NaN elements, and NaNs in
    corresponding locations.  False otherwise. It is assumed that left and right
    are NumPy arrays of the same dtype. The behavior of this function
    (particularly with respect to NaNs) is not defined if the dtypes are
    different.

    Parameters
    ----------
    left, right : ndarrays

    Returns
    -------
    b : bool
        Returns True if the arrays are equivalent.

    Examples
    --------
    >>> array_equivalent(np.array([1, 2, nan]), np.array([1, 2, nan]))
    True
    >>> array_equivalent(np.array([1, nan, 2]), np.array([1, 2, nan]))
    False
    """
    if left.shape != right.shape: return False
    # NaNs occur only in object arrays, float or complex arrays.
    if not issubclass(left.dtype.type, (np.floating, np.complexfloating)):
        return np.array_equal(left, right)
    return  ((left == right) | (np.isnan(left) & np.isnan(right))).all()

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
    nonna = values_to_mask[-na_mask]

    mask = None
    for x in nonna:
        if mask is None:
            mask = arr == x

            # if x is a string and mask is not, then we get a scalar
            # return value, which is not good
            if not isinstance(mask, np.ndarray):
                m = mask
                mask = np.empty(arr.shape, dtype=np.bool)
                mask.fill(m)
        else:
            mask = mask | (arr == x)

    if na_mask.any():
        if mask is None:
            mask = isnull(arr)
        else:
            mask = mask | isnull(arr)

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
    func(arr, indexer, out, fill_value)
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
    dtype = arr.dtype
    na = np.nan

    if is_timedelta64_dtype(arr) or is_datetime64_dtype(arr):
        dtype = 'timedelta64[ns]'
        arr = arr.view('i8')
        na = tslib.iNaT
    elif issubclass(dtype.type, np.integer):
        dtype = np.float64
    elif issubclass(dtype.type, np.bool_):
        dtype = np.object_

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
        if dtype == 'timedelta64[ns]':
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

    return out_arr


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
                r = bool(r)
            elif dtype.kind == 'f':
                r = float(r)
            elif dtype.kind == 'i':
                r = int(r)
        except:
            pass

        return r

    return np.array([conv(r, dtype) for r, dtype in zip(result, dtypes)])


def _infer_dtype_from_scalar(val):
    """ interpret the dtype from a scalar, upcast floats and ints
        return the new value and the dtype """

    dtype = np.object_

    # a 1-element ndarray
    if isinstance(val, pa.Array):
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


def _maybe_cast_scalar(dtype, value):
    """ if we a scalar value and are casting to a dtype that needs nan -> NaT
    conversion
    """
    if np.isscalar(value) and dtype in _DATELIKE_DTYPES and isnull(value):
        return tslib.iNaT
    return value


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
    if issubclass(np.dtype(dtype).type, compat.string_types):
        dtype = np.object_

    return dtype, fill_value


def _maybe_upcast_putmask(result, mask, other, dtype=None, change=None):
    """ a safe version of put mask that (potentially upcasts the result
    return the result
    if change is not None, then MUTATE the change (and change the dtype)
    return a changed flag
    """

    if mask.any():

        other = _maybe_cast_scalar(result.dtype, other)

        def changeit():

            # try to directly set by expanding our array to full
            # length of the boolean
            try:
                om = other[mask]
                om_at = om.astype(result.dtype)
                if (om == om_at).all():
                    new_other = result.values.copy()
                    new_other[mask] = om_at
                    result[:] = new_other
                    return result, False
            except:
                pass

            # we are forced to change the dtype of the result as the input
            # isn't compatible
            r, fill_value = _maybe_upcast(
                result, fill_value=other, dtype=dtype, copy=True)
            np.putmask(r, mask, other)

            # we need to actually change the dtype here
            if change is not None:

                # if we are trying to do something unsafe
                # like put a bigger dtype in a smaller one, use the smaller one
                # pragma: no cover
                if change.dtype.itemsize < r.dtype.itemsize:
                    raise AssertionError(
                        "cannot change dtype of input to smaller size")
                change.dtype = r.dtype
                change[:] = r

            return r, True

        # we want to decide whether putmask will work
        # if we have nans in the False portion of our mask then we need to
        # upcast (possibily) otherwise we DON't want to upcast (e.g. if we are
        # have values, say integers in the success portion then its ok to not
        # upcast)
        new_dtype, fill_value = _maybe_promote(result.dtype, other)
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
            np.putmask(result, mask, other)
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

    if np.isscalar(result) or not len(result):
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

        # don't allow upcasts here
        if dtype.kind == result.dtype.kind:
            if result.dtype.itemsize <= dtype.itemsize:
                return result

        if issubclass(dtype.type, np.floating):
            return result.astype(dtype)
        elif dtype == np.bool_ or issubclass(dtype.type, np.integer):

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


def _fill_zeros(result, y, fill):
    """ if we have an integer value (or array in y)
    and we have 0's, fill them with the fill,
    return the result
    """

    if fill is not None:
        if not isinstance(y, np.ndarray):
            dtype, value = _infer_dtype_from_scalar(y)
            y = pa.empty(result.shape, dtype=dtype)
            y.fill(value)

        if is_integer_dtype(y):

            mask = y.ravel() == 0
            if mask.any():
                shape = result.shape
                result, changed = _maybe_upcast_putmask(
                    result.ravel(), mask, fill)
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


def pad_1d(values, limit=None, mask=None):

    dtype = values.dtype.name
    _method = None
    if is_float_dtype(values):
        _method = getattr(algos, 'pad_inplace_%s' % dtype, None)
    elif is_datetime64_dtype(values):
        _method = _pad_1d_datetime
    elif is_integer_dtype(values):
        values = _ensure_float64(values)
        _method = algos.pad_inplace_float64
    elif values.dtype == np.object_:
        _method = algos.pad_inplace_object

    if _method is None:
        raise ValueError('Invalid dtype for pad_1d [%s]' % dtype)

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)
    _method(values, mask, limit=limit)
    return values


def backfill_1d(values, limit=None, mask=None):

    dtype = values.dtype.name
    _method = None
    if is_float_dtype(values):
        _method = getattr(algos, 'backfill_inplace_%s' % dtype, None)
    elif is_datetime64_dtype(values):
        _method = _backfill_1d_datetime
    elif is_integer_dtype(values):
        values = _ensure_float64(values)
        _method = algos.backfill_inplace_float64
    elif values.dtype == np.object_:
        _method = algos.backfill_inplace_object

    if _method is None:
        raise ValueError('Invalid dtype for backfill_1d [%s]' % dtype)

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)

    _method(values, mask, limit=limit)
    return values


def pad_2d(values, limit=None, mask=None):

    dtype = values.dtype.name
    _method = None
    if is_float_dtype(values):
        _method = getattr(algos, 'pad_2d_inplace_%s' % dtype, None)
    elif is_datetime64_dtype(values):
        _method = _pad_2d_datetime
    elif is_integer_dtype(values):
        values = _ensure_float64(values)
        _method = algos.pad_2d_inplace_float64
    elif values.dtype == np.object_:
        _method = algos.pad_2d_inplace_object

    if _method is None:
        raise ValueError('Invalid dtype for pad_2d [%s]' % dtype)

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)

    if np.all(values.shape):
        _method(values, mask, limit=limit)
    else:
        # for test coverage
        pass
    return values


def backfill_2d(values, limit=None, mask=None):

    dtype = values.dtype.name
    _method = None
    if is_float_dtype(values):
        _method = getattr(algos, 'backfill_2d_inplace_%s' % dtype, None)
    elif is_datetime64_dtype(values):
        _method = _backfill_2d_datetime
    elif is_integer_dtype(values):
        values = _ensure_float64(values)
        _method = algos.backfill_2d_inplace_float64
    elif values.dtype == np.object_:
        _method = algos.backfill_2d_inplace_object

    if _method is None:
        raise ValueError('Invalid dtype for backfill_2d [%s]' % dtype)

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)

    if np.all(values.shape):
        _method(values, mask, limit=limit)
    else:
        # for test coverage
        pass
    return values


def _clean_interp_method(method, order=None, **kwargs):
    valid = ['linear', 'time', 'values', 'nearest', 'zero', 'slinear',
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
                   fill_value=None, bounds_error=False, **kwargs):
    """
    Logic for the 1-d interpolation.  The result should be 1-d, inputs
    xvalues and yvalues will each be 1-d arrays of the same length.

    Bounds_error is currently hardcoded to False since non-scipy ones don't
    take it as an argumnet.
    """
    # Treat the original, non-scipy methods first.

    invalid = isnull(yvalues)
    valid = ~invalid

    valid_y = yvalues[valid]
    valid_x = xvalues[valid]
    new_x = xvalues[invalid]

    if method == 'time':
        if not getattr(xvalues, 'is_all_dates', None):
        # if not issubclass(xvalues.dtype.type, np.datetime64):
            raise ValueError('time-weighted interpolation only works '
                             'on Series or DataFrames with a '
                             'DatetimeIndex')
        method = 'values'

    def _interp_limit(invalid, limit):
        """mask off values that won't be filled since they exceed the limit"""
        all_nans = np.where(invalid)[0]
        violate = [invalid[x:x + limit + 1] for x in all_nans]
        violate = np.array([x.all() & (x.size > limit) for x in violate])
        return all_nans[violate] + limit

    xvalues = getattr(xvalues, 'values', xvalues)
    yvalues = getattr(yvalues, 'values', yvalues)

    if limit:
        violate_limit = _interp_limit(invalid, limit)
    if valid.any():
        firstIndex = valid.argmax()
        valid = valid[firstIndex:]
        invalid = invalid[firstIndex:]
        result = yvalues.copy()
        if valid.all():
            return yvalues
    else:
        # have to call np.array(xvalues) since xvalues could be an Index
        # which cant be mutated
        result = np.empty_like(np.array(xvalues), dtype=np.float64)
        result.fill(np.nan)
        return result

    if method in ['linear', 'time', 'values']:
        if method in ('values', 'index'):
            inds = np.asarray(xvalues)
            # hack for DatetimeIndex, #1646
            if issubclass(inds.dtype.type, np.datetime64):
                inds = inds.view(pa.int64)

            if inds.dtype == np.object_:
                inds = lib.maybe_convert_objects(inds)
        else:
            inds = xvalues

        inds = inds[firstIndex:]

        result[firstIndex:][invalid] = np.interp(inds[invalid], inds[valid],
                                                 yvalues[firstIndex:][valid])

        if limit:
            result[violate_limit] = np.nan
        return result

    sp_methods = ['nearest', 'zero', 'slinear', 'quadratic', 'cubic',
                  'barycentric', 'krogh', 'spline', 'polynomial',
                  'piecewise_polynomial', 'pchip']
    if method in sp_methods:
        new_x = new_x[firstIndex:]
        xvalues = xvalues[firstIndex:]

        result[firstIndex:][invalid] = _interpolate_scipy_wrapper(
            valid_x, valid_y, new_x, method=method, fill_value=fill_value,
            bounds_error=bounds_error, **kwargs)
        if limit:
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
        x, new_x = x.values.astype('i8'), new_x.astype('i8')

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
        terp = interpolate.UnivariateSpline(x, y, k=order)
        new_y = terp(new_x)
    else:
        method = alt_methods[method]
        new_y = method(x, y, new_x)
    return new_y


def interpolate_2d(values, method='pad', axis=0, limit=None, fill_value=None):
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
        values = transf(pad_2d(transf(values), limit=limit, mask=mask))
    else:
        values = transf(backfill_2d(transf(values), limit=limit, mask=mask))

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


def _maybe_box(indexer, values, obj, key):

    # if we have multiples coming back, box em
    if isinstance(values, np.ndarray):
        return obj[indexer.get_loc(key)]

    # return the value
    return values


def _values_from_object(o):
    """ return my values or the object if we are say an ndarray """
    f = getattr(o, 'get_values', None)
    if f is not None:
        o = f()
    return o


def _possibly_convert_objects(values, convert_dates=True,
                              convert_numeric=True,
                              convert_timedeltas=True):
    """ if we have an object dtype, try to coerce dates and/or numbers """

    # if we have passed in a list or scalar
    if isinstance(values, (list, tuple)):
        values = np.array(values, dtype=np.object_)
    if not hasattr(values, 'dtype'):
        values = np.array([values], dtype=np.object_)

    # convert dates
    if convert_dates and values.dtype == np.object_:

        # we take an aggressive stance and convert to datetime64[ns]
        if convert_dates == 'coerce':
            new_values = _possibly_cast_to_datetime(
                values, 'M8[ns]', coerce=True)

            # if we are all nans then leave me alone
            if not isnull(new_values).all():
                values = new_values

        else:
            values = lib.maybe_convert_objects(
                values, convert_datetime=convert_dates)

    # convert timedeltas
    if convert_timedeltas and values.dtype == np.object_:

        if convert_timedeltas == 'coerce':
            from pandas.tseries.timedeltas import \
                 _possibly_cast_to_timedelta
            values = _possibly_cast_to_timedelta(values, coerce=True)

            # if we are all nans then leave me alone
            if not isnull(new_values).all():
                values = new_values

        else:
            values = lib.maybe_convert_objects(
                values, convert_timedelta=convert_timedeltas)

    # convert to numeric
    if values.dtype == np.object_:
        if convert_numeric:
            try:
                new_values = lib.maybe_convert_numeric(
                    values, set(), coerce_numeric=True)

                # if we are all nans then leave me alone
                if not isnull(new_values).all():
                    values = new_values

            except:
                pass
        else:

            # soft-conversion
            values = lib.maybe_convert_objects(values)

    return values


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
        if hasattr(values, 'values'):
            values = values.values
        values = lib.maybe_convert_objects(values)

    return values


def _possibly_cast_to_datetime(value, dtype, coerce=False):
    """ try to cast the array/value to a datetimelike dtype, converting float
    nan to iNaT
    """

    if dtype is not None:
        if isinstance(dtype, compat.string_types):
            dtype = np.dtype(dtype)

        is_datetime64 = is_datetime64_dtype(dtype)
        is_timedelta64 = is_timedelta64_dtype(dtype)

        if is_datetime64 or is_timedelta64:

            # force the dtype if needed
            if is_datetime64 and dtype != _NS_DTYPE:
                if dtype.name == 'datetime64[ns]':
                    dtype = _NS_DTYPE
                else:
                    raise TypeError(
                        "cannot convert datetimelike to dtype [%s]" % dtype)
            elif is_timedelta64 and dtype != _TD_DTYPE:
                if dtype.name == 'timedelta64[ns]':
                    dtype = _TD_DTYPE
                else:
                    raise TypeError(
                        "cannot convert timedeltalike to dtype [%s]" % dtype)

            if np.isscalar(value):
                if value == tslib.iNaT or isnull(value):
                    value = tslib.iNaT
            else:
                value = np.array(value)

                # have a scalar array-like (e.g. NaT)
                if value.ndim == 0:
                    value = tslib.iNaT

                # we have an array of datetime or timedeltas & nulls
                elif np.prod(value.shape) and value.dtype != dtype:
                    try:
                        if is_datetime64:
                            from pandas.tseries.tools import to_datetime
                            value = to_datetime(value, coerce=coerce).values
                        elif is_timedelta64:
                            from pandas.tseries.timedeltas import \
                                _possibly_cast_to_timedelta
                            value = _possibly_cast_to_timedelta(value, coerce='compat')
                    except:
                        pass

    else:

        is_array = isinstance(value, np.ndarray)

        # catch a datetime/timedelta that is not of ns variety
        # and no coercion specified
        if (is_array and value.dtype.kind in ['M','m']):
            dtype = value.dtype

            if dtype.kind == 'M' and dtype != _NS_DTYPE:
                try:
                    value = tslib.array_to_datetime(value)
                except:
                    raise

            elif dtype.kind == 'm' and dtype != _TD_DTYPE:
                from pandas.tseries.timedeltas import \
                     _possibly_cast_to_timedelta
                value = _possibly_cast_to_timedelta(value, coerce='compat')

        # only do this if we have an array and the dtype of the array is not
        # setup already we are not an integer/object, so don't bother with this
        # conversion
        elif (is_array and not (
            issubclass(value.dtype.type, np.integer) or
            value.dtype == np.object_)):
            pass

        else:
            # we might have a array (or single object) that is datetime like,
            # and no dtype is passed don't change the value unless we find a
            # datetime set
            v = value
            if not is_list_like(v):
                v = [v]
            if len(v):
                inferred_type = lib.infer_dtype(v)
                if inferred_type in ['datetime', 'datetime64']:
                    try:
                        value = tslib.array_to_datetime(np.array(v))
                    except:
                        pass
                elif inferred_type in ['timedelta', 'timedelta64']:
                    from pandas.tseries.timedeltas import \
                        _possibly_cast_to_timedelta
                    value = _possibly_cast_to_timedelta(value, coerce='compat')

    return value


def _is_bool_indexer(key):
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
    result = values.view(Int64Index)
    result.name = None
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


def rands(n):
    """Generates a random alphanumeric string of length *n*"""
    from random import Random
    import string
    return ''.join(Random().sample(string.ascii_letters + string.digits, n))


def adjoin(space, *lists):
    """
    Glues together two sets of strings using the amount of space requested.
    The idea is to prettify.
    """
    out_lines = []
    newLists = []
    lengths = [max(map(len, x)) + space for x in lists[:-1]]

    # not the last one
    lengths.append(max(map(len, lists[-1])))

    maxLen = max(map(len, lists))
    for i, lst in enumerate(lists):
        nl = [x.ljust(lengths[i]) for x in lst]
        nl.extend([' ' * lengths[i]] * (maxLen - len(lst)))
        newLists.append(nl)
    toJoin = zip(*newLists)
    for lines in toJoin:
        out_lines.append(_join_unicode(lines))
    return _join_unicode(out_lines, sep='\n')


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
    >>> iterpairs([1, 2, 3, 4])
    [(1, 2), (2, 3), (3, 4)
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


def _shift_indexer(N, periods):
    # small reusable utility
    indexer = np.zeros(N, dtype=int)

    if periods > 0:
        indexer[periods:] = np.arange(N - periods)
    else:
        indexer[:periods] = np.arange(-periods, N)

    return indexer


def _asarray_tuplesafe(values, dtype=None):
    from pandas.core.index import Index

    if not isinstance(values, (list, tuple, np.ndarray)):
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


def is_bool(obj):
    return isinstance(obj, (bool, np.bool_))


def is_integer(obj):
    return isinstance(obj, (numbers.Integral, np.integer))


def is_float(obj):
    return isinstance(obj, (float, np.floating))


def is_complex(obj):
    return isinstance(obj, (numbers.Complex, np.complexfloating))


def is_iterator(obj):
    # python 3 generators have __next__ instead of next
    return hasattr(obj, 'next') or hasattr(obj, '__next__')


def is_number(obj):
    return isinstance(obj, (numbers.Number, np.number))


def is_integer_dtype(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        tipo = arr_or_dtype.type
    else:
        tipo = arr_or_dtype.dtype.type
    return (issubclass(tipo, np.integer) and not
            (issubclass(tipo, np.datetime64) or
             issubclass(tipo, np.timedelta64)))


def _is_int_or_datetime_dtype(arr_or_dtype):
    # also timedelta64
    if isinstance(arr_or_dtype, np.dtype):
        tipo = arr_or_dtype.type
    else:
        tipo = arr_or_dtype.dtype.type
    return issubclass(tipo, np.integer)


def is_datetime64_dtype(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        tipo = arr_or_dtype.type
    elif isinstance(arr_or_dtype, type):
        tipo = np.dtype(arr_or_dtype).type
    else:
        tipo = arr_or_dtype.dtype.type
    return issubclass(tipo, np.datetime64)


def is_datetime64_ns_dtype(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        tipo = arr_or_dtype
    elif isinstance(arr_or_dtype, type):
        tipo = np.dtype(arr_or_dtype)
    else:
        tipo = arr_or_dtype.dtype
    return tipo == _NS_DTYPE


def is_timedelta64_dtype(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        tipo = arr_or_dtype.type
    elif isinstance(arr_or_dtype, type):
        tipo = np.dtype(arr_or_dtype).type
    else:
        tipo = arr_or_dtype.dtype.type
    return issubclass(tipo, np.timedelta64)


def needs_i8_conversion(arr_or_dtype):
    return (is_datetime64_dtype(arr_or_dtype) or
            is_timedelta64_dtype(arr_or_dtype))


def is_numeric_dtype(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        tipo = arr_or_dtype.type
    else:
        tipo = arr_or_dtype.dtype.type
    return (issubclass(tipo, (np.number, np.bool_))
            and not issubclass(tipo, (np.datetime64, np.timedelta64)))

def is_float_dtype(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        tipo = arr_or_dtype.type
    else:
        tipo = arr_or_dtype.dtype.type
    return issubclass(tipo, np.floating)


def is_complex_dtype(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        tipo = arr_or_dtype.type
    else:
        tipo = arr_or_dtype.dtype.type
    return issubclass(tipo, np.complexfloating)


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


def _is_sequence(x):
    try:
        iter(x)
        len(x)  # it has a length
        return not isinstance(x, compat.string_and_binary_types)
    except (TypeError, AttributeError):
        return False


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
        dtype = np.dtype(dtype)

    if is_datetime64_dtype(arr):
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
            return arr.astype(object)

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
    elif issubclass(dtype.type, compat.string_types):
        return lib.astype_str(arr.ravel()).reshape(arr.shape)

    if copy:
        return arr.astype(dtype)
    return arr.view(dtype)


def _clean_fill_method(method):
    if method is None:
        return None
    method = method.lower()
    if method == 'ffill':
        method = 'pad'
    if method == 'bfill':
        method = 'backfill'
    if method not in ['pad', 'backfill']:
        msg = ('Invalid fill method. Expecting pad (ffill) or backfill '
               '(bfill). Got %s' % method)
        raise ValueError(msg)
    return method


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
        if compat.PY3_2:
            # gzip and bz2 don't work with TextIOWrapper in 3.2
            encoding = encoding or get_option('display.encoding')
            f = StringIO(f.read().decode(encoding))
        elif compat.PY3:
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


def _concat_compat(to_concat, axis=0):
    # filter empty arrays
    nonempty = [x for x in to_concat if x.shape[axis] > 0]

    # If all arrays are empty, there's nothing to convert, just short-cut to
    # the concatenation, #3121.
    #
    # Creating an empty array directly is tempting, but the winnings would be
    # marginal given that it would still require shape & dtype calculation and
    # np.concatenate which has them both implemented is compiled.
    if nonempty:
        is_datetime64 = [x.dtype == _NS_DTYPE for x in nonempty]
        if all(is_datetime64):
            # work around NumPy 1.6 bug
            new_values = np.concatenate([x.view(np.int64) for x in nonempty],
                                        axis=axis)
            return new_values.view(_NS_DTYPE)
        elif any(is_datetime64):
            to_concat = [_to_pydatetime(x) for x in nonempty]

    return np.concatenate(to_concat, axis=axis)


def _to_pydatetime(x):
    if x.dtype == _NS_DTYPE:
        shape = x.shape
        x = tslib.ints_to_pydatetime(x.view(np.int64).ravel())
        x = x.reshape(shape)

    return x


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


def _pprint_seq(seq, _nest_lvl=0, **kwds):
    """
    internal. pprinter for iterables. you should probably use pprint_thing()
    rather then calling this directly.

    bounds length of printed sequence, depending on options
    """
    if isinstance(seq, set):
        fmt = u("set([%s])")
    else:
        fmt = u("[%s]") if hasattr(seq, '__setitem__') else u("(%s)")

    nitems = get_option("max_seq_items") or len(seq)

    s = iter(seq)
    r = []
    for i in range(min(nitems, len(seq))):  # handle sets, no slicing
        r.append(pprint_thing(next(s), _nest_lvl + 1, **kwds))
    body = ", ".join(r)

    if nitems < len(seq):
        body += ", ..."
    elif isinstance(seq, tuple) and len(seq) == 1:
        body += ','

    return fmt % body


def _pprint_dict(seq, _nest_lvl=0, **kwds):
    """
    internal. pprinter for iterables. you should probably use pprint_thing()
    rather then calling this directly.
    """
    fmt = u("{%s}")
    pairs = []

    pfmt = u("%s: %s")

    nitems = get_option("max_seq_items") or len(seq)

    for k, v in list(seq.items())[:nitems]:
        pairs.append(pfmt % (pprint_thing(k, _nest_lvl + 1, **kwds),
                             pprint_thing(v, _nest_lvl + 1, **kwds)))

    if nitems < len(seq):
        return fmt % (", ".join(pairs) + ", ...")
    else:
        return fmt % ", ".join(pairs)


def pprint_thing(thing, _nest_lvl=0, escape_chars=None, default_escapes=False,
                 quote_strings=False):
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
        result = _pprint_dict(thing, _nest_lvl, quote_strings=True)
    elif _is_sequence(thing) and _nest_lvl < \
            get_option("display.pprint_nest_depth"):
        result = _pprint_seq(thing, _nest_lvl, escape_chars=escape_chars,
                             quote_strings=quote_strings)
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


def load(path):  # TODO remove in 0.13
    """
    Load pickled pandas object (or any other pickled object) from the specified
    file path

    Warning: Loading pickled data received from untrusted sources can be
    unsafe. See: http://docs.python.org/2.7/library/pickle.html

    Parameters
    ----------
    path : string
        File path

    Returns
    -------
    unpickled : type of object stored in file
    """
    import warnings
    warnings.warn("load is deprecated, use read_pickle", FutureWarning)
    from pandas.io.pickle import read_pickle
    return read_pickle(path)


def save(obj, path):  # TODO remove in 0.13
    """
    Pickle (serialize) object to input file path

    Parameters
    ----------
    obj : any object
    path : string
        File path
    """
    import warnings
    warnings.warn("save is deprecated, use obj.to_pickle", FutureWarning)
    from pandas.io.pickle import to_pickle
    return to_pickle(obj, path)


def _maybe_match_name(a, b):
    a_name = getattr(a, 'name', None)
    b_name = getattr(b, 'name', None)
    if a_name == b_name:
        return a_name
    return None
