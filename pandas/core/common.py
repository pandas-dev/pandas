"""
Misc tools for implementing data structures
"""

import itertools
import re
from datetime import datetime

from numpy.lib.format import read_array, write_array
import numpy as np

import pandas.algos as algos
import pandas.lib as lib
import pandas.tslib as tslib

from pandas.util import py3compat
import codecs
import csv

from pandas.util.py3compat import StringIO, BytesIO

from pandas.core.config import get_option
from pandas.core import array as pa

# XXX: HACK for NumPy 1.5.1 to suppress warnings
try:
    np.seterr(all='ignore')
    # np.set_printoptions(suppress=True)
except Exception:  # pragma: no cover
    pass


class PandasError(Exception):
    pass


class AmbiguousIndexError(PandasError, KeyError):
    pass


_POSSIBLY_CAST_DTYPES = set([ np.dtype(t) for t in ['M8[ns]','m8[ns]','O','int8','uint8','int16','uint16','int32','uint32','int64','uint64'] ])
_NS_DTYPE = np.dtype('M8[ns]')
_TD_DTYPE = np.dtype('m8[ns]')
_INT64_DTYPE = np.dtype(np.int64)
_DATELIKE_DTYPES = set([ np.dtype(t) for t in ['M8[ns]','m8[ns]'] ])

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

    from pandas.core.generic import PandasContainer
    if isinstance(obj, np.ndarray):
        return _isnull_ndarraylike(obj)
    elif isinstance(obj, PandasContainer):
        # TODO: optimize for DataFrame, etc.
        return obj.apply(isnull)
    elif isinstance(obj, list) or hasattr(obj, '__array__'):
        return _isnull_ndarraylike(obj)
    else:
        return obj is None


def _isnull_old(obj):
    '''
    Detect missing values. Treat None, NaN, INF, -INF as null.

    Parameters
    ----------
    arr: ndarray or object value

    Returns
    -------
    boolean ndarray or boolean
    '''
    if lib.isscalar(obj):
        return lib.checknull_old(obj)

    from pandas.core.generic import PandasContainer
    if isinstance(obj, np.ndarray):
        return _isnull_ndarraylike_old(obj)
    elif isinstance(obj, PandasContainer):
        # TODO: optimize for DataFrame, etc.
        return obj.apply(_isnull_old)
    elif isinstance(obj, list) or hasattr(obj, '__array__'):
        return _isnull_ndarraylike_old(obj)
    else:
        return obj is None

_isnull = _isnull_new

def _use_inf_as_null(key):
    '''Option change callback for null/inf behaviour
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
    '''
    flag = get_option(key)
    if flag:
        globals()['_isnull'] = _isnull_old
    else:
        globals()['_isnull'] = _isnull_new


def _isnull_ndarraylike(obj):
    from pandas import Series
    values = np.asarray(obj)

    if values.dtype.kind in ('O', 'S', 'U'):
        # Working around NumPy ticket 1542
        shape = values.shape

        if values.dtype.kind in ('S', 'U'):
            result = np.zeros(values.shape, dtype=bool)
        else:
            result = np.empty(shape, dtype=bool)
            vec = lib.isnullobj(values.ravel())
            result[:] = vec.reshape(shape)

        if isinstance(obj, Series):
            result = Series(result, index=obj.index, copy=False)
    elif values.dtype == np.dtype('M8[ns]'):
        # this is the NaT pattern
        result = values.view('i8') == tslib.iNaT
    elif values.dtype == np.dtype('m8[ns]'):
        # this is the NaT pattern
        result = values.view('i8') == tslib.iNaT
    else:
        # -np.isfinite(obj)
        result = np.isnan(obj)
    return result


def _isnull_ndarraylike_old(obj):
    from pandas import Series
    values = np.asarray(obj)

    if values.dtype.kind in ('O', 'S', 'U'):
        # Working around NumPy ticket 1542
        shape = values.shape

        if values.dtype.kind in ('S', 'U'):
            result = np.zeros(values.shape, dtype=bool)
        else:
            result = np.empty(shape, dtype=bool)
            vec = lib.isnullobj_old(values.ravel())
            result[:] = vec.reshape(shape)

        if isinstance(obj, Series):
            result = Series(result, index=obj.index, copy=False)
    elif values.dtype == np.dtype('M8[ns]'):
        # this is the NaT pattern
        result = values.view('i8') == tslib.iNaT
    else:
        result = -np.isfinite(obj)
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
            if not isinstance(mask,np.ndarray):
                m = mask
                mask = np.empty(arr.shape,dtype=np.bool)
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
        u = row_idx[i]
        for j in range(len(col_idx)):
            v = col_idx[j]
            out[i, j] = arr[u, v]


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
    ('datetime64[ns]','datetime64[ns]'):
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
    ('datetime64[ns]','datetime64[ns]'):
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
    ('datetime64[ns]','datetime64[ns]'):
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
    ('datetime64[ns]','datetime64[ns]'):
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
                        raise Exception('Incompatible type for fill_value')
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
                        raise Exception('Incompatible type for fill_value')
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

            result = res-lag
            result[mask] = na
            out_arr[res_indexer] = result
        else:
            out_arr[res_indexer] = arr[res_indexer] - arr[lag_indexer]

    return out_arr


def _infer_dtype_from_scalar(val):
    """ interpret the dtype from a scalar, upcast floats and ints
        return the new value and the dtype """

    dtype = np.object_

    # a 1-element ndarray
    if isinstance(val, pa.Array):
        if val.ndim != 0:
            raise ValueError("invalid ndarray passed to _infer_dtype_from_scalar")

        dtype = val.dtype
        val   = val.item()

    elif isinstance(val, basestring):

        # If we create an empty array using a string to infer
        # the dtype, NumPy will only allocate one character per entry
        # so this is kind of bad. Alternately we could use np.repeat
        # instead of np.empty (but then you still don't want things
        # coming out as np.str_!

        dtype = np.object_

    elif isinstance(val, np.datetime64):
        # ugly hacklet
        val   = lib.Timestamp(val).value
        dtype = np.dtype('M8[ns]')

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
    """ if we a scalar value and are casting to a dtype that needs nan -> NaT conversion """
    if np.isscalar(value) and dtype in _DATELIKE_DTYPES and isnull(value):
        return tslib.iNaT
    return value

def _maybe_promote(dtype, fill_value=np.nan):

    # if we passed an array here, determine the fill value by dtype
    if isinstance(fill_value,np.ndarray):
        if issubclass(fill_value.dtype.type, (np.datetime64,np.timedelta64)):
            fill_value = tslib.iNaT
        else:

            # we need to change to object type as our
            # fill_value is of object type
            if fill_value.dtype == np.object_:
                dtype = np.dtype(np.object_)
            fill_value = np.nan

    # returns tuple of (dtype, fill_value)
    if issubclass(dtype.type, (np.datetime64,np.timedelta64)):
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
                    # the proper thing to do here would probably be to upcast to
                    # object (but numpy 1.6.1 doesn't do this properly)
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
    if issubclass(np.dtype(dtype).type, basestring):
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

            # we are forced to change the dtype of the result as the input isn't compatible
            r, fill_value = _maybe_upcast(result, fill_value=other, dtype=dtype, copy=True)
            np.putmask(r, mask, other)

            # we need to actually change the dtype here
            if change is not None:

                # if we are trying to do something unsafe
                # like put a bigger dtype in a smaller one, use the smaller one
                if change.dtype.itemsize < r.dtype.itemsize:
                    raise Exception("cannot change dtype of input to smaller size")
                change.dtype = r.dtype
                change[:] = r

            return r, True

        # we want to decide whether putmask will work
        # if we have nans in the False portion of our mask then we need to upcast (possibily)
        # otherwise we DON't want to upcast (e.g. if we are have values, say integers in
        # the success portion then its ok to not upcast)
        new_dtype, fill_value = _maybe_promote(result.dtype,other)
        if new_dtype != result.dtype:

            # we have a scalar or len 0 ndarray
            # and its nan and we are changing some values
            if np.isscalar(other) or (isinstance(other,np.ndarray) and other.ndim < 1):
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

def _maybe_upcast_indexer(result, indexer, other, dtype=None):
    """ a safe version of setitem that (potentially upcasts the result
        return the result and a changed flag
        """

    other = _maybe_cast_scalar(result.dtype, other)
    original_dtype = result.dtype
    def changeit():
        # our type is wrong here, need to upcast
        r, fill_value = _maybe_upcast(result, fill_value=other, dtype=dtype, copy=True)
        try:
            r[indexer] = other
        except:

            # if we hit this then we still have an incompatible type
            r[indexer] = fill_value

        # if we have changed to floats, might want to cast back if we can
        r = _possibly_downcast_to_dtype(r,original_dtype)
        return r, True

    new_dtype, fill_value = _maybe_promote(original_dtype,other)
    if new_dtype != result.dtype:
        return changeit()

    try:
        result[indexer] = other
    except:
        return changeit()

    return result, False

def _maybe_upcast(values, fill_value=np.nan, dtype=None, copy=False):
    """ provide explicty type promotion and coercion

        Parameters
        ----------
        values : the ndarray that we want to maybe upcast
        fill_value : what we want to fill with
        dtype : if None, then use the dtype of the values, else coerce to this type
        copy : if True always make a copy even if no upcast is required """

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
        or could be an astype of float64->float32 """

    if not isinstance(result, np.ndarray):
        return result

    try:
        if issubclass(dtype.type,np.floating):
            return result.astype(dtype)
        elif dtype == np.bool_ or issubclass(dtype.type,np.integer):
            if issubclass(result.dtype.type, np.number) and notnull(result).all():
                new_result = result.astype(dtype)
                if (new_result == result).all():
                    return new_result
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
        return the result """

    if fill is not None:
        if not isinstance(y, np.ndarray):
            dtype, value = _infer_dtype_from_scalar(y)
            y = pa.empty(result.shape,dtype=dtype)
            y.fill(value)

        if is_integer_dtype(y):

            mask = y.ravel() == 0
            if mask.any():
                shape = result.shape
                result, changed = _maybe_upcast_putmask(result.ravel(),mask,fill)
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

    dtype   = values.dtype.name
    _method = None
    if is_float_dtype(values):
        _method = getattr(algos,'pad_inplace_%s' % dtype,None)
    elif is_datetime64_dtype(values):
        _method = _pad_1d_datetime
    elif values.dtype == np.object_:
        _method = algos.pad_inplace_object

    if _method is None:
        raise ValueError('Invalid dtype for pad_1d [%s]' % dtype)

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)
    _method(values, mask, limit=limit)


def backfill_1d(values, limit=None, mask=None):

    dtype   = values.dtype.name
    _method = None
    if is_float_dtype(values):
        _method = getattr(algos,'backfill_inplace_%s' % dtype,None)
    elif is_datetime64_dtype(values):
        _method = _backfill_1d_datetime
    elif values.dtype == np.object_:
        _method = algos.backfill_inplace_object

    if _method is None:
        raise ValueError('Invalid dtype for backfill_1d [%s]' % dtype)

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)

    _method(values, mask, limit=limit)


def pad_2d(values, limit=None, mask=None):

    dtype   = values.dtype.name
    _method = None
    if is_float_dtype(values):
        _method = getattr(algos,'pad_2d_inplace_%s' % dtype,None)
    elif is_datetime64_dtype(values):
        _method = _pad_2d_datetime
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


def backfill_2d(values, limit=None, mask=None):

    dtype   = values.dtype.name
    _method = None
    if is_float_dtype(values):
        _method = getattr(algos,'backfill_2d_inplace_%s' % dtype,None)
    elif is_datetime64_dtype(values):
        _method = _backfill_2d_datetime
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


def _consensus_name_attr(objs):
    name = objs[0].name
    for obj in objs[1:]:
        if obj.name != name:
            return None
    return name

#----------------------------------------------------------------------
# Lots of little utilities


def _possibly_convert_objects(values, convert_dates=True, convert_numeric=True):
    """ if we have an object dtype, try to coerce dates and/or numers """

    # convert dates
    if convert_dates and values.dtype == np.object_:

        # we take an aggressive stance and convert to datetime64[ns]
        if convert_dates == 'coerce':
            new_values = _possibly_cast_to_datetime(values, 'M8[ns]', coerce = True)

            # if we are all nans then leave me alone
            if not isnull(new_values).all():
                values = new_values

        else:
            values = lib.maybe_convert_objects(values, convert_datetime=convert_dates)

    # convert to numeric
    if convert_numeric and values.dtype == np.object_:
        try:
            new_values = lib.maybe_convert_numeric(values,set(),coerce_numeric=True)

            # if we are all nans then leave me alone
            if not isnull(new_values).all():
                values = new_values

        except:
            pass

    return values

def _possibly_castable(arr):
    return arr.dtype not in _POSSIBLY_CAST_DTYPES

def _possibly_convert_platform(values):
    """ try to do platform conversion, allow ndarray or list here """

    if isinstance(values, (list,tuple)):
        values = lib.list_to_object_array(values)
    if getattr(values,'dtype',None) == np.object_:
        values = lib.maybe_convert_objects(values)

    return values

def _possibly_cast_to_timedelta(value, coerce=True):
    """ try to cast to timedelta64, if already a timedeltalike, then make
        sure that we are [ns] (as numpy 1.6.2 is very buggy in this regards,
        don't force the conversion unless coerce is True """

    # deal with numpy not being able to handle certain timedelta operations
    if isinstance(value,np.ndarray) and value.dtype.kind == 'm':
        if value.dtype != 'timedelta64[ns]':
            value = value.astype('timedelta64[ns]')
        return value

    # we don't have a timedelta, but we want to try to convert to one (but don't force it)
    if coerce:
        new_value = tslib.array_to_timedelta64(value.astype(object), coerce=False)
        if new_value.dtype == 'i8':
            value = np.array(new_value,dtype='timedelta64[ns]')

    return value

def _possibly_cast_to_datetime(value, dtype, coerce = False):
    """ try to cast the array/value to a datetimelike dtype, converting float nan to iNaT """

    if dtype is not None:
        if isinstance(dtype, basestring):
            dtype = np.dtype(dtype)

        is_datetime64  = is_datetime64_dtype(dtype)
        is_timedelta64 = is_timedelta64_dtype(dtype)

        if is_datetime64 or is_timedelta64:

            # force the dtype if needed
            if is_datetime64 and dtype != _NS_DTYPE:
                  raise TypeError("cannot convert datetimelike to dtype [%s]" % dtype)
            elif is_timedelta64 and dtype != _TD_DTYPE:
                raise TypeError("cannot convert timedeltalike to dtype [%s]" % dtype)

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
                            value = _possibly_cast_to_timedelta(value)
                    except:
                        pass

    else:

        # only do this if we have an array and the dtype of the array is not setup already
        # we are not an integer/object, so don't bother with this conversion
        if isinstance(value, np.ndarray) and not (issubclass(value.dtype.type, np.integer) or value.dtype == np.object_):
            pass

        else:
            # we might have a array (or single object) that is datetime like, and no dtype is passed
            # don't change the value unless we find a datetime set
            v = value
            if not is_list_like(v):
                v = [ v ]
            if len(v):
                inferred_type = lib.infer_dtype(v)
                if inferred_type in ['datetime','datetime64']:
                    try:
                        value = tslib.array_to_datetime(np.array(v))
                    except:
                        pass
                elif inferred_type in ['timedelta','timedelta64']:
                    value = _possibly_cast_to_timedelta(value)

    return value


def _is_bool_indexer(key):
    if isinstance(key, np.ndarray) and key.dtype == np.object_:
        key = np.asarray(key)

        if not lib.is_bool_array(key):
            if isnull(key).any():
                raise ValueError('cannot index with vector containing '
                                 'NA / NaN values')
            return False
        return True
    elif isinstance(key, np.ndarray) and key.dtype == np.bool_:
        return True
    elif isinstance(key, list):
        try:
            return np.asarray(key).dtype == np.bool_
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


def _mut_exclusive(arg1, arg2):
    if arg1 is not None and arg2 is not None:
        raise Exception('mutually exclusive arguments')
    elif arg1 is not None:
        return arg1
    else:
        return arg2


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
        sep = unicode(sep)
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

    return itertools.izip(seq_it, seq_it_next)


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
    result = 1L
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

    if issubclass(result.dtype.type, basestring):
        result = np.asarray(values, dtype=object)

    if result.ndim == 2:
        if isinstance(values, list):
            return lib.list_to_object_array(values)
        else:
            # Making a 1D array that safely contains tuples is a bit tricky
            # in numpy, leading to the following
            result = np.empty(len(values), dtype=object)
            result[:] = values

    return result


def _index_labels_to_array(labels):
    if isinstance(labels, (basestring, tuple)):
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
    return isinstance(obj, (int, long, np.integer))


def is_float(obj):
    return isinstance(obj, (float, np.floating))


def is_complex(obj):
    return isinstance(obj, (complex, np.complexfloating))


def is_iterator(obj):
    # python 3 generators have __next__ instead of next
    return hasattr(obj, 'next') or hasattr(obj, '__next__')


def is_number(obj):
    return isinstance(obj, (np.number, int, long, float, complex))


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


def is_timedelta64_dtype(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        tipo = arr_or_dtype.type
    elif isinstance(arr_or_dtype, type):
        tipo = np.dtype(arr_or_dtype).type
    else:
        tipo = arr_or_dtype.dtype.type
    return issubclass(tipo, np.timedelta64)

def needs_i8_conversion(arr_or_dtype):
    return is_datetime64_dtype(arr_or_dtype) or is_timedelta64_dtype(arr_or_dtype)

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
    return hasattr(arg, '__iter__') and not isinstance(arg, basestring)

def _is_sequence(x):
    try:
        iter(x)
        len(x) # it has a length
        return not isinstance(x, basestring) and True
    except Exception:
        return False

_ensure_float64 = algos.ensure_float64
_ensure_float32 = algos.ensure_float32
_ensure_int64 = algos.ensure_int64
_ensure_int32 = algos.ensure_int32
_ensure_int16 = algos.ensure_int16
_ensure_int8 = algos.ensure_int8
_ensure_platform_int = algos.ensure_platform_int
_ensure_object = algos.ensure_object


def _astype_nansafe(arr, dtype, copy = True):
    """ return a view if copy is False """
    if not isinstance(dtype, np.dtype):
        dtype = np.dtype(dtype)

    if is_datetime64_dtype(arr):
        if dtype == object:
            return tslib.ints_to_pydatetime(arr.view(np.int64))
        elif dtype == np.int64:
            return arr.view(dtype)
        elif dtype != _NS_DTYPE:
            raise TypeError("cannot astype a datetimelike from [%s] to [%s]" % (arr.dtype,dtype))
        return arr.astype(_NS_DTYPE)
    elif is_timedelta64_dtype(arr):
        if dtype == np.int64:
            return arr.view(dtype)
        elif dtype == object:
            return arr.astype(object)

        # in py3, timedelta64[ns] are int64
        elif (py3compat.PY3 and dtype not in [_INT64_DTYPE,_TD_DTYPE]) or (not py3compat.PY3 and dtype != _TD_DTYPE):
            raise TypeError("cannot astype a timedelta from [%s] to [%s]" % (arr.dtype,dtype))
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


def _clean_fill_method(method):
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
        return self.reader.next().encode("utf-8")


def _get_handle(path, mode, encoding=None, compression=None):
    if compression is not None:
        if encoding is not None:
            raise ValueError('encoding + compression not yet supported')

        if compression == 'gzip':
            import gzip
            return gzip.GzipFile(path, 'rb')
        elif compression == 'bz2':
            import bz2
            return bz2.BZ2File(path, 'rb')
        else:
            raise ValueError('Unrecognized compression type: %s' %
                             compression)

    if py3compat.PY3:  # pragma: no cover
        if encoding:
            f = open(path, mode, encoding=encoding)
        else:
            f = open(path, mode, errors='replace')
    else:
        f = open(path, mode)
    return f

if py3compat.PY3:  # pragma: no cover
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
            row = self.reader.next()
            return [unicode(s, "utf-8") for s in row]

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
    to_concat = [x for x in to_concat if x.shape[axis] > 0]

    # return the empty np array, if nothing to concatenate, #3121
    if not to_concat: return np.array([], dtype=object)

    is_datetime64 = [x.dtype == _NS_DTYPE for x in to_concat]
    if all(is_datetime64):
        # work around NumPy 1.6 bug
        new_values = np.concatenate([x.view(np.int64) for x in to_concat],
                                    axis=axis)
        return new_values.view(_NS_DTYPE)
    elif any(is_datetime64):
        to_concat = [_to_pydatetime(x) for x in to_concat]

    return np.concatenate(to_concat, axis=axis)


def _to_pydatetime(x):
    if x.dtype == _NS_DTYPE:
        shape = x.shape
        x = tslib.ints_to_pydatetime(x.view(np.int64).ravel())
        x = x.reshape(shape)

    return x

def _where_compat(mask, arr1, arr2):
    if arr1.dtype == _NS_DTYPE and arr2.dtype == _NS_DTYPE:
        new_vals = np.where(mask, arr1.view(np.int64), arr2.view(np.int64))
        return new_vals.view(_NS_DTYPE)

    import pandas.tslib as tslib
    if arr1.dtype == _NS_DTYPE:
        arr1 = tslib.ints_to_pydatetime(arr1.view(np.int64))
    if arr2.dtype == _NS_DTYPE:
        arr2 = tslib.ints_to_pydatetime(arr2.view(np.int64))

    return np.where(mask, arr1, arr2)

def sentinal_factory():
    class Sentinal(object):
        pass

    return Sentinal()

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
        front_end = (ip.config.get('KernelApp',{}).get('parent_appname',"") or
                         ip.config.get('IPKernelApp',{}).get('parent_appname',""))
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
        front_end = (ip.config.get('KernelApp',{}).get('parent_appname',"") or
                         ip.config.get('IPKernelApp',{}).get('parent_appname',""))
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
    if isinstance(seq,set):
        fmt = u"set([%s])"
    else:
        fmt = u"[%s]" if hasattr(seq, '__setitem__') else u"(%s)"

    nitems = get_option("max_seq_items") or len(seq)

    s = iter(seq)
    r = []
    for i in range(min(nitems,len(seq))): # handle sets, no slicing
        r.append(pprint_thing(next(s), _nest_lvl + 1, **kwds))
    body = ", ".join(r)

    if nitems < len(seq):
        body+= ", ..."
    elif isinstance(seq,tuple) and len(seq) == 1:
        body += ','

    return fmt % body


def _pprint_dict(seq, _nest_lvl=0,**kwds):
    """
    internal. pprinter for iterables. you should probably use pprint_thing()
    rather then calling this directly.
    """
    fmt = u"{%s}"
    pairs = []

    pfmt = u"%s: %s"

    nitems = get_option("max_seq_items") or len(seq)

    for k, v in seq.items()[:nitems]:
        pairs.append(pfmt % (pprint_thing(k,_nest_lvl+1,**kwds),
                             pprint_thing(v,_nest_lvl+1,**kwds)))

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
    def as_escaped_unicode(thing,escape_chars=escape_chars):
        # Unicode is fine, else we try to decode using utf-8 and 'replace'
        # if that's not it either, we have no way of knowing and the user
        #should deal with it himself.

        try:
            result = unicode(thing)  # we should try this first
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
            escape_chars = escape_chars.keys()
        else:
            escape_chars = escape_chars or tuple()
        for c in escape_chars:
            result = result.replace(c, translate[c])

        return unicode(result)

    if (py3compat.PY3 and hasattr(thing, '__next__')) or \
            hasattr(thing, 'next'):
        return unicode(thing)
    elif (isinstance(thing, dict) and
          _nest_lvl < get_option("display.pprint_nest_depth")):
        result = _pprint_dict(thing, _nest_lvl,quote_strings=True)
    elif _is_sequence(thing) and _nest_lvl < \
            get_option("display.pprint_nest_depth"):
        result = _pprint_seq(thing, _nest_lvl, escape_chars=escape_chars,
                             quote_strings=quote_strings)
    elif isinstance(thing,basestring) and quote_strings:
        if py3compat.PY3:
            fmt = "'%s'"
        else:
            fmt = "u'%s'"
        result = fmt % as_escaped_unicode(thing)
    else:
        result = as_escaped_unicode(thing)

    return unicode(result)  # always unicode


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

    Warning: Loading pickled data received from untrusted sources can be unsafe.
    See: http://docs.python.org/2.7/library/pickle.html

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
    '''
    Pickle (serialize) object to input file path

    Parameters
    ----------
    obj : any object
    path : string
        File path
    '''
    import warnings
    warnings.warn("save is deprecated, use obj.to_pickle", FutureWarning)
    from pandas.io.pickle import to_pickle
    return to_pickle(obj, path)
