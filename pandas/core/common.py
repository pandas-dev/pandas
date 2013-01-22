"""
Misc tools for implementing data structures
"""
try:
    import cPickle as pickle
except ImportError:  # pragma: no cover
    import pickle

import itertools

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


def isnull(obj):
    '''
    Detect missing values (NaN in numeric arrays, None/NaN in object arrays)

    Parameters
    ----------
    arr: ndarray or object value

    Returns
    -------
    boolean ndarray or boolean
    '''
    return _isnull(obj)


def _isnull_new(obj):
    if lib.isscalar(obj):
        return lib.checknull(obj)

    from pandas.core.generic import PandasObject
    if isinstance(obj, np.ndarray):
        return _isnull_ndarraylike(obj)
    elif isinstance(obj, PandasObject):
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

    from pandas.core.generic import PandasObject
    if isinstance(obj, np.ndarray):
        return _isnull_ndarraylike_old(obj)
    elif isinstance(obj, PandasObject):
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
    elif issubclass(values.dtype.type, np.timedelta64):
        # -np.isfinite(values.view('i8'))
        result = np.ones(values.shape, dtype=bool)
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
    '''
    Replacement for numpy.isfinite / -numpy.isnan which is suitable
    for use on object arrays.

    Parameters
    ----------
    arr: ndarray or object value

    Returns
    -------
    boolean ndarray or boolean
    '''
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


def _view_wrapper(f, wrap_dtype, na_override=None):
    def wrapper(arr, indexer, out, fill_value=np.nan):
        if na_override is not None and np.isnan(fill_value):
            fill_value = na_override
        view = arr.view(wrap_dtype)
        outview = out.view(wrap_dtype)
        f(view, indexer, outview, fill_value=fill_value)
    return wrapper


_take1d_dict = {
    'float64': algos.take_1d_float64,
    'int32': algos.take_1d_int32,
    'int64': algos.take_1d_int64,
    'object': algos.take_1d_object,
    'bool': _view_wrapper(algos.take_1d_bool, np.uint8),
    'datetime64[ns]': _view_wrapper(algos.take_1d_int64, np.int64,
                                    na_override=tslib.iNaT),
}

_take2d_axis0_dict = {
    'float64': algos.take_2d_axis0_float64,
    'int32': algos.take_2d_axis0_int32,
    'int64': algos.take_2d_axis0_int64,
    'object': algos.take_2d_axis0_object,
    'bool': _view_wrapper(algos.take_2d_axis0_bool, np.uint8),
    'datetime64[ns]': _view_wrapper(algos.take_2d_axis0_int64, np.int64,
                                    na_override=tslib.iNaT),
}

_take2d_axis1_dict = {
    'float64': algos.take_2d_axis1_float64,
    'int32': algos.take_2d_axis1_int32,
    'int64': algos.take_2d_axis1_int64,
    'object': algos.take_2d_axis1_object,
    'bool': _view_wrapper(algos.take_2d_axis1_bool, np.uint8),
    'datetime64[ns]': _view_wrapper(algos.take_2d_axis1_int64, np.int64,
                                    na_override=tslib.iNaT),
}

_take2d_multi_dict = {
    'float64': algos.take_2d_multi_float64,
    'int32': algos.take_2d_multi_int32,
    'int64': algos.take_2d_multi_int64,
    'object': algos.take_2d_multi_object,
    'bool': _view_wrapper(algos.take_2d_multi_bool, np.uint8),
    'datetime64[ns]': _view_wrapper(algos.take_2d_multi_int64, np.int64,
                                    na_override=tslib.iNaT),
}


def _get_take2d_function(dtype_str, axis=0):
    if axis == 0:
        return _take2d_axis0_dict[dtype_str]
    elif axis == 1:
        return _take2d_axis1_dict[dtype_str]
    elif axis == 'multi':
        return _take2d_multi_dict[dtype_str]
    else:  # pragma: no cover
        raise ValueError('bad axis: %s' % axis)


def take_1d(arr, indexer, out=None, fill_value=np.nan):
    """
    Specialized Cython take which sets NaN values in one pass
    """
    dtype_str = arr.dtype.name

    n = len(indexer)

    indexer = _ensure_int64(indexer)

    out_passed = out is not None
    take_f = _take1d_dict.get(dtype_str)

    if dtype_str in ('int32', 'int64', 'bool'):
        try:
            if out is None:
                out = np.empty(n, dtype=arr.dtype)
            take_f(arr, _ensure_int64(indexer), out=out, fill_value=fill_value)
        except ValueError:
            mask = indexer == -1
            if len(arr) == 0:
                if not out_passed:
                    out = np.empty(n, dtype=arr.dtype)
            else:
                out = ndtake(arr, indexer, out=out)
            if mask.any():
                if out_passed:
                    raise Exception('out with dtype %s does not support NA' %
                                    out.dtype)
                out = _maybe_upcast(out)
                np.putmask(out, mask, fill_value)
    elif dtype_str in ('float64', 'object', 'datetime64[ns]'):
        if out is None:
            out = np.empty(n, dtype=arr.dtype)
        take_f(arr, _ensure_int64(indexer), out=out, fill_value=fill_value)
    else:
        out = ndtake(arr, indexer, out=out)
        mask = indexer == -1
        if mask.any():
            if out_passed:
                raise Exception('out with dtype %s does not support NA' %
                                out.dtype)
            out = _maybe_upcast(out)
            np.putmask(out, mask, fill_value)

    return out


def take_2d_multi(arr, row_idx, col_idx, fill_value=np.nan, out=None):

    dtype_str = arr.dtype.name

    out_shape = len(row_idx), len(col_idx)

    if dtype_str in ('int32', 'int64', 'bool'):
        row_mask = row_idx == -1
        col_mask = col_idx == -1
        needs_masking = row_mask.any() or col_mask.any()

        if needs_masking:
            return take_2d_multi(_maybe_upcast(arr), row_idx, col_idx,
                                 fill_value=fill_value, out=out)
        else:
            if out is None:
                out = np.empty(out_shape, dtype=arr.dtype)
            take_f = _get_take2d_function(dtype_str, axis='multi')
            take_f(arr, _ensure_int64(row_idx),
                   _ensure_int64(col_idx), out=out,
                   fill_value=fill_value)
            return out
    elif dtype_str in ('float64', 'object', 'datetime64[ns]'):
        if out is None:
            out = np.empty(out_shape, dtype=arr.dtype)
        take_f = _get_take2d_function(dtype_str, axis='multi')
        take_f(arr, _ensure_int64(row_idx), _ensure_int64(col_idx), out=out,
               fill_value=fill_value)
        return out
    else:
        if out is not None:
            raise ValueError('Cannot pass out in this case')

        return take_2d(take_2d(arr, row_idx, axis=0, fill_value=fill_value),
                       col_idx, axis=1, fill_value=fill_value)


def take_2d(arr, indexer, out=None, mask=None, needs_masking=None, axis=0,
            fill_value=np.nan):
    """
    Specialized Cython take which sets NaN values in one pass
    """
    dtype_str = arr.dtype.name

    out_shape = list(arr.shape)
    out_shape[axis] = len(indexer)
    out_shape = tuple(out_shape)

    if not isinstance(indexer, np.ndarray):
        indexer = np.array(indexer, dtype=np.int64)

    if dtype_str in ('int32', 'int64', 'bool'):
        if mask is None:
            mask = indexer == -1
            needs_masking = mask.any()

        if needs_masking:
            # upcasting may be required
            result = ndtake(arr, indexer, axis=axis, out=out)
            result = _maybe_mask(result, mask, needs_masking, axis=axis,
                                 out_passed=out is not None,
                                 fill_value=fill_value)
            return result
        else:
            if out is None:
                out = np.empty(out_shape, dtype=arr.dtype)
            take_f = _get_take2d_function(dtype_str, axis=axis)
            take_f(arr, _ensure_int64(indexer), out=out, fill_value=fill_value)
            return out
    elif dtype_str in ('float64', 'object', 'datetime64[ns]'):
        if out is None:
            out = np.empty(out_shape, dtype=arr.dtype)
        take_f = _get_take2d_function(dtype_str, axis=axis)
        take_f(arr, _ensure_int64(indexer), out=out, fill_value=fill_value)
        return out
    else:
        if mask is None:
            mask = indexer == -1
            needs_masking = mask.any()

        # GH #486
        if out is not None and arr.dtype != out.dtype:
            arr = arr.astype(out.dtype)

        result = ndtake(arr, indexer, axis=axis, out=out)
        result = _maybe_mask(result, mask, needs_masking, axis=axis,
                             out_passed=out is not None,
                             fill_value=fill_value)
        return result


def ndtake(arr, indexer, axis=0, out=None):
    return arr.take(_ensure_platform_int(indexer), axis=axis, out=out)


def mask_out_axis(arr, mask, axis, fill_value=np.nan):
    indexer = [slice(None)] * arr.ndim
    indexer[axis] = mask

    arr[tuple(indexer)] = fill_value

_diff_special = {
    'float64': algos.diff_2d_float64,
    'int64': algos.diff_2d_int64,
    'int32': algos.diff_2d_int32
}


def diff(arr, n, axis=0):
    n = int(n)
    dtype = arr.dtype
    if issubclass(dtype.type, np.integer):
        dtype = np.float64
    elif issubclass(dtype.type, np.bool_):
        dtype = np.object_

    out_arr = np.empty(arr.shape, dtype=dtype)

    na_indexer = [slice(None)] * arr.ndim
    na_indexer[axis] = slice(None, n) if n >= 0 else slice(n, None)
    out_arr[tuple(na_indexer)] = np.nan

    if arr.ndim == 2 and arr.dtype.name in _diff_special:
        f = _diff_special[arr.dtype.name]
        f(arr, out_arr, n, axis)
    else:
        res_indexer = [slice(None)] * arr.ndim
        res_indexer[axis] = slice(n, None) if n >= 0 else slice(None, n)
        res_indexer = tuple(res_indexer)

        lag_indexer = [slice(None)] * arr.ndim
        lag_indexer[axis] = slice(None, -n) if n >= 0 else slice(-n, None)
        lag_indexer = tuple(lag_indexer)

        out_arr[res_indexer] = arr[res_indexer] - arr[lag_indexer]

    return out_arr


def take_fast(arr, indexer, mask, needs_masking, axis=0, out=None,
              fill_value=np.nan):
    if arr.ndim == 2:
        return take_2d(arr, indexer, out=out, mask=mask,
                       needs_masking=needs_masking,
                       axis=axis, fill_value=fill_value)
    indexer = _ensure_platform_int(indexer)
    result = ndtake(arr, indexer, axis=axis, out=out)
    result = _maybe_mask(result, mask, needs_masking, axis=axis,
                         out_passed=out is not None, fill_value=fill_value)
    return result


def _maybe_mask(result, mask, needs_masking, axis=0, out_passed=False,
                fill_value=np.nan):
    if needs_masking:
        if out_passed and _need_upcast(result):
            raise Exception('incompatible type for NAs')
        else:
            # a bit spaghettified
            result = _maybe_upcast(result)
            mask_out_axis(result, mask, axis, fill_value)
    return result


def _maybe_upcast(values):
    if issubclass(values.dtype.type, np.integer):
        values = values.astype(float)
    elif issubclass(values.dtype.type, np.bool_):
        values = values.astype(object)

    return values


def _need_upcast(values):
    if issubclass(values.dtype.type, (np.integer, np.bool_)):
        return True
    return False


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
    if is_float_dtype(values):
        _method = algos.pad_inplace_float64
    elif is_datetime64_dtype(values):
        _method = _pad_1d_datetime
    elif values.dtype == np.object_:
        _method = algos.pad_inplace_object
    else:  # pragma: no cover
        raise ValueError('Invalid dtype for padding')

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)
    _method(values, mask, limit=limit)


def backfill_1d(values, limit=None, mask=None):
    if is_float_dtype(values):
        _method = algos.backfill_inplace_float64
    elif is_datetime64_dtype(values):
        _method = _backfill_1d_datetime
    elif values.dtype == np.object_:
        _method = algos.backfill_inplace_object
    else:  # pragma: no cover
        raise ValueError('Invalid dtype for padding')

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)

    _method(values, mask, limit=limit)


def pad_2d(values, limit=None, mask=None):
    if is_float_dtype(values):
        _method = algos.pad_2d_inplace_float64
    elif is_datetime64_dtype(values):
        _method = _pad_2d_datetime
    elif values.dtype == np.object_:
        _method = algos.pad_2d_inplace_object
    else:  # pragma: no cover
        raise ValueError('Invalid dtype for padding')

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)

    if np.all(values.shape):
        _method(values, mask, limit=limit)
    else:
        # for test coverage
        pass


def backfill_2d(values, limit=None, mask=None):
    if is_float_dtype(values):
        _method = algos.backfill_2d_inplace_float64
    elif is_datetime64_dtype(values):
        _method = _backfill_2d_datetime
    elif values.dtype == np.object_:
        _method = algos.backfill_2d_inplace_object
    else:  # pragma: no cover
        raise ValueError('Invalid dtype for padding')

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


def _possibly_cast_to_datetime(value, dtype):
    """ try to cast the array/value to a datetimelike dtype, converting float nan to iNaT """

    if dtype == 'M8[ns]':
        if np.isscalar(value):
            if value == tslib.iNaT or isnull(value):
                value = tslib.iNaT
        else:
            value = np.array(value)

            # have a scalar array-like (e.g. NaT)
            if value.ndim == 0:
                value = tslib.iNaT

            # we have an array of datetime & nulls
            elif np.prod(value.shape):
                try:
                    value = tslib.array_to_datetime(value)
                except:
                    pass

    return value


def _infer_dtype(value):
    if isinstance(value, (float, np.floating)):
        return np.float_
    elif isinstance(value, (bool, np.bool_)):
        return np.bool_
    elif isinstance(value, (int, np.integer)):
        return np.int_
    else:
        return np.object_


def _possibly_cast_item(obj, item, dtype):
    chunk = obj[item]

    if chunk.values.dtype != dtype:
        if dtype in (np.object_, np.bool_):
            obj[item] = chunk.astype(np.object_)
        elif not issubclass(dtype, (np.integer, np.bool_)):  # pragma: no cover
            raise ValueError("Unexpected dtype encountered: %s" % dtype)


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

    Example
    -------
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


def is_integer(obj):
    return isinstance(obj, (int, long, np.integer))


def is_float(obj):
    return isinstance(obj, (float, np.floating))


def is_iterator(obj):
    # python 3 generators have __next__ instead of next
    return hasattr(obj, 'next') or hasattr(obj, '__next__')


def is_number(obj):
    return isinstance(obj, (np.number, int, long, float))


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
    else:
        tipo = arr_or_dtype.dtype.type
    return issubclass(tipo, np.datetime64)


def is_float_dtype(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        tipo = arr_or_dtype.type
    else:
        tipo = arr_or_dtype.dtype.type
    return issubclass(tipo, np.floating)


def is_list_like(arg):
    return hasattr(arg, '__iter__') and not isinstance(arg, basestring)


def _is_sequence(x):
    try:
        iter(x)
        return not isinstance(x, basestring) and True
    except Exception:
        return False

_ensure_float64 = algos.ensure_float64
_ensure_int64 = algos.ensure_int64
_ensure_int32 = algos.ensure_int32
_ensure_platform_int = algos.ensure_platform_int
_ensure_object = algos.ensure_object


def _astype_nansafe(arr, dtype):
    if not isinstance(dtype, np.dtype):
        dtype = np.dtype(dtype)

    if issubclass(arr.dtype.type, np.datetime64):
        if dtype == object:
            return tslib.ints_to_pydatetime(arr.view(np.int64))
    elif (np.issubdtype(arr.dtype, np.floating) and
          np.issubdtype(dtype, np.integer)):

        if np.isnan(arr).any():
            raise ValueError('Cannot convert NA to integer')
    elif arr.dtype == np.object_ and np.issubdtype(dtype.type, np.integer):
        # work around NumPy brokenness, #1987
        return lib.astype_intsafe(arr.ravel(), dtype).reshape(arr.shape)

    return arr.astype(dtype)


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


def save(obj, path):
    """
    Pickle (serialize) object to input file path

    Parameters
    ----------
    obj : any object
    path : string
        File path
    """
    f = open(path, 'wb')
    try:
        pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
    finally:
        f.close()


def load(path):
    """
    Load pickled pandas object (or any other pickled object) from the specified
    file path

    Parameters
    ----------
    path : string
        File path

    Returns
    -------
    unpickled : type of object stored in file
    """
    f = open(path, 'rb')
    try:
        return pickle.load(f)
    finally:
        f.close()


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


_NS_DTYPE = np.dtype('M8[ns]')


def _concat_compat(to_concat, axis=0):
    # filter empty arrays
    to_concat = [x for x in to_concat if x.shape[axis] > 0]

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
        if ip.config['KernelApp']['parent_appname'] == 'ipython-qtconsole':
            return True
    except:
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
    """
    fmt = u"[%s]" if hasattr(seq, '__setitem__') else u"(%s)"
    return fmt % ", ".join(pprint_thing(e, _nest_lvl + 1, **kwds) for e in seq)


def _pprint_dict(seq, _nest_lvl=0):
    """
    internal. pprinter for iterables. you should probably use pprint_thing()
    rather then calling this directly.
    """
    fmt = u"{%s}"
    pairs = []

    pfmt = u"%s: %s"
    for k, v in seq.items():
        pairs.append(pfmt % (repr(k), repr(v)))
    return fmt % ", ".join(pairs)


def pprint_thing(thing, _nest_lvl=0, escape_chars=None, default_escapes=False):
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

    if thing is None:
        result = ''
    elif (py3compat.PY3 and hasattr(thing, '__next__')) or \
            hasattr(thing, 'next'):
        return unicode(thing)
    elif (isinstance(thing, dict) and
          _nest_lvl < get_option("display.pprint_nest_depth")):
        result = _pprint_dict(thing, _nest_lvl)
    elif _is_sequence(thing) and _nest_lvl < \
            get_option("display.pprint_nest_depth"):
        result = _pprint_seq(thing, _nest_lvl, escape_chars=escape_chars)
    else:
        # when used internally in the package, everything
        # should be unicode text. However as an aid to transition
        # we also accept utf8 encoded strings,
        # if that's not it either, we have no way of knowing,
        # and the user should deal with it himself.
        # we resort to utf-8 with replacing errors, rather then throwing
        # an exception.

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
