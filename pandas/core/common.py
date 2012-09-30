"""
Misc tools for implementing data structures
"""
try:
    import cPickle as pickle
except ImportError:  # pragma: no cover
    import pickle

import itertools

try:
    next
except NameError:  # pragma: no cover
    # Python < 2.6
    def next(x):
        return x.next()

from numpy.lib.format import read_array, write_array
import numpy as np

import pandas._algos as _algos
import pandas.lib as lib
from pandas.util import py3compat
import codecs
import csv

from pandas.util.py3compat import StringIO, BytesIO

# XXX: HACK for NumPy 1.5.1 to suppress warnings
try:
    np.seterr(all='ignore')
    # np.set_printoptions(suppress=True)
except Exception: # pragma: no cover
    pass

class PandasError(Exception):
    pass

class AmbiguousIndexError(PandasError, KeyError):
    pass

def isnull(obj):
    '''
    Replacement for numpy.isnan / -numpy.isfinite which is suitable
    for use on object arrays.

    Parameters
    ----------
    arr: ndarray or object value

    Returns
    -------
    boolean ndarray or boolean
    '''
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
        result = values.view('i8') == lib.iNaT
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
    'float64' : _algos.take_1d_float64,
    'int32' : _algos.take_1d_int32,
    'int64' : _algos.take_1d_int64,
    'object' : _algos.take_1d_object,
    'bool' : _view_wrapper(_algos.take_1d_bool, np.uint8),
    'datetime64[ns]' : _view_wrapper(_algos.take_1d_int64, np.int64,
                                     na_override=lib.iNaT),
}

_take2d_axis0_dict = {
    'float64' : _algos.take_2d_axis0_float64,
    'int32' : _algos.take_2d_axis0_int32,
    'int64' : _algos.take_2d_axis0_int64,
    'object' : _algos.take_2d_axis0_object,
    'bool' : _view_wrapper(_algos.take_2d_axis0_bool, np.uint8),
    'datetime64[ns]' : _view_wrapper(_algos.take_2d_axis0_int64, np.int64,
                                     na_override=lib.iNaT),
}

_take2d_axis1_dict = {
    'float64' : _algos.take_2d_axis1_float64,
    'int32' : _algos.take_2d_axis1_int32,
    'int64' : _algos.take_2d_axis1_int64,
    'object' : _algos.take_2d_axis1_object,
    'bool' : _view_wrapper(_algos.take_2d_axis1_bool, np.uint8),
    'datetime64[ns]' : _view_wrapper(_algos.take_2d_axis1_int64, np.int64,
                                     na_override=lib.iNaT),
}

_take2d_multi_dict = {
    'float64' : _algos.take_2d_multi_float64,
    'int32' : _algos.take_2d_multi_int32,
    'int64' : _algos.take_2d_multi_int64,
    'object' : _algos.take_2d_multi_object,
    'bool' : _view_wrapper(_algos.take_2d_multi_bool, np.uint8),
    'datetime64[ns]' : _view_wrapper(_algos.take_2d_multi_int64, np.int64,
                                     na_override=lib.iNaT),
}

def _get_take2d_function(dtype_str, axis=0):
    if axis == 0:
        return _take2d_axis0_dict[dtype_str]
    elif axis == 1:
        return _take2d_axis1_dict[dtype_str]
    elif axis == 'multi':
        return _take2d_multi_dict[dtype_str]
    else: # pragma: no cover
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
        col_mask=  col_idx == -1
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

_pad_1d_datetime = _interp_wrapper(_algos.pad_inplace_int64, np.int64)
_pad_2d_datetime = _interp_wrapper(_algos.pad_2d_inplace_int64, np.int64)
_backfill_1d_datetime = _interp_wrapper(_algos.backfill_inplace_int64, np.int64)
_backfill_2d_datetime = _interp_wrapper(_algos.backfill_2d_inplace_int64, np.int64)

def pad_1d(values, limit=None, mask=None):
    if is_float_dtype(values):
        _method = _algos.pad_inplace_float64
    elif is_datetime64_dtype(values):
        _method = _pad_1d_datetime
    elif values.dtype == np.object_:
        _method = _algos.pad_inplace_object
    else: # pragma: no cover
        raise ValueError('Invalid dtype for padding')

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)
    _method(values, mask, limit=limit)

def backfill_1d(values, limit=None, mask=None):
    if is_float_dtype(values):
        _method = _algos.backfill_inplace_float64
    elif is_datetime64_dtype(values):
        _method = _backfill_1d_datetime
    elif values.dtype == np.object_:
        _method = _algos.backfill_inplace_object
    else: # pragma: no cover
        raise ValueError('Invalid dtype for padding')

    if mask is None:
        mask = isnull(values)
    mask = mask.view(np.uint8)

    _method(values, mask, limit=limit)

def pad_2d(values, limit=None, mask=None):
    if is_float_dtype(values):
        _method = _algos.pad_2d_inplace_float64
    elif is_datetime64_dtype(values):
        _method = _pad_2d_datetime
    elif values.dtype == np.object_:
        _method = _algos.pad_2d_inplace_object
    else: # pragma: no cover
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
        _method = _algos.backfill_2d_inplace_float64
    elif is_datetime64_dtype(values):
        _method = _backfill_2d_datetime
    elif values.dtype == np.object_:
        _method = _algos.backfill_2d_inplace_object
    else: # pragma: no cover
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
        elif not issubclass(dtype, (np.integer, np.bool_)): # pragma: no cover
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
        except TypeError: # pragma: no cover
            return False

    return False

def _default_index(n):
    from pandas.core.index import Index
    return Index(np.arange(n))

def ensure_float(arr):
    if issubclass(arr.dtype.type, np.integer):
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
    return ''.join(Random().sample(string.ascii_letters+string.digits, n))

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
    _ = next(seq_it_next)

    return itertools.izip(seq_it, seq_it_next)

def indent(string, spaces=4):
    dent = ' ' * spaces
    return '\n'.join([dent + x for x in string.split('\n')])

def banner(message):
    """
    Return 80-char width message declaration with = bars on top and bottom.
    """
    bar = '=' * 80
    return '%s\n%s\n%s' % (bar, message, bar)

class groupby(dict):
    """
    A simple groupby different from the one in itertools.

    Does not require the sequence elements to be sorted by keys,
    however it is slower.
    """
    def __init__(self, seq, key=lambda x:x):
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
        except TypeError: # non-iterable
            labels = [labels]

    labels = _asarray_tuplesafe(labels)

    return labels

def _stringify(col, encoding='UTF8'):
    # unicode workaround
    try:
        return unicode(col)
    except UnicodeError:
        try:
            if isinstance(col, str):
                return col.decode(encoding)
        except UnicodeError:
            pass
        return console_encode(col)

def _stringify_seq(values):
    if any(isinstance(x, unicode) for x in values):
        return [_stringify(x) for x in values]
    return [str(x) for x in values]

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

def is_integer_dtype(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        tipo = arr_or_dtype.type
    else:
        tipo = arr_or_dtype.dtype.type
    return (issubclass(tipo, np.integer) and not
            issubclass(tipo, np.datetime64))

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

_ensure_float64 = _algos.ensure_float64
_ensure_int64 = _algos.ensure_int64
_ensure_int32 = _algos.ensure_int32
_ensure_platform_int = _algos.ensure_platform_int
_ensure_object = _algos.ensure_object


def _astype_nansafe(arr, dtype):
    if not isinstance(dtype, np.dtype):
        dtype = np.dtype(dtype)

    if issubclass(arr.dtype.type, np.datetime64):
        if dtype == object:
            return lib.ints_to_pydatetime(arr.view(np.int64))
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
        msg = ('Invalid fill method. Expecting pad (ffill) or backfill (bfill).'
               ' Got %s' % method)
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

def console_encode(value):
    if py3compat.PY3 or not isinstance(value, unicode):
        return value

    try:
        import sys
        return value.encode(sys.stdin.encoding or 'utf-8', 'replace')
    except (AttributeError, TypeError):
        return value.encode('ascii', 'replace')

class UTF8Recoder:
    """
    Iterator that reads an encoded stream and reencodes the input to UTF-8
    """
    def __init__(self, f, encoding):
        self.reader = codecs.getreader(encoding)(f)

    def __iter__(self):
        return self

    def next(self):
        return self.reader.next().encode("utf-8")

def _get_handle(path, mode, encoding=None):
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

        def writerow(self, row):
            row = [x if isinstance(x, basestring) else str(x) for x in row]
            self.writer.writerow([s.encode("utf-8") for s in row])
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

    if all(x.dtype == _NS_DTYPE for x in to_concat):
        # work around NumPy 1.6 bug
        new_values = np.concatenate([x.view(np.int64) for x in to_concat],
                                    axis=axis)
        return new_values.view(_NS_DTYPE)
    else:
        return np.concatenate(to_concat, axis=axis)
