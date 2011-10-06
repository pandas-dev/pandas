"""
Misc tools for implementing data structures
"""
try:
    from io import BytesIO
except ImportError:   # Python < 2.6
    from cStringIO import StringIO as BytesIO
import itertools

from numpy.lib.format import read_array, write_array
import numpy as np

import pandas._tseries as lib

# XXX: HACK for NumPy 1.5.1 to suppress warnings
try:
    np.seterr(all='ignore')
except Exception: # pragma: no cover
    pass

class PandasError(Exception):
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
    if np.isscalar(obj) or obj is None:
        return lib.checknull(obj)

    from pandas.core.generic import PandasObject
    from pandas import Series
    if isinstance(obj, np.ndarray):
        if obj.dtype.kind in ('O', 'S'):
            # Working around NumPy ticket 1542
            shape = obj.shape
            result = np.empty(shape, dtype=bool)
            vec = lib.isnullobj(obj.ravel())
            result[:] = vec.reshape(shape)

            if isinstance(obj, Series):
                result = Series(result, index=obj.index, copy=False)
        else:
            result = -np.isfinite(obj)
        return result
    elif isinstance(obj, PandasObject):
        # TODO: optimize for DataFrame, etc.
        return obj.apply(isnull)
    else:
        return obj is None

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

def _pickle_array(arr):
    arr = arr.view(np.ndarray)

    buf = BytesIO()
    write_array(buf, arr)

    return buf.getvalue()

def _unpickle_array(bytes):
    arr = read_array(BytesIO(bytes))
    return arr

def _take_1d_bool(arr, indexer, out):
    view = arr.view(np.uint8)
    outview = out.view(np.uint8)
    lib.take_1d_bool(view, indexer, outview)

def _take_2d_axis0_bool(arr, indexer, out):
    view = arr.view(np.uint8)
    outview = out.view(np.uint8)
    lib.take_2d_axis0_bool(view, indexer, outview)

def _take_2d_axis1_bool(arr, indexer, out):
    view = arr.view(np.uint8)
    outview = out.view(np.uint8)
    lib.take_2d_axis1_bool(view, indexer, outview)

_take1d_dict = {
    'float64' : lib.take_1d_float64,
    'int32' : lib.take_1d_int32,
    'int64' : lib.take_1d_int64,
    'object' : lib.take_1d_object,
    'bool' : _take_1d_bool
}

_take2d_axis0_dict = {
    'float64' : lib.take_2d_axis0_float64,
    'int32' : lib.take_2d_axis0_int32,
    'int64' : lib.take_2d_axis0_int64,
    'object' : lib.take_2d_axis0_object,
    'bool' : _take_2d_axis0_bool
}

_take2d_axis1_dict = {
    'float64' : lib.take_2d_axis1_float64,
    'int32' : lib.take_2d_axis1_int32,
    'int64' : lib.take_2d_axis1_int64,
    'object' : lib.take_2d_axis1_object,
    'bool' : _take_2d_axis1_bool
}

def _get_take2d_function(dtype_str, axis=0):
    if axis == 0:
        return _take2d_axis0_dict[dtype_str]
    else:
        return _take2d_axis1_dict[dtype_str]

def take_1d(arr, indexer, out=None):
    """
    Specialized Cython take which sets NaN values in one pass
    """
    dtype_str = arr.dtype.name

    n = len(indexer)

    if not isinstance(indexer, np.ndarray):
        # Cython methods expects 32-bit integers
        indexer = np.array(indexer, dtype=np.int32)

    out_passed = out is not None

    if dtype_str in ('int32', 'int64', 'bool'):
        try:
            if out is None:
                out = np.empty(n, dtype=arr.dtype)
            take_f = _take1d_dict[dtype_str]
            take_f(arr, indexer, out=out)
        except ValueError:
            mask = indexer == -1
            out = arr.take(indexer, out=out)
            if mask.any():
                if out_passed:
                    raise Exception('out with dtype %s does not support NA' %
                                    out.dtype)
                out = _maybe_upcast(out)
                np.putmask(out, mask, np.nan)
    elif dtype_str in ('float64', 'object'):
        if out is None:
            out = np.empty(n, dtype=arr.dtype)
        take_f = _take1d_dict[dtype_str]
        take_f(arr, indexer, out=out)
    else:
        out = arr.take(indexer, out=out)
        mask = indexer == -1
        if mask.any():
            if out_passed:
                raise Exception('out with dtype %s does not support NA' %
                                out.dtype)
            out = _maybe_upcast(out)
            np.putmask(out, mask, np.nan)

    return out

def take_2d(arr, indexer, out=None, mask=None, needs_masking=None, axis=0):
    """
    Specialized Cython take which sets NaN values in one pass
    """
    dtype_str = arr.dtype.name

    out_shape = list(arr.shape)
    out_shape[axis] = len(indexer)
    out_shape = tuple(out_shape)

    if not isinstance(indexer, np.ndarray):
        # Cython methods expects 32-bit integers
        indexer = np.array(indexer, dtype=np.int32)

    if dtype_str in ('int32', 'int64', 'bool'):
        if mask is None:
            mask = indexer == -1
            needs_masking = mask.any()

        if needs_masking:
            # upcasting may be required
            result = arr.take(indexer, axis=axis, out=out)
            result = _maybe_mask(result, mask, needs_masking, axis=axis,
                                 out_passed=out is not None)
            return result
        else:
            if out is None:
                out = np.empty(out_shape, dtype=arr.dtype)
            take_f = _get_take2d_function(dtype_str, axis=axis)
            take_f(arr, indexer, out=out)
            return out
    elif dtype_str in ('float64', 'object'):
        if out is None:
            out = np.empty(out_shape, dtype=arr.dtype)
        take_f = _get_take2d_function(dtype_str, axis=axis)
        take_f(arr, indexer, out=out)
        return out
    else:
        if mask is None:
            mask = indexer == -1
            needs_masking = mask.any()

        result = arr.take(indexer, axis=axis, out=out)
        result = _maybe_mask(result, mask, needs_masking, axis=axis,
                             out_passed=out is not None)
        return result

def null_out_axis(arr, mask, axis):
    indexer = [slice(None)] * arr.ndim
    indexer[axis] = mask

    arr[tuple(indexer)] = np.NaN

def take_fast(arr, indexer, mask, needs_masking, axis=0, out=None):
    if arr.ndim == 2:
        return take_2d(arr, indexer, out=out, mask=mask,
                       needs_masking=needs_masking,
                       axis=axis)

    result = arr.take(indexer, axis=axis, out=out)
    result = _maybe_mask(result, mask, needs_masking, axis=axis,
                         out_passed=out is not None)
    return result

def _maybe_mask(result, mask, needs_masking, axis=0, out_passed=False):
    if needs_masking:
        if out_passed and _need_upcast(result):
            raise Exception('incompatible type for NAs')
        else:
            # a bit spaghettified
            result = _maybe_upcast(result)
            null_out_axis(result, mask, axis)
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

#-------------------------------------------------------------------------------
# Lots of little utilities

def _infer_dtype(value):
    if isinstance(value, (float, np.floating)):
        return float
    elif isinstance(value, (bool, np.bool_)):
        return bool
    elif isinstance(value, (int, np.integer)):
        return int
    else:
        return object

def _is_bool_indexer(key):
    if isinstance(key, np.ndarray) and key.dtype == np.object_:
        mask = isnull(key)
        if mask.any():
            raise ValueError('cannot index with vector containing '
                             'NA / NaN values')
        return set([True, False]).issubset(set(key))
    elif isinstance(key, np.ndarray) and key.dtype == np.bool_:
        return True
    elif isinstance(key, list):
        try:
            return np.asarray(key).dtype == np.bool_
        except TypeError: # pragma: no cover
            return False

    return False

def _default_index(n):
    from pandas.core.index import NULL_INDEX
    if n == 0:
        return NULL_INDEX
    else:
        return np.arange(n)

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

def set_printoptions(precision=None, column_space=None):
    """
    Alter default behavior of DataFrame.toString

    precision : int
        Floating point output precision
    column_space : int
        Default space for DataFrame columns, defaults to 12
    """
    global _float_format, _column_space
    if precision is not None:
        float_format = '%.' + '%d' % precision + 'g'
        _float_format = lambda x: float_format % x
    if column_space is not None:
        _column_space = column_space

_float_format = lambda x: '%.4g' % x
_column_space = 12

def _pfixed(s, space, nanRep=None, float_format=None):
    if isinstance(s, float):
        if nanRep is not None and isnull(s):
            if np.isnan(s):
                s = nanRep
            return (' %s' % s).ljust(space)

        if float_format:
            formatted = float_format(s)
        else:
            is_neg = s < 0
            formatted = _float_format(np.abs(s))

            if is_neg:
                formatted = '-' + formatted
            else:
                formatted = ' ' + formatted

        return formatted.ljust(space)
    else:
        stringified = _stringify(s)
        return (' %s' % stringified)[:space].ljust(space)

def _stringify(col):
    # unicode workaround
    if isinstance(col, tuple):
        return str(col)
    else:
        return '%s' % col

def _format(s, nanRep=None, float_format=None):
    if isinstance(s, float):
        if nanRep is not None and isnull(s):
            if np.isnan(s):
                s = nanRep
            return ' %s' % s

        if float_format:
            formatted = float_format(s)
        else:
            is_neg = s < 0
            formatted = _float_format(np.abs(s))

            if is_neg:
                formatted = '-' + formatted
            else:
                formatted = ' ' + formatted

        return formatted
    else:
        return ' %s' % _stringify(s)

#-------------------------------------------------------------------------------
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
    outLines = []
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
        outLines.append(''.join(lines))
    return '\n'.join(outLines)

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
    _ = seq_it_next.next()

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
    except AttributeError:  # Python 3
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

def _asarray_tuplesafe(values):
    if not isinstance(values, (list, np.ndarray)):
        values = list(values)

    result = np.asarray(values)

    if issubclass(result.dtype.type, basestring):
        result = np.asarray(values, dtype=object)

    if result.ndim == 2:
        result = np.empty(len(values), dtype=object)
        result[:] = values

    return result
