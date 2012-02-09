"""
Misc tools for implementing data structures
"""
try:
    import cPickle as pickle
except ImportError:  # pragma: no cover
    import pickle

try:
    from io import BytesIO
except ImportError:  # pragma: no cover
    # Python < 2.6
    from cStringIO import StringIO as BytesIO
import itertools

from cStringIO import StringIO

from numpy.lib.format import read_array, write_array
import numpy as np

import pandas._tseries as lib
from pandas.util import py3compat
import codecs
import csv

# XXX: HACK for NumPy 1.5.1 to suppress warnings
try:
    np.seterr(all='ignore')
    np.set_printoptions(suppress=True)
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
            if len(arr) == 0:
                if not out_passed:
                    out = np.empty(n, dtype=arr.dtype)
            else:
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

        # GH #486
        if out is not None and arr.dtype != out.dtype:
            arr = arr.astype(out.dtype)

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
    from pandas.core.index import NULL_INDEX, Index
    if n == 0:
        return NULL_INDEX
    else:
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
    if not isinstance(values, (list, tuple, np.ndarray)):
        values = list(values)

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

def _stringify(col):
    # unicode workaround
    return unicode(col)

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
    return issubclass(tipo, np.integer)

def is_float_dtype(arr_or_dtype):
    if isinstance(arr_or_dtype, np.dtype):
        tipo = arr_or_dtype.type
    else:
        tipo = arr_or_dtype.dtype.type
    return issubclass(tipo, np.floating)


def _ensure_float64(arr):
    if arr.dtype != np.float64:
        arr = arr.astype(np.float64)
    return arr

def _ensure_int64(arr):
    if arr.dtype != np.int64:
        arr = arr.astype(np.int64)
    return arr

def _ensure_platform_int(labels):
    if labels.dtype != np.int_:  # pragma: no cover
        labels = labels.astype(np.int_)
    return labels

def _ensure_int32(arr):
    if arr.dtype != np.int32:
        arr = arr.astype(np.int32)
    return arr

def _ensure_object(arr):
    if arr.dtype != np.object_:
        arr = arr.astype('O')
    return arr


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
        return value.encode(sys.stdin.encoding, 'replace')
    except (AttributeError, TypeError):
        return value.encode('ascii', 'replace')

def csv_encode(value, encoding='UTF-8'):
    if py3compat.PY3 or not isinstance(value, unicode):
        return value

    return value.encode(encoding, 'replace')

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
        f = open(path, mode, encoding=encoding)
    else:
        f = open(path, mode)
    return f

if py3compat.PY3:  # pragma: no cover
    def UnicodeReader(f, dialect=csv.excel, encoding="utf-8", **kwds):
        # ignore encoding
        return csv.reader(f, dialect=csv.excel, **kwds)

    def UnicodeWriter(f, dialect=csv.excel, encoding="utf-8", **kwds):
        return csv.writer(f, dialect=csv.excel, **kwds)
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

        def writerows(self, rows):
            for row in rows:
                self.writerow(row)
