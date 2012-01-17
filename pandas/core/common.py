"""
Misc tools for implementing data structures
"""
import cPickle
try:
    from io import BytesIO
except ImportError:  # pragma: no cover
    # Python < 2.6
    from cStringIO import StringIO as BytesIO
import itertools

from numpy.lib.format import read_array, write_array
import numpy as np

import decimal
import math

import pandas._tseries as lib

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

#-------------------------------------------------------------------------------
# Global formatting options

def set_printoptions(precision=None, column_space=None, max_rows=None,
                     max_columns=None, colheader_justify='right'):
    """
    Alter default behavior of DataFrame.toString

    precision : int
        Floating point output precision (number of significant digits)
    column_space : int
        Default space for DataFrame columns, defaults to 12
    max_rows : int
    max_columns : int
        max_rows and max_columns are used in __repr__() methods to decide if
        to_string() or info() is used to render an object to a string.
        Either one, or both can be set to 0 (experimental). Pandas will figure
        out how big the terminal is and will not display more rows or/and
        columns that can fit on it.
    """
    if precision is not None:
        print_config.precision = precision
    if column_space is not None:
        print_config.column_space = column_space
    if max_rows is not None:
        print_config.max_rows = max_rows
    if max_columns is not None:
        print_config.max_columns = max_columns
    if colheader_justify is not None:
        print_config.colheader_justify = colheader_justify

def reset_printoptions():
    print_config.reset()

class EngFormatter(object):
    """
    Formats float values according to engineering format.

    Based on matplotlib.ticker.EngFormatter
    """

    # The SI engineering prefixes
    ENG_PREFIXES = {
        -24: "y",
        -21: "z",
        -18: "a",
        -15: "f",
        -12: "p",
         -9: "n",
         -6: "u",
         -3: "m",
          0: "",
          3: "k",
          6: "M",
          9: "G",
         12: "T",
         15: "P",
         18: "E",
         21: "Z",
         24: "Y"
      }

    def __init__(self, accuracy=None, use_eng_prefix=False):
        self.accuracy = accuracy
        self.use_eng_prefix = use_eng_prefix

    def __call__(self, num):
        """ Formats a number in engineering notation, appending a letter
        representing the power of 1000 of the original number. Some examples:

        >>> format_eng(0)       # for self.accuracy = 0
        ' 0'

        >>> format_eng(1000000) # for self.accuracy = 1,
                                #     self.use_eng_prefix = True
        ' 1.0M'

        >>> format_eng("-1e-6") # for self.accuracy = 2
                                #     self.use_eng_prefix = False
        '-1.00E-06'

        @param num: the value to represent
        @type num: either a numeric value or a string that can be converted to
                   a numeric value (as per decimal.Decimal constructor)

        @return: engineering formatted string
        """

        dnum = decimal.Decimal(str(num))

        sign = 1

        if dnum < 0:  # pragma: no cover
            sign = -1
            dnum = -dnum

        if dnum != 0:
            pow10 = decimal.Decimal(int(math.floor(dnum.log10()/3)*3))
        else:
            pow10 = decimal.Decimal(0)

        pow10 = pow10.min(max(self.ENG_PREFIXES.keys()))
        pow10 = pow10.max(min(self.ENG_PREFIXES.keys()))
        int_pow10 = int(pow10)

        if self.use_eng_prefix:
            prefix = self.ENG_PREFIXES[int_pow10]
        else:
            if int_pow10 < 0:
                prefix = 'E-%02d' % (-int_pow10)
            else:
                prefix = 'E+%02d' % int_pow10

        mant = sign*dnum/(10**pow10)

        if self.accuracy is None:  # pragma: no cover
            format_str = u"% g%s"
        else:
            format_str = (u"%% .%if%%s" % self.accuracy )

        formatted = format_str % (mant, prefix)

        return formatted #.strip()

def set_eng_float_format(precision=None, accuracy=3, use_eng_prefix=False):
    """
    Alter default behavior on how float is formatted in DataFrame.
    Format float in engineering format. By accuracy, we mean the number of
    decimal digits after the floating point.

    See also EngFormatter.
    """
    if precision is not None: # pragma: no cover
        import warnings
        warnings.warn("'precision' parameter in set_eng_float_format is "
                      "being renamed to 'accuracy'" , FutureWarning)
        accuracy = precision

    print_config.float_format = EngFormatter(accuracy, use_eng_prefix)
    print_config.column_space = max(12, accuracy + 9)

def _stringify(col):
    # unicode workaround
    if isinstance(col, tuple):
        return str(col)
    else:
        return '%s' % col

def _float_format_default(v, width=None):
    """
    Take a float and its formatted representation and if it needs extra space
    to fit the width, reformat it to that width.
    """

    fmt_str   = '%% .%dg' % print_config.precision
    formatted = fmt_str % v

    if width is None:
        return formatted

    extra_spc = width - len(formatted)

    if extra_spc <= 0:
        return formatted

    if 'e' in formatted:
        # we have something like 1e13 or 1.23e13
        base, exp = formatted.split('e')

        if '.' in base:
            # expand fraction by extra space
            whole, frac = base.split('.')
            fmt_str = '%%.%df' % (len(frac) + extra_spc)
            frac = fmt_str % float("0.%s" % frac)
            base = whole + frac[1:]
        else:
            if extra_spc > 1:
                # enough room for fraction
                fmt_str = '%% .%df' % (extra_spc - 1)
                base = fmt_str % float(base)
            else:
                # enough room for decimal point only
                base += '.'

        return base + 'e' + exp
    else:
        # we have something like 123 or 123.456
        if '.' in formatted:
            # expand fraction by extra space
            wholel, fracl = map(len, formatted.split("."))
            fmt_str = '%% .%df' % (fracl + extra_spc)
        else:
            if extra_spc > 1:
                # enough room for fraction
                fmt_str = '%% .%df' % (extra_spc - 1)
            else:
                # enough room for decimal point only
                fmt_str = '% d.'

        return fmt_str % v

def _format(s, dtype, space=None, na_rep=None, float_format=None,
            col_width=None):
    def _just_help(x):
        if space is None:
            return x
        return x[:space].ljust(space)

    def _make_float_format(x):
        if na_rep is not None and isnull(x):
            if np.isnan(x):
                x = ' ' + na_rep
            return _just_help('%s' % x)

        if float_format:
            formatted = float_format(x)
        elif print_config.float_format:
            formatted = print_config.float_format(x)
        else:
            formatted = _float_format_default(x, col_width)

        return _just_help(formatted)

    def _make_int_format(x):
        return _just_help('% d' % x)

    if is_float_dtype(dtype):
        return _make_float_format(s)
    elif is_integer_dtype(dtype):
        return _make_int_format(s)
    else:
        if na_rep is not None and lib.checknull(s):
            return na_rep
        else:
            # object dtype
            return _just_help('%s' % _stringify(s))

class _GlobalPrintConfig(object):
    def __init__(self):
        self.precision = 4
        self.float_format = None
        self.column_space = 12
        self.max_rows = 200
        self.max_columns = 0
        self.colheader_justify = 'right'

    def reset(self):
        self.__init__()

print_config = _GlobalPrintConfig()

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

def _maybe_make_list(obj):
    if obj is not None and not isinstance(obj, (tuple, list)):
        return [obj]
    return obj

def is_integer(obj):
    return isinstance(obj, (int, long, np.integer))

def is_float(obj):
    return isinstance(obj, (float, np.floating))

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
        cPickle.dump(obj, f, protocol=cPickle.HIGHEST_PROTOCOL)
    finally:
        f.close()


def load(path):
    """
    Load pickled pandas object (or any other pickled object) from the specified
    file path

    Parameters
    ----------
p    path : string
        File path

    Returns
    -------
    unpickled : type of object stored in file
    """
    f = open(path, 'rb')
    try:
        return cPickle.load(f)
    finally:
        f.close()


