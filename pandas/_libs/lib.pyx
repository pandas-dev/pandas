# cython: profile=False
cimport numpy as np
cimport cython
import numpy as np
import sys

cdef bint PY3 = (sys.version_info[0] >= 3)

from numpy cimport *

# initialize numpy
np.import_array()
np.import_ufunc()

from libc.stdlib cimport malloc, free

from cpython cimport (Py_INCREF, PyTuple_SET_ITEM,
                      PyList_Check, PyFloat_Check,
                      PyString_Check,
                      PyBytes_Check,
                      PyUnicode_Check,
                      PyTuple_New,
                      PyObject_RichCompareBool,
                      PyBytes_GET_SIZE,
                      PyUnicode_GET_SIZE,
                      PyObject)

try:
    from cpython cimport PyString_GET_SIZE
except ImportError:
    from cpython cimport PyUnicode_GET_SIZE as PyString_GET_SIZE

cdef extern from "Python.h":
    Py_ssize_t PY_SSIZE_T_MAX

cdef extern from "compat_helper.h":

    cdef int slice_get_indices(
        PyObject* s, Py_ssize_t length,
        Py_ssize_t *start, Py_ssize_t *stop, Py_ssize_t *step,
        Py_ssize_t *slicelength) except -1

cimport cpython

isnan = np.isnan
cdef double NaN = <double> np.NaN
cdef double nan = NaN

from cpython.datetime cimport (PyDateTime_Check, PyDate_Check,
                               PyTime_Check, PyDelta_Check,
                               PyDateTime_IMPORT)
PyDateTime_IMPORT

from tslibs.np_datetime cimport get_timedelta64_value, get_datetime64_value

from tslib import NaT, Timestamp, Timedelta, array_to_datetime
from interval import Interval
from missing cimport checknull

cdef int64_t NPY_NAT = util.get_nat()

cimport util
from util cimport is_array, _checknull

from libc.math cimport sqrt, fabs


def values_from_object(object o):
    """ return my values or the object if we are say an ndarray """
    cdef f

    f = getattr(o, 'get_values', None)
    if f is not None:
        o = f()

    return o


@cython.wraparound(False)
@cython.boundscheck(False)
def memory_usage_of_objects(ndarray[object, ndim=1] arr):
    """ return the memory usage of an object array in bytes,
    does not include the actual bytes of the pointers """
    cdef Py_ssize_t i, n
    cdef int64_t s = 0

    n = len(arr)
    for i from 0 <= i < n:
        s += arr[i].__sizeof__()
    return s


# ----------------------------------------------------------------------


cpdef bint isscalar(object val):
    """
    Return True if given value is scalar.

    This includes:
    - numpy array scalar (e.g. np.int64)
    - Python builtin numerics
    - Python builtin byte arrays and strings
    - None
    - instances of datetime.datetime
    - instances of datetime.timedelta
    - Period
    - instances of decimal.Decimal
    - Interval

    """

    return (np.PyArray_IsAnyScalar(val)
            # As of numpy-1.9, PyArray_IsAnyScalar misses bytearrays on Py3.
            or PyBytes_Check(val)
            # We differ from numpy (as of 1.10), which claims that None is
            # not scalar in np.isscalar().
            or val is None
            or PyDate_Check(val)
            or PyDelta_Check(val)
            or PyTime_Check(val)
            or util.is_period_object(val)
            or is_decimal(val)
            or is_interval(val))


def item_from_zerodim(object val):
    """
    If the value is a zerodim array, return the item it contains.

    Examples
    --------
    >>> item_from_zerodim(1)
    1
    >>> item_from_zerodim('foobar')
    'foobar'
    >>> item_from_zerodim(np.array(1))
    1
    >>> item_from_zerodim(np.array([1]))
    array([1])

    """
    return util.unbox_if_zerodim(val)


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef ndarray[object] list_to_object_array(list obj):
    """
    Convert list to object ndarray. Seriously can\'t believe
    I had to write this function.
    """
    cdef:
        Py_ssize_t i, n = len(obj)
        ndarray[object] arr = np.empty(n, dtype=object)

    for i in range(n):
        arr[i] = obj[i]

    return arr


@cython.wraparound(False)
@cython.boundscheck(False)
def fast_unique(ndarray[object] values):
    cdef:
        Py_ssize_t i, n = len(values)
        list uniques = []
        dict table = {}
        object val, stub = 0

    for i from 0 <= i < n:
        val = values[i]
        if val not in table:
            table[val] = stub
            uniques.append(val)
    try:
        uniques.sort()
    except Exception:
        pass

    return uniques


@cython.wraparound(False)
@cython.boundscheck(False)
def fast_unique_multiple(list arrays):
    cdef:
        ndarray[object] buf
        Py_ssize_t k = len(arrays)
        Py_ssize_t i, j, n
        list uniques = []
        dict table = {}
        object val, stub = 0

    for i from 0 <= i < k:
        buf = arrays[i]
        n = len(buf)
        for j from 0 <= j < n:
            val = buf[j]
            if val not in table:
                table[val] = stub
                uniques.append(val)
    try:
        uniques.sort()
    except Exception:
        pass

    return uniques


@cython.wraparound(False)
@cython.boundscheck(False)
def fast_unique_multiple_list(list lists):
    cdef:
        list buf
        Py_ssize_t k = len(lists)
        Py_ssize_t i, j, n
        list uniques = []
        dict table = {}
        object val, stub = 0

    for i from 0 <= i < k:
        buf = lists[i]
        n = len(buf)
        for j from 0 <= j < n:
            val = buf[j]
            if val not in table:
                table[val] = stub
                uniques.append(val)
    try:
        uniques.sort()
    except Exception:
        pass

    return uniques


@cython.wraparound(False)
@cython.boundscheck(False)
def fast_unique_multiple_list_gen(object gen, bint sort=True):
    """
    Generate a list of unique values from a generator of lists.

    Parameters
    ----------
    gen : generator object
        A generator of lists from which the unique list is created
    sort : boolean
        Whether or not to sort the resulting unique list

    Returns
    -------
    unique_list : list of unique values
    """
    cdef:
        list buf
        Py_ssize_t j, n
        list uniques = []
        dict table = {}
        object val, stub = 0

    for buf in gen:
        n = len(buf)
        for j from 0 <= j < n:
            val = buf[j]
            if val not in table:
                table[val] = stub
                uniques.append(val)
    if sort:
        try:
            uniques.sort()
        except Exception:
            pass

    return uniques


@cython.wraparound(False)
@cython.boundscheck(False)
def dicts_to_array(list dicts, list columns):
    cdef:
        Py_ssize_t i, j, k, n
        ndarray[object, ndim=2] result
        dict row
        object col, onan = np.nan

    k = len(columns)
    n = len(dicts)

    result = np.empty((n, k), dtype='O')

    for i in range(n):
        row = dicts[i]
        for j in range(k):
            col = columns[j]
            if col in row:
                result[i, j] = row[col]
            else:
                result[i, j] = onan

    return result


def fast_zip(list ndarrays):
    """
    For zipping multiple ndarrays into an ndarray of tuples
    """
    cdef:
        Py_ssize_t i, j, k, n
        ndarray[object] result
        flatiter it
        object val, tup

    k = len(ndarrays)
    n = len(ndarrays[0])

    result = np.empty(n, dtype=object)

    # initialize tuples on first pass
    arr = ndarrays[0]
    it = <flatiter> PyArray_IterNew(arr)
    for i in range(n):
        val = PyArray_GETITEM(arr, PyArray_ITER_DATA(it))
        tup = PyTuple_New(k)

        PyTuple_SET_ITEM(tup, 0, val)
        Py_INCREF(val)
        result[i] = tup
        PyArray_ITER_NEXT(it)

    for j in range(1, k):
        arr = ndarrays[j]
        it = <flatiter> PyArray_IterNew(arr)
        if len(arr) != n:
            raise ValueError('all arrays must be same length')

        for i in range(n):
            val = PyArray_GETITEM(arr, PyArray_ITER_DATA(it))
            PyTuple_SET_ITEM(result[i], j, val)
            Py_INCREF(val)
            PyArray_ITER_NEXT(it)

    return result


def get_reverse_indexer(ndarray[int64_t] indexer, Py_ssize_t length):
    """
    Reverse indexing operation.

    Given `indexer`, make `indexer_inv` of it, such that::

        indexer_inv[indexer[x]] = x

    .. note:: If indexer is not unique, only first occurrence is accounted.

    """

    cdef:
        Py_ssize_t i, n = len(indexer)
        ndarray[int64_t] rev_indexer
        int64_t idx

    rev_indexer = np.empty(length, dtype=np.int64)
    rev_indexer.fill(-1)
    for i in range(n):
        idx = indexer[i]
        if idx != -1:
            rev_indexer[idx] = i

    return rev_indexer


def has_infs_f4(ndarray[float32_t] arr):
    cdef:
        Py_ssize_t i, n = len(arr)
        float32_t inf, neginf, val

    inf = np.inf
    neginf = -inf

    for i in range(n):
        val = arr[i]
        if val == inf or val == neginf:
            return True
    return False


def has_infs_f8(ndarray[float64_t] arr):
    cdef:
        Py_ssize_t i, n = len(arr)
        float64_t inf, neginf, val

    inf = np.inf
    neginf = -inf

    for i in range(n):
        val = arr[i]
        if val == inf or val == neginf:
            return True
    return False


def convert_timestamps(ndarray values):
    cdef:
        object val, f, result
        dict cache = {}
        Py_ssize_t i, n = len(values)
        ndarray[object] out

    # for HDFStore, a bit temporary but...

    from datetime import datetime
    f = datetime.fromtimestamp

    out = np.empty(n, dtype='O')

    for i in range(n):
        val = util.get_value_1d(values, i)
        if val in cache:
            out[i] = cache[val]
        else:
            cache[val] = out[i] = f(val)

    return out


def maybe_indices_to_slice(ndarray[int64_t] indices, int max_len):
    cdef:
        Py_ssize_t i, n = len(indices)
        int k, vstart, vlast, v

    if n == 0:
        return slice(0, 0)

    vstart = indices[0]
    if vstart < 0 or max_len <= vstart:
        return indices

    if n == 1:
        return slice(vstart, vstart + 1)

    vlast = indices[n - 1]
    if vlast < 0 or max_len <= vlast:
        return indices

    k = indices[1] - indices[0]
    if k == 0:
        return indices
    else:
        for i in range(2, n):
            v = indices[i]
            if v - indices[i - 1] != k:
                return indices

        if k > 0:
            return slice(vstart, vlast + 1, k)
        else:
            if vlast == 0:
                return slice(vstart, None, k)
            else:
                return slice(vstart, vlast - 1, k)


def maybe_booleans_to_slice(ndarray[uint8_t] mask):
    cdef:
        Py_ssize_t i, n = len(mask)
        Py_ssize_t start, end
        bint started = 0, finished = 0

    for i in range(n):
        if mask[i]:
            if finished:
                return mask.view(np.bool_)
            if not started:
                started = 1
                start = i
        else:
            if finished:
                continue

            if started:
                end = i
                finished = 1

    if not started:
        return slice(0, 0)
    if not finished:
        return slice(start, None)
    else:
        return slice(start, end)


@cython.wraparound(False)
@cython.boundscheck(False)
def scalar_compare(ndarray[object] values, object val, object op):
    import operator
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[uint8_t, cast=True] result
        bint isnull_val
        int flag
        object x

    if op is operator.lt:
        flag = cpython.Py_LT
    elif op is operator.le:
        flag = cpython.Py_LE
    elif op is operator.gt:
        flag = cpython.Py_GT
    elif op is operator.ge:
        flag = cpython.Py_GE
    elif op is operator.eq:
        flag = cpython.Py_EQ
    elif op is operator.ne:
        flag = cpython.Py_NE
    else:
        raise ValueError('Unrecognized operator')

    result = np.empty(n, dtype=bool).view(np.uint8)
    isnull_val = checknull(val)

    if flag == cpython.Py_NE:
        for i in range(n):
            x = values[i]
            if checknull(x):
                result[i] = True
            elif isnull_val:
                result[i] = True
            else:
                try:
                    result[i] = cpython.PyObject_RichCompareBool(x, val, flag)
                except (TypeError):
                    result[i] = True
    elif flag == cpython.Py_EQ:
        for i in range(n):
            x = values[i]
            if checknull(x):
                result[i] = False
            elif isnull_val:
                result[i] = False
            else:
                try:
                    result[i] = cpython.PyObject_RichCompareBool(x, val, flag)
                except (TypeError):
                    result[i] = False

    else:
        for i in range(n):
            x = values[i]
            if checknull(x):
                result[i] = False
            elif isnull_val:
                result[i] = False
            else:
                result[i] = cpython.PyObject_RichCompareBool(x, val, flag)

    return result.view(bool)


@cython.wraparound(False)
@cython.boundscheck(False)
cpdef bint array_equivalent_object(object[:] left, object[:] right):
    """ perform an element by element comparion on 1-d object arrays
        taking into account nan positions """
    cdef:
        Py_ssize_t i, n = left.shape[0]
        object x, y

    for i in range(n):
        x = left[i]
        y = right[i]

        # we are either not equal or both nan
        # I think None == None will be true here
        if not (PyObject_RichCompareBool(x, y, cpython.Py_EQ) or
                _checknull(x) and _checknull(y)):
            return False
    return True


@cython.wraparound(False)
@cython.boundscheck(False)
def vec_compare(ndarray[object] left, ndarray[object] right, object op):
    import operator
    cdef:
        Py_ssize_t i, n = len(left)
        ndarray[uint8_t, cast=True] result
        int flag

    if n != len(right):
        raise ValueError('Arrays were different lengths: %d vs %d'
                         % (n, len(right)))

    if op is operator.lt:
        flag = cpython.Py_LT
    elif op is operator.le:
        flag = cpython.Py_LE
    elif op is operator.gt:
        flag = cpython.Py_GT
    elif op is operator.ge:
        flag = cpython.Py_GE
    elif op is operator.eq:
        flag = cpython.Py_EQ
    elif op is operator.ne:
        flag = cpython.Py_NE
    else:
        raise ValueError('Unrecognized operator')

    result = np.empty(n, dtype=bool).view(np.uint8)

    if flag == cpython.Py_NE:
        for i in range(n):
            x = left[i]
            y = right[i]

            if checknull(x) or checknull(y):
                result[i] = True
            else:
                result[i] = cpython.PyObject_RichCompareBool(x, y, flag)
    else:
        for i in range(n):
            x = left[i]
            y = right[i]

            if checknull(x) or checknull(y):
                result[i] = False
            else:
                result[i] = cpython.PyObject_RichCompareBool(x, y, flag)

    return result.view(bool)


@cython.wraparound(False)
@cython.boundscheck(False)
def scalar_binop(ndarray[object] values, object val, object op):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[object] result
        object x

    result = np.empty(n, dtype=object)
    if _checknull(val):
        result.fill(val)
        return result

    for i in range(n):
        x = values[i]
        if _checknull(x):
            result[i] = x
        else:
            result[i] = op(x, val)

    return maybe_convert_bool(result)


@cython.wraparound(False)
@cython.boundscheck(False)
def vec_binop(ndarray[object] left, ndarray[object] right, object op):
    cdef:
        Py_ssize_t i, n = len(left)
        ndarray[object] result

    if n != len(right):
        raise ValueError('Arrays were different lengths: %d vs %d'
                         % (n, len(right)))

    result = np.empty(n, dtype=object)

    for i in range(n):
        x = left[i]
        y = right[i]
        try:
            result[i] = op(x, y)
        except TypeError:
            if _checknull(x):
                result[i] = x
            elif _checknull(y):
                result[i] = y
            else:
                raise

    return maybe_convert_bool(result)


def astype_intsafe(ndarray[object] arr, new_dtype):
    cdef:
        Py_ssize_t i, n = len(arr)
        object v
        bint is_datelike
        ndarray result

    # on 32-bit, 1.6.2 numpy M8[ns] is a subdtype of integer, which is weird
    is_datelike = new_dtype in ['M8[ns]', 'm8[ns]']

    result = np.empty(n, dtype=new_dtype)
    for i in range(n):
        v = arr[i]
        if is_datelike and checknull(v):
            result[i] = NPY_NAT
        else:
            # we can use the unsafe version because we know `result` is mutable
            # since it was created from `np.empty`
            util.set_value_at_unsafe(result, i, v)

    return result


cpdef ndarray[object] astype_unicode(ndarray arr):
    cdef:
        Py_ssize_t i, n = arr.size
        ndarray[object] result = np.empty(n, dtype=object)

    for i in range(n):
        # we can use the unsafe version because we know `result` is mutable
        # since it was created from `np.empty`
        util.set_value_at_unsafe(result, i, unicode(arr[i]))

    return result


cpdef ndarray[object] astype_str(ndarray arr):
    cdef:
        Py_ssize_t i, n = arr.size
        ndarray[object] result = np.empty(n, dtype=object)

    for i in range(n):
        # we can use the unsafe version because we know `result` is mutable
        # since it was created from `np.empty`
        util.set_value_at_unsafe(result, i, str(arr[i]))

    return result


def clean_index_list(list obj):
    """
    Utility used in pandas.core.index._ensure_index
    """
    cdef:
        Py_ssize_t i, n = len(obj)
        object v
        bint all_arrays = 1

    for i in range(n):
        v = obj[i]
        if not (PyList_Check(v) or np.PyArray_Check(v) or hasattr(v, '_data')):
            all_arrays = 0
            break

    if all_arrays:
        return obj, all_arrays

    # don't force numpy coerce with nan's
    inferred = infer_dtype(obj)
    if inferred in ['string', 'bytes', 'unicode',
                    'mixed', 'mixed-integer']:
        return np.asarray(obj, dtype=object), 0
    elif inferred in ['integer']:

        # TODO: we infer an integer but it *could* be a unint64
        try:
            return np.asarray(obj, dtype='int64'), 0
        except OverflowError:
            return np.asarray(obj, dtype='object'), 0

    return np.asarray(obj), 0


ctypedef fused pandas_string:
    str
    unicode
    bytes


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef Py_ssize_t max_len_string_array(pandas_string[:] arr):
    """ return the maximum size of elements in a 1-dim string array """
    cdef:
        Py_ssize_t i, m = 0, l = 0, length = arr.shape[0]
        pandas_string v

    for i in range(length):
        v = arr[i]
        if PyString_Check(v):
            l = PyString_GET_SIZE(v)
        elif PyBytes_Check(v):
            l = PyBytes_GET_SIZE(v)
        elif PyUnicode_Check(v):
            l = PyUnicode_GET_SIZE(v)

        if l > m:
            m = l

    return m


@cython.boundscheck(False)
@cython.wraparound(False)
def string_array_replace_from_nan_rep(
        ndarray[object, ndim=1] arr, object nan_rep,
        object replace=None):
    """
    Replace the values in the array with 'replacement' if
    they are 'nan_rep'. Return the same array.
    """

    cdef int length = arr.shape[0], i = 0
    if replace is None:
        replace = np.nan

    for i from 0 <= i < length:
        if arr[i] == nan_rep:
            arr[i] = replace

    return arr


@cython.boundscheck(False)
@cython.wraparound(False)
def convert_json_to_lines(object arr):
    """
    replace comma separated json with line feeds, paying special attention
    to quotes & brackets
    """
    cdef:
        Py_ssize_t i = 0, num_open_brackets_seen = 0, length
        bint in_quotes = 0, is_escaping = 0
        ndarray[uint8_t] narr
        unsigned char v, comma, left_bracket, right_brack, newline

    newline = ord('\n')
    comma = ord(',')
    left_bracket = ord('{')
    right_bracket = ord('}')
    quote = ord('"')
    backslash = ord('\\')

    narr = np.frombuffer(arr.encode('utf-8'), dtype='u1').copy()
    length = narr.shape[0]
    for i in range(length):
        v = narr[i]
        if v == quote and i > 0 and not is_escaping:
            in_quotes = ~in_quotes
        if v == backslash or is_escaping:
            is_escaping = ~is_escaping
        if v == comma:  # commas that should be \n
            if num_open_brackets_seen == 0 and not in_quotes:
                narr[i] = newline
        elif v == left_bracket:
            if not in_quotes:
                num_open_brackets_seen += 1
        elif v == right_bracket:
            if not in_quotes:
                num_open_brackets_seen -= 1

    return narr.tostring().decode('utf-8')


@cython.boundscheck(False)
@cython.wraparound(False)
def write_csv_rows(list data, ndarray data_index,
                   int nlevels, ndarray cols, object writer):

    cdef int N, j, i, ncols
    cdef list rows
    cdef object val

    # In crude testing, N>100 yields little marginal improvement
    N=100

    # pre-allocate rows
    ncols = len(cols)
    rows = [[None] * (nlevels + ncols) for x in range(N)]

    j = -1
    if nlevels == 1:
        for j in range(len(data_index)):
            row = rows[j % N]
            row[0] = data_index[j]
            for i in range(ncols):
                row[1 + i] = data[i][j]

            if j >= N - 1 and j % N == N - 1:
                writer.writerows(rows)
    elif nlevels > 1:
        for j in range(len(data_index)):
            row = rows[j % N]
            row[:nlevels] = list(data_index[j])
            for i in range(ncols):
                row[nlevels + i] = data[i][j]

            if j >= N - 1 and j % N == N - 1:
                writer.writerows(rows)
    else:
        for j in range(len(data_index)):
            row = rows[j % N]
            for i in range(ncols):
                row[i] = data[i][j]

            if j >= N - 1 and j % N == N - 1:
                writer.writerows(rows)

    if j >= 0 and (j < N - 1 or (j % N) != N - 1):
        writer.writerows(rows[:((j + 1) % N)])


# ------------------------------------------------------------------------------
# Groupby-related functions

@cython.wraparound(False)
@cython.boundscheck(False)
def is_lexsorted(list list_of_arrays):
    cdef:
        int i
        Py_ssize_t n, nlevels
        int64_t k, cur, pre
        ndarray arr

    nlevels = len(list_of_arrays)
    n = len(list_of_arrays[0])

    cdef int64_t **vecs = <int64_t**> malloc(nlevels * sizeof(int64_t*))
    for i from 0 <= i < nlevels:
        arr = list_of_arrays[i]
        vecs[i] = <int64_t *> arr.data

    # Assume uniqueness??
    for i from 1 <= i < n:
        for k from 0 <= k < nlevels:
            cur = vecs[k][i]
            pre = vecs[k][i - 1]
            if cur == pre:
                continue
            elif cur > pre:
                break
            else:
                return False
    free(vecs)
    return True


# TODO: could do even better if we know something about the data. eg, index has
# 1-min data, binner has 5-min data, then  bins are just strides in index. This
# is a general, O(max(len(values), len(binner))) method.
@cython.boundscheck(False)
@cython.wraparound(False)
def generate_bins_dt64(ndarray[int64_t] values, ndarray[int64_t] binner,
                       object closed='left', bint hasnans=0):
    """
    Int64 (datetime64) version of generic python version in groupby.py
    """
    cdef:
        Py_ssize_t lenidx, lenbin, i, j, bc, vc
        ndarray[int64_t] bins
        int64_t l_bin, r_bin, nat_count
        bint right_closed = closed == 'right'

    nat_count = 0
    if hasnans:
        mask = values == iNaT
        nat_count = np.sum(mask)
        values = values[~mask]

    lenidx = len(values)
    lenbin = len(binner)

    if lenidx <= 0 or lenbin <= 0:
        raise ValueError("Invalid length for values or for binner")

    # check binner fits data
    if values[0] < binner[0]:
        raise ValueError("Values falls before first bin")

    if values[lenidx - 1] > binner[lenbin - 1]:
        raise ValueError("Values falls after last bin")

    bins = np.empty(lenbin - 1, dtype=np.int64)

    j = 0  # index into values
    bc = 0  # bin count

    # linear scan
    if right_closed:
        for i in range(0, lenbin - 1):
            r_bin = binner[i + 1]
            # count values in current bin, advance to next bin
            while j < lenidx and values[j] <= r_bin:
                j += 1
            bins[bc] = j
            bc += 1
    else:
        for i in range(0, lenbin - 1):
            r_bin = binner[i + 1]
            # count values in current bin, advance to next bin
            while j < lenidx and values[j] < r_bin:
                j += 1
            bins[bc] = j
            bc += 1

    if nat_count > 0:
        # shift bins by the number of NaT
        bins = bins + nat_count
        bins = np.insert(bins, 0, nat_count)

    return bins


@cython.boundscheck(False)
@cython.wraparound(False)
def row_bool_subset(ndarray[float64_t, ndim=2] values,
                    ndarray[uint8_t, cast=True] mask):
    cdef:
        Py_ssize_t i, j, n, k, pos = 0
        ndarray[float64_t, ndim=2] out

    n, k = (<object> values).shape
    assert(n == len(mask))

    out = np.empty((mask.sum(), k), dtype=np.float64)

    for i in range(n):
        if mask[i]:
            for j in range(k):
                out[pos, j] = values[i, j]
            pos += 1

    return out


@cython.boundscheck(False)
@cython.wraparound(False)
def row_bool_subset_object(ndarray[object, ndim=2] values,
                           ndarray[uint8_t, cast=True] mask):
    cdef:
        Py_ssize_t i, j, n, k, pos = 0
        ndarray[object, ndim=2] out

    n, k = (<object> values).shape
    assert(n == len(mask))

    out = np.empty((mask.sum(), k), dtype=object)

    for i in range(n):
        if mask[i]:
            for j in range(k):
                out[pos, j] = values[i, j]
            pos += 1

    return out


@cython.boundscheck(False)
@cython.wraparound(False)
def get_level_sorter(ndarray[int64_t, ndim=1] label,
                     ndarray[int64_t, ndim=1] starts):
    """
    argsort for a single level of a multi-index, keeping the order of higher
    levels unchanged. `starts` points to starts of same-key indices w.r.t
    to leading levels; equivalent to:
        np.hstack([label[starts[i]:starts[i+1]].argsort(kind='mergesort')
            + starts[i] for i in range(len(starts) - 1)])
    """
    cdef:
        int64_t l, r
        Py_ssize_t i
        ndarray[int64_t, ndim=1] out = np.empty(len(label), dtype=np.int64)

    for i in range(len(starts) - 1):
        l, r = starts[i], starts[i + 1]
        out[l:r] = l + label[l:r].argsort(kind='mergesort')

    return out


@cython.boundscheck(False)
@cython.wraparound(False)
def count_level_2d(ndarray[uint8_t, ndim=2, cast=True] mask,
                   ndarray[int64_t, ndim=1] labels,
                   Py_ssize_t max_bin,
                   int axis):
    cdef:
        Py_ssize_t i, j, k, n
        ndarray[int64_t, ndim=2] counts

    assert(axis == 0 or axis == 1)
    n, k = (<object> mask).shape

    if axis == 0:
        counts = np.zeros((max_bin, k), dtype='i8')
        with nogil:
            for i from 0 <= i < n:
                for j from 0 <= j < k:
                    counts[labels[i], j] += mask[i, j]

    else:  # axis == 1
        counts = np.zeros((n, max_bin), dtype='i8')
        with nogil:
            for i from 0 <= i < n:
                for j from 0 <= j < k:
                    counts[i, labels[j]] += mask[i, j]

    return counts


def generate_slices(ndarray[int64_t] labels, Py_ssize_t ngroups):
    cdef:
        Py_ssize_t i, group_size, n, start
        int64_t lab
        object slobj
        ndarray[int64_t] starts, ends

    n = len(labels)

    starts = np.zeros(ngroups, dtype=np.int64)
    ends = np.zeros(ngroups, dtype=np.int64)

    start = 0
    group_size = 0
    for i in range(n):
        lab = labels[i]
        if lab < 0:
            start += 1
        else:
            group_size += 1
            if i == n - 1 or lab != labels[i + 1]:
                starts[lab] = start
                ends[lab] = start + group_size
                start += group_size
                group_size = 0

    return starts, ends


def indices_fast(object index, ndarray[int64_t] labels, list keys,
                 list sorted_labels):
    cdef:
        Py_ssize_t i, j, k, lab, cur, start, n = len(labels)
        dict result = {}
        object tup

    k = len(keys)

    if n == 0:
        return result

    start = 0
    cur = labels[0]
    for i in range(1, n):
        lab = labels[i]

        if lab != cur:
            if lab != -1:
                tup = PyTuple_New(k)
                for j in range(k):
                    val = util.get_value_at(keys[j],
                                            sorted_labels[j][i - 1])
                    PyTuple_SET_ITEM(tup, j, val)
                    Py_INCREF(val)

                result[tup] = index[start:i]
            start = i
        cur = lab

    tup = PyTuple_New(k)
    for j in range(k):
        val = util.get_value_at(keys[j],
                                sorted_labels[j][n - 1])
        PyTuple_SET_ITEM(tup, j, val)
        Py_INCREF(val)
    result[tup] = index[start:]

    return result


@cython.boundscheck(False)
@cython.wraparound(False)
def get_blkno_indexers(int64_t[:] blknos, bint group=True):
    """
    Enumerate contiguous runs of integers in ndarray.

    Iterate over elements of `blknos` yielding ``(blkno, slice(start, stop))``
    pairs for each contiguous run found.

    If `group` is True and there is more than one run for a certain blkno,
    ``(blkno, array)`` with an array containing positions of all elements equal
    to blkno.

    Returns
    -------
    iter : iterator of (int, slice or array)

    """
    # There's blkno in this function's name because it's used in block &
    # blockno handling.
    cdef:
        int64_t cur_blkno
        Py_ssize_t i, start, stop, n, diff

        object blkno
        list group_order
        dict group_slices
        int64_t[:] res_view

    n = blknos.shape[0]

    if n == 0:
        return

    start = 0
    cur_blkno = blknos[start]

    if group == False:
        for i in range(1, n):
            if blknos[i] != cur_blkno:
                yield cur_blkno, slice(start, i)

                start = i
                cur_blkno = blknos[i]

        yield cur_blkno, slice(start, n)
    else:
        group_order = []
        group_dict = {}

        for i in range(1, n):
            if blknos[i] != cur_blkno:
                if cur_blkno not in group_dict:
                    group_order.append(cur_blkno)
                    group_dict[cur_blkno] = [(start, i)]
                else:
                    group_dict[cur_blkno].append((start, i))

                start = i
                cur_blkno = blknos[i]

        if cur_blkno not in group_dict:
            group_order.append(cur_blkno)
            group_dict[cur_blkno] = [(start, n)]
        else:
            group_dict[cur_blkno].append((start, n))

        for blkno in group_order:
            slices = group_dict[blkno]
            if len(slices) == 1:
                yield blkno, slice(slices[0][0], slices[0][1])
            else:
                tot_len = sum(stop - start for start, stop in slices)
                result = np.empty(tot_len, dtype=np.int64)
                res_view = result

                i = 0
                for start, stop in slices:
                    for diff in range(start, stop):
                        res_view[i] = diff
                        i += 1

                yield blkno, result


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef slice indexer_as_slice(int64_t[:] vals):
    cdef:
        Py_ssize_t i, n, start, stop
        int64_t d

    if vals is None:
        raise TypeError("vals must be ndarray")

    n = vals.shape[0]

    if n == 0 or vals[0] < 0:
        return None

    if n == 1:
        return slice(vals[0], vals[0] + 1, 1)

    if vals[1] < 0:
        return None

    # n > 2
    d = vals[1] - vals[0]

    if d == 0:
        return None

    for i in range(2, n):
        if vals[i] < 0 or vals[i] - vals[i - 1] != d:
            return None

    start = vals[0]
    stop = start + n * d
    if stop < 0 and d < 0:
        return slice(start, None, d)
    else:
        return slice(start, stop, d)


cpdef slice_canonize(slice s):
    """
    Convert slice to canonical bounded form.
    """
    cdef:
        Py_ssize_t start = 0, stop = 0, step = 1, length

    if s.step is None:
        step = 1
    else:
        step = <Py_ssize_t>s.step
        if step == 0:
            raise ValueError("slice step cannot be zero")

    if step > 0:
        if s.stop is None:
            raise ValueError("unbounded slice")

        stop = <Py_ssize_t>s.stop
        if s.start is None:
            start = 0
        else:
            start = <Py_ssize_t>s.start
            if start > stop:
                start = stop
    elif step < 0:
        if s.start is None:
            raise ValueError("unbounded slice")

        start = <Py_ssize_t>s.start
        if s.stop is None:
            stop = -1
        else:
            stop = <Py_ssize_t>s.stop
            if stop > start:
                stop = start

    if start < 0 or (stop < 0 and s.stop is not None):
        raise ValueError("unbounded slice")

    if stop < 0:
        return slice(start, None, step)
    else:
        return slice(start, stop, step)


cpdef slice_get_indices_ex(slice slc, Py_ssize_t objlen=PY_SSIZE_T_MAX):
    """
    Get (start, stop, step, length) tuple for a slice.

    If `objlen` is not specified, slice must be bounded, otherwise the result
    will be wrong.

    """
    cdef:
        Py_ssize_t start, stop, step, length

    if slc is None:
        raise TypeError("slc should be a slice")

    slice_get_indices(<PyObject *>slc, objlen,
                      &start, &stop, &step, &length)

    return start, stop, step, length


cpdef Py_ssize_t slice_len(
        slice slc, Py_ssize_t objlen=PY_SSIZE_T_MAX) except -1:
    """
    Get length of a bounded slice.

    The slice must not have any "open" bounds that would create dependency on
    container size, i.e.:
    - if ``s.step is None or s.step > 0``, ``s.stop`` is not ``None``
    - if ``s.step < 0``, ``s.start`` is not ``None``

    Otherwise, the result is unreliable.

    """
    cdef:
        Py_ssize_t start, stop, step, length

    if slc is None:
        raise TypeError("slc must be slice")

    slice_get_indices(<PyObject *>slc, objlen,
                      &start, &stop, &step, &length)

    return length


def slice_getitem(slice slc not None, ind):
    cdef:
        Py_ssize_t s_start, s_stop, s_step, s_len
        Py_ssize_t ind_start, ind_stop, ind_step, ind_len

    s_start, s_stop, s_step, s_len = slice_get_indices_ex(slc)

    if isinstance(ind, slice):
        ind_start, ind_stop, ind_step, ind_len = slice_get_indices_ex(ind,
                                                                      s_len)

        if ind_step > 0 and ind_len == s_len:
            # short-cut for no-op slice
            if ind_len == s_len:
                return slc

        if ind_step < 0:
            s_start = s_stop - s_step
            ind_step = -ind_step

        s_step *= ind_step
        s_stop = s_start + ind_stop * s_step
        s_start = s_start + ind_start * s_step

        if s_step < 0 and s_stop < 0:
            return slice(s_start, None, s_step)
        else:
            return slice(s_start, s_stop, s_step)

    else:
        return np.arange(s_start, s_stop, s_step, dtype=np.int64)[ind]


cdef class BlockPlacement:
    # __slots__ = '_as_slice', '_as_array', '_len'
    cdef slice _as_slice
    cdef object _as_array

    cdef bint _has_slice, _has_array, _is_known_slice_like

    def __init__(self, val):
        cdef slice slc

        self._has_slice = False
        self._has_array = False

        if isinstance(val, slice):
            slc = slice_canonize(val)

            if slc.start != slc.stop:
                self._as_slice = slc
                self._has_slice = True
            else:
                arr = np.empty(0, dtype=np.int64)
                self._as_array = arr
                self._has_array = True
        else:
            # Cython memoryview interface requires ndarray to be writeable.
            arr = np.require(val, dtype=np.int64, requirements='W')
            assert arr.ndim == 1
            self._as_array = arr
            self._has_array = True

    def __str__(self):
        cdef slice s = self._ensure_has_slice()
        if s is not None:
            v = self._as_slice
        else:
            v = self._as_array

        return '%s(%r)' % (self.__class__.__name__, v)

    __repr__ = __str__

    def __len__(self):
        cdef slice s = self._ensure_has_slice()
        if s is not None:
            return slice_len(s)
        else:
            return len(self._as_array)

    def __iter__(self):
        cdef slice s = self._ensure_has_slice()
        cdef Py_ssize_t start, stop, step, _
        if s is not None:
            start, stop, step, _ = slice_get_indices_ex(s)
            return iter(range(start, stop, step))
        else:
            return iter(self._as_array)

    @property
    def as_slice(self):
        cdef slice s = self._ensure_has_slice()
        if s is None:
            raise TypeError('Not slice-like')
        else:
            return s

    @property
    def indexer(self):
        cdef slice s = self._ensure_has_slice()
        if s is not None:
            return s
        else:
            return self._as_array

    def isin(self, arr):
        from pandas.core.index import Int64Index
        return Int64Index(self.as_array, copy=False).isin(arr)

    @property
    def as_array(self):
        cdef Py_ssize_t start, stop, end, _
        if not self._has_array:
            start, stop, step, _ = slice_get_indices_ex(self._as_slice)
            self._as_array = np.arange(start, stop, step,
                                       dtype=np.int64)
            self._has_array = True
        return self._as_array

    @property
    def is_slice_like(self):
        cdef slice s = self._ensure_has_slice()
        return s is not None

    def __getitem__(self, loc):
        cdef slice s = self._ensure_has_slice()
        if s is not None:
            val = slice_getitem(s, loc)
        else:
            val = self._as_array[loc]

        if not isinstance(val, slice) and val.ndim == 0:
            return val

        return BlockPlacement(val)

    def delete(self, loc):
        return BlockPlacement(np.delete(self.as_array, loc, axis=0))

    def append(self, others):
        if len(others) == 0:
            return self

        return BlockPlacement(np.concatenate([self.as_array] +
                                             [o.as_array for o in others]))

    cdef iadd(self, other):
        cdef slice s = self._ensure_has_slice()
        cdef Py_ssize_t other_int, start, stop, step, l

        if isinstance(other, int) and s is not None:
            other_int = <Py_ssize_t>other

            if other_int == 0:
                return self

            start, stop, step, l = slice_get_indices_ex(s)
            start += other_int
            stop += other_int

            if ((step > 0 and start < 0) or
                    (step < 0 and stop < step)):
                raise ValueError("iadd causes length change")

            if stop < 0:
                self._as_slice = slice(start, None, step)
            else:
                self._as_slice = slice(start, stop, step)

            self._has_array = False
            self._as_array = None
        else:
            newarr = self.as_array + other
            if (newarr < 0).any():
                raise ValueError("iadd causes length change")

            self._as_array = newarr
            self._has_array = True
            self._has_slice = False
            self._as_slice = None

        return self

    cdef BlockPlacement copy(self):
        cdef slice s = self._ensure_has_slice()
        if s is not None:
            return BlockPlacement(s)
        else:
            return BlockPlacement(self._as_array)

    def add(self, other):
        return self.copy().iadd(other)

    def sub(self, other):
        return self.add(-other)

    cdef slice _ensure_has_slice(self):
        if not self._has_slice:
            self._as_slice = indexer_as_slice(self._as_array)
            self._has_slice = True
        return self._as_slice


include "reduce.pyx"
include "inference.pyx"
