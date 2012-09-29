cimport numpy as np
cimport cython
import numpy as np

from numpy cimport *
from numpy cimport NPY_INT32 as NPY_int32
from numpy cimport NPY_INT64 as NPY_int64
from numpy cimport NPY_FLOAT32 as NPY_float32
from numpy cimport NPY_FLOAT64 as NPY_float64

int32 = np.dtype(np.int32)
int64 = np.dtype(np.int64)
float32 = np.dtype(np.float32)
float64 = np.dtype(np.float64)

cdef np.int32_t MINint32 = np.iinfo(np.int32).min
cdef np.int64_t MINint64 = np.iinfo(np.int64).min
cdef np.float32_t MINfloat32 = np.NINF
cdef np.float64_t MINfloat64 = np.NINF

cdef np.int32_t MAXint32 = np.iinfo(np.int32).max
cdef np.int64_t MAXint64 = np.iinfo(np.int64).max
cdef np.float32_t MAXfloat32 = np.inf
cdef np.float64_t MAXfloat64 = np.inf


cdef extern from "numpy/arrayobject.h":
    cdef enum NPY_TYPES:
        NPY_intp "NPY_INTP"

from cpython cimport (PyDict_New, PyDict_GetItem, PyDict_SetItem,
                      PyDict_Contains, PyDict_Keys,
                      Py_INCREF, PyTuple_SET_ITEM,
                      PyTuple_SetItem,
                      PyTuple_New)
from cpython cimport PyFloat_Check
cimport cpython

isnan = np.isnan
cdef double NaN = <double> np.NaN
cdef double nan = NaN
cdef double NAN = nan

from datetime import datetime as pydatetime

# this is our datetime.pxd
from datetime cimport *

cdef int64_t NPY_NAT = util.get_nat()

from khash cimport *

cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

ctypedef unsigned char UChar

cimport util
from util cimport is_array, _checknull, _checknan

cdef extern from "math.h":
    double sqrt(double x)
    double fabs(double)

# import datetime C API
PyDateTime_IMPORT

# initialize numpy
import_array()
import_ufunc()

cpdef map_indices_list(list index):
    '''
    Produce a dict mapping the values of the input array to their respective
    locations.

    Example:
        array(['hi', 'there']) --> {'hi' : 0 , 'there' : 1}

    Better to do this with Cython because of the enormous speed boost.
    '''
    cdef Py_ssize_t i, length
    cdef dict result = {}

    length = len(index)

    for i from 0 <= i < length:
        result[index[i]] = i

    return result


from libc.stdlib cimport malloc, free

def ismember(ndarray arr, set values):
    '''
    Checks whether

    Parameters
    ----------
    arr : ndarray
    values : set

    Returns
    -------
    ismember : ndarray (boolean dtype)
    '''
    cdef:
        Py_ssize_t i, n
        ndarray[uint8_t] result
        object val

    n = len(arr)
    result = np.empty(n, dtype=np.uint8)
    for i in range(n):
        val = util.get_value_at(arr, i)
        if val in values:
            result[i] = 1
        else:
            result[i] = 0

    return result.view(np.bool_)

#----------------------------------------------------------------------
# datetime / io related

cdef int _EPOCH_ORD = 719163

from datetime import date as pydate

cdef inline int64_t gmtime(object date):
    cdef int y, m, d, h, mn, s, days

    y = PyDateTime_GET_YEAR(date)
    m = PyDateTime_GET_MONTH(date)
    d = PyDateTime_GET_DAY(date)
    h = PyDateTime_DATE_GET_HOUR(date)
    mn = PyDateTime_DATE_GET_MINUTE(date)
    s = PyDateTime_DATE_GET_SECOND(date)

    days = pydate(y, m, 1).toordinal() - _EPOCH_ORD + d - 1
    return ((<int64_t> (((days * 24 + h) * 60 + mn))) * 60 + s) * 1000

cpdef object to_datetime(int64_t timestamp):
    return pydatetime.utcfromtimestamp(timestamp / 1000.0)

cpdef object to_timestamp(object dt):
    return gmtime(dt)

def array_to_timestamp(ndarray[object, ndim=1] arr):
    cdef int i, n
    cdef ndarray[int64_t, ndim=1] result

    n = len(arr)
    result = np.empty(n, dtype=np.int64)

    for i from 0 <= i < n:
        result[i] = gmtime(arr[i])

    return result

def time64_to_datetime(ndarray[int64_t, ndim=1] arr):
    cdef int i, n
    cdef ndarray[object, ndim=1] result

    n = len(arr)
    result = np.empty(n, dtype=object)

    for i from 0 <= i < n:
        result[i] = to_datetime(arr[i])

    return result

#----------------------------------------------------------------------
# isnull / notnull related

cdef double INF = <double> np.inf
cdef double NEGINF = -INF

cpdef checknull(object val):
    if util.is_float_object(val) or util.is_complex_object(val):
        return val != val or val == INF or val == NEGINF
    elif util.is_datetime64_object(val):
        return get_datetime64_value(val) == NPY_NAT
    elif isinstance(val, _NaT):
        return True
    elif is_array(val):
        return False
    else:
        return util._checknull(val)

def isscalar(object val):
    return np.isscalar(val) or val is None or isinstance(val, _Timestamp)


@cython.wraparound(False)
@cython.boundscheck(False)
def isnullobj(ndarray[object] arr):
    cdef Py_ssize_t i, n
    cdef object val
    cdef ndarray[uint8_t] result

    n = len(arr)
    result = np.zeros(n, dtype=np.uint8)
    for i from 0 <= i < n:
        result[i] = util._checknull(arr[i])
    return result.view(np.bool_)


@cython.wraparound(False)
@cython.boundscheck(False)
def isnullobj2d(ndarray[object, ndim=2] arr):
    cdef Py_ssize_t i, j, n, m
    cdef object val
    cdef ndarray[uint8_t, ndim=2] result

    n, m = (<object> arr).shape
    result = np.zeros((n, m), dtype=np.uint8)
    for i from 0 <= i < n:
        for j from 0 <= j < m:
            val = arr[i, j]
            if checknull(val):
                result[i, j] = 1
    return result.view(np.bool_)

def list_to_object_array(list obj):
    '''
    Convert list to object ndarray. Seriously can't believe I had to write this
    function
    '''
    cdef:
        Py_ssize_t i, n
        ndarray[object] arr

    n = len(obj)
    arr = np.empty(n, dtype=object)

    for i from 0 <= i < n:
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
def fast_unique_multiple_list_gen(object gen):
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
    '''
    For zipping multiple ndarrays into an ndarray of tuples
    '''
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

def maybe_indices_to_slice(ndarray[int64_t] indices):
    cdef:
        Py_ssize_t i, n = len(indices)

    if n == 0:
        return indices

    for i in range(1, n):
        if indices[i] - indices[i - 1] != 1:
            return indices
    return slice(indices[0], indices[n - 1] + 1)


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

    if flag == cpython.Py_NE:
        for i in range(n):
            x = values[i]
            if _checknull(x):
                result[i] = True
            else:
                result[i] = cpython.PyObject_RichCompareBool(x, val, flag)
    else:
        for i in range(n):
            x = values[i]
            if _checknull(x):
                result[i] = False
            else:
                result[i] = cpython.PyObject_RichCompareBool(x, val, flag)

    return result.view(bool)

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

            if _checknull(x) or _checknull(y):
                result[i] = True
            else:
                result[i] = cpython.PyObject_RichCompareBool(x, y, flag)
    else:
        for i in range(n):
            x = left[i]
            y = right[i]

            if _checknull(x) or _checknull(y):
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

    for i in range(n):
        x = values[i]
        if util._checknull(x):
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
            if util._checknull(x):
                result[i] = x
            elif util._checknull(y):
                result[i] = y
            else:
                raise

    return maybe_convert_bool(result)


def value_count_int64(ndarray[int64_t] values):
    cdef:
        Py_ssize_t i, n = len(values)
        kh_int64_t *table
        int ret = 0
        list uniques = []

    table = kh_init_int64()
    kh_resize_int64(table, n)

    for i in range(n):
        val = values[i]
        k = kh_get_int64(table, val)
        if k != table.n_buckets:
            table.vals[k] += 1
        else:
            k = kh_put_int64(table, val, &ret)
            table.vals[k] = 1

    # for (k = kh_begin(h); k != kh_end(h); ++k)
    # 	if (kh_exist(h, k)) kh_value(h, k) = 1;
    i = 0
    result_keys = np.empty(table.n_occupied, dtype=np.int64)
    result_counts = np.zeros(table.n_occupied, dtype=np.int64)
    for k in range(table.n_buckets):
        if kh_exist_int64(table, k):
            result_keys[i] = table.keys[k]
            result_counts[i] = table.vals[k]
            i += 1
    kh_destroy_int64(table)

    return result_keys, result_counts

def astype_intsafe(ndarray[object] arr, new_dtype):
    cdef:
        Py_ssize_t i, n = len(arr)
        ndarray result

    result = np.empty(n, dtype=new_dtype)
    for i in range(n):
        util.set_value_at(result, i, arr[i])

    return result

include "hashtable.pyx"
include "datetime.pyx"
include "skiplist.pyx"
include "groupby.pyx"
include "moments.pyx"
include "reindex.pyx"
include "reduce.pyx"
include "stats.pyx"
include "properties.pyx"
include "inference.pyx"
include "join.pyx"
include "engines.pyx"
