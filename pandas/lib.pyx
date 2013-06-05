cimport numpy as np
cimport cython
import numpy as np

from numpy cimport *


cdef extern from "numpy/arrayobject.h":
    cdef enum NPY_TYPES:
        NPY_intp "NPY_INTP"

from cpython cimport (PyDict_New, PyDict_GetItem, PyDict_SetItem,
                      PyDict_Contains, PyDict_Keys,
                      Py_INCREF, PyTuple_SET_ITEM,
                      PyList_Check, PyFloat_Check,
                      PyString_Check,
		      PyBytes_Check,
                      PyTuple_SetItem,
                      PyTuple_New,
                      PyObject_SetAttrString)

cimport cpython

isnan = np.isnan
cdef double NaN = <double> np.NaN
cdef double nan = NaN
cdef double NAN = nan

from datetime import datetime as pydatetime

# this is our tseries.pxd
from datetime cimport *

from tslib cimport convert_to_tsobject
import tslib
from tslib import NaT, Timestamp, repr_timedelta64

cdef int64_t NPY_NAT = util.get_nat()

ctypedef unsigned char UChar

cimport util
from util cimport is_array, _checknull, _checknan

cdef extern from "headers/stdint.h":
    enum: UINT8_MAX
    enum: INT64_MAX
    enum: INT64_MIN


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

cdef inline int64_t get_timedelta64_value(val):
    return val.view('i8')

#----------------------------------------------------------------------
# isnull / notnull related

cdef double INF = <double> np.inf
cdef double NEGINF = -INF


cpdef checknull(object val):
    if util.is_float_object(val) or util.is_complex_object(val):
        return val != val # and val != INF and val != NEGINF
    elif util.is_datetime64_object(val):
        return get_datetime64_value(val) == NPY_NAT
    elif val is NaT:
        return True
    elif util.is_timedelta64_object(val):
        return get_timedelta64_value(val) == NPY_NAT
    elif is_array(val):
        return False
    else:
        return util._checknull(val)

cpdef checknull_old(object val):
    if util.is_float_object(val) or util.is_complex_object(val):
        return val != val or val == INF or val == NEGINF
    elif util.is_datetime64_object(val):
        return get_datetime64_value(val) == NPY_NAT
    elif val is NaT:
        return True
    elif util.is_timedelta64_object(val):
        return get_timedelta64_value(val) == NPY_NAT
    elif is_array(val):
        return False
    else:
        return util._checknull(val)

def isscalar(object val):
    return np.isscalar(val) or val is None or PyDateTime_Check(val)


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
def isnullobj_old(ndarray[object] arr):
    cdef Py_ssize_t i, n
    cdef object val
    cdef ndarray[uint8_t] result

    n = len(arr)
    result = np.zeros(n, dtype=np.uint8)
    for i from 0 <= i < n:
        result[i] = util._checknull_old(arr[i])
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

@cython.wraparound(False)
@cython.boundscheck(False)
def isnullobj_old(ndarray[object] arr):
    cdef Py_ssize_t i, n
    cdef object val
    cdef ndarray[uint8_t] result

    n = len(arr)
    result = np.zeros(n, dtype=np.uint8)
    for i from 0 <= i < n:
        result[i] = util._checknull_old(arr[i])
    return result.view(np.bool_)


@cython.wraparound(False)
@cython.boundscheck(False)
def isnullobj2d_old(ndarray[object, ndim=2] arr):
    cdef Py_ssize_t i, j, n, m
    cdef object val
    cdef ndarray[uint8_t, ndim=2] result

    n, m = (<object> arr).shape
    result = np.zeros((n, m), dtype=np.uint8)
    for i from 0 <= i < n:
        for j from 0 <= j < m:
            val = arr[i, j]
            if checknull_old(val):
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


def astype_intsafe(ndarray[object] arr, new_dtype):
    cdef:
        Py_ssize_t i, n = len(arr)
        ndarray result

    result = np.empty(n, dtype=new_dtype)
    for i in range(n):
        util.set_value_at(result, i, arr[i])

    return result

def clean_index_list(list obj):
    '''
    Utility used in pandas.core.index._ensure_index
    '''
    cdef:
        ndarray[object] converted
        Py_ssize_t i, n = len(obj)
        object v
        bint all_arrays = 1

    for i in range(n):
        v = obj[i]
        if not (PyList_Check(v) or np.PyArray_Check(v)):
            all_arrays = 0
            break

    if all_arrays:
        return obj, all_arrays

    converted = np.empty(n, dtype=object)
    for i in range(n):
        v = obj[i]
        if PyList_Check(v) or np.PyArray_Check(v):
            converted[i] = tuple(v)
        else:
            converted[i] = v

    return maybe_convert_objects(converted), 0

@cython.boundscheck(False)
@cython.wraparound(False)
def max_len_string_array(ndarray[object, ndim=1] arr):
    """ return the maximum size of elements in a 1-dim string array """
    cdef:
        int i, m, l
        length = arr.shape[0]
        object v

    m = 0
    for i from 0 <= i < length:
        v = arr[i]
        if PyString_Check(v) or PyBytes_Check(v):
            l = len(v)

            if l > m:
                m = l

    return m

@cython.boundscheck(False)
@cython.wraparound(False)
def string_array_replace_from_nan_rep(ndarray[object, ndim=1] arr, object nan_rep, object replace = None):
    """ replace the values in the array with replacement if they are nan_rep; return the same array """

    cdef int length = arr.shape[0], i = 0
    if replace is None:
        replace = np.nan

    for i from 0 <= i < length:
        if arr[i] == nan_rep:
            arr[i] = replace

    return arr

@cython.boundscheck(False)
@cython.wraparound(False)
def write_csv_rows(list data, list data_index, int nlevels, list cols, object writer):

    cdef int N, j, i, ncols
    cdef list rows
    cdef object val

    # In crude testing, N>100 yields little marginal improvement
    N=100

    # pre-allocate  rows
    ncols = len(cols)
    rows = [[None]*(nlevels+ncols) for x in range(N)]

    j = -1
    if nlevels == 1:
        for j in range(len(data_index)):
            row = rows[j % N]
            row[0] = data_index[j]
            for i in range(ncols):
                row[1+i] = data[i][j]

            if j >= N-1 and j % N == N-1:
                writer.writerows(rows)
    elif nlevels > 1:
        for j in range(len(data_index)):
            row = rows[j % N]
            row[:nlevels] = list(data_index[j])
            for i in range(ncols):
                row[nlevels+i] = data[i][j]

            if j >= N-1 and j % N == N-1:
                writer.writerows(rows)
    else:
        for j in range(len(data_index)):
            row = rows[j % N]
            for i in range(ncols):
                row[i] = data[i][j]

            if j >= N-1 and j % N == N-1:
                writer.writerows(rows)

    if  j >= 0 and (j < N-1 or (j % N) != N-1 ):
        writer.writerows(rows[:((j+1) % N)])


@cython.boundscheck(False)
@cython.wraparound(False)
def create_hdf_rows_2d(ndarray indexer0,
                       object dtype,
                       ndarray[np.uint8_t, ndim=1] mask,
                       ndarray[np.uint8_t, ndim=1] searchable,	 
                       list values):
    """ return a list of objects ready to be converted to rec-array format """

    cdef:
        int i, l, b, n_indexer0, n_blocks, tup_size
        ndarray result
        tuple tup
        object v

    n_indexer0 = indexer0.shape[0]
    n_blocks   = len(values)
    tup_size   = n_blocks+1

    result = np.empty(n_indexer0,dtype=dtype)
    l = 0
    for i in range(n_indexer0):

        if not mask[i]:
         
            tup = PyTuple_New(tup_size)

            v  = indexer0[i]
            PyTuple_SET_ITEM(tup, 0, v)
            Py_INCREF(v)

            for b in range(n_blocks):

                v = values[b][i]
                if searchable[b]:
                    v = v[0]
        
                PyTuple_SET_ITEM(tup, b+1, v)
                Py_INCREF(v)

            result[l] = tup
            l += 1

    return result[0:l]

@cython.boundscheck(False)
@cython.wraparound(False)
def create_hdf_rows_3d(ndarray indexer0, ndarray indexer1,
                       object dtype,
                       ndarray[np.uint8_t, ndim=2] mask, 
                       ndarray[np.uint8_t, ndim=1] searchable,	 
                       list values):
    """ return a list of objects ready to be converted to rec-array format """

    cdef:
        int i, j, l, b, n_indexer0, n_indexer1, n_blocks, tup_size
        tuple tup
        object v
        ndarray result

    n_indexer0 = indexer0.shape[0]
    n_indexer1 = indexer1.shape[0]
    n_blocks   = len(values)
    tup_size   = n_blocks+2
    result = np.empty(n_indexer0*n_indexer1,dtype=dtype)
    l = 0
    for i from 0 <= i < n_indexer0:

        for j from 0 <= j < n_indexer1:

            if not mask[i, j]:

                tup = PyTuple_New(tup_size)

                v = indexer0[i]
                PyTuple_SET_ITEM(tup, 0, v)
                Py_INCREF(v)
                v = indexer1[j]
                PyTuple_SET_ITEM(tup, 1, v)
                Py_INCREF(v)

                for b from 0 <= b < n_blocks:

                    v   = values[b][i, j]
                    if searchable[b]:
                        v = v[0]

                    PyTuple_SET_ITEM(tup, b+2, v)
                    Py_INCREF(v)

                result[l] = tup
                l += 1

    return result[0:l]

@cython.boundscheck(False)
@cython.wraparound(False)
def create_hdf_rows_4d(ndarray indexer0, ndarray indexer1, ndarray indexer2,
                       object dtype,
                       ndarray[np.uint8_t, ndim=3] mask, 
                       ndarray[np.uint8_t, ndim=1] searchable,	 
                       list values):
    """ return a list of objects ready to be converted to rec-array format """

    cdef:
        int i, j, k, l, b, n_indexer0, n_indexer1, n_indexer2, n_blocks, tup_size
        tuple tup
        object v
        ndarray result

    n_indexer0 = indexer0.shape[0]
    n_indexer1 = indexer1.shape[0]
    n_indexer2 = indexer2.shape[0]
    n_blocks   = len(values)
    tup_size   = n_blocks+3
    result = np.empty(n_indexer0*n_indexer1*n_indexer2,dtype=dtype)
    l = 0
    for i from 0 <= i < n_indexer0:

        for j from 0 <= j < n_indexer1:

            for k from 0 <= k < n_indexer2:

                if not mask[i, j, k]:

                    tup = PyTuple_New(tup_size)

                    v = indexer0[i]
                    PyTuple_SET_ITEM(tup, 0, v)
                    Py_INCREF(v)
                    v = indexer1[j]
                    PyTuple_SET_ITEM(tup, 1, v)
                    Py_INCREF(v)
                    v = indexer2[k]
                    PyTuple_SET_ITEM(tup, 2, v)
                    Py_INCREF(v)

                    for b from 0 <= b < n_blocks:

                        v   = values[b][i, j, k]
                        if searchable[b]:
                            v = v[0]
                        PyTuple_SET_ITEM(tup, b+3, v)
                        Py_INCREF(v)

                    result[l] = tup
                    l += 1

    return result[0:l]

#-------------------------------------------------------------------------------
# Groupby-related functions

@cython.boundscheck(False)
def arrmap(ndarray[object] index, object func):
    cdef int length = index.shape[0]
    cdef int i = 0

    cdef ndarray[object] result = np.empty(length, dtype=np.object_)

    for i from 0 <= i < length:
        result[i] = func(index[i])

    return result

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
        # vecs[i] = <int64_t *> (<ndarray> list_of_arrays[i]).data

        arr = list_of_arrays[i]
        vecs[i] = <int64_t *> arr.data
    # assume uniqueness??

    for i from 1 <= i < n:
        for k from 0 <= k < nlevels:
            cur = vecs[k][i]
            pre = vecs[k][i-1]
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
                       object closed='left'):
    """
    Int64 (datetime64) version of generic python version in groupby.py
    """
    cdef:
        Py_ssize_t lenidx, lenbin, i, j, bc, vc
        ndarray[int64_t] bins
        int64_t l_bin, r_bin
        bint right_closed = closed == 'right'

    lenidx = len(values)
    lenbin = len(binner)

    if lenidx <= 0 or lenbin <= 0:
        raise ValueError("Invalid length for values or for binner")

    # check binner fits data
    if values[0] < binner[0]:
        raise ValueError("Values falls before first bin")

    if values[lenidx-1] > binner[lenbin-1]:
        raise ValueError("Values falls after last bin")

    bins   = np.empty(lenbin - 1, dtype=np.int64)

    j  = 0 # index into values
    bc = 0 # bin count

    # linear scan
    for i in range(0, lenbin - 1):
        l_bin = binner[i]
        r_bin = binner[i+1]

        # count values in current bin, advance to next bin
        while j < lenidx and (values[j] < r_bin or
                              (right_closed and values[j] == r_bin)):
            j += 1

        bins[bc] = j
        bc += 1

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


def group_count(ndarray[int64_t] values, Py_ssize_t size):
    cdef:
        Py_ssize_t i, n = len(values)
        ndarray[int64_t] counts

    counts = np.zeros(size, dtype=np.int64)
    for i in range(n):
        counts[values[i]] += 1
    return counts

def lookup_values(ndarray[object] values, dict mapping):
    cdef:
        Py_ssize_t i, n = len(values)

    result = np.empty(n, dtype='O')
    for i in range(n):
        result[i] = mapping[values[i]]
    return maybe_convert_objects(result)


def count_level_1d(ndarray[uint8_t, cast=True] mask,
                   ndarray[int64_t] labels, Py_ssize_t max_bin):
    cdef:
        Py_ssize_t i, n
        ndarray[int64_t] counts

    counts = np.zeros(max_bin, dtype='i8')

    n = len(mask)

    for i from 0 <= i < n:
        if mask[i]:
            counts[labels[i]] += 1

    return counts


def count_level_2d(ndarray[uint8_t, ndim=2, cast=True] mask,
                   ndarray[int64_t] labels, Py_ssize_t max_bin):
    cdef:
        Py_ssize_t i, j, k, n
        ndarray[int64_t, ndim=2] counts

    n, k = (<object> mask).shape
    counts = np.zeros((max_bin, k), dtype='i8')

    for i from 0 <= i < n:
        for j from 0 <= j < k:
            if mask[i, j]:
                counts[labels[i], j] += 1

    return counts

cdef class _PandasNull:

    def __richcmp__(_PandasNull self, object other, int op):
        if op == 2: # ==
            return isinstance(other, _PandasNull)
        elif op == 3: # !=
            return not isinstance(other, _PandasNull)
        else:
            return False

    def __hash__(self):
        return 0

pandas_null = _PandasNull()

def fast_zip_fillna(list ndarrays, fill_value=pandas_null):
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

        if val != val:
            val = fill_value

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
            if val != val:
                val = fill_value

            PyTuple_SET_ITEM(result[i], j, val)
            Py_INCREF(val)
            PyArray_ITER_NEXT(it)

    return result

def duplicated(ndarray[object] values, take_last=False):
    cdef:
        Py_ssize_t i, n
        dict seen = {}
        object row

    n = len(values)
    cdef ndarray[uint8_t] result = np.zeros(n, dtype=np.uint8)

    if take_last:
        for i from n > i >= 0:
            row = values[i]

            if row in seen:
                result[i] = 1
            else:
                seen[row] = None
                result[i] = 0
    else:
        for i from 0 <= i < n:
            row = values[i]
            if row in seen:
                result[i] = 1
            else:
                seen[row] = None
                result[i] = 0

    return result.view(np.bool_)

def generate_slices(ndarray[int64_t] labels, Py_ssize_t ngroups):
    cdef:
        Py_ssize_t i, group_size, n, lab, start
        object slobj
        ndarray[int64_t] starts

    n = len(labels)

    starts = np.zeros(ngroups, dtype=np.int64)
    ends = np.zeros(ngroups, dtype=np.int64)

    start = 0
    group_size = 0
    for i in range(n):
        group_size += 1
        lab = labels[i]
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
                                            sorted_labels[j][i-1])
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

include "reduce.pyx"
include "properties.pyx"
include "inference.pyx"
