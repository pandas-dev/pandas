# -*- coding: utf-8 -*-
# cython: profile=False
cimport cython
from cython cimport Py_ssize_t

import numpy as np
cimport numpy as cnp
from numpy cimport (ndarray, PyArray_NDIM, PyArray_GETITEM,
                    PyArray_ITER_DATA, PyArray_ITER_NEXT, PyArray_IterNew,
                    flatiter, NPY_OBJECT,
                    int64_t,
                    float32_t, float64_t,
                    uint8_t, uint64_t,
                    complex128_t)
cnp.import_array()

from cpython cimport (Py_INCREF, PyTuple_SET_ITEM,
                      PyList_Check, PyFloat_Check,
                      PyString_Check,
                      PyBytes_Check,
                      PyUnicode_Check,
                      PyTuple_New,
                      Py_EQ,
                      PyObject_RichCompareBool)

from cpython.datetime cimport (PyDateTime_Check, PyDate_Check,
                               PyTime_Check, PyDelta_Check,
                               PyDateTime_IMPORT)
PyDateTime_IMPORT

from tslib import NaT, array_to_datetime
from missing cimport checknull


cimport util
cdef int64_t NPY_NAT = util.get_nat()
from util cimport is_array, is_nan


def values_from_object(object obj):
    """ return my values or the object if we are say an ndarray """
    cdef func  # TODO: Does declaring this without a type accomplish anything?

    func = getattr(obj, 'get_values', None)
    if func is not None:
        obj = func()

    return obj


@cython.wraparound(False)
@cython.boundscheck(False)
def memory_usage_of_objects(object[:] arr):
    """ return the memory usage of an object array in bytes,
    does not include the actual bytes of the pointers """
    cdef:
        Py_ssize_t i, n
        int64_t size = 0

    n = len(arr)
    for i in range(n):
        size += arr[i].__sizeof__()
    return size


# ----------------------------------------------------------------------


cpdef bint is_scalar(object val):
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
    - DateOffset

    """

    return (cnp.PyArray_IsAnyScalar(val)
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
            or is_interval(val)
            or is_offset(val))


def item_from_zerodim(object val):
    """
    If the value is a zerodim array, return the item it contains.

    Parameters
    ----------
    val : object

    Returns
    -------
    result : object

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
    if cnp.PyArray_IsZeroDim(val):
        return cnp.PyArray_ToScalar(cnp.PyArray_DATA(val), val)
    return val


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

    for i in range(n):
        buf = arrays[i]
        n = len(buf)
        for j in range(n):
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
def fast_unique_multiple_list(list lists, bint sort=True):
    cdef:
        list buf
        Py_ssize_t k = len(lists)
        Py_ssize_t i, j, n
        list uniques = []
        dict table = {}
        object val, stub = 0

    for i in range(k):
        buf = lists[i]
        n = len(buf)
        for j in range(n):
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
        for j in range(n):
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
        if not (PyObject_RichCompareBool(x, y, Py_EQ) or
                (x is None or is_nan(x)) and (y is None or is_nan(y))):
            return False
    return True


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
    Utility used in pandas.core.index.ensure_index
    """
    cdef:
        Py_ssize_t i, n = len(obj)
        object v
        bint all_arrays = 1

    for i in range(n):
        v = obj[i]
        if not (PyList_Check(v) or util.is_array(v) or hasattr(v, '_data')):
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


# ------------------------------------------------------------------------------
# Groupby-related functions

# TODO: could do even better if we know something about the data. eg, index has
# 1-min data, binner has 5-min data, then bins are just strides in index. This
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
            for i in range(n):
                for j in range(k):
                    counts[labels[i], j] += mask[i, j]

    else:  # axis == 1
        counts = np.zeros((n, max_bin), dtype='i8')
        with nogil:
            for i in range(n):
                for j in range(k):
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


include "src/inference.pyx"
