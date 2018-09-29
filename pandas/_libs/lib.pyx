# -*- coding: utf-8 -*-
from decimal import Decimal
import sys

import cython
from cython import Py_ssize_t

from cpython cimport (Py_INCREF, PyTuple_SET_ITEM,
                      PyTuple_New,
                      Py_EQ,
                      PyObject_RichCompareBool)

from cpython.datetime cimport (PyDateTime_Check, PyDate_Check,
                               PyTime_Check, PyDelta_Check,
                               PyDateTime_IMPORT)
PyDateTime_IMPORT


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

cdef extern from "numpy/arrayobject.h":
    # cython's numpy.dtype specification is incorrect, which leads to
    # errors in issubclass(self.dtype.type, np.bool_), so we directly
    # include the correct version
    # https://github.com/cython/cython/issues/2022

    ctypedef class numpy.dtype [object PyArray_Descr]:
        # Use PyDataType_* macros when possible, however there are no macros
        # for accessing some of the fields, so some are defined. Please
        # ask on cython-dev if you need more.
        cdef int type_num
        cdef int itemsize "elsize"
        cdef char byteorder
        cdef object fields
        cdef tuple names


cdef extern from "src/parse_helper.h":
    int floatify(object, double *result, int *maybe_int) except -1

cimport util
from util cimport (is_nan,
                   UINT8_MAX, UINT64_MAX, INT64_MAX, INT64_MIN)

from tslib import array_to_datetime
from tslibs.nattype import NaT
from tslibs.conversion cimport convert_to_tsobject
from tslibs.timedeltas cimport convert_to_timedelta64
from tslibs.timezones cimport get_timezone, tz_compare

from missing cimport (checknull,
                      is_null_datetime64, is_null_timedelta64, is_null_period)


# constants that will be compared to potentially arbitrarily large
# python int
cdef object oINT64_MAX = <int64_t>INT64_MAX
cdef object oINT64_MIN = <int64_t>INT64_MIN
cdef object oUINT64_MAX = <uint64_t>UINT64_MAX

cdef int64_t NPY_NAT = util.get_nat()
iNaT = util.get_nat()

cdef bint PY2 = sys.version_info[0] == 2
cdef double nan = <double>np.NaN


def values_from_object(object obj):
    """ return my values or the object if we are say an ndarray """
    cdef func  # TODO: Does declaring this without a type accomplish anything?

    func = getattr(obj, 'get_values', None)
    if func is not None:
        obj = func()

    return obj


@cython.wraparound(False)
@cython.boundscheck(False)
def memory_usage_of_objects(arr: object[:]) -> int64_t:
    """ return the memory usage of an object array in bytes,
    does not include the actual bytes of the pointers """
    i: Py_ssize_t
    n: Py_ssize_t
    size: int64_t

    size = 0
    n = len(arr)
    for i in range(n):
        size += arr[i].__sizeof__()
    return size


# ----------------------------------------------------------------------


def is_scalar(val: object) -> bint:
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
            or isinstance(val, bytes)
            # We differ from numpy (as of 1.10), which claims that None is
            # not scalar in np.isscalar().
            or val is None
            or PyDate_Check(val)
            or PyDelta_Check(val)
            or PyTime_Check(val)
            or util.is_period_object(val)
            or is_decimal(val)
            or is_interval(val)
            or util.is_offset_object(val))


def item_from_zerodim(val: object) -> object:
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

    for i in range(k):
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


def has_infs_f4(ndarray[float32_t] arr) -> bint:
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


def has_infs_f8(ndarray[float64_t] arr) -> bint:
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
def array_equivalent_object(left: object[:], right: object[:]) -> bint:
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
            result[i] = v

    return result


def astype_unicode(arr: ndarray) -> ndarray[object]:
    cdef:
        Py_ssize_t i, n = arr.size
        ndarray[object] result = np.empty(n, dtype=object)

    for i in range(n):
        result[i] = unicode(arr[i])

    return result


def astype_str(arr: ndarray) -> ndarray[object]:
    cdef:
        Py_ssize_t i, n = arr.size
        ndarray[object] result = np.empty(n, dtype=object)

    for i in range(n):
        result[i] = str(arr[i])

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
        if not (isinstance(v, list) or
                util.is_array(v) or hasattr(v, '_data')):
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
    assert (n == len(mask))

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
    assert (n == len(mask))

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

    assert (axis == 0 or axis == 1)
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


# core.common import for fast inference checks

def is_float(obj: object) -> bint:
    return util.is_float_object(obj)


def is_integer(obj: object) -> bint:
    return util.is_integer_object(obj)


def is_bool(obj: object) -> bint:
    return util.is_bool_object(obj)


def is_complex(obj: object) -> bint:
    return util.is_complex_object(obj)


cpdef bint is_decimal(object obj):
    return isinstance(obj, Decimal)


cpdef bint is_interval(object obj):
    return getattr(obj, '_typ', '_typ') == 'interval'


def is_period(val: object) -> bint:
    """ Return a boolean if this is a Period object """
    return util.is_period_object(val)


_TYPE_MAP = {
    'categorical': 'categorical',
    'category': 'categorical',
    'int8': 'integer',
    'int16': 'integer',
    'int32': 'integer',
    'int64': 'integer',
    'i': 'integer',
    'uint8': 'integer',
    'uint16': 'integer',
    'uint32': 'integer',
    'uint64': 'integer',
    'u': 'integer',
    'float32': 'floating',
    'float64': 'floating',
    'f': 'floating',
    'complex128': 'complex',
    'c': 'complex',
    'string': 'string' if PY2 else 'bytes',
    'S': 'string' if PY2 else 'bytes',
    'unicode': 'unicode' if PY2 else 'string',
    'U': 'unicode' if PY2 else 'string',
    'bool': 'boolean',
    'b': 'boolean',
    'datetime64[ns]': 'datetime64',
    'M': 'datetime64',
    'timedelta64[ns]': 'timedelta64',
    'm': 'timedelta64',
}

# types only exist on certain platform
try:
    np.float128
    _TYPE_MAP['float128'] = 'floating'
except AttributeError:
    pass
try:
    np.complex256
    _TYPE_MAP['complex256'] = 'complex'
except AttributeError:
    pass
try:
    np.float16
    _TYPE_MAP['float16'] = 'floating'
except AttributeError:
    pass


cdef class Seen(object):
    """
    Class for keeping track of the types of elements
    encountered when trying to perform type conversions.
    """

    cdef:
        bint int_             # seen_int
        bint bool_            # seen_bool
        bint null_            # seen_null
        bint uint_            # seen_uint (unsigned integer)
        bint sint_            # seen_sint (signed integer)
        bint float_           # seen_float
        bint object_          # seen_object
        bint complex_         # seen_complex
        bint datetime_        # seen_datetime
        bint coerce_numeric   # coerce data to numeric
        bint timedelta_       # seen_timedelta
        bint datetimetz_      # seen_datetimetz

    def __cinit__(self, bint coerce_numeric=0):
        """
        Initialize a Seen instance.

        Parameters
        ----------
        coerce_numeric : bint, default 0
            Whether or not to force conversion to a numeric data type if
            initial methods to convert to numeric fail.
        """
        self.int_ = 0
        self.bool_ = 0
        self.null_ = 0
        self.uint_ = 0
        self.sint_ = 0
        self.float_ = 0
        self.object_ = 0
        self.complex_ = 0
        self.datetime_ = 0
        self.timedelta_ = 0
        self.datetimetz_ = 0
        self.coerce_numeric = coerce_numeric

    cdef inline bint check_uint64_conflict(self) except -1:
        """
        Check whether we can safely convert a uint64 array to a numeric dtype.

        There are two cases when conversion to numeric dtype with a uint64
        array is not safe (and will therefore not be performed)

        1) A NaN element is encountered.

           uint64 cannot be safely cast to float64 due to truncation issues
           at the extreme ends of the range.

        2) A negative number is encountered.

           There is no numerical dtype that can hold both negative numbers
           and numbers greater than INT64_MAX. Hence, at least one number
           will be improperly cast if we convert to a numeric dtype.

        Returns
        -------
        return_values : bool
            Whether or not we should return the original input array to avoid
            data truncation.

        Raises
        ------
        ValueError : uint64 elements were detected, and at least one of the
                     two conflict cases was also detected. However, we are
                     trying to force conversion to a numeric dtype.
        """
        return (self.uint_ and (self.null_ or self.sint_)
                and not self.coerce_numeric)

    cdef inline saw_null(self):
        """
        Set flags indicating that a null value was encountered.
        """
        self.null_ = 1
        self.float_ = 1

    cdef saw_int(self, object val):
        """
        Set flags indicating that an integer value was encountered.

        In addition to setting a flag that an integer was seen, we
        also set two flags depending on the type of integer seen:

        1) sint_ : a negative (signed) number in the
                   range of [-2**63, 0) was encountered
        2) uint_ : a positive number in the range of
                   [2**63, 2**64) was encountered

        Parameters
        ----------
        val : Python int
            Value with which to set the flags.
        """
        self.int_ = 1
        self.sint_ = self.sint_ or (oINT64_MIN <= val < 0)
        self.uint_ = self.uint_ or (oINT64_MAX < val <= oUINT64_MAX)

    @property
    def numeric_(self):
        return self.complex_ or self.float_ or self.int_

    @property
    def is_bool(self):
        return not (self.datetime_ or self.numeric_ or self.timedelta_)

    @property
    def is_float_or_complex(self):
        return not (self.bool_ or self.datetime_ or self.timedelta_)


cdef _try_infer_map(v):
    """ if its in our map, just return the dtype """
    cdef:
        object attr, val
    for attr in ['name', 'kind', 'base']:
        val = getattr(v.dtype, attr)
        if val in _TYPE_MAP:
            return _TYPE_MAP[val]
    return None


def infer_dtype(object value, bint skipna=False):
    """
    Efficiently infer the type of a passed val, or list-like
    array of values. Return a string describing the type.

    Parameters
    ----------
    value : scalar, list, ndarray, or pandas type
    skipna : bool, default False
        Ignore NaN values when inferring the type. The default of ``False``
        will be deprecated in a later version of pandas.

        .. versionadded:: 0.21.0

    Returns
    -------
    string describing the common type of the input data.
    Results can include:

    - string
    - unicode
    - bytes
    - floating
    - integer
    - mixed-integer
    - mixed-integer-float
    - decimal
    - complex
    - categorical
    - boolean
    - datetime64
    - datetime
    - date
    - timedelta64
    - timedelta
    - time
    - period
    - mixed

    Raises
    ------
    TypeError if ndarray-like but cannot infer the dtype

    Notes
    -----
    - 'mixed' is the catchall for anything that is not otherwise
      specialized
    - 'mixed-integer-float' are floats and integers
    - 'mixed-integer' are integers mixed with non-integers

    Examples
    --------
    >>> infer_dtype(['foo', 'bar'])
    'string'

    >>> infer_dtype(['a', np.nan, 'b'], skipna=True)
    'string'

    >>> infer_dtype(['a', np.nan, 'b'], skipna=False)
    'mixed'

    >>> infer_dtype([b'foo', b'bar'])
    'bytes'

    >>> infer_dtype([1, 2, 3])
    'integer'

    >>> infer_dtype([1, 2, 3.5])
    'mixed-integer-float'

    >>> infer_dtype([1.0, 2.0, 3.5])
    'floating'

    >>> infer_dtype(['a', 1])
    'mixed-integer'

    >>> infer_dtype([Decimal(1), Decimal(2.0)])
    'decimal'

    >>> infer_dtype([True, False])
    'boolean'

    >>> infer_dtype([True, False, np.nan])
    'mixed'

    >>> infer_dtype([pd.Timestamp('20130101')])
    'datetime'

    >>> infer_dtype([datetime.date(2013, 1, 1)])
    'date'

    >>> infer_dtype([np.datetime64('2013-01-01')])
    'datetime64'

    >>> infer_dtype([datetime.timedelta(0, 1, 1)])
    'timedelta'

    >>> infer_dtype(pd.Series(list('aabc')).astype('category'))
    'categorical'
    """
    cdef:
        Py_ssize_t i, n
        object val
        ndarray values
        bint seen_pdnat = False
        bint seen_val = False

    if util.is_array(value):
        values = value
    elif hasattr(value, 'dtype'):
        # this will handle ndarray-like
        # e.g. categoricals
        try:
            values = getattr(value, '_values', getattr(value, 'values', value))
        except:
            value = _try_infer_map(value)
            if value is not None:
                return value

            # its ndarray like but we can't handle
            raise ValueError("cannot infer type for {typ}"
                             .format(typ=type(value)))

    else:
        if not isinstance(value, list):
            value = list(value)
        from pandas.core.dtypes.cast import (
            construct_1d_object_array_from_listlike)
        values = construct_1d_object_array_from_listlike(value)

    values = getattr(values, 'values', values)
    val = _try_infer_map(values)
    if val is not None:
        return val

    if values.dtype != np.object_:
        values = values.astype('O')

    # make contiguous
    values = values.ravel()

    n = len(values)
    if n == 0:
        return 'empty'

    # try to use a valid value
    for i in range(n):
        val = values[i]

        # do not use is_nul_datetimelike to keep
        # np.datetime64('nat') and np.timedelta64('nat')
        if val is None or util.is_nan(val):
            pass
        elif val is NaT:
            seen_pdnat = True
        else:
            seen_val = True
            break

    # if all values are nan/NaT
    if seen_val is False and seen_pdnat is True:
        return 'datetime'
        # float/object nan is handled in latter logic

    if util.is_datetime64_object(val):
        if is_datetime64_array(values):
            return 'datetime64'
        elif is_timedelta_or_timedelta64_array(values):
            return 'timedelta'

    elif is_timedelta(val):
        if is_timedelta_or_timedelta64_array(values):
            return 'timedelta'

    elif util.is_integer_object(val):
        # a timedelta will show true here as well
        if is_timedelta(val):
            if is_timedelta_or_timedelta64_array(values):
                return 'timedelta'

        if is_integer_array(values):
            return 'integer'
        elif is_integer_float_array(values):
            return 'mixed-integer-float'
        elif is_timedelta_or_timedelta64_array(values):
            return 'timedelta'
        return 'mixed-integer'

    elif PyDateTime_Check(val):
        if is_datetime_array(values):
            return 'datetime'

    elif PyDate_Check(val):
        if is_date_array(values, skipna=skipna):
            return 'date'

    elif PyTime_Check(val):
        if is_time_array(values, skipna=skipna):
            return 'time'

    elif is_decimal(val):
        return 'decimal'

    elif util.is_float_object(val):
        if is_float_array(values):
            return 'floating'
        elif is_integer_float_array(values):
            return 'mixed-integer-float'

    elif util.is_bool_object(val):
        if is_bool_array(values, skipna=skipna):
            return 'boolean'

    elif isinstance(val, str):
        if is_string_array(values, skipna=skipna):
            return 'string'

    elif isinstance(val, unicode):
        if is_unicode_array(values, skipna=skipna):
            return 'unicode'

    elif isinstance(val, bytes):
        if is_bytes_array(values, skipna=skipna):
            return 'bytes'

    elif util.is_period_object(val):
        if is_period_array(values):
            return 'period'

    elif is_interval(val):
        if is_interval_array(values):
            return 'interval'

    for i in range(n):
        val = values[i]
        if (util.is_integer_object(val) and
                not util.is_timedelta64_object(val) and
                not util.is_datetime64_object(val)):
            return 'mixed-integer'

    return 'mixed'


def infer_datetimelike_array(arr: object) -> object:
    """
    infer if we have a datetime or timedelta array
    - date: we have *only* date and maybe strings, nulls
    - datetime: we have *only* datetimes and maybe strings, nulls
    - timedelta: we have *only* timedeltas and maybe strings, nulls
    - nat: we do not have *any* date, datetimes or timedeltas, but do have
      at least a NaT
    - mixed: other objects (strings, a mix of tz-aware and tz-naive, or
                            actual objects)

    Parameters
    ----------
    arr : object array

    Returns
    -------
    string: {datetime, timedelta, date, nat, mixed}
    """

    cdef:
        Py_ssize_t i, n = len(arr)
        bint seen_timedelta = 0, seen_date = 0, seen_datetime = 0
        bint seen_tz_aware = 0, seen_tz_naive = 0
        bint seen_nat = 0
        list objs = []
        object v

    for i in range(n):
        v = arr[i]
        if util.is_string_object(v):
            objs.append(v)

            if len(objs) == 3:
                break

        elif v is None or util.is_nan(v):
            # nan or None
            pass
        elif v is NaT:
            seen_nat = 1
        elif PyDateTime_Check(v):
            # datetime
            seen_datetime = 1

            # disambiguate between tz-naive and tz-aware
            if v.tzinfo is None:
                seen_tz_naive = 1
            else:
                seen_tz_aware = 1

            if seen_tz_naive and seen_tz_aware:
                return 'mixed'
        elif util.is_datetime64_object(v):
            # np.datetime64
            seen_datetime = 1
        elif PyDate_Check(v):
            seen_date = 1
        elif is_timedelta(v) or util.is_timedelta64_object(v):
            # timedelta, or timedelta64
            seen_timedelta = 1
        else:
            return 'mixed'

    if seen_date and not (seen_datetime or seen_timedelta):
        return 'date'
    elif seen_datetime and not seen_timedelta:
        return 'datetime'
    elif seen_timedelta and not seen_datetime:
        return 'timedelta'
    elif seen_nat:
        return 'nat'

    # short-circuit by trying to
    # actually convert these strings
    # this is for performance as we don't need to try
    # convert *every* string array
    if len(objs):
        try:
            array_to_datetime(objs, errors='raise')
            return 'datetime'
        except:
            pass

        # we are *not* going to infer from strings
        # for timedelta as too much ambiguity

    return 'mixed'


cdef inline bint is_timedelta(object o):
    return PyDelta_Check(o) or util.is_timedelta64_object(o)


cdef class Validator:

    cdef:
        Py_ssize_t n
        dtype dtype
        bint skipna

    def __cinit__(self, Py_ssize_t n, dtype dtype=np.dtype(np.object_),
                  bint skipna=False):
        self.n = n
        self.dtype = dtype
        self.skipna = skipna

    cdef bint validate(self, ndarray values) except -1:
        if not self.n:
            return False

        if self.is_array_typed():
            return True
        elif self.dtype.type_num == NPY_OBJECT:
            if self.skipna:
                return self._validate_skipna(values)
            else:
                return self._validate(values)
        else:
            return False

    @cython.wraparound(False)
    @cython.boundscheck(False)
    cdef bint _validate(self, ndarray values) except -1:
        cdef:
            Py_ssize_t i
            Py_ssize_t n = self.n

        for i in range(n):
            if not self.is_valid(values[i]):
                return False

        return self.finalize_validate()

    @cython.wraparound(False)
    @cython.boundscheck(False)
    cdef bint _validate_skipna(self, ndarray values) except -1:
        cdef:
            Py_ssize_t i
            Py_ssize_t n = self.n

        for i in range(n):
            if not self.is_valid_skipna(values[i]):
                return False

        return self.finalize_validate_skipna()

    cdef bint is_valid(self, object value) except -1:
        return self.is_value_typed(value)

    cdef bint is_valid_skipna(self, object value) except -1:
        return self.is_valid(value) or self.is_valid_null(value)

    cdef bint is_value_typed(self, object value) except -1:
        raise NotImplementedError(
            '{typ} child class must define is_value_typed'
            .format(typ=type(self).__name__))

    cdef bint is_valid_null(self, object value) except -1:
        return value is None or util.is_nan(value)

    cdef bint is_array_typed(self) except -1:
        return False

    cdef inline bint finalize_validate(self):
        return True

    cdef bint finalize_validate_skipna(self):
        # TODO(phillipc): Remove the existing validate methods and replace them
        # with the skipna versions upon full deprecation of skipna=False
        return True


cdef class BoolValidator(Validator):
    cdef inline bint is_value_typed(self, object value) except -1:
        return util.is_bool_object(value)

    cdef inline bint is_array_typed(self) except -1:
        return issubclass(self.dtype.type, np.bool_)


cpdef bint is_bool_array(ndarray values, bint skipna=False):
    cdef:
        BoolValidator validator = BoolValidator(len(values),
                                                values.dtype,
                                                skipna=skipna)
    return validator.validate(values)


cdef class IntegerValidator(Validator):
    cdef inline bint is_value_typed(self, object value) except -1:
        return util.is_integer_object(value)

    cdef inline bint is_array_typed(self) except -1:
        return issubclass(self.dtype.type, np.integer)


cpdef bint is_integer_array(ndarray values):
    cdef:
        IntegerValidator validator = IntegerValidator(len(values),
                                                      values.dtype)
    return validator.validate(values)


cdef class IntegerFloatValidator(Validator):
    cdef inline bint is_value_typed(self, object value) except -1:
        return util.is_integer_object(value) or util.is_float_object(value)

    cdef inline bint is_array_typed(self) except -1:
        return issubclass(self.dtype.type, np.integer)


cdef bint is_integer_float_array(ndarray values):
    cdef:
        IntegerFloatValidator validator = IntegerFloatValidator(len(values),
                                                                values.dtype)
    return validator.validate(values)


cdef class FloatValidator(Validator):
    cdef inline bint is_value_typed(self, object value) except -1:
        return util.is_float_object(value)

    cdef inline bint is_array_typed(self) except -1:
        return issubclass(self.dtype.type, np.floating)


cpdef bint is_float_array(ndarray values):
    cdef:
        FloatValidator validator = FloatValidator(len(values), values.dtype)
    return validator.validate(values)


cdef class StringValidator(Validator):
    cdef inline bint is_value_typed(self, object value) except -1:
        return isinstance(value, str)

    cdef inline bint is_array_typed(self) except -1:
        return issubclass(self.dtype.type, np.str_)


cpdef bint is_string_array(ndarray values, bint skipna=False):
    cdef:
        StringValidator validator = StringValidator(len(values),
                                                    values.dtype,
                                                    skipna=skipna)
    return validator.validate(values)


cdef class UnicodeValidator(Validator):
    cdef inline bint is_value_typed(self, object value) except -1:
        return isinstance(value, unicode)

    cdef inline bint is_array_typed(self) except -1:
        return issubclass(self.dtype.type, np.unicode_)


cdef bint is_unicode_array(ndarray values, bint skipna=False):
    cdef:
        UnicodeValidator validator = UnicodeValidator(len(values),
                                                      values.dtype,
                                                      skipna=skipna)
    return validator.validate(values)


cdef class BytesValidator(Validator):
    cdef inline bint is_value_typed(self, object value) except -1:
        return isinstance(value, bytes)

    cdef inline bint is_array_typed(self) except -1:
        return issubclass(self.dtype.type, np.bytes_)


cdef bint is_bytes_array(ndarray values, bint skipna=False):
    cdef:
        BytesValidator validator = BytesValidator(len(values), values.dtype,
                                                  skipna=skipna)
    return validator.validate(values)


cdef class TemporalValidator(Validator):
    cdef:
        Py_ssize_t generic_null_count

    def __cinit__(self, Py_ssize_t n, dtype dtype=np.dtype(np.object_),
                  bint skipna=False):
        self.n = n
        self.dtype = dtype
        self.skipna = skipna
        self.generic_null_count = 0

    cdef inline bint is_valid(self, object value) except -1:
        return self.is_value_typed(value) or self.is_valid_null(value)

    cdef bint is_valid_null(self, object value) except -1:
        raise NotImplementedError(
            '{typ} child class must define is_valid_null'
            .format(typ=type(self).__name__))

    cdef inline bint is_valid_skipna(self, object value) except -1:
        cdef:
            bint is_typed_null = self.is_valid_null(value)
            bint is_generic_null = value is None or util.is_nan(value)
        self.generic_null_count += is_typed_null and is_generic_null
        return self.is_value_typed(value) or is_typed_null or is_generic_null

    cdef inline bint finalize_validate_skipna(self):
        return self.generic_null_count != self.n


cdef class DatetimeValidator(TemporalValidator):
    cdef bint is_value_typed(self, object value) except -1:
        return PyDateTime_Check(value)

    cdef inline bint is_valid_null(self, object value) except -1:
        return is_null_datetime64(value)


cpdef bint is_datetime_array(ndarray values):
    cdef:
        DatetimeValidator validator = DatetimeValidator(len(values),
                                                        skipna=True)
    return validator.validate(values)


cdef class Datetime64Validator(DatetimeValidator):
    cdef inline bint is_value_typed(self, object value) except -1:
        return util.is_datetime64_object(value)


cpdef bint is_datetime64_array(ndarray values):
    cdef:
        Datetime64Validator validator = Datetime64Validator(len(values),
                                                            skipna=True)
    return validator.validate(values)


def is_datetime_with_singletz_array(values: ndarray) -> bint:
    """
    Check values have the same tzinfo attribute.
    Doesn't check values are datetime-like types.
    """
    cdef:
        Py_ssize_t i, j, n = len(values)
        object base_val, base_tz, val, tz

    if n == 0:
        return False

    for i in range(n):
        base_val = values[i]
        if base_val is not NaT:
            base_tz = get_timezone(getattr(base_val, 'tzinfo', None))

            for j in range(i, n):
                val = values[j]
                if val is not NaT:
                    tz = getattr(val, 'tzinfo', None)
                    if not tz_compare(base_tz, tz):
                        return False
            break

    return True


cdef class TimedeltaValidator(TemporalValidator):
    cdef bint is_value_typed(self, object value) except -1:
        return PyDelta_Check(value)

    cdef inline bint is_valid_null(self, object value) except -1:
        return is_null_timedelta64(value)


# TODO: Not used outside of tests; remove?
def is_timedelta_array(values: ndarray) -> bint:
    cdef:
        TimedeltaValidator validator = TimedeltaValidator(len(values),
                                                          skipna=True)
    return validator.validate(values)


cdef class Timedelta64Validator(TimedeltaValidator):
    cdef inline bint is_value_typed(self, object value) except -1:
        return util.is_timedelta64_object(value)


# TODO: Not used outside of tests; remove?
def is_timedelta64_array(values: ndarray) -> bint:
    cdef:
        Timedelta64Validator validator = Timedelta64Validator(len(values),
                                                              skipna=True)
    return validator.validate(values)


cdef class AnyTimedeltaValidator(TimedeltaValidator):
    cdef inline bint is_value_typed(self, object value) except -1:
        return is_timedelta(value)


cpdef bint is_timedelta_or_timedelta64_array(ndarray values):
    """ infer with timedeltas and/or nat/none """
    cdef:
        AnyTimedeltaValidator validator = AnyTimedeltaValidator(len(values),
                                                                skipna=True)
    return validator.validate(values)


cdef class DateValidator(Validator):
    cdef inline bint is_value_typed(self, object value) except -1:
        return PyDate_Check(value)


cpdef bint is_date_array(ndarray values, bint skipna=False):
    cdef:
        DateValidator validator = DateValidator(len(values), skipna=skipna)
    return validator.validate(values)


cdef class TimeValidator(Validator):
    cdef inline bint is_value_typed(self, object value) except -1:
        return PyTime_Check(value)


cpdef bint is_time_array(ndarray values, bint skipna=False):
    cdef:
        TimeValidator validator = TimeValidator(len(values), skipna=skipna)
    return validator.validate(values)


cdef class PeriodValidator(TemporalValidator):
    cdef inline bint is_value_typed(self, object value) except -1:
        return util.is_period_object(value)

    cdef inline bint is_valid_null(self, object value) except -1:
        return is_null_period(value)


cpdef bint is_period_array(ndarray values):
    cdef:
        PeriodValidator validator = PeriodValidator(len(values), skipna=True)
    return validator.validate(values)


cdef class IntervalValidator(Validator):
    cdef inline bint is_value_typed(self, object value) except -1:
        return is_interval(value)


cpdef bint is_interval_array(ndarray values):
    cdef:
        IntervalValidator validator = IntervalValidator(len(values),
                                                        skipna=True)
    return validator.validate(values)


@cython.boundscheck(False)
@cython.wraparound(False)
def maybe_convert_numeric(ndarray[object] values, set na_values,
                          bint convert_empty=True, bint coerce_numeric=False):
    """
    Convert object array to a numeric array if possible.

    Parameters
    ----------
    values : ndarray
        Array of object elements to convert.
    na_values : set
        Set of values that should be interpreted as NaN.
    convert_empty : bool, default True
        If an empty array-like object is encountered, whether to interpret
        that element as NaN or not. If set to False, a ValueError will be
        raised if such an element is encountered and 'coerce_numeric' is False.
    coerce_numeric : bool, default False
        If initial attempts to convert to numeric have failed, whether to
        force conversion to numeric via alternative methods or by setting the
        element to NaN. Otherwise, an Exception will be raised when such an
        element is encountered.

        This boolean also has an impact on how conversion behaves when a
        numeric array has no suitable numerical dtype to return (i.e. uint64,
        int32, uint8). If set to False, the original object array will be
        returned. Otherwise, a ValueError will be raised.

    Returns
    -------
    numeric_array : array of converted object values to numerical ones
    """

    if len(values) == 0:
        return np.array([], dtype='i8')

    # fastpath for ints - try to convert all based on first value
    cdef:
        object val = values[0]

    if util.is_integer_object(val):
        try:
            maybe_ints = values.astype('i8')
            if (maybe_ints == values).all():
                return maybe_ints
        except (ValueError, OverflowError, TypeError):
            pass

    # otherwise, iterate and do full infererence
    cdef:
        int status, maybe_int
        Py_ssize_t i, n = values.size
        Seen seen = Seen(coerce_numeric)
        ndarray[float64_t] floats = np.empty(n, dtype='f8')
        ndarray[complex128_t] complexes = np.empty(n, dtype='c16')
        ndarray[int64_t] ints = np.empty(n, dtype='i8')
        ndarray[uint64_t] uints = np.empty(n, dtype='u8')
        ndarray[uint8_t] bools = np.empty(n, dtype='u1')
        float64_t fval

    for i in range(n):
        val = values[i]

        if val.__hash__ is not None and val in na_values:
            seen.saw_null()
            floats[i] = complexes[i] = nan
        elif util.is_float_object(val):
            fval = val
            if fval != fval:
                seen.null_ = True

            floats[i] = complexes[i] = fval
            seen.float_ = True
        elif util.is_integer_object(val):
            floats[i] = complexes[i] = val

            val = int(val)
            seen.saw_int(val)

            if val >= 0:
                if val <= oUINT64_MAX:
                    uints[i] = val
                else:
                    seen.float_ = True

            if val <= oINT64_MAX:
                ints[i] = val

            if seen.sint_ and seen.uint_:
                seen.float_ = True

        elif util.is_bool_object(val):
            floats[i] = uints[i] = ints[i] = bools[i] = val
            seen.bool_ = True
        elif val is None:
            seen.saw_null()
            floats[i] = complexes[i] = nan
        elif hasattr(val, '__len__') and len(val) == 0:
            if convert_empty or seen.coerce_numeric:
                seen.saw_null()
                floats[i] = complexes[i] = nan
            else:
                raise ValueError('Empty string encountered')
        elif util.is_complex_object(val):
            complexes[i] = val
            seen.complex_ = True
        elif is_decimal(val):
            floats[i] = complexes[i] = val
            seen.float_ = True
        else:
            try:
                status = floatify(val, &fval, &maybe_int)

                if fval in na_values:
                    seen.saw_null()
                    floats[i] = complexes[i] = nan
                else:
                    if fval != fval:
                        seen.null_ = True

                    floats[i] = fval

                if maybe_int:
                    as_int = int(val)

                    if as_int in na_values:
                        seen.saw_null()
                    else:
                        seen.saw_int(as_int)

                    if not (seen.float_ or as_int in na_values):
                        if as_int < oINT64_MIN or as_int > oUINT64_MAX:
                            raise ValueError('Integer out of range.')

                        if as_int >= 0:
                            uints[i] = as_int
                        if as_int <= oINT64_MAX:
                            ints[i] = as_int

                    seen.float_ = seen.float_ or (seen.uint_ and seen.sint_)
                else:
                    seen.float_ = True
            except (TypeError, ValueError) as e:
                if not seen.coerce_numeric:
                    raise type(e)(str(e) + ' at position {pos}'.format(pos=i))
                elif "uint64" in str(e):  # Exception from check functions.
                    raise
                seen.saw_null()
                floats[i] = nan

    if seen.check_uint64_conflict():
        return values

    if seen.complex_:
        return complexes
    elif seen.float_:
        return floats
    elif seen.int_:
        if seen.uint_:
            return uints
        else:
            return ints
    elif seen.bool_:
        return bools.view(np.bool_)
    elif seen.uint_:
        return uints
    return ints


@cython.boundscheck(False)
@cython.wraparound(False)
def maybe_convert_objects(ndarray[object] objects, bint try_float=0,
                          bint safe=0, bint convert_datetime=0,
                          bint convert_timedelta=0):
    """
    Type inference function-- convert object array to proper dtype
    """
    cdef:
        Py_ssize_t i, n
        ndarray[float64_t] floats
        ndarray[complex128_t] complexes
        ndarray[int64_t] ints
        ndarray[uint64_t] uints
        ndarray[uint8_t] bools
        ndarray[int64_t] idatetimes
        ndarray[int64_t] itimedeltas
        Seen seen = Seen()
        object val
        float64_t fval, fnan

    n = len(objects)

    floats = np.empty(n, dtype='f8')
    complexes = np.empty(n, dtype='c16')
    ints = np.empty(n, dtype='i8')
    uints = np.empty(n, dtype='u8')
    bools = np.empty(n, dtype=np.uint8)

    if convert_datetime:
        datetimes = np.empty(n, dtype='M8[ns]')
        idatetimes = datetimes.view(np.int64)

    if convert_timedelta:
        timedeltas = np.empty(n, dtype='m8[ns]')
        itimedeltas = timedeltas.view(np.int64)

    fnan = np.nan

    for i in range(n):
        val = objects[i]

        if val is None:
            seen.null_ = 1
            floats[i] = complexes[i] = fnan
        elif val is NaT:
            if convert_datetime:
                idatetimes[i] = iNaT
                seen.datetime_ = 1
            if convert_timedelta:
                itimedeltas[i] = iNaT
                seen.timedelta_ = 1
            if not (convert_datetime or convert_timedelta):
                seen.object_ = 1
        elif util.is_bool_object(val):
            seen.bool_ = 1
            bools[i] = val
        elif util.is_float_object(val):
            floats[i] = complexes[i] = val
            seen.float_ = 1
        elif util.is_datetime64_object(val):
            if convert_datetime:
                idatetimes[i] = convert_to_tsobject(
                    val, None, None, 0, 0).value
                seen.datetime_ = 1
            else:
                seen.object_ = 1
                break
        elif is_timedelta(val):
            if convert_timedelta:
                itimedeltas[i] = convert_to_timedelta64(val, 'ns')
                seen.timedelta_ = 1
            else:
                seen.object_ = 1
                break
        elif util.is_integer_object(val):
            seen.int_ = 1
            floats[i] = <float64_t> val
            complexes[i] = <double complex> val
            if not seen.null_:
                seen.saw_int(int(val))

                if ((seen.uint_ and seen.sint_) or
                        val > oUINT64_MAX or val < oINT64_MIN):
                    seen.object_ = 1
                    break

                if seen.uint_:
                    uints[i] = val
                elif seen.sint_:
                    ints[i] = val
                else:
                    uints[i] = val
                    ints[i] = val

        elif util.is_complex_object(val):
            complexes[i] = val
            seen.complex_ = 1
        elif PyDateTime_Check(val) or util.is_datetime64_object(val):

            # if we have an tz's attached then return the objects
            if convert_datetime:
                if getattr(val, 'tzinfo', None) is not None:
                    seen.datetimetz_ = 1
                    break
                else:
                    seen.datetime_ = 1
                    idatetimes[i] = convert_to_tsobject(
                        val, None, None, 0, 0).value
            else:
                seen.object_ = 1
                break
        elif try_float and not util.is_string_object(val):
            # this will convert Decimal objects
            try:
                floats[i] = float(val)
                complexes[i] = complex(val)
                seen.float_ = 1
            except Exception:
                seen.object_ = 1
                break
        else:
            seen.object_ = 1
            break

    # we try to coerce datetime w/tz but must all have the same tz
    if seen.datetimetz_:
        if len({getattr(val, 'tzinfo', None) for val in objects}) == 1:
            from pandas import DatetimeIndex
            return DatetimeIndex(objects)
        seen.object_ = 1

    if not seen.object_:
        if not safe:
            if seen.null_:
                if seen.is_float_or_complex:
                    if seen.complex_:
                        return complexes
                    elif seen.float_ or seen.int_:
                        return floats
            else:
                if not seen.bool_:
                    if seen.datetime_:
                        if not seen.numeric_:
                            return datetimes
                    elif seen.timedelta_:
                        if not seen.numeric_:
                            return timedeltas
                    else:
                        if seen.complex_:
                            return complexes
                        elif seen.float_:
                            return floats
                        elif seen.int_:
                            if seen.uint_:
                                return uints
                            else:
                                return ints
                elif seen.is_bool:
                    return bools.view(np.bool_)

        else:
            # don't cast int to float, etc.
            if seen.null_:
                if seen.is_float_or_complex:
                    if seen.complex_:
                        if not seen.int_:
                            return complexes
                    elif seen.float_:
                        if not seen.int_:
                            return floats
            else:
                if not seen.bool_:
                    if seen.datetime_:
                        if not seen.numeric_:
                            return datetimes
                    elif seen.timedelta_:
                        if not seen.numeric_:
                            return timedeltas
                    else:
                        if seen.complex_:
                            if not seen.int_:
                                return complexes
                        elif seen.float_:
                            if not seen.int_:
                                return floats
                        elif seen.int_:
                            if seen.uint_:
                                return uints
                            else:
                                return ints
                elif seen.is_bool:
                    return bools.view(np.bool_)

    return objects


def map_infer_mask(ndarray arr, object f, ndarray[uint8_t] mask,
                   bint convert=1):
    """
    Substitute for np.vectorize with pandas-friendly dtype inference

    Parameters
    ----------
    arr : ndarray
    f : function

    Returns
    -------
    mapped : ndarray
    """
    cdef:
        Py_ssize_t i, n
        ndarray[object] result
        object val

    n = len(arr)
    result = np.empty(n, dtype=object)
    for i in range(n):
        if mask[i]:
            val = util.get_value_at(arr, i)
        else:
            val = f(util.get_value_at(arr, i))

            if util.is_array(val) and PyArray_NDIM(val) == 0:
                # unbox 0-dim arrays, GH#690
                # TODO: is there a faster way to unbox?
                #   item_from_zerodim?
                val = val.item()

        result[i] = val

    if convert:
        return maybe_convert_objects(result,
                                     try_float=0,
                                     convert_datetime=0,
                                     convert_timedelta=0)

    return result


def map_infer(ndarray arr, object f, bint convert=1):
    """
    Substitute for np.vectorize with pandas-friendly dtype inference

    Parameters
    ----------
    arr : ndarray
    f : function

    Returns
    -------
    mapped : ndarray
    """
    cdef:
        Py_ssize_t i, n
        ndarray[object] result
        object val

    n = len(arr)
    result = np.empty(n, dtype=object)
    for i in range(n):
        val = f(util.get_value_at(arr, i))

        if util.is_array(val) and PyArray_NDIM(val) == 0:
            # unbox 0-dim arrays, GH#690
            # TODO: is there a faster way to unbox?
            #   item_from_zerodim?
            val = val.item()

        result[i] = val

    if convert:
        return maybe_convert_objects(result,
                                     try_float=0,
                                     convert_datetime=0,
                                     convert_timedelta=0)

    return result


def to_object_array(list rows, int min_width=0):
    """
    Convert a list of lists into an object array.

    Parameters
    ----------
    rows : 2-d array (N, K)
        A list of lists to be converted into an array
    min_width : int
        The minimum width of the object array. If a list
        in `rows` contains fewer than `width` elements,
        the remaining elements in the corresponding row
        will all be `NaN`.

    Returns
    -------
    obj_array : numpy array of the object dtype
    """
    cdef:
        Py_ssize_t i, j, n, k, tmp
        ndarray[object, ndim=2] result
        list row

    n = len(rows)

    k = min_width
    for i in range(n):
        tmp = len(rows[i])
        if tmp > k:
            k = tmp

    result = np.empty((n, k), dtype=object)

    for i in range(n):
        row = rows[i]

        for j in range(len(row)):
            result[i, j] = row[j]

    return result


def tuples_to_object_array(ndarray[object] tuples):
    cdef:
        Py_ssize_t i, j, n, k, tmp
        ndarray[object, ndim=2] result
        tuple tup

    n = len(tuples)
    k = len(tuples[0])
    result = np.empty((n, k), dtype=object)
    for i in range(n):
        tup = tuples[i]
        for j in range(k):
            result[i, j] = tup[j]

    return result


def to_object_array_tuples(list rows):
    cdef:
        Py_ssize_t i, j, n, k, tmp
        ndarray[object, ndim=2] result
        tuple row

    n = len(rows)

    k = 0
    for i in range(n):
        tmp = len(rows[i])
        if tmp > k:
            k = tmp

    result = np.empty((n, k), dtype=object)

    try:
        for i in range(n):
            row = rows[i]
            for j in range(len(row)):
                result[i, j] = row[j]
    except Exception:
        # upcast any subclasses to tuple
        for i in range(n):
            row = tuple(rows[i])
            for j in range(len(row)):
                result[i, j] = row[j]

    return result


def fast_multiget(dict mapping, ndarray keys, default=np.nan):
    cdef:
        Py_ssize_t i, n = len(keys)
        object val
        ndarray[object] output = np.empty(n, dtype='O')

    if n == 0:
        # kludge, for Series
        return np.empty(0, dtype='f8')

    keys = getattr(keys, 'values', keys)

    for i in range(n):
        val = keys[i]
        if val in mapping:
            output[i] = mapping[val]
        else:
            output[i] = default

    return maybe_convert_objects(output)
