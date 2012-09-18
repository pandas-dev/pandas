from numpy cimport ndarray

from numpy cimport float64_t, int32_t, int64_t, uint8_t
cimport cython

cimport numpy as cnp

cnp.import_array()
cnp.import_ufunc()

cimport util

import numpy as np

import _algos

# include "hashtable.pyx"

cdef extern from "datetime.h":
    bint PyDateTime_Check(object o)
    void PyDateTime_IMPORT()

PyDateTime_IMPORT

cdef extern from "Python.h":
    int PySlice_Check(object)

#     int PyList_Check(object)
#     int PyTuple_Check(object)

cdef inline is_definitely_invalid_key(object val):
    if PyTuple_Check(val):
        try:
            hash(val)
        except TypeError:
            return True

    return (PySlice_Check(val) or cnp.PyArray_Check(val)
            or PyList_Check(val))

def get_value_at(ndarray arr, object loc):
    return util.get_value_at(arr, loc)

def set_value_at(ndarray arr, object loc, object val):
    return util.set_value_at(arr, loc, val)


# Don't populate hash tables in monotonic indexes larger than this
_SIZE_CUTOFF = 1000000


cdef class IndexEngine:

    cdef readonly:
        object vgetter
        HashTable mapping
        bint over_size_threshold

    cdef:
        bint unique, monotonic
        bint initialized, monotonic_check, unique_check

    def __init__(self, vgetter, n):
        self.vgetter = vgetter

        self.over_size_threshold = n >= _SIZE_CUTOFF

        self.initialized = 0
        self.monotonic_check = 0

        self.unique = 0
        self.monotonic = 0

    def __contains__(self, object val):
        self._ensure_mapping_populated()
        hash(val)
        return val in self.mapping

    cpdef get_value(self, ndarray arr, object key):
        '''
        arr : 1-dimensional ndarray
        '''
        cdef:
            object loc
            void* data_ptr

        loc = self.get_loc(key)
        if PySlice_Check(loc) or cnp.PyArray_Check(loc):
            return arr[loc]
        else:
            if arr.descr.type_num == NPY_DATETIME:
                return Timestamp(util.get_value_at(arr, loc))
            return util.get_value_at(arr, loc)

    cpdef set_value(self, ndarray arr, object key, object value):
        '''
        arr : 1-dimensional ndarray
        '''
        cdef:
            object loc
            void* data_ptr

        loc = self.get_loc(key)
        value = convert_scalar(arr, value)

        if PySlice_Check(loc) or cnp.PyArray_Check(loc):
            arr[loc] = value
        else:
            util.set_value_at(arr, loc, value)

    cpdef get_loc(self, object val):
        if is_definitely_invalid_key(val):
            raise TypeError

        if self.over_size_threshold and self.is_monotonic:
            if not self.is_unique:
                return self._get_loc_duplicates(val)
            values = self._get_index_values()
            loc = _bin_search(values, val) # .searchsorted(val, side='left')
            if util.get_value_at(values, loc) != val:
                raise KeyError(val)
            return loc

        self._ensure_mapping_populated()
        if not self.unique:
            return self._get_loc_duplicates(val)

        try:
            return self.mapping.get_item(val)
        except TypeError:
            self._check_type(val)
            raise KeyError(val)

    cdef inline _get_loc_duplicates(self, object val):
        cdef:
            Py_ssize_t diff

        if self.is_monotonic:
            values = self._get_index_values()
            left = values.searchsorted(val, side='left')
            right = values.searchsorted(val, side='right')

            diff = right - left
            if diff == 0:
                raise KeyError(val)
            elif diff == 1:
                return left
            else:
                return slice(left, right)
        else:
            return self._maybe_get_bool_indexer(val)

    cdef _maybe_get_bool_indexer(self, object val):
        cdef:
            ndarray[uint8_t] indexer
            ndarray[object] values
            int count = 0
            Py_ssize_t i, n
            int last_true

        values = self._get_index_values()
        n = len(values)

        result = np.empty(n, dtype=bool)
        indexer = result.view(np.uint8)

        for i in range(n):
            if values[i] == val:
                count += 1
                indexer[i] = 1
                last_true = i
            else:
                indexer[i] = 0

        if count == 0:
            raise KeyError(val)
        if count == 1:
            return last_true

        return result

    property is_unique:

        def __get__(self):
            if not self.unique_check:
                self._do_unique_check()

            return self.unique == 1

    property is_monotonic:

        def __get__(self):
            if not self.monotonic_check:
                self._do_monotonic_check()

            return self.monotonic == 1

    cdef inline _do_monotonic_check(self):
        try:
            values = self._get_index_values()
            self.monotonic, unique = self._call_monotonic(values)

            if unique is not None:
                self.unique = unique
                self.unique_check = 1

        except TypeError:
            self.monotonic = 0
        self.monotonic_check = 1

    cdef _get_index_values(self):
        return self.vgetter()

    cdef inline _do_unique_check(self):
        self._ensure_mapping_populated()

    def _call_monotonic(self, values):
        raise NotImplementedError

    cdef _make_hash_table(self, n):
        raise NotImplementedError

    cdef inline _check_type(self, object val):
        hash(val)

    cdef inline _ensure_mapping_populated(self):
        if not self.initialized:
            self.initialize()

    cdef initialize(self):
        values = self._get_index_values()

        self.mapping = self._make_hash_table(len(values))
        self.mapping.map_locations(values)

        if len(self.mapping) == len(values):
            self.unique = 1
            self.unique_check = 1

        self.initialized = 1

    def clear_mapping(self):
        self.mapping = None
        self.initialized = 0

    def get_indexer(self, values):
        self._ensure_mapping_populated()
        return self.mapping.lookup(values)



# @cache_readonly
# def _monotonicity_check(self):
#     try:
#         f = self._algos['is_monotonic']
#         # wrong buffer type raises ValueError
#         return f(self.values)
#     except TypeError:
#         return False, None



cdef class Int64Engine(IndexEngine):

    # cdef Int64HashTable mapping

    cdef _make_hash_table(self, n):
        return Int64HashTable(n)

    def _call_monotonic(self, values):
        return _algos.is_monotonic_int64(values)

    def get_pad_indexer(self, other, limit=None):
        return _algos.pad_int64(self._get_index_values(), other,
                                  limit=limit)

    def get_backfill_indexer(self, other, limit=None):
        return _algos.backfill_int64(self._get_index_values(), other,
                                       limit=limit)

    cdef _maybe_get_bool_indexer(self, object val):
        cdef:
            ndarray[uint8_t, cast=True] indexer
            ndarray[int64_t] values
            int count = 0
            Py_ssize_t i, n
            int64_t ival
            int last_true

        if not util.is_integer_object(val):
            raise KeyError(val)

        ival = val

        values = self._get_index_values()
        n = len(values)

        result = np.empty(n, dtype=bool)
        indexer = result.view(np.uint8)

        for i in range(n):
            if values[i] == val:
                count += 1
                indexer[i] = 1
                last_true = i
            else:
                indexer[i] = 0

        if count == 0:
            raise KeyError(val)
        if count == 1:
            return last_true

        return result

cdef class Float64Engine(IndexEngine):

    # cdef Float64HashTable mapping

    cdef _make_hash_table(self, n):
        return Float64HashTable(n)

    def _call_monotonic(self, values):
        return _algos.is_monotonic_float64(values)

    def get_pad_indexer(self, other, limit=None):
        return _algos.pad_float64(self._get_index_values(), other,
                                    limit=limit)

    def get_backfill_indexer(self, other, limit=None):
        return _algos.backfill_float64(self._get_index_values(), other,
                                         limit=limit)


cdef Py_ssize_t _bin_search(ndarray values, object val):
    cdef:
        Py_ssize_t mid, lo = 0, hi = len(values) - 1
        object pval

    if hi >= 0 and val > util.get_value_at(values, hi):
        return len(values)

    while lo < hi:
        mid = (lo + hi) // 2
        pval = util.get_value_at(values, mid)
        if val < pval:
            hi = mid
        elif val > pval:
            lo = mid + 1
        else:
            while mid > 0 and val == util.get_value_at(values, mid - 1):
                mid -= 1
            return mid

    if val <= util.get_value_at(values, mid):
        return mid
    else:
        return mid + 1

_pad_functions = {
    'object' : _algos.pad_object,
    'int64' : _algos.pad_int64,
    'float64' : _algos.pad_float64
}

_backfill_functions = {
    'object': _algos.backfill_object,
    'int64': _algos.backfill_int64,
    'float64': _algos.backfill_float64
}

cdef class ObjectEngine(IndexEngine):

    # cdef PyObjectHashTable mapping

    cdef _make_hash_table(self, n):
        return PyObjectHashTable(n)

    def _call_monotonic(self, values):
        return _algos.is_monotonic_object(values)

    def get_pad_indexer(self, other, limit=None):
        return _algos.pad_object(self._get_index_values(), other,
                                   limit=limit)

    def get_backfill_indexer(self, other, limit=None):
        return _algos.backfill_object(self._get_index_values(), other,
                                        limit=limit)


cdef class DatetimeEngine(Int64Engine):

    def __contains__(self, object val):
        if self.over_size_threshold and self.is_monotonic:
            if not self.is_unique:
                return self._get_loc_duplicates(val)
            values = self._get_index_values()
            conv = _to_i8(val)
            loc = values.searchsorted(conv, side='left')
            return util.get_value_at(values, loc) == conv

        self._ensure_mapping_populated()
        return _to_i8(val) in self.mapping

    cdef _get_index_values(self):
        return self.vgetter().view('i8')

    def _call_monotonic(self, values):
        return _algos.is_monotonic_int64(values)

    cpdef get_loc(self, object val):
        if is_definitely_invalid_key(val):
            raise TypeError

        # Welcome to the spaghetti factory

        if self.over_size_threshold and self.is_monotonic:
            if not self.is_unique:
                val = _to_i8(val)
                return self._get_loc_duplicates(val)
            values = self._get_index_values()
            conv = _to_i8(val)
            loc = values.searchsorted(conv, side='left')
            if loc == len(values) or util.get_value_at(values, loc) != conv:
                raise KeyError(val)
            return loc

        self._ensure_mapping_populated()
        if not self.unique:
            val = _to_i8(val)
            return self._get_loc_duplicates(val)

        try:
            return self.mapping.get_item(val.value)
        except KeyError:
            raise KeyError(val)
        except AttributeError:
            pass

        try:
            val = _to_i8(val)
            return self.mapping.get_item(val)
        except TypeError:
            self._date_check_type(val)
            raise KeyError(val)

    cdef inline _date_check_type(self, object val):
        hash(val)
        if not util.is_integer_object(val):
            raise KeyError(val)

    def get_indexer(self, values):
        self._ensure_mapping_populated()
        if values.dtype != 'M8[ns]':
            return np.repeat(-1, len(values)).astype('i4')
        values = np.asarray(values).view('i8')
        return self.mapping.lookup(values)

    def get_pad_indexer(self, other, limit=None):
        if other.dtype != 'M8[ns]':
            return np.repeat(-1, len(other)).astype('i4')
        other = np.asarray(other).view('i8')
        return _algos.pad_int64(self._get_index_values(), other,
                                limit=limit)

    def get_backfill_indexer(self, other, limit=None):
        if other.dtype != 'M8[ns]':
            return np.repeat(-1, len(other)).astype('i4')
        other = np.asarray(other).view('i8')
        return _algos.backfill_int64(self._get_index_values(), other,
                                     limit=limit)


cpdef convert_scalar(ndarray arr, object value):
    if arr.descr.type_num == NPY_DATETIME:
        if isinstance(value, _Timestamp):
            return (<_Timestamp> value).value
        elif value is None or value != value:
            return iNaT
        else:
            return Timestamp(value).value

    if issubclass(arr.dtype.type, (np.integer, np.bool_)):
        if util.is_float_object(value) and value != value:
            raise ValueError('Cannot assign nan to integer series')

    return value

cdef inline _to_i8(object val):
    cdef pandas_datetimestruct dts
    try:
        return val.value
    except AttributeError:
        if util.is_datetime64_object(val):
            return get_datetime64_value(val)
        elif PyDateTime_Check(val):
            return _pydatetime_to_dts(val, &dts)
        return val


# ctypedef fused idxvalue_t:
#     object
#     int
#     float64_t
#     int32_t
#     int64_t

# @cython.boundscheck(False)
# @cython.wraparound(False)
# def is_monotonic(ndarray[idxvalue_t] arr):
#     '''
#     Returns
#     -------
#     is_monotonic, is_unique
#     '''
#     cdef:
#         Py_ssize_t i, n
#         idxvalue_t prev, cur
#         bint is_unique = 1

#     n = len(arr)

#     if n < 2:
#         return True, True

#     prev = arr[0]
#     for i in range(1, n):
#         cur = arr[i]
#         if cur < prev:
#             return False, None
#         elif cur == prev:
#             is_unique = 0
#         prev = cur
#     return True, is_unique


# @cython.wraparound(False)
# @cython.boundscheck(False)
# def groupby_index(ndarray[idxvalue_t] index, ndarray labels):
#     cdef dict result = {}
#     cdef Py_ssize_t i, length
#     cdef list members
#     cdef object idx, key

#     length = len(index)

#     for i in range(length):
#         key = util.get_value_1d(labels, i)

#         if util._checknull(key):
#             continue

#         idx = index[i]
#         if key in result:
#             members = result[key]
#             members.append(idx)
#         else:
#             result[key] = [idx]

#     return result
