# cython: profile=False

from numpy cimport ndarray

from numpy cimport (float64_t, int32_t, int64_t, uint8_t,
                    NPY_DATETIME, NPY_TIMEDELTA)
cimport cython

cimport numpy as cnp

cnp.import_array()
cnp.import_ufunc()

cimport util

import numpy as np

cimport tslib
from hashtable cimport *
from pandas import algos, tslib, hashtable as _hash
from pandas.tslib import Timestamp, Timedelta

from datetime cimport (get_datetime64_value, _pydatetime_to_dts,
                       pandas_datetimestruct)

from cpython cimport PyTuple_Check, PyList_Check

cdef extern from "datetime.h":
    bint PyDateTime_Check(object o)
    void PyDateTime_IMPORT()

cdef int64_t iNaT = util.get_nat()

try:
    from dateutil.tz import tzutc as _du_utc
    import pytz
    UTC = pytz.utc
    have_pytz = True
except ImportError:
    have_pytz = False

PyDateTime_IMPORT

cdef extern from "Python.h":
    int PySlice_Check(object)


cdef inline is_definitely_invalid_key(object val):
    if PyTuple_Check(val):
        try:
            hash(val)
        except TypeError:
            return True

    # we have a _data, means we are a NDFrame
    return (PySlice_Check(val) or cnp.PyArray_Check(val)
            or PyList_Check(val) or hasattr(val, '_data'))


def get_value_at(ndarray arr, object loc):
    if arr.descr.type_num == NPY_DATETIME:
        return Timestamp(util.get_value_at(arr, loc))
    elif arr.descr.type_num == NPY_TIMEDELTA:
        return Timedelta(util.get_value_at(arr, loc))
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
        bint unique, monotonic_inc, monotonic_dec
        bint need_monotonic_check, need_unique_check

    def __init__(self, vgetter, n):
        self.vgetter = vgetter

        self.over_size_threshold = n >= _SIZE_CUTOFF
        self.clear_mapping()

    def __contains__(self, object val):
        self._ensure_mapping_populated()
        hash(val)
        return val in self.mapping

    cpdef get_value(self, ndarray arr, object key, object tz=None):
        """
        arr : 1-dimensional ndarray
        """
        cdef:
            object loc
            void* data_ptr

        loc = self.get_loc(key)
        if PySlice_Check(loc) or cnp.PyArray_Check(loc):
            return arr[loc]
        else:
            if arr.descr.type_num == NPY_DATETIME:
                return Timestamp(util.get_value_at(arr, loc), tz=tz)
            elif arr.descr.type_num == NPY_TIMEDELTA:
                return Timedelta(util.get_value_at(arr, loc))
            return util.get_value_at(arr, loc)

    cpdef set_value(self, ndarray arr, object key, object value):
        """
        arr : 1-dimensional ndarray
        """
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
            raise TypeError("'{val}' is an invalid key".format(val=val))

        if self.over_size_threshold and self.is_monotonic_increasing:
            if not self.is_unique:
                return self._get_loc_duplicates(val)
            values = self._get_index_values()
            loc = _bin_search(values, val) # .searchsorted(val, side='left')
            if loc >= len(values):
                raise KeyError(val)
            if util.get_value_at(values, loc) != val:
                raise KeyError(val)
            return loc

        self._ensure_mapping_populated()
        if not self.unique:
            return self._get_loc_duplicates(val)

        self._check_type(val)

        try:
            return self.mapping.get_item(val)
        except TypeError:
            raise KeyError(val)

    cdef inline _get_loc_duplicates(self, object val):
        cdef:
            Py_ssize_t diff

        if self.is_monotonic_increasing:
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
            if self.need_unique_check:
                self._do_unique_check()

            return self.unique == 1

    cdef inline _do_unique_check(self):

        # this de-facto the same
        self._ensure_mapping_populated()

    property is_monotonic_increasing:

        def __get__(self):
            if self.need_monotonic_check:
                self._do_monotonic_check()

            return self.monotonic_inc == 1

    property is_monotonic_decreasing:

        def __get__(self):
            if self.need_monotonic_check:
                self._do_monotonic_check()

            return self.monotonic_dec == 1

    cdef inline _do_monotonic_check(self):
        cdef object is_unique
        try:
            values = self._get_index_values()
            self.monotonic_inc, self.monotonic_dec, is_unique = \
                self._call_monotonic(values)
        except TypeError:
            self.monotonic_inc = 0
            self.monotonic_dec = 0
            is_unique = 0

        self.need_monotonic_check = 0

        # we can only be sure of uniqueness if is_unique=1
        if is_unique:
            self.unique = 1
            self.need_unique_check = 0

    cdef _get_index_values(self):
        return self.vgetter()

    def _call_monotonic(self, values):
        raise NotImplementedError

    cdef _make_hash_table(self, n):
        raise NotImplementedError

    cdef _check_type(self, object val):
        hash(val)

    property is_mapping_populated:

        def __get__(self):
            return self.mapping is not None

    cdef inline _ensure_mapping_populated(self):
        # this populates the mapping
        # if its not already populated
        # also satisfies the need_unique_check

        if not self.is_mapping_populated:

            values = self._get_index_values()

            self.mapping = self._make_hash_table(len(values))
            self.mapping.map_locations(values)

            if len(self.mapping) == len(values):
                self.unique = 1

        self.need_unique_check = 0

    def clear_mapping(self):
        self.mapping = None
        self.need_monotonic_check = 1
        self.need_unique_check = 1

        self.unique = 0
        self.monotonic_inc = 0
        self.monotonic_dec = 0

    def get_indexer(self, values):
        self._ensure_mapping_populated()
        return self.mapping.lookup(values)

    def get_indexer_non_unique(self, targets):
        """ return an indexer suitable for takng from a non unique index
            return the labels in the same order ast the target
            and a missing indexer into the targets (which correspond
            to the -1 indicies in the results """

        cdef:
            ndarray values, x
            ndarray[int64_t] result, missing
            set stargets
            dict d = {}
            object val
            int count = 0, count_missing = 0
            Py_ssize_t i, j, n, n_t, n_alloc

        self._ensure_mapping_populated()
        values = self._get_index_values()
        stargets = set(targets)
        n = len(values)
        n_t = len(targets)
        if n > 10000:
            n_alloc = 10000
        else:
            n_alloc = n

        result = np.empty(n_alloc, dtype=np.int64)
        missing = np.empty(n_t, dtype=np.int64)

        # form the set of the results (like ismember)
        members = np.empty(n, dtype=np.uint8)
        for i in range(n):
            val = util.get_value_1d(values, i)
            if val in stargets:
                if val not in d:
                    d[val] = []
                d[val].append(i)

        for i in range(n_t):

            val = util.get_value_1d(targets, i)

            # found
            if val in d:
                for j in d[val]:

                    # realloc if needed
                    if count >= n_alloc:
                        n_alloc += 10000
                        result = np.resize(result, n_alloc)

                    result[count] = j
                    count += 1

            # value not found
            else:

                if count >= n_alloc:
                    n_alloc += 10000
                    result = np.resize(result, n_alloc)
                result[count] = -1
                count += 1
                missing[count_missing] = i
                count_missing += 1

        return result[0:count], missing[0:count_missing]

cdef class Int64Engine(IndexEngine):

    cdef _get_index_values(self):
        return algos.ensure_int64(self.vgetter())

    cdef _make_hash_table(self, n):
        return _hash.Int64HashTable(n)

    def _call_monotonic(self, values):
        return algos.is_monotonic_int64(values, timelike=False)

    def get_pad_indexer(self, other, limit=None):
        return algos.pad_int64(self._get_index_values(), other,
                               limit=limit)

    def get_backfill_indexer(self, other, limit=None):
        return algos.backfill_int64(self._get_index_values(), other,
                                    limit=limit)

    cdef _check_type(self, object val):
        hash(val)
        if util.is_bool_object(val):
            raise KeyError(val)
        elif util.is_float_object(val):
            raise KeyError(val)

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

    cdef _make_hash_table(self, n):
        return _hash.Float64HashTable(n)

    cdef _get_index_values(self):
        return algos.ensure_float64(self.vgetter())

    cdef _maybe_get_bool_indexer(self, object val):
        cdef:
            ndarray[uint8_t] indexer
            ndarray[float64_t] values
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

    def _call_monotonic(self, values):
        return algos.is_monotonic_float64(values, timelike=False)

    def get_pad_indexer(self, other, limit=None):
        return algos.pad_float64(self._get_index_values(), other,
                                    limit=limit)

    def get_backfill_indexer(self, other, limit=None):
        return algos.backfill_float64(self._get_index_values(), other,
                                         limit=limit)


cdef Py_ssize_t _bin_search(ndarray values, object val) except -1:
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
    'object': algos.pad_object,
    'int64': algos.pad_int64,
    'float64': algos.pad_float64
}

_backfill_functions = {
    'object': algos.backfill_object,
    'int64': algos.backfill_int64,
    'float64': algos.backfill_float64
}

cdef class ObjectEngine(IndexEngine):

    cdef _make_hash_table(self, n):
        return _hash.PyObjectHashTable(n)

    def _call_monotonic(self, values):
        return algos.is_monotonic_object(values, timelike=False)

    def get_pad_indexer(self, other, limit=None):
        return algos.pad_object(self._get_index_values(), other,
                                   limit=limit)

    def get_backfill_indexer(self, other, limit=None):
        return algos.backfill_object(self._get_index_values(), other,
                                        limit=limit)


cdef class DatetimeEngine(Int64Engine):

    cdef _get_box_dtype(self):
        return 'M8[ns]'

    def __contains__(self, object val):
        if self.over_size_threshold and self.is_monotonic_increasing:
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
        return algos.is_monotonic_int64(values, timelike=True)

    cpdef get_loc(self, object val):
        if is_definitely_invalid_key(val):
            raise TypeError

        # Welcome to the spaghetti factory
        if self.over_size_threshold and self.is_monotonic_increasing:
            if not self.is_unique:
                val = _to_i8(val)
                return self._get_loc_duplicates(val)
            values = self._get_index_values()

            try:
                conv = _to_i8(val)
                loc = values.searchsorted(conv, side='left')
            except TypeError:
                self._date_check_type(val)
                raise KeyError(val)

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
        if values.dtype != self._get_box_dtype():
            return np.repeat(-1, len(values)).astype('i4')
        values = np.asarray(values).view('i8')
        return self.mapping.lookup(values)

    def get_pad_indexer(self, other, limit=None):
        if other.dtype != self._get_box_dtype():
            return np.repeat(-1, len(other)).astype('i4')
        other = np.asarray(other).view('i8')
        return algos.pad_int64(self._get_index_values(), other,
                                limit=limit)

    def get_backfill_indexer(self, other, limit=None):
        if other.dtype != self._get_box_dtype():
            return np.repeat(-1, len(other)).astype('i4')
        other = np.asarray(other).view('i8')
        return algos.backfill_int64(self._get_index_values(), other,
                                     limit=limit)


cdef class TimedeltaEngine(DatetimeEngine):

    cdef _get_box_dtype(self):
        return 'm8[ns]'

cpdef convert_scalar(ndarray arr, object value):
    if arr.descr.type_num == NPY_DATETIME:
        if isinstance(value, np.ndarray):
            pass
        elif isinstance(value, Timestamp):
            return value.value
        elif value is None or value != value:
            return iNaT
        else:
            return Timestamp(value).value
    elif arr.descr.type_num == NPY_TIMEDELTA:
        if isinstance(value, np.ndarray):
            pass
        elif isinstance(value, Timedelta):
            return value.value
        elif value is None or value != value:
            return iNaT
        else:
            return Timedelta(value).value

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
            tzinfo = getattr(val, 'tzinfo', None)
            # Save the original date value so we can get the utcoffset from it.
            ival = _pydatetime_to_dts(val, &dts)
            if tzinfo is not None and not _is_utc(tzinfo):
                offset = tslib._get_utcoffset(tzinfo, val)
                ival -= tslib._delta_to_nanoseconds(offset)
            return ival
        return val

cdef inline bint _is_utc(object tz):
    return tz is UTC or isinstance(tz, _du_utc)
