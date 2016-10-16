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
    bint PyTime_Check(object o)
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
        bint avoid_hashtable

    cdef:
        # Metadata indicating properties of our index labels. These are lazily
        # computed and used to dispatch to more efficient algorithms
        # (e.g. binary search instead of linear) when possible.
        bint unique, monotonic_inc, monotonic_dec

        # Flags indicating whether we've populated a hash table, whether we've
        # checked our index labels for uniqueness, and whether we've checked
        # our index labels for monotonicity.
        bint hashtable_populated, uniqueness_checked, monotonicity_checked

    def __init__(self, vgetter, n, avoid_hashtable=None):
        self.vgetter = vgetter

        if avoid_hashtable is not None:
            self.avoid_hashtable = avoid_hashtable
        else:
            self.avoid_hashtable = (n >= _SIZE_CUTOFF)

        self.unique = 0
        self.monotonic_inc = 0
        self.monotonic_dec = 0

        self.hashtable_populated = 0
        self.uniqueness_checked = 0
        self.monotonicity_checked = 0

    def _call_monotonic(self, ndarray values):
        """Check whether ``values`` are monotonic."""
        raise NotImplementedError('_call_monotonic')

    cdef _make_hash_table(self, ndarray values):
        """Make a hash table specific to our value type."""
        raise NotImplementedError('_make_hash_table')

    def get_pad_indexer(self, other, limit=None):
        raise NotImplementedError('get_pad_indexer')

    def get_backfill_indexer(self, other, limit=None):
        raise NotImplementedError('get_backfill_indexer')

    cpdef get_value(self, ndarray arr, object key, object tz=None):
        """
        Look up value(s) from ``arr`` at the index(es) of ``key``.

        Roughly equivalent to ``arr[self.get_loc(key)]``, with special handling
        for datetime types.

        Parameters
        ----------
        arr : ndarray
            Array from which to look up values.
        key : object
            Index key to use to find values in ``arr``.
        tz : object, optional
            Timezone to associate with ``key``.  Ignored unless ``arr`` is of
            datetime dtype.
        """
        cdef:
            object loc

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
        Set ``value`` into ``arr`` at the index(es) of ``key``.

        Roughly equivalent to ``arr[self.get_loc(key)] = value``.

        Parameters
        ----------
        arr : ndarray
            Array from which to look up values.
        key : object
            Index key to use to find values in ``arr``.
        """
        # XXX: Why does get_value take a tz but this method doesn't?
        cdef:
            object loc

        loc = self.get_loc(key)
        value = convert_scalar(arr, value)

        if PySlice_Check(loc) or cnp.PyArray_Check(loc):
            arr[loc] = value
        else:
            util.set_value_at(arr, loc, value)

    def __contains__(self, object val):
        hash(val)  # Force a TypeError on non-hashable input.
        try:
            self.get_loc(val)
            return True
        except KeyError:
            return False

    cpdef get_loc(self, object val):
        """
        Compute an indexer for a single key.

        May return an integer, a slice, or a boolean array, depending on the
        structure of the index values.

        - If the index values are unique, returns an integer.

        - If the index values are non-unique but monotonically-increasing,
          returns an integer or a slice.

        - If the index values are non-unique and non-monotonically-increasing,
          returns a boolean array.

        Raises a KeyError if `val` does not appear in the index.
        """
        if is_definitely_invalid_key(val):
            raise TypeError("'{val}' is an invalid key".format(val=val))

        self._check_type(val)

        if self.avoid_hashtable and self.is_monotonic_increasing:
            # XXX: This branch is duplicated with an identical branch below
            # because the first access of is_monotonic_increasing can set
            # self.is_unique without requiring the creation of a hashtable over
            # our values.
            if not self.is_unique:
                return self._get_loc_duplicates(val)
            return self._get_loc_binsearch_scalar(val)

        elif not self.is_unique:
            return self._get_loc_duplicates(val)

        try:
            self._ensure_hashtable_populated()
            return self.mapping.get_item(val)
        except TypeError:
            raise KeyError(val)

    cdef _get_loc_binsearch_scalar(self, object val):
        values = self._get_index_values()
        loc = _bin_search(values, val) # .searchsorted(val, side='left')
        if loc >= len(values):
            raise KeyError(val)
        if util.get_value_at(values, loc) != val:
            raise KeyError(val)
        return loc

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
            if not self.uniqueness_checked:
                if self.avoid_hashtable:
                    # Create a table to determine uniqueness, but don't store.
                    values = self._get_index_values()
                    tmp_table = self._make_hash_table(values)
                    self.unique = (len(tmp_table) == len(values))
                    self.uniqueness_checked = 1
                else:
                    self._ensure_hashtable_populated()
                    assert self.uniqueness_checked, \
                        "Failed to set uniqueness in populate_hashtable!"

            # Either we had already checked uniqueness, or one of the two
            # branches above must have performed the check.
            return self.unique == 1

    property is_monotonic_increasing:

        def __get__(self):
            if not self.monotonicity_checked:
                self._do_monotonic_check()
            return self.monotonic_inc == 1

    property is_monotonic_decreasing:

        def __get__(self):
            if not self.monotonicity_checked:
                self._do_monotonic_check()
            return self.monotonic_dec == 1

    cdef inline _do_monotonic_check(self):
        cdef object is_unique  # This is either a bint or None.
        try:
            values = self._get_index_values()
            self.monotonic_inc, self.monotonic_dec, is_unique = \
                self._call_monotonic(values)
        except TypeError:
            self.monotonic_inc = 0
            self.monotonic_dec = 0
            is_unique = None

        self.monotonicity_checked = 1

        # _call_monotonic returns None for is_unique if uniqueness could not be
        # determined from the monotonicity check.
        if is_unique is not None:
            self.unique = is_unique
            self.uniqueness_checked = 1

    cdef _get_index_values(self):
        return self.vgetter()

    cdef _check_type(self, object val):
        hash(val)

    cdef inline _ensure_hashtable_populated(self):
        if not self.hashtable_populated:
            values = self._get_index_values()

            self.mapping = self._make_hash_table(values)
            self.hashtable_populated = 1

            self.unique = (len(self.mapping) == len(values))
            self.uniqueness_checked = 1

    def clear_mapping(self):
        self.mapping = None
        self.hashtable_populated = 0
        self.monotonicity_checked = 0
        self.uniqueness_checked = 0

        self.unique = 0
        self.monotonic_inc = 0
        self.monotonic_dec = 0

    def get_indexer(self, values):
        self._ensure_hashtable_populated()
        return self.mapping.lookup(values)

    def get_indexer_non_unique(self, targets):
        """ return an indexer suitable for takng from a non unique index
            return the labels in the same order ast the target
            and a missing indexer into the targets (which correspond
            to the -1 indicies in the results """

        cdef:
            ndarray values
            ndarray[int64_t] result, missing
            set stargets
            dict d = {}
            object val
            int count = 0, count_missing = 0
            Py_ssize_t i, j, n, n_t, n_alloc

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

    cdef _make_hash_table(self, ndarray values):
        t = _hash.Int64HashTable(len(values))
        t.map_locations(values)
        return t

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
        if not util.is_integer_object(val):
            raise KeyError(val)

    cdef _maybe_get_bool_indexer(self, object val):
        cdef:
            ndarray[uint8_t, cast=True] indexer
            ndarray[int64_t] values
            int count = 0
            Py_ssize_t i, n
            int64_t ival
            int last_true

        self._check_type(val)

        ival = val

        values = self._get_index_values()
        n = len(values)

        result = np.empty(n, dtype=bool)
        indexer = result.view(np.uint8)

        for i in range(n):
            if values[i] == ival:
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

    cdef _make_hash_table(self, ndarray values):
        t = _hash.Float64HashTable(len(values))
        t.map_locations(values)
        return t

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


cdef class ObjectEngine(IndexEngine):

    cdef _make_hash_table(self, ndarray values):
        t = _hash.PyObjectHashTable(len(values))
        t.map_locations(values)
        return t

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

    cdef _get_index_values(self):
        return self.vgetter().view('i8')

    def _call_monotonic(self, values):
        return algos.is_monotonic_int64(values, timelike=True)

    cpdef get_loc(self, object val):
        if PyTime_Check(val):  # TODO: Document this.
            raise KeyError(val)
        return super(DatetimeEngine, self).get_loc(_to_i8(val))

    cdef inline _check_type(self, object val):
        hash(val)
        if not util.is_integer_object(val):
            raise KeyError(val)

    def get_indexer(self, values):
        self._ensure_hashtable_populated()
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
