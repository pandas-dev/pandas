import warnings

import numpy as np
cimport numpy as cnp
from numpy cimport (
    float32_t,
    float64_t,
    int8_t,
    int16_t,
    int32_t,
    int64_t,
    intp_t,
    ndarray,
    uint8_t,
    uint16_t,
    uint32_t,
    uint64_t,
)
cnp.import_array()


cimport pandas._libs.util as util

from pandas._libs.tslibs import Period
from pandas._libs.tslibs.nattype cimport c_NaT as NaT
from pandas._libs.tslibs.c_timestamp cimport _Timestamp

from pandas._libs.hashtable cimport HashTable

from pandas._libs import algos, hashtable as _hash
from pandas._libs.tslibs import Timedelta, period as periodlib
from pandas._libs.missing import checknull


cdef inline bint is_definitely_invalid_key(object val):
    try:
        hash(val)
    except TypeError:
        return True
    return False


# Don't populate hash tables in monotonic indexes larger than this
_SIZE_CUTOFF = 1_000_000


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

    def __contains__(self, val: object) -> bool:
        # We assume before we get here:
        #  - val is hashable
        self._ensure_mapping_populated()
        return val in self.mapping

    cpdef get_loc(self, object val):
        cdef:
            Py_ssize_t loc

        if is_definitely_invalid_key(val):
            raise TypeError(f"'{val}' is an invalid key")

        if self.over_size_threshold and self.is_monotonic_increasing:
            if not self.is_unique:
                return self._get_loc_duplicates(val)
            values = self._get_index_values()

            self._check_type(val)
            loc = _bin_search(values, val)  # .searchsorted(val, side='left')
            if loc >= len(values):
                raise KeyError(val)
            if values[loc] != val:
                raise KeyError(val)
            return loc

        self._ensure_mapping_populated()
        if not self.unique:
            return self._get_loc_duplicates(val)

        self._check_type(val)

        try:
            return self.mapping.get_item(val)
        except (TypeError, ValueError):
            raise KeyError(val)

    cdef inline _get_loc_duplicates(self, object val):
        cdef:
            Py_ssize_t diff

        if self.is_monotonic_increasing:
            values = self._get_index_values()
            try:
                left = values.searchsorted(val, side='left')
                right = values.searchsorted(val, side='right')
            except TypeError:
                # e.g. GH#29189 get_loc(None) with a Float64Index
                raise KeyError(val)

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
            ndarray[uint8_t, ndim=1, cast=True] indexer

        indexer = self._get_index_values() == val
        return self._unpack_bool_indexer(indexer, val)

    cdef _unpack_bool_indexer(self,
                              ndarray[uint8_t, ndim=1, cast=True] indexer,
                              object val):
        cdef:
            ndarray[intp_t, ndim=1] found
            int count

        found = np.where(indexer)[0]
        count = len(found)

        if count > 1:
            return indexer
        if count == 1:
            return int(found[0])

        raise KeyError(val)

    def sizeof(self, deep: bool = False) -> int:
        """ return the sizeof our mapping """
        if not self.is_mapping_populated:
            return 0
        return self.mapping.sizeof(deep=deep)

    def __sizeof__(self) -> int:
        return self.sizeof()

    @property
    def is_unique(self) -> bool:
        if self.need_unique_check:
            self._do_unique_check()

        return self.unique == 1

    cdef inline _do_unique_check(self):

        # this de-facto the same
        self._ensure_mapping_populated()

    @property
    def is_monotonic_increasing(self) -> bool:
        if self.need_monotonic_check:
            self._do_monotonic_check()

        return self.monotonic_inc == 1

    @property
    def is_monotonic_decreasing(self) -> bool:
        if self.need_monotonic_check:
            self._do_monotonic_check()

        return self.monotonic_dec == 1

    cdef inline _do_monotonic_check(self):
        cdef:
            bint is_unique
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

    cdef _call_monotonic(self, values):
        return algos.is_monotonic(values, timelike=False)

    def get_backfill_indexer(self, other: np.ndarray, limit=None) -> np.ndarray:
        return algos.backfill(self._get_index_values(), other, limit=limit)

    def get_pad_indexer(self, other: np.ndarray, limit=None) -> np.ndarray:
        return algos.pad(self._get_index_values(), other, limit=limit)

    cdef _make_hash_table(self, Py_ssize_t n):
        raise NotImplementedError

    cdef _check_type(self, object val):
        hash(val)

    @property
    def is_mapping_populated(self) -> bool:
        return self.mapping is not None

    cdef inline _ensure_mapping_populated(self):
        # this populates the mapping
        # if its not already populated
        # also satisfies the need_unique_check

        if not self.is_mapping_populated:

            values = self._get_index_values()
            self.mapping = self._make_hash_table(len(values))
            self._call_map_locations(values)

            if len(self.mapping) == len(values):
                self.unique = 1

        self.need_unique_check = 0

    cdef void _call_map_locations(self, values):
        self.mapping.map_locations(values)

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
        """
        Return an indexer suitable for taking from a non unique index
        return the labels in the same order ast the target
        and a missing indexer into the targets (which correspond
        to the -1 indices in the results
        """
        cdef:
            ndarray values, x
            ndarray[int64_t] result, missing
            set stargets, remaining_stargets
            dict d = {}
            object val
            int count = 0, count_missing = 0
            Py_ssize_t i, j, n, n_t, n_alloc

        self._ensure_mapping_populated()
        values = np.array(self._get_index_values(), copy=False)
        stargets = set(targets)
        n = len(values)
        n_t = len(targets)
        if n > 10_000:
            n_alloc = 10_000
        else:
            n_alloc = n

        result = np.empty(n_alloc, dtype=np.int64)
        missing = np.empty(n_t, dtype=np.int64)

        # map each starget to its position in the index
        if stargets and len(stargets) < 5 and self.is_monotonic_increasing:
            # if there are few enough stargets and the index is monotonically
            # increasing, then use binary search for each starget
            remaining_stargets = set()
            for starget in stargets:
                try:
                    start = values.searchsorted(starget, side='left')
                    end = values.searchsorted(starget, side='right')
                except TypeError:  # e.g. if we tried to search for string in int array
                    remaining_stargets.add(starget)
                else:
                    if start != end:
                        d[starget] = list(range(start, end))

            stargets = remaining_stargets

        if stargets:
            # otherwise, map by iterating through all items in the index
            for i in range(n):
                val = values[i]
                if val in stargets:
                    if val not in d:
                        d[val] = []
                    d[val].append(i)

        for i in range(n_t):
            val = targets[i]

            # found
            if val in d:
                for j in d[val]:

                    # realloc if needed
                    if count >= n_alloc:
                        n_alloc += 10_000
                        result = np.resize(result, n_alloc)

                    result[count] = j
                    count += 1

            # value not found
            else:

                if count >= n_alloc:
                    n_alloc += 10_000
                    result = np.resize(result, n_alloc)
                result[count] = -1
                count += 1
                missing[count_missing] = i
                count_missing += 1

        return result[0:count], missing[0:count_missing]


cdef Py_ssize_t _bin_search(ndarray values, object val) except -1:
    cdef:
        Py_ssize_t mid = 0, lo = 0, hi = len(values) - 1
        object pval

    if hi == 0 or (hi > 0 and val > values[hi]):
        return len(values)

    while lo < hi:
        mid = (lo + hi) // 2
        pval = values[mid]
        if val < pval:
            hi = mid
        elif val > pval:
            lo = mid + 1
        else:
            while mid > 0 and val == values[mid - 1]:
                mid -= 1
            return mid

    if val <= values[mid]:
        return mid
    else:
        return mid + 1


cdef class ObjectEngine(IndexEngine):
    """
    Index Engine for use with object-dtype Index, namely the base class Index.
    """
    cdef _make_hash_table(self, Py_ssize_t n):
        return _hash.PyObjectHashTable(n)


cdef class DatetimeEngine(Int64Engine):

    cdef str _get_box_dtype(self):
        return 'M8[ns]'

    cdef int64_t _unbox_scalar(self, scalar) except? -1:
        # NB: caller is responsible for ensuring tzawareness compat
        #  before we get here
        if not (isinstance(scalar, _Timestamp) or scalar is NaT):
            raise TypeError(scalar)
        return scalar.value

    def __contains__(self, val: object) -> bool:
        # We assume before we get here:
        #  - val is hashable
        cdef:
            int64_t loc, conv

        conv = self._unbox_scalar(val)
        if self.over_size_threshold and self.is_monotonic_increasing:
            if not self.is_unique:
                return self._get_loc_duplicates(conv)
            values = self._get_index_values()
            loc = values.searchsorted(conv, side='left')
            return values[loc] == conv

        self._ensure_mapping_populated()
        return conv in self.mapping

    cdef _get_index_values(self):
        return self.vgetter().view('i8')

    cdef _call_monotonic(self, values):
        return algos.is_monotonic(values, timelike=True)

    cpdef get_loc(self, object val):
        # NB: the caller is responsible for ensuring that we are called
        #  with either a Timestamp or NaT (Timedelta or NaT for TimedeltaEngine)

        cdef:
            int64_t loc
        if is_definitely_invalid_key(val):
            raise TypeError(f"'{val}' is an invalid key")

        try:
            conv = self._unbox_scalar(val)
        except TypeError:
            raise KeyError(val)

        # Welcome to the spaghetti factory
        if self.over_size_threshold and self.is_monotonic_increasing:
            if not self.is_unique:
                return self._get_loc_duplicates(conv)
            values = self._get_index_values()

            loc = values.searchsorted(conv, side='left')

            if loc == len(values) or values[loc] != conv:
                raise KeyError(val)
            return loc

        self._ensure_mapping_populated()
        if not self.unique:
            return self._get_loc_duplicates(conv)

        try:
            return self.mapping.get_item(conv)
        except KeyError:
            raise KeyError(val)

    def get_indexer(self, values):
        self._ensure_mapping_populated()
        if values.dtype != self._get_box_dtype():
            return np.repeat(-1, len(values)).astype('i4')
        values = np.asarray(values).view('i8')
        return self.mapping.lookup(values)

    def get_pad_indexer(self, other: np.ndarray, limit=None) -> np.ndarray:
        if other.dtype != self._get_box_dtype():
            return np.repeat(-1, len(other)).astype('i4')
        other = np.asarray(other).view('i8')
        return algos.pad(self._get_index_values(), other, limit=limit)

    def get_backfill_indexer(self, other: np.ndarray, limit=None) -> np.ndarray:
        if other.dtype != self._get_box_dtype():
            return np.repeat(-1, len(other)).astype('i4')
        other = np.asarray(other).view('i8')
        return algos.backfill(self._get_index_values(), other, limit=limit)


cdef class TimedeltaEngine(DatetimeEngine):

    cdef str _get_box_dtype(self):
        return 'm8[ns]'

    cdef int64_t _unbox_scalar(self, scalar) except? -1:
        if not (isinstance(scalar, Timedelta) or scalar is NaT):
            raise TypeError(scalar)
        return scalar.value


cdef class PeriodEngine(Int64Engine):

    cdef int64_t _unbox_scalar(self, scalar) except? -1:
        if scalar is NaT:
            return scalar.value
        if isinstance(scalar, Period):
            # NB: we assume that we have the correct freq here.
            return scalar.ordinal
        raise TypeError(scalar)

    cpdef get_loc(self, object val):
        # NB: the caller is responsible for ensuring that we are called
        #  with either a Period or NaT
        cdef:
            int64_t conv

        try:
            conv = self._unbox_scalar(val)
        except TypeError:
            raise KeyError(val)

        return Int64Engine.get_loc(self, conv)

    cdef _get_index_values(self):
        return super(PeriodEngine, self).vgetter().view("i8")

    cdef _call_monotonic(self, values):
        return algos.is_monotonic(values, timelike=True)

    def get_indexer(self, values):
        cdef:
            ndarray[int64_t, ndim=1] ordinals

        super(PeriodEngine, self)._ensure_mapping_populated()

        freq = super(PeriodEngine, self).vgetter().freq
        ordinals = periodlib.extract_ordinals(values, freq)

        return self.mapping.lookup(ordinals)

    def get_pad_indexer(self, other: np.ndarray, limit=None) -> np.ndarray:
        freq = super(PeriodEngine, self).vgetter().freq
        ordinal = periodlib.extract_ordinals(other, freq)

        return algos.pad(self._get_index_values(),
                         np.asarray(ordinal), limit=limit)

    def get_backfill_indexer(self, other: np.ndarray, limit=None) -> np.ndarray:
        freq = super(PeriodEngine, self).vgetter().freq
        ordinal = periodlib.extract_ordinals(other, freq)

        return algos.backfill(self._get_index_values(),
                              np.asarray(ordinal), limit=limit)

    def get_indexer_non_unique(self, targets):
        freq = super(PeriodEngine, self).vgetter().freq
        ordinal = periodlib.extract_ordinals(targets, freq)
        ordinal_array = np.asarray(ordinal)

        return super(PeriodEngine, self).get_indexer_non_unique(ordinal_array)


cdef class BaseMultiIndexCodesEngine:
    """
    Base class for MultiIndexUIntEngine and MultiIndexPyIntEngine, which
    represent each label in a MultiIndex as an integer, by juxtaposing the bits
    encoding each level, with appropriate offsets.

    For instance: if 3 levels have respectively 3, 6 and 1 possible values,
    then their labels can be represented using respectively 2, 3 and 1 bits,
    as follows:
     _ _ _ _____ _ __ __ __
    |0|0|0| ... |0| 0|a1|a0| -> offset 0 (first level)
     — — — ————— — —— —— ——
    |0|0|0| ... |0|b2|b1|b0| -> offset 2 (bits required for first level)
     — — — ————— — —— —— ——
    |0|0|0| ... |0| 0| 0|c0| -> offset 5 (bits required for first two levels)
     ‾ ‾ ‾ ‾‾‾‾‾ ‾ ‾‾ ‾‾ ‾‾
    and the resulting unsigned integer representation will be:
     _ _ _ _____ _ __ __ __ __ __ __
    |0|0|0| ... |0|c0|b2|b1|b0|a1|a0|
     ‾ ‾ ‾ ‾‾‾‾‾ ‾ ‾‾ ‾‾ ‾‾ ‾‾ ‾‾ ‾‾

    Offsets are calculated at initialization, labels are transformed by method
    _codes_to_ints.

    Keys are located by first locating each component against the respective
    level, then locating (the integer representation of) codes.
    """
    def __init__(self, object levels, object labels,
                 ndarray[uint64_t, ndim=1] offsets):
        """
        Parameters
        ----------
        levels : list-like of numpy arrays
            Levels of the MultiIndex.
        labels : list-like of numpy arrays of integer dtype
            Labels of the MultiIndex.
        offsets : numpy array of uint64 dtype
            Pre-calculated offsets, one for each level of the index.
        """
        self.levels = levels
        self.offsets = offsets

        # Transform labels in a single array, and add 1 so that we are working
        # with positive integers (-1 for NaN becomes 0):
        codes = (np.array(labels, dtype='int64').T + 1).astype('uint64',
                                                               copy=False)

        # Map each codes combination in the index to an integer unambiguously
        # (no collisions possible), based on the "offsets", which describe the
        # number of bits to switch labels for each level:
        lab_ints = self._codes_to_ints(codes)

        # Initialize underlying index (e.g. libindex.UInt64Engine) with
        # integers representing labels: we will use its get_loc and get_indexer
        self._base.__init__(self, lambda: lab_ints, len(lab_ints))

    def _codes_to_ints(self, codes):
        raise NotImplementedError("Implemented by subclass")

    def _extract_level_codes(self, object target):
        """
        Map the requested list of (tuple) keys to their integer representations
        for searching in the underlying integer index.

        Parameters
        ----------
        target : list-like of keys
            Each key is a tuple, with a label for each level of the index.

        Returns
        ------
        int_keys : 1-dimensional array of dtype uint64 or object
            Integers representing one combination each
        """
        level_codes = [lev.get_indexer(codes, method=method) + 1 for lev, codes
                       in zip(self.levels, zip(*target))]
        return self._codes_to_ints(np.array(level_codes, dtype='uint64').T)
        print('\n\n<DEBUGGING>\n\n')
        print('extracting level codes from:\n{}\nwith:\n{}\nwith method: {}\n\n'.format(
              str(self.levels),
              str(zip(*target)),
              str(method)))
        level_codes = np.array([
            lev.get_indexer(codes, method=method) for lev, codes
            in zip(self.levels, zip(*target))
        ], dtype='uint64').T + 1

        print('this gives us:\n{}'.format(str(level_codes)))
        # idk if truthy/falsy stuff works w/cython...
        # also this entire block should basically be moved into its own helper function
        # or something like that
        if method is not None:
            print('method is {}, so we need any fills at level n to override '
                  'all fills at level n + k, for all k > 0'.format(str(method)))
            level_codes_no_fill = np.array([
                lev.get_indexer(codes) for lev, codes
                in zip(self.levels, zip(*target))
            ], dtype='uint64').T + 1
            print('without {}-style filling, the level_codes are:\n{}'.format(
                  str(method),
                  str(level_codes_no_fill)))


            # TODO: add arithmetic for (backfill-only) "incremementing" when we hit
            # a "0" (i.e. NaN, i.e. too large value) at level i (i.e. incremementing
            # levels 1, ..., i-1).  this necessarily involves "place-value arithmetic"
            # w.r.t. map(len, self.levels), i.e. the max value we can have at each
            # level, after which we have to set it to 0 (i.e. "1") and then "carry"
            #
            # then, when we hit a case where we need to "carry" past level 0, we need
            # to set the whole row to -1
            #
            # it will LIKELY be the case that we need to pass that past this function,
            # but, for the moment, let's see if we can just do something with level
            # codes
            # 
            # my HOPE is that the result will be something like:
            # sum(1 << offset for offset in offsets) (or maybe offsets[:-1], not sure
            # yet....)
            #
            # let's impl this now!
            # eventually let's see if this can be pulled out into a helper function
            if method == 'backfill':
                print('we\'re backfilling, so we want to find values which we might '
                      'need to be bumped to the next value in one or more levels or '
                      'removed entirely')
                for i, row in enumerate(level_codes):
                    print('examining row: {}'.format(str(row)))
                    need_to_carry = False
                    # go from right to left for place-value arithmetic
                    for j in range(len(row) - 1, -1, -1):
                        print('looking at row[{}], i.e. {}'.format(str(j), str(row[j])))
                        # this is AFTER backfilling, i.e. this means value was too large
                        # subtract 1 bc all values here have 1 added to them
                        max_val = len(self.levels[j])
                        print('the highest value you can have in this row is: {}'.format(str(max_val)))
                        if row[j] == 0 or (level_codes_no_fill[i][j] == max_val):
                            need_to_carry = True
                            print('row[{}], i.e. {}, is already too large (max val is {}), '
                                  'or will be after incrementing, so we need to continue '
                                  'carrying'.format(str(j), str(row[j]), str(max_val)))
                        elif (row[j] == max_val):
                            print('this row is at the max value, but since it was backfilled '
                                  'up to it, we dont actually need to do anything')
                            need_to_carry = False

                        if need_to_carry:
                            print('we need to carry')
                            new_val = row[j] + 1 if row[j] == level_codes_no_fill[i][j] else new_val
                            print('new value would be {}'.format(new_val))
                            if new_val > max_val or row[j] == 0:
                                print('which it too big (if new_val is 1 and we\'re here, that\'s because '
                                      'it was previously 0 *AFTER* backfilling), so we will see what '
                                      'happens now')
                                # if we're at the first row, then we're done, this whole thing
                                # is too big
                                if j == 0:
                                    print('we\'re at the end, i.e. level 0, so we gotta just set the '
                                          'whole thing to NaN, which we\'ll in turn need to clean up '
                                          'later to prevent it from being backfilled to the smallest '
                                          'value')
                                    for k in range(len(row)):
                                        row[k] = 0
                                        print('the whole row is now: {}'.format(str(level_codes[i])))
                                # we need to keep carrying, but will not necessarily be a problem
                                else:
                                    print('we\'re not at the end, so there\'s still hope -- setting this '
                                          'to 1, i.e. min value, since that\'s what it should be if the "next", '
                                          'i.e. next one to the *left*, can be bumped')
                                    row[j] = 1
                            # done carrying, for now
                            else:
                                print('new value, {}, is legit, so we can stop carrying for now'.format(
                                      str(new_val)))
                                row[j] = new_val
                                need_to_carry = False
                                print('`need_to_carry` is now False, continuing')

            # MOTIVATION:
            # this is basically just backfilling/forward-filling a tuple into a list
            # of tuples.  the ordering of tuples (i.e. x_i vs. y_i, with ties broken
            # by y_{i+1], y_{i+1}, etc.) is such that we want to explicitly avoid
            # considering values in levels i+1, ..., n (for an n-level MultiIndex)
            # when the value at level i is different w/and w/o filling
            # 
            # - when backfilling, it'll bump the value at level i, so we want to
            # decrease the values at levels i+1, ..., n, to their min value (i.e. "1"),
            # since we're using an arithmetic here where 1 is 0, more or less (TODO: add
            # caveats, formalizations, etc., to this).  this will work since, for all
            # possible values x_1, ..., x_{i-1}, x_{i+1}, ..., x_n, it's necessarily
            # the case that (x_1, ..., x_n) < (x_1, ..., x_{i-1}, x_{i}', 1, ..., 1),
            # for all x_{i}' > x_i
            #
            # - when forward-filling (aka "padding"), it'll drop the value at level i,
            # and so we want to increase the values at levels i+1, ..., n, to their
            # max possible values, i.e. map(len, self.levels[i+1:])
            #
            #
            # TODO: see if this can be replaced w/higher-level functions w/o
            # sacrificing performance.  presently unclear if that can be
            # accomplished
            # als TODO: clean this up
            for i, row in enumerate(level_codes):
                # THIS is where i can apply the algorithm described earlier
                for j, level in enumerate(row):
                  # if it was filled from the prev/next value, then everything after
                  # that should either be the min val, if backfill, or the max, if
                  # using pad.  i think.  lemme mull this one over a bit more
                  if level_codes_no_fill[i][j] == 0 and level_codes[i][j] >= 1:
                      print('without filling, level_codes[{}][{}], curently {}, would be {}'.format(
                            str(i),
                            str(j),
                            str(level_codes[i][j]),
                            '0'))
                      for k in range(j + 1, len(row)):
                          old_val = row[k]
                          row[k] = 1 if method == 'backfill' else len(self.levels[k])
                          print('replaced level_codes[{}][{}], previously {}, with {}'.format(
                                str(i),
                                str(k),
                                str(old_val),
                                str(row[k])))

        print('our cleaned-up level codes are now:\n{},\nwhich we\'ll pass to self._codes_to_ints()'.format(str(level_codes)))
        int_reprs = self._codes_to_ints(level_codes)
        print('which, using the index backend\'s int representations, is:\n{}'.format(str(int_reprs)))
        return int_reprs
        #return self._codes_to_ints(level_codes)

    def get_indexer(self, object target, object method=None,
                    object limit=None):
        lab_ints = self._extract_level_codes(target, method=method)

        print('extracting level codes from {} on {} (w/method = {}) returned:\n{}'.format(
              str(self),
              str(target),
              str(method),
              str(lab_ints)))
        # All methods (exact, backfill, pad) directly map to the respective
        # methods of the underlying (integers) index...
        if method is not None:
            # but underlying backfill and pad methods require index and keys
            # to be sorted. The index already is (checked in
            # Index._get_fill_indexer), sort (integer representations of) keys:
            order = np.argsort(lab_ints)
            lab_ints = lab_ints[order]
            indexer = (getattr(self._base, f'get_{method}_indexer')
                       (self, lab_ints, limit=limit))

            # TODO: completely replace this with correctly-generalized fix
            # with pre-processes by "adding 1" when backfilling when it finds
            # 0s (i.e. -1 (i.e. NaN repr) + 1)
            #
            # HACK: because the backfilled representation of NAs at any level l
            # will always be 2^(offset_l), as a result of adding 1 to -1 and
            # then shifting by offset_l bits, when backfilling, we need to
            # prevent too-large values from being treated as too small values
            # the fix for this, therefore, is to change all too-small values to
            # 1 larger than the max of the entity being backfilled
            if method == 'backfill':
              # inherently necessary bc NaNs are stored as a value which would
              # otherwise be backfilled to zero
              print('method is backfill, doing some cleanup, let\'s see if it\'s sufficient... ')
              print('lab_ints (after re-ordering) is: {}'.format(str(lab_ints)))
              print('indexer is: {}'.format(str(indexer)))
              #indexer = np.array([
              #    idx if lab_ints[i] != nan_val else -1
              #    for i, idx in enumerate(indexer)
              #], dtype=indexer.dtype)
              # any too-large values will be backfilled into the minimum index value
              for i in range(len(indexer)):
                  if lab_ints[i] == 0:
                      print('lab_ints[{}] = {}, but was (incorrectly) backfilled to {}, re-setting to -1'.format(
                            str(i),
                            str(lab_ints[i]),
                            str(indexer[i])))
                      indexer[i] = -1
              print('after cleanup, indexer is: {}'.format(str(indexer)))

            # TODO: try using argsort more to accomplish this.  does not matter for now, though
            new_indexer = [0] * len(indexer)
            for i, idx in enumerate(order):
              new_indexer[idx] = indexer[i]
            print('after fixing order, indexer is: {}'.format(str(new_indexer)))
            return new_indexer
        else:
            indexer = self._base.get_indexer(self, lab_ints)

        return indexer

    def get_loc(self, object key):
        if is_definitely_invalid_key(key):
            raise TypeError(f"'{key}' is an invalid key")
        if not isinstance(key, tuple):
            raise KeyError(key)
        try:
            indices = [0 if checknull(v) else lev.get_loc(v) + 1
                       for lev, v in zip(self.levels, key)]
        except KeyError:
            raise KeyError(key)

        # Transform indices into single integer:
        lab_int = self._codes_to_ints(np.array(indices, dtype='uint64'))

        return self._base.get_loc(self, lab_int)

    def get_indexer_non_unique(self, object target):
        # This needs to be overridden just because the default one works on
        # target._values, and target can be itself a MultiIndex.

        lab_ints = self._extract_level_codes(target)
        indexer = self._base.get_indexer_non_unique(self, lab_ints)

        return indexer

    def __contains__(self, val: object) -> bool:
        # We assume before we get here:
        #  - val is hashable
        # Default __contains__ looks in the underlying mapping, which in this
        # case only contains integer representations.
        try:
            self.get_loc(val)
            return True
        except (KeyError, TypeError, ValueError):
            return False


# Generated from template.
include "index_class_helper.pxi"
