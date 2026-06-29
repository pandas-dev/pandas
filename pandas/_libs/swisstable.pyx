# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
"""
Cython wrapper for Swiss Table hash maps

This provides Python-accessible interfaces to the C Swiss table implementation
for all pandas numeric types.

Available types:
- SwissInt64Map: int64_t -> size_t
- SwissUInt64Map: uint64_t -> size_t
- SwissInt32Map: int32_t -> size_t
- SwissUInt32Map: uint32_t -> size_t
- SwissInt16Map: int16_t -> size_t
- SwissUInt16Map: uint16_t -> size_t
- SwissInt8Map: int8_t -> size_t
- SwissUInt8Map: uint8_t -> size_t
- SwissFloat64Map: double -> size_t (with NaN handling)
- SwissFloat32Map: float -> size_t (with NaN handling)
- SwissComplex64Map: complex64 -> size_t (with NaN handling)
- SwissComplex128Map: complex128 -> size_t (with NaN handling)

Integration Methods (matching khash HashTable API):
- map_locations(): Build table mapping values to their array positions
- lookup(): Look up array values, returning positions or -1 for misses
- ismember_*(): Standalone functions for membership testing (isin)
- unique(): Get unique values with optional inverse mapping
- factorize(): Get unique values and labels (for categorical encoding)
- value_count_*(): Standalone functions for counting value occurrences
- duplicated_*(): Standalone functions for finding duplicated values
"""

cimport cython
from cpython.object cimport PyObject
from cpython.ref cimport Py_INCREF, Py_DECREF
from libc.math cimport isnan
from libc.stddef cimport size_t
from libc.stdint cimport (
    int8_t,
    int16_t,
    int32_t,
    int64_t,
    uint8_t,
    uint16_t,
    uint32_t,
    uint64_t,
)

import numpy as np

cimport numpy as cnp
from numpy cimport (
    intp_t,
    ndarray,
    uint8_t as np_uint8_t,
    PyArray_DATA,
)

from pandas._libs.missing cimport (
    checknull,
)
# intp_t matches ssize_t in the C header
ctypedef double float64_t
ctypedef float float32_t

cnp.import_array()

include "swisstable_class_helper.pxi"


# =============================================================================
# Float64 Map (with NaN handling)
# =============================================================================
cdef class SwissFloat64Map(HashTable):
    """
    Swiss Table hash map: float64 -> size_t

    NaN handling: All NaN values are considered equal.
    Zero handling: +0.0 and -0.0 are considered equal.
    """

    def __cinit__(self, size_t size_hint=0):
        self.uses_mask = False
        if size_hint > 0:
            self.table.reserve(size_hint)

    def insert(self, double key, size_t value):
        cdef int ret = self.table.insert(key, value)
        return ret == 1

    def get(self, double key, default=None):
        cdef size_t val
        if self.table.get(key, &val):
            return val
        return default

    def __getitem__(self, double key):
        cdef size_t val
        if self.table.get(key, &val):
            return val
        raise KeyError(key)

    def __setitem__(self, double key, size_t value):
        cdef int ret = self.table.insert(key, value)
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    cpdef get_item(self, double key):
        cdef size_t val
        if self.table.get(key, &val):
            return val
        raise KeyError(key)

    cpdef set_item(self, double key, Py_ssize_t value):
        cdef int ret = self.table.insert(key, value)
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    cpdef get_na(self):
        """Extracts the position of na_value from the hashtable.

        Returns
        -------
        The position of the last na value.
        """
        raise NotImplementedError("not implement for get_na")

    cpdef set_na(self, Py_ssize_t val):
        # Caller is responsible for checking for pd.NA
        raise NotImplementedError("not implement for set_na")

    def __contains__(self, double key):
        return self.table.contains(key)

    def __len__(self):
        return self.table.size()

    def __iter__(self):
        cdef size_t index = 0
        while index < self.table.capacity():
            index = self.table.iter_next(index)
            if index < self.table.capacity():
                yield self.table.key_at(index)
                index += 1

    def items(self):
        cdef size_t index = 0
        while index < self.table.capacity():
            index = self.table.iter_next(index)
            if index < self.table.capacity():
                yield (self.table.key_at(index),
                       self.table.val_at(index))
                index += 1

    @property
    def size(self):
        return self.table.size()

    @property
    def capacity(self):
        return self.table.capacity()

    def __repr__(self):
        return f"SwissFloat64Map(size={self.size}, capacity={self.capacity})"

    def map_locations(self, const double[:] values, const uint8_t[:] mask = None) -> None:
        """Build table mapping values to their array positions."""
        cdef:
            Py_ssize_t i, n = len(values)
            double val
            int ret = 0

        if self.uses_mask and mask is not None:
            raise NotImplementedError("mask not implement for swisstable")  # pragma: no cover

        with nogil:
            if n > 0 and values.strides[0] == sizeof(double):
                ret = self.table.map_locations(&values[0], <size_t>n)
            else:
                for i in range(n):
                    val = values[i]
                    ret = self.table.insert(val, <size_t>i)
                    if ret == -1:
                        break

        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    def lookup(self, const double[:] values) -> ndarray:
        """Look up array values in table, returning positions (-1 if not found)."""
        cdef:
            Py_ssize_t i, n = len(values)
            double val
            size_t loc
            intp_t[::1] locs = np.empty(n, dtype=np.intp)

        with nogil:
            if n > 0 and values.strides[0] == sizeof(double):
                self.table.lookup_batch(&values[0], <size_t>n, &locs[0])
            else:
                for i in range(n):
                    val = values[i]
                    if self.table.get(val, &loc):
                        locs[i] = <intp_t>loc
                    else:
                        locs[i] = -1

        return np.asarray(locs)

    def unique(self, const double[:] values, bint return_inverse=False,
               object mask=None):
        """Calculate unique values and optionally their inverse mapping."""
        cdef:
            Py_ssize_t i, n = len(values), count = 0
            double val
            size_t idx
            int ret = 0
            double[::1] uniques_arr = np.empty(n, dtype=np.float64)
            intp_t[::1] labels
            bint uses_mask = mask is not None
            bint seen_na = False
            const uint8_t[:] mask_view
            np_uint8_t[::1] result_mask_arr

        if not uses_mask and n > 0 and values.strides[0] == sizeof(double):
            if return_inverse:
                labels = np.empty(n, dtype=np.intp)
                with nogil:
                    count = self.table.unique_with_inverse(
                        &values[0], <size_t>n,
                        &uniques_arr[0], &labels[0]
                    )
                if count == -1:
                    raise MemoryError("Failed to insert into Swiss table")
                uniques = np.asarray(uniques_arr)[:count].copy()
                return uniques, np.asarray(labels)
            else:
                with nogil:
                    count = self.table.unique_batch(&values[0], <size_t>n, &uniques_arr[0])
                if count == -1:
                    raise MemoryError("Failed to insert into Swiss table")
                uniques = np.asarray(uniques_arr)[:count].copy()
                return uniques

        if uses_mask:
            mask_view = mask if mask.dtype == np.uint8 else mask.view(np.uint8)
            result_mask_arr = np.empty(n, dtype=np.uint8)

        if return_inverse:
            labels = np.empty(n, dtype=np.intp)
            with nogil:
                for i in range(n):
                    if uses_mask and mask_view[i]:
                        if not seen_na:
                            uniques_arr[count] = values[i]
                            result_mask_arr[count] = 1
                            labels[i] = count
                            seen_na = True
                            count += 1
                        else:
                            for idx in range(<size_t>count):
                                if result_mask_arr[idx] == 1:
                                    labels[i] = <intp_t>idx
                                    break
                    else:
                        val = values[i]
                        ret = self.table.insert_if_absent(val, <size_t>count, &idx)  # noqa: E501
                        if ret == 1:
                            uniques_arr[count] = val
                            if uses_mask:
                                result_mask_arr[count] = 0
                            labels[i] = count
                            count += 1
                        elif ret == 0:
                            labels[i] = <intp_t>idx
                        else:
                            break

            if ret == -1:
                raise MemoryError("Failed to insert into Swiss table")

            uniques = np.asarray(uniques_arr)[:count].copy()
            if uses_mask:
                result_mask = np.asarray(result_mask_arr)[:count].copy()
                return uniques, np.asarray(labels), result_mask
            return uniques, np.asarray(labels)
        else:
            with nogil:
                for i in range(n):
                    if uses_mask and mask_view[i]:
                        if not seen_na:
                            uniques_arr[count] = values[i]
                            result_mask_arr[count] = 1
                            seen_na = True
                            count += 1
                    else:
                        val = values[i]
                        ret = self.table.insert_key_only(val)
                        if ret == 1:
                            uniques_arr[count] = val
                            if uses_mask:
                                result_mask_arr[count] = 0
                            count += 1
                        elif ret == -1:
                            break

            if ret == -1:
                raise MemoryError("Failed to insert into Swiss table")

            uniques = np.asarray(uniques_arr)[:count].copy()
            if uses_mask:
                return uniques, np.asarray(result_mask_arr)[:count].copy()
            return uniques

    def factorize(self, const double[:] values, Py_ssize_t na_sentinel=-1,
                  object na_value=None, object mask=None, bint ignore_na=True):
        """Calculate unique values and labels for categorical encoding (NaN-aware)."""
        cdef:
            Py_ssize_t i, n = len(values), count = 0, na_idx = -1
            double val
            size_t idx
            int ret = 0
            double[::1] uniques_arr = np.empty(n, dtype=np.float64)
            intp_t[::1] labels = np.empty(n, dtype=np.intp)
            bint uses_mask = mask is not None
            bint seen_na = False
            const uint8_t[:] mask_view
            double na_val = na_value if na_value is not None else NaNTraits[double].NaN()

        # Fast path: no mask - use C batch function for maximum performance
        if not uses_mask and n > 0 and values.strides[0] == sizeof(double):
            with nogil:
                if ignore_na:
                    count = self.table.factorize_batch(
                        &values[0], <size_t>n,
                        &uniques_arr[0], &labels[0],
                        na_sentinel, na_val
                    )
                else:
                    count = self.table.unique_with_inverse(
                        &values[0], <size_t>n,
                        &uniques_arr[0], &labels[0]
                    )
            if count == -1:
                raise MemoryError("Failed to insert into Swiss table")
            uniques = np.asarray(uniques_arr)[:count].copy()
            return uniques, np.asarray(labels)

        if uses_mask:
            mask_view = mask if mask.dtype == np.uint8 else mask.view(np.uint8)

        # Pre-reserve capacity to avoid resize checks in hot loop
        if not self.table.reserve(<size_t>n):
            raise MemoryError("Failed to reserve capacity in Swiss table")

        with nogil:
            for i in range(n):
                val = values[i]
                # Check for NA: mask takes precedence, then NaN
                if uses_mask and mask_view[i]:
                    if ignore_na:
                        labels[i] = na_sentinel
                    else:
                        if not seen_na:
                            uniques_arr[count] = val
                            na_idx = count
                            labels[i] = count
                            seen_na = True
                            count += 1
                        else:
                            labels[i] = na_idx
                elif ignore_na and NaNTraits[double].AreEqual(val, na_val):
                    labels[i] = na_sentinel
                else:
                    if not self.table.get(val, &idx):
                        ret = self.table.insert(val, <size_t>count)
                        if ret == -1:
                            break
                        uniques_arr[count] = val
                        labels[i] = count
                        count += 1
                    else:
                        labels[i] = <intp_t>idx

        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

        uniques = np.asarray(uniques_arr)[:count].copy()
        return uniques, np.asarray(labels)

    def get_labels(self, const double[:] values, ndarray uniques, bint uniques_view_exits,
        Py_ssize_t count_prior=0, Py_ssize_t na_sentinel=-1,
        object na_value=None, object mask=None):
        """Calculate unique values and labels for categorical encoding."""
        cdef:
            Py_ssize_t i, n = len(values), count = count_prior
            double val
            size_t idx
            int ret = 0
            double[::1] uniques_arr = uniques
            intp_t[::1] labels = np.empty(n, dtype=np.intp)
            bint uses_na_value = na_value is not None
            bint uses_mask = mask is not None
            bint seen_na = False
            const uint8_t[:] mask_view
            double na_val = na_value if uses_na_value else NaNTraits[double].NaN()

        # Fast path: no mask - use C batch function for maximum performance
        if not uses_mask and (count_prior + n <= uniques_arr.shape[0]) and (
            n > 0 and values.strides[0] == sizeof(double)):
            with nogil:
                count = self.table.factorize_batch(
                    &values[0], <size_t>n,
                    &uniques_arr[0], &labels[0],
                    na_sentinel, na_val, count_prior
                )
            if count == -1:
                raise MemoryError("Failed to insert into Swiss table")
            return count, np.asarray(labels)

        if uses_mask:
            mask_view = mask if mask.dtype == np.uint8 else mask.view(np.uint8)

        with nogil:
            for i in range(n):
                val = values[i]
                if (uses_mask and mask_view[i]) or NaNTraits[double].AreEqual(val, na_val):
                    labels[i] = na_sentinel
                else:
                    if not self.table.get(val, &idx):
                        ret = self.table.insert(val, <size_t>count)
                        if ret == -1:
                            break
                        if count >= uniques_arr.shape[0]:
                            with gil:
                                uniques.resize(count * 2, refcheck=uniques_view_exits)
                                uniques_arr = uniques
                        uniques_arr[count] = val
                        labels[i] = count
                        count += 1
                    else:
                        labels[i] = <intp_t>idx

        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

        return count, np.asarray(labels)


def ismember_float64(const double[:] arr, const double[:] values) -> ndarray:
    """Return boolean array indicating membership of arr elements in values."""
    cdef:
        Py_ssize_t n_arr = len(arr), n_values = len(values)
        SwissTable[float64_t, size_t] table
        int ret = 0
        np_uint8_t[::1] result_view

    # Build the lookup table from values
    table.reserve(<size_t>n_values)

    # Pre-allocate result array
    result = np.empty(n_arr, dtype=np.uint8)
    result_view = result

    # Build set from values (key-only, no value storage needed)
    with nogil:
        ret = table.build_set(&values[0], <size_t>n_values)
    if ret == -1:
        raise MemoryError("Failed to insert into Swiss table")

    # Check membership using batch contains
    with nogil:
        table.contains_batch(&arr[0], <size_t>n_arr, &result_view[0])

    return result.view(np.bool_)


def value_count_float64(const double[:] values, bint dropna=True, const uint8_t[:] mask=None):
    """Count occurrences of each unique value (NaN-aware)."""
    cdef:
        Py_ssize_t i, n = len(values), n_unique = 0, na_count = 0, na_add = 0
        SwissTable[float64_t, size_t] table
        float64_t val
        size_t idx
        int ret = 0
        float64_t[::1] keys_arr = np.empty(n, dtype=np.float64)
        size_t[::1] indices_arr = np.empty(n, dtype=np.uintp)
        int64_t[::1] counts_arr
        bint uses_mask = mask is not None
        bint isna_entry = False

    table.reserve(<size_t>n)

    # Use batch value_count which stores indices for direct access
    if not uses_mask and n > 0 and values.strides[0] == sizeof(double):
        with nogil:
            n_unique = table.value_count_batch(
                &values[0], <size_t>n,
                &keys_arr[0], &indices_arr[0], dropna
            )
    else:
        with nogil:
            for i in range(n):
                val = values[i]
                if uses_mask:
                    if mask[i]:
                        if not dropna:
                            na_count += 1
                        continue
                else:
                    if dropna and NaNTraits[float64_t].IsNaN(val):
                        continue

                ret = table.increment(val, &idx)
                if ret == 1:
                    keys_arr[n_unique] = val
                    indices_arr[n_unique] = idx
                    n_unique += 1
                elif ret == -1:
                    break

    if ret == -1:
        raise MemoryError("Failed to insert into Swiss table")

    if na_count > 0:
        na_add = 1

    # Extract counts using stored indices (direct access, no lookup)
    counts_arr = np.empty(n_unique + na_add, dtype=np.int64)
    with nogil:
        for i in range(n_unique):
            counts_arr[i] = <int64_t>table.vals()[indices_arr[i]]

    if na_count > 0:
        keys_arr[n_unique] = val
        counts_arr[n_unique] = na_count

    keys = np.asarray(keys_arr)[:n_unique + na_add].copy()
    counts = np.asarray(counts_arr)
    return keys, counts, na_count


def duplicated_float64(const double[:] values, object keep="first",
                       const uint8_t[:] mask=None):
    """
    Return boolean array indicating duplicated values (NaN-aware).

    Parameters
    ----------
    values : ndarray[float64]
        Array to check for duplicates
    keep : {'first', 'last', False}, default 'first'
        - 'first': Mark duplicates as True except for first occurrence
        - 'last': Mark duplicates as True except for last occurrence
        - False: Mark all duplicates as True
    mask : ndarray[uint8], optional
        If provided, mask[i] == 1 indicates that values[i] is NA/missing.
        NA values are treated as a single group for duplicate detection.
        Mask takes precedence over NaN check in values.

    Returns
    -------
    ndarray[bool]
        Boolean array indicating duplicates
    """
    cdef:
        Py_ssize_t i, n = len(values), first_na = -1
        SwissTable[float64_t, size_t] table
        double val
        size_t idx
        int ret = 0
        np_uint8_t[::1] result
        bint keep_first = keep == "first"
        bint keep_last = keep == "last"
        bint _keep_none = keep is False
        bint uses_mask = mask is not None
        bint seen_na = False
        bint seen_multiple_na = False

    if keep not in ("first", "last", False):
        raise ValueError('keep must be either "first", "last" or False')

    table.reserve(<size_t>n)

    if not uses_mask and n > 0 and values.strides[0] == sizeof(double):
        result = np.empty(n, dtype=np.uint8)
        with nogil:
            if keep_first or keep_last:
                ret = table.duplicated_keep_batch(&values[0], n, keep_first, &result[0])
            else:  # keep == False
                ret = table.duplicated_false_batch(&values[0], n, &result[0])
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

        return np.asarray(result).view(np.bool_)

    result = np.zeros(n, dtype=np.uint8)

    if keep_first:
        # Use insert_key_only: no value storage needed
        with nogil:
            for i in range(n):
                # When mask is provided, use mask only. Otherwise check NaN.
                if uses_mask and mask[i]:
                    if seen_na:
                        result[i] = 1  # Duplicate NA
                    else:
                        seen_na = True
                        result[i] = 0  # First NA
                else:
                    val = values[i]
                    ret = table.insert_key_only(val)
                    if ret == 0:
                        result[i] = 1  # Duplicate (NaN handled by hash table)
                    elif ret == -1:
                        break
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    elif keep_last:
        # Use insert_key_only: no value storage needed
        with nogil:
            for i in range(n - 1, -1, -1):
                # When mask is provided, use mask only. Otherwise check NaN.
                if uses_mask and mask[i]:
                    if seen_na:
                        result[i] = 1  # Duplicate NA
                    else:
                        seen_na = True
                        result[i] = 0  # First NA (from end)
                else:
                    val = values[i]
                    ret = table.insert_key_only(val)
                    if ret == 0:
                        result[i] = 1  # Duplicate (NaN handled by hash table)
                    elif ret == -1:
                        break
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    else:  # keep == False
        # Need insert_if_absent: must store index to mark first occurrence
        with nogil:
            for i in range(n):
                # When mask is provided, use mask only. Otherwise check NaN.
                if uses_mask and mask[i]:
                    # Handle NA values - mark all as duplicates if multiple
                    if not seen_na:
                        first_na = i
                        seen_na = True
                        result[i] = 0
                    elif not seen_multiple_na:
                        result[i] = 1
                        result[first_na] = 1
                        seen_multiple_na = True
                    else:
                        result[i] = 1
                else:
                    val = values[i]
                    ret = table.insert_if_absent(val, <size_t>i, &idx)  # noqa: E501
                    if ret == 0:
                        result[idx] = 1
                        result[i] = 1
                    elif ret == -1:
                        break
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    return np.asarray(result).view(np.bool_)


# =============================================================================
# Float32 Map (with NaN handling)
# =============================================================================
cdef class SwissFloat32Map(HashTable):
    """
    Swiss Table hash map: float32 -> size_t

    NaN handling: All NaN values are considered equal.
    Zero handling: +0.0f and -0.0f are considered equal.
    """

    def __cinit__(self, size_t size_hint=0):
        self.uses_mask = False
        if size_hint > 0:
            self.table.reserve(size_hint)

    def insert(self, float key, size_t value):
        cdef int ret = self.table.insert(key, value)
        return ret == 1

    def get(self, float key, default=None):
        cdef size_t val
        if self.table.get(key, &val):
            return val
        return default

    def __getitem__(self, float key):
        cdef size_t val
        if self.table.get(key, &val):
            return val
        raise KeyError(key)

    def __setitem__(self, float key, size_t value):
        cdef int ret = self.table.insert(key, value)

    cpdef get_item(self, float key):
        cdef size_t val
        if self.table.get(key, &val):
            return val
        raise KeyError(key)

    cpdef set_item(self, float key, Py_ssize_t value):
        cdef int ret = self.table.insert(key, value)

    cpdef get_na(self):
        """Extracts the position of na_value from the hashtable.

        Returns
        -------
        The position of the last na value.
        """
        raise NotImplementedError("not implement for get_na")

    cpdef set_na(self, Py_ssize_t val):
        # Caller is responsible for checking for pd.NA
        raise NotImplementedError("not implement for set_na")

    def __contains__(self, float key):
        return self.table.contains(key)

    def __len__(self):
        return self.table.size()

    def __iter__(self):
        cdef size_t index = 0
        while index < self.table.capacity():
            index = self.table.iter_next(index)
            if index < self.table.capacity():
                yield self.table.key_at(index)
                index += 1

    def items(self):
        cdef size_t index = 0
        while index < self.table.capacity():
            index = self.table.iter_next(index)
            if index < self.table.capacity():
                yield (self.table.key_at(index),
                       self.table.val_at(index))
                index += 1

    @property
    def size(self):
        return self.table.size()

    @property
    def capacity(self):
        return self.table.capacity()

    def __repr__(self):
        return f"SwissFloat32Map(size={self.size}, capacity={self.capacity})"

    def map_locations(self, const float[:] values, const uint8_t[:] mask = None) -> None:
        """Build table mapping values to their array positions."""
        cdef:
            Py_ssize_t i, n = len(values)
            float val
            int ret = 0

        if self.uses_mask and mask is not None:
            raise NotImplementedError("mask not implement for swisstable")  # pragma: no cover

        with nogil:
            if n > 0 and values.strides[0] == sizeof(float):
                ret = self.table.map_locations(&values[0], <size_t>n)
            else:
                for i in range(n):
                    val = values[i]
                    ret = self.table.insert(val, <size_t>i)
                    if ret == -1:
                        break

        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    def lookup(self, const float[:] values) -> ndarray:
        """Look up array values in table, returning positions (-1 if not found)."""
        cdef:
            Py_ssize_t i, n = len(values)
            float val
            size_t loc
            intp_t[::1] locs = np.empty(n, dtype=np.intp)

        with nogil:
            if n > 0 and values.strides[0] == sizeof(float):
                self.table.lookup_batch(&values[0], <size_t>n, &locs[0])
            else:
                for i in range(n):
                    val = values[i]
                    if self.table.get(val, &loc):
                        locs[i] = <intp_t>loc
                    else:
                        locs[i] = -1

        return np.asarray(locs)

    def unique(self, const float[:] values, bint return_inverse=False,
               object mask=None):
        """Calculate unique values and optionally their inverse mapping."""
        cdef:
            Py_ssize_t i, n = len(values), count = 0
            float val
            size_t idx
            int ret = 0
            float[::1] uniques_arr = np.empty(n, dtype=np.float32)
            intp_t[::1] labels
            bint uses_mask = mask is not None
            bint seen_na = False
            const uint8_t[:] mask_view
            np_uint8_t[::1] result_mask_arr

        if not uses_mask and n > 0 and values.strides[0] == sizeof(float):
            if return_inverse:
                labels = np.empty(n, dtype=np.intp)
                with nogil:
                    count = self.table.unique_with_inverse(
                        &values[0], <size_t>n,
                        &uniques_arr[0], &labels[0]
                    )
                if count == -1:
                    raise MemoryError("Failed to insert into Swiss table")
                uniques = np.asarray(uniques_arr)[:count].copy()
                return uniques, np.asarray(labels)
            else:
                with nogil:
                    count = self.table.unique_batch(&values[0], <size_t>n, &uniques_arr[0])
                if count == -1:
                    raise MemoryError("Failed to insert into Swiss table")
                uniques = np.asarray(uniques_arr)[:count].copy()
                return uniques

        if uses_mask:
            mask_view = mask if mask.dtype == np.uint8 else mask.view(np.uint8)
            result_mask_arr = np.empty(n, dtype=np.uint8)

        if return_inverse:
            labels = np.empty(n, dtype=np.intp)
            with nogil:
                for i in range(n):
                    if uses_mask and mask_view[i]:
                        if not seen_na:
                            uniques_arr[count] = values[i]
                            result_mask_arr[count] = 1
                            labels[i] = count
                            seen_na = True
                            count += 1
                        else:
                            for idx in range(<size_t>count):
                                if result_mask_arr[idx] == 1:
                                    labels[i] = <intp_t>idx
                                    break
                    else:
                        val = values[i]
                        ret = self.table.insert_if_absent(val, <size_t>count, &idx)  # noqa: E501
                        if ret == 1:
                            uniques_arr[count] = val
                            if uses_mask:
                                result_mask_arr[count] = 0
                            labels[i] = count
                            count += 1
                        elif ret == 0:
                            labels[i] = <intp_t>idx
                        else:
                            break

            if ret == -1:
                raise MemoryError("Failed to insert into Swiss table")

            uniques = np.asarray(uniques_arr)[:count].copy()
            if uses_mask:
                result_mask = np.asarray(result_mask_arr)[:count].copy()
                return uniques, np.asarray(labels), result_mask
            return uniques, np.asarray(labels)
        else:
            with nogil:
                for i in range(n):
                    if uses_mask and mask_view[i]:
                        if not seen_na:
                            uniques_arr[count] = values[i]
                            result_mask_arr[count] = 1
                            seen_na = True
                            count += 1
                    else:
                        val = values[i]
                        ret = self.table.insert_key_only(val)
                        if ret == 1:
                            uniques_arr[count] = val
                            if uses_mask:
                                result_mask_arr[count] = 0
                            count += 1
                        elif ret == -1:
                            break

            if ret == -1:
                raise MemoryError("Failed to insert into Swiss table")

            uniques = np.asarray(uniques_arr)[:count].copy()
            if uses_mask:
                return uniques, np.asarray(result_mask_arr)[:count].copy()
            return uniques

    def factorize(self, const float[:] values, Py_ssize_t na_sentinel=-1,
                  object na_value=None, object mask=None, bint ignore_na=True):
        """Calculate unique values and labels for categorical encoding (NaN-aware)."""
        cdef:
            Py_ssize_t i, n = len(values), count = 0, na_idx = -1
            float val
            size_t idx
            int ret = 0
            float[::1] uniques_arr = np.empty(n, dtype=np.float32)
            intp_t[::1] labels = np.empty(n, dtype=np.intp)
            bint uses_mask = mask is not None
            bint seen_na = False
            const uint8_t[:] mask_view
            float na_val = na_value if na_value is not None else NaNTraits[float].NaN()

        # Fast path: no mask - use C batch function for maximum performance
        if not uses_mask and n > 0 and values.strides[0] == sizeof(float):
            with nogil:
                if ignore_na:
                    count = self.table.factorize_batch(
                        &values[0], <size_t>n,
                        &uniques_arr[0], &labels[0],
                        na_sentinel, na_val
                    )
                else:
                    count = self.table.unique_with_inverse(
                        &values[0], <size_t>n,
                        &uniques_arr[0], &labels[0]
                    )
            if count == -1:
                raise MemoryError("Failed to insert into Swiss table")
            uniques = np.asarray(uniques_arr)[:count].copy()
            return uniques, np.asarray(labels)

        if uses_mask:
            mask_view = mask if mask.dtype == np.uint8 else mask.view(np.uint8)

        # Pre-reserve capacity to avoid resize checks in hot loop
        if not self.table.reserve(<size_t>n):
            raise MemoryError("Failed to reserve capacity in Swiss table")

        with nogil:
            for i in range(n):
                val = values[i]
                # Check for NA: mask takes precedence, then NaN
                if uses_mask and mask_view[i]:
                    if ignore_na:
                        labels[i] = na_sentinel
                    else:
                        if not seen_na:
                            uniques_arr[count] = val
                            na_idx = count
                            labels[i] = count
                            seen_na = True
                            count += 1
                        else:
                            labels[i] = na_idx
                elif ignore_na and NaNTraits[float].AreEqual(val, na_val):
                    labels[i] = na_sentinel
                else:
                    if not self.table.get(val, &idx):
                        ret = self.table.insert(val, <size_t>count)
                        if ret == -1:
                            break
                        uniques_arr[count] = val
                        labels[i] = count
                        count += 1
                    else:
                        labels[i] = <intp_t>idx

        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

        uniques = np.asarray(uniques_arr)[:count].copy()
        return uniques, np.asarray(labels)

    def get_labels(self, const float[:] values, ndarray uniques, bint uniques_view_exits,
        Py_ssize_t count_prior=0, Py_ssize_t na_sentinel=-1,
        object na_value=None, object mask=None):
        """Calculate unique values and labels for categorical encoding."""
        cdef:
            Py_ssize_t i, n = len(values), count = count_prior
            float val
            size_t idx
            int ret = 0
            float[::1] uniques_arr = uniques
            intp_t[::1] labels = np.empty(n, dtype=np.intp)
            bint uses_na_value = na_value is not None
            bint uses_mask = mask is not None
            bint seen_na = False
            const uint8_t[:] mask_view
            float na_val = na_value if uses_na_value else NaNTraits[float].NaN()

        # Fast path: no mask - use C batch function for maximum performance
        if not uses_mask and (count_prior + n <= uniques_arr.shape[0]) and (
            n > 0 and values.strides[0] == sizeof(float)):
            with nogil:
                count = self.table.factorize_batch(
                    &values[0], <size_t>n,
                    &uniques_arr[0], &labels[0],
                    na_sentinel, na_val, count_prior
                )
            if count == -1:
                raise MemoryError("Failed to insert into Swiss table")
            return count, np.asarray(labels)

        if uses_mask:
            mask_view = mask if mask.dtype == np.uint8 else mask.view(np.uint8)

        with nogil:
            for i in range(n):
                val = values[i]
                if (uses_mask and mask_view[i]) or NaNTraits[float].AreEqual(val, na_val):
                    labels[i] = na_sentinel
                else:
                    if not self.table.get(val, &idx):
                        ret = self.table.insert(val, <size_t>count)
                        if ret == -1:
                            break
                        if count >= uniques_arr.shape[0]:
                            with gil:
                                uniques.resize(count * 2, refcheck=uniques_view_exits)
                                uniques_arr = uniques
                        uniques_arr[count] = val
                        labels[i] = count
                        count += 1
                    else:
                        labels[i] = <intp_t>idx

        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

        return count, np.asarray(labels)


def ismember_float32(const float[:] arr, const float[:] values) -> ndarray:
    """Return boolean array indicating membership of arr elements in values."""
    cdef:
        Py_ssize_t n_arr = len(arr), n_values = len(values)
        SwissTable[float32_t, size_t] table
        int ret = 0
        np_uint8_t[::1] result_view

    # Build the lookup table from values
    table.reserve(<size_t>n_values)

    # Pre-allocate result array
    result = np.empty(n_arr, dtype=np.uint8)
    result_view = result

    # Build set from values (key-only, no value storage needed)
    with nogil:
        ret = table.build_set(&values[0], <size_t>n_values)
    if ret == -1:
        raise MemoryError("Failed to insert into Swiss table")

    # Check membership using batch contains
    with nogil:
        table.contains_batch(&arr[0], <size_t>n_arr, &result_view[0])

    return result.view(np.bool_)


def value_count_float32(const float[:] values, bint dropna=True, const uint8_t[:] mask=None):
    """Count occurrences of each unique value (NaN-aware)."""
    cdef:
        Py_ssize_t i, n = len(values), n_unique = 0, na_count = 0, na_add = 0
        SwissTable[float32_t, size_t] table
        float32_t val
        size_t idx
        int ret = 0
        float32_t[::1] keys_arr = np.empty(n, dtype=np.float32)
        size_t[::1] indices_arr = np.empty(n, dtype=np.uintp)
        int64_t[::1] counts_arr
        bint uses_mask = mask is not None
        bint isna_entry = False

    table.reserve(<size_t>n)

    # Use batch value_count which stores indices for direct access
    if not uses_mask and n > 0 and values.strides[0] == sizeof(float):
        with nogil:
            n_unique = table.value_count_batch(
                &values[0], <size_t>n,
                &keys_arr[0], &indices_arr[0], dropna
            )
    else:
        with nogil:
            for i in range(n):
                val = values[i]
                if uses_mask:
                    if mask[i]:
                        if not dropna:
                            na_count += 1
                        continue
                else:
                    if dropna and NaNTraits[float32_t].IsNaN(val):
                        continue

                ret = table.increment(val, &idx)
                if ret == 1:
                    keys_arr[n_unique] = val
                    indices_arr[n_unique] = idx
                    n_unique += 1
                elif ret == -1:
                    break

    if ret == -1:
        raise MemoryError("Failed to insert into Swiss table")

    if na_count > 0:
        na_add = 1

    # Extract counts using stored indices (direct access, no lookup)
    counts_arr = np.empty(n_unique + na_add, dtype=np.int64)
    with nogil:
        for i in range(n_unique):
            counts_arr[i] = <int64_t>table.vals()[indices_arr[i]]

    if na_count > 0:
        keys_arr[n_unique] = val
        counts_arr[n_unique] = na_count

    keys = np.asarray(keys_arr)[:n_unique + na_add].copy()
    counts = np.asarray(counts_arr)
    return keys, counts, na_count


def duplicated_float32(const float[:] values, object keep="first",
                       const uint8_t[:] mask=None):
    """
    Return boolean array indicating duplicated values (NaN-aware).

    Parameters
    ----------
    values : ndarray[float32]
        Array to check for duplicates
    keep : {'first', 'last', False}, default 'first'
        - 'first': Mark duplicates as True except for first occurrence
        - 'last': Mark duplicates as True except for last occurrence
        - False: Mark all duplicates as True
    mask : ndarray[uint8], optional
        If provided, mask[i] == 1 indicates that values[i] is NA/missing.
        NA values are treated as a single group for duplicate detection.
        Mask takes precedence over NaN check in values.

    Returns
    -------
    ndarray[bool]
        Boolean array indicating duplicates
    """
    cdef:
        Py_ssize_t i, n = len(values), first_na = -1
        SwissTable[float32_t, size_t] table
        float val
        size_t idx
        int ret = 0
        np_uint8_t[::1] result
        bint keep_first = keep == "first"
        bint keep_last = keep == "last"
        bint keep_none = keep is False
        bint uses_mask = mask is not None
        bint seen_na = False
        bint seen_multiple_na = False

    if not (keep_first or keep_last or keep_none):
        raise ValueError('keep must be either "first", "last" or False')

    table.reserve(<size_t>n)

    if not uses_mask and n > 0 and values.strides[0] == sizeof(float):
        result = np.empty(n, dtype=np.uint8)
        with nogil:
            if keep_first or keep_last:
                ret = table.duplicated_keep_batch(&values[0], n, keep_first, &result[0])
            else:  # keep == False
                ret = table.duplicated_false_batch(&values[0], n, &result[0])
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

        return np.asarray(result).view(np.bool_)

    result = np.zeros(n, dtype=np.uint8)

    if keep_first:
        # Use insert_key_only: no value storage needed
        with nogil:
            for i in range(n):
                # When mask is provided, use mask only. Otherwise check NaN.
                if uses_mask and mask[i]:
                    if seen_na:
                        result[i] = 1  # Duplicate NA
                    else:
                        seen_na = True
                        result[i] = 0  # First NA
                else:
                    val = values[i]
                    ret = table.insert_key_only(val)
                    if ret == 0:
                        result[i] = 1  # Duplicate (NaN handled by hash table)
                    elif ret == -1:
                        break
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    elif keep_last:
        # Use insert_key_only: no value storage needed
        with nogil:
            for i in range(n - 1, -1, -1):
                # When mask is provided, use mask only. Otherwise check NaN.
                if uses_mask and mask[i]:
                    if seen_na:
                        result[i] = 1  # Duplicate NA
                    else:
                        seen_na = True
                        result[i] = 0  # First NA (from end)
                else:
                    val = values[i]
                    ret = table.insert_key_only(val)
                    if ret == 0:
                        result[i] = 1  # Duplicate (NaN handled by hash table)
                    elif ret == -1:
                        break
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    else:  # keep == False
        # Need insert_if_absent: must store index to mark first occurrence
        with nogil:
            for i in range(n):
                # When mask is provided, use mask only. Otherwise check NaN.
                if uses_mask and mask[i]:
                    # Handle NA values - mark all as duplicates if multiple
                    if not seen_na:
                        first_na = i
                        seen_na = True
                        result[i] = 0
                    elif not seen_multiple_na:
                        result[i] = 1
                        result[first_na] = 1
                        seen_multiple_na = True
                    else:
                        result[i] = 1
                else:
                    val = values[i]
                    ret = table.insert_if_absent(val, <size_t>i, &idx)  # noqa: E501
                    if ret == 0:
                        result[idx] = 1
                        result[i] = 1
                    elif ret == -1:
                        break
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    return np.asarray(result).view(np.bool_)


# =============================================================================
# Complex64 Map
# =============================================================================
cdef class SwissComplex64Map(HashTable):
    """
    Swiss Table hash map: complex64 -> size_t

    NaN handling: Complex numbers with NaN components are considered equal.
    """

    def __cinit__(self, size_t size_hint=0):
        self.uses_mask = False
        if size_hint > 0:
            self.table.reserve(size_hint)

    cdef swiss_complex64_t _to_c_complex(self, object key):
        cdef swiss_complex64_t c_key
        c_key.real = key.real
        c_key.imag = key.imag
        return c_key

    def insert(self, key, size_t value):
        cdef swiss_complex64_t c_key = self._to_c_complex(key)
        cdef int ret = self.table.insert(c_key, value)
        return ret == 1

    def get(self, key, default=None):
        cdef swiss_complex64_t c_key = self._to_c_complex(key)
        cdef size_t val
        if self.table.get(c_key, &val):
            return val
        return default

    def __getitem__(self, key):
        cdef swiss_complex64_t c_key = self._to_c_complex(key)
        cdef size_t val
        if self.table.get(c_key, &val):
            return val
        raise KeyError(key)

    def __setitem__(self, key, size_t value):
        cdef swiss_complex64_t c_key = self._to_c_complex(key)
        cdef int ret = self.table.insert(c_key, value)

    cpdef get_item(self, complex64_t key):
        cdef swiss_complex64_t c_key
        c_key.real = key.real
        c_key.imag = key.imag
        cdef size_t val
        if self.table.get(c_key, &val):
            return val
        raise KeyError(key)

    cpdef set_item(self, complex64_t key, Py_ssize_t value):
        cdef swiss_complex64_t c_key
        c_key.real = key.real
        c_key.imag = key.imag
        cdef int ret = self.table.insert(c_key, value)

    cpdef get_na(self):
        """Extracts the position of na_value from the hashtable.

        Returns
        -------
        The position of the last na value.
        """
        raise NotImplementedError("not implement for get_na")

    cpdef set_na(self, Py_ssize_t val):
        # Caller is responsible for checking for pd.NA
        raise NotImplementedError("not implement for set_na")

    def __contains__(self, key):
        cdef swiss_complex64_t c_key = self._to_c_complex(key)
        return self.table.contains(c_key)

    def __len__(self):
        return self.table.size()

    def __iter__(self):
        cdef size_t index = 0
        cdef swiss_complex64_t c_key
        while index < self.table.capacity():
            index = self.table.iter_next(index)
            if index < self.table.capacity():
                c_key = self.table.key_at(index)
                yield complex(c_key.real, c_key.imag)
                index += 1

    def items(self):
        cdef size_t index = 0
        cdef swiss_complex64_t c_key
        while index < self.table.capacity():
            index = self.table.iter_next(index)
            if index < self.table.capacity():
                c_key = self.table.key_at(index)
                yield (complex(c_key.real, c_key.imag),
                       self.table.val_at(index))
                index += 1

    @property
    def size(self):
        return self.table.size()

    @property
    def capacity(self):
        return self.table.capacity()

    def __repr__(self):
        return f"SwissComplex64Map(size={self.size}, capacity={self.capacity})"

    def map_locations(self, const float complex[:] values, const uint8_t[:] mask = None) -> None:
        """Build table mapping values to their array positions."""
        cdef:
            Py_ssize_t i, n = len(values)
            swiss_complex64_t c_key
            int ret = 0

        if self.uses_mask and mask is not None:
            raise NotImplementedError("mask not implement for swisstable")  # pragma: no cover

        with nogil:
            if n > 0 and values.strides[0] == sizeof(swiss_complex64_t):
                ret = self.table.map_locations(<swiss_complex64_t*>&values[0], <size_t>n)
            else:
                for i in range(n):
                    c_key.real = values[i].real
                    c_key.imag = values[i].imag
                    ret = self.table.insert(c_key, <size_t>i)
                    if ret == -1:
                        break

        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    def lookup(self, const float complex[:] values) -> ndarray:
        """Look up array values in table, returning positions (-1 if not found)."""
        cdef:
            Py_ssize_t i, n = len(values)
            swiss_complex64_t c_key
            size_t loc
            intp_t[::1] locs = np.empty(n, dtype=np.intp)

        with nogil:
            if n > 0 and values.strides[0] == sizeof(swiss_complex64_t):
                self.table.lookup_batch(<swiss_complex64_t*>&values[0], <size_t>n, &locs[0])
            else:
                for i in range(n):
                    c_key.real = values[i].real
                    c_key.imag = values[i].imag
                    if self.table.get(c_key, &loc):
                        locs[i] = <intp_t>loc
                    else:
                        locs[i] = -1

        return np.asarray(locs)

    def unique(self, const float complex[:] values, bint return_inverse=False,
               object mask=None):
        """Calculate unique values and optionally their inverse mapping."""
        cdef:
            Py_ssize_t i, n = len(values), count = 0
            swiss_complex64_t c_key
            size_t idx
            int ret = 0
            float complex[::1] uniques_arr = np.empty(n, dtype=np.complex64)
            intp_t[::1] labels
            bint uses_mask = mask is not None
            bint seen_na = False
            const uint8_t[:] mask_view
            np_uint8_t[::1] result_mask_arr

        if not uses_mask and n > 0 and values.strides[0] == sizeof(swiss_complex64_t):
            if return_inverse:
                labels = np.empty(n, dtype=np.intp)
                with nogil:
                    count = self.table.unique_with_inverse(
                        <swiss_complex64_t*>&values[0], <size_t>n,
                        <swiss_complex64_t*>&uniques_arr[0], &labels[0]
                    )
                if count == -1:
                    raise MemoryError("Failed to insert into Swiss table")
                uniques = np.asarray(uniques_arr)[:count].copy()
                return uniques, np.asarray(labels)
            else:
                with nogil:
                    count = self.table.unique_batch(<swiss_complex64_t*>&values[0], <size_t>n,
                        <swiss_complex64_t*>&uniques_arr[0])
                if count == -1:
                    raise MemoryError("Failed to insert into Swiss table")
                uniques = np.asarray(uniques_arr)[:count].copy()
                return uniques

        if uses_mask:
            mask_view = mask if mask.dtype == np.uint8 else mask.view(np.uint8)
            result_mask_arr = np.empty(n, dtype=np.uint8)

        if return_inverse:
            labels = np.empty(n, dtype=np.intp)
            with nogil:
                for i in range(n):
                    if uses_mask and mask_view[i]:
                        if not seen_na:
                            uniques_arr[count] = values[i]
                            result_mask_arr[count] = 1
                            labels[i] = count
                            seen_na = True
                            count += 1
                        else:
                            for idx in range(<size_t>count):
                                if result_mask_arr[idx] == 1:
                                    labels[i] = <intp_t>idx
                                    break
                    else:
                        c_key.real = values[i].real
                        c_key.imag = values[i].imag
                        ret = self.table.insert_if_absent(c_key, <size_t>count, &idx)  # noqa: E501
                        if ret == 1:
                            uniques_arr[count] = values[i]
                            if uses_mask:
                                result_mask_arr[count] = 0
                            labels[i] = count
                            count += 1
                        elif ret == 0:
                            labels[i] = <intp_t>idx
                        else:
                            break

            if ret == -1:
                raise MemoryError("Failed to insert into Swiss table")

            uniques = np.asarray(uniques_arr)[:count].copy()
            if uses_mask:
                result_mask = np.asarray(result_mask_arr)[:count].copy()
                return uniques, np.asarray(labels), result_mask
            return uniques, np.asarray(labels)
        else:
            with nogil:
                for i in range(n):
                    if uses_mask and mask_view[i]:
                        if not seen_na:
                            uniques_arr[count] = values[i]
                            result_mask_arr[count] = 1
                            seen_na = True
                            count += 1
                    else:
                        c_key.real = values[i].real
                        c_key.imag = values[i].imag
                        ret = self.table.insert_key_only(c_key)
                        if ret == 1:
                            uniques_arr[count] = values[i]
                            if uses_mask:
                                result_mask_arr[count] = 0
                            count += 1
                        elif ret == -1:
                            break

            if ret == -1:
                raise MemoryError("Failed to insert into Swiss table")

            uniques = np.asarray(uniques_arr)[:count].copy()
            if uses_mask:
                return uniques, np.asarray(result_mask_arr)[:count].copy()
            return uniques

    def factorize(self, const float complex[:] values, Py_ssize_t na_sentinel=-1,
                  object na_value=None, object mask=None, bint ignore_na=True):
        """Calculate unique values and labels for categorical encoding (NaN-aware)."""
        cdef:
            Py_ssize_t i, n = len(values), count = 0, na_idx = -1
            swiss_complex64_t c_key
            size_t idx
            int ret = 0
            float complex[::1] uniques_arr = np.empty(n, dtype=np.complex64)
            intp_t[::1] labels = np.empty(n, dtype=np.intp)
            bint uses_mask = mask is not None
            bint seen_na = False
            const uint8_t[:] mask_view
            float real_part, imag_part
            swiss_complex64_t na_val = na_value if na_value is not None else NaNTraits[swiss_complex64_t].NaN()

        # Fast path: no mask - use C batch function for maximum performance
        if not uses_mask and n > 0 and values.strides[0] == sizeof(swiss_complex64_t):
            with nogil:
                if ignore_na:
                    count = self.table.factorize_batch(
                        <swiss_complex64_t*>&values[0], <size_t>n,
                        <swiss_complex64_t*>&uniques_arr[0], &labels[0],
                        na_sentinel, na_val
                    )
                else:
                    count = self.table.unique_with_inverse(
                        <swiss_complex64_t*>&values[0], <size_t>n,
                        <swiss_complex64_t*>&uniques_arr[0], &labels[0]
                    )
            if count == -1:
                raise MemoryError("Failed to insert into Swiss table")
            uniques = np.asarray(uniques_arr)[:count].copy()
            return uniques, np.asarray(labels)

        if uses_mask:
            mask_view = mask if mask.dtype == np.uint8 else mask.view(np.uint8)

        # Pre-reserve capacity to avoid resize checks in hot loop
        self.table.reserve(<size_t>n)

        with nogil:
            for i in range(n):
                c_key.real = values[i].real
                c_key.imag = values[i].imag
                # Check for NA: mask takes precedence, then NaN in components
                if uses_mask and mask_view[i]:
                    if ignore_na:
                        labels[i] = na_sentinel
                    else:
                        if not seen_na:
                            uniques_arr[count] = values[i]
                            na_idx = count
                            labels[i] = count
                            seen_na = True
                            count += 1
                        else:
                            labels[i] = na_idx
                elif ignore_na and NaNTraits[swiss_complex64_t].AreEqual(c_key, na_val):
                    labels[i] = na_sentinel
                else:
                    if not self.table.get(c_key, &idx):
                        ret = self.table.insert(c_key, <size_t>count)
                        if ret == -1:
                            break
                        uniques_arr[count] = values[i]
                        labels[i] = count
                        count += 1
                    else:
                        labels[i] = <intp_t>idx

        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

        uniques = np.asarray(uniques_arr)[:count].copy()
        return uniques, np.asarray(labels)

    def get_labels(self, const float complex[:] values, ndarray uniques, bint uniques_view_exits,
        Py_ssize_t count_prior=0, Py_ssize_t na_sentinel=-1,
        object na_value=None, object mask=None):
        """Calculate unique values and labels for categorical encoding."""
        cdef:
            Py_ssize_t i, n = len(values), count = count_prior
            swiss_complex64_t c_key
            size_t idx
            int ret = 0
            float complex[::1] uniques_arr = uniques
            intp_t[::1] labels = np.empty(n, dtype=np.intp)
            bint uses_na_value = na_value is not None
            bint uses_mask = mask is not None
            bint seen_na = False
            const uint8_t[:] mask_view
            swiss_complex64_t na_val = na_value if uses_na_value else NaNTraits[swiss_complex64_t].NaN()

        # Fast path: no mask - use C batch function for maximum performance
        if not uses_mask and (count_prior + n <= uniques_arr.shape[0]) and (
            n > 0 and values.strides[0] == sizeof(swiss_complex64_t)):
            with nogil:
                count = self.table.factorize_batch(
                    <swiss_complex64_t*>&values[0], <size_t>n,
                    <swiss_complex64_t*>&uniques_arr[0], &labels[0],
                    na_sentinel, na_val, count_prior
                )
            if count == -1:
                raise MemoryError("Failed to insert into Swiss table")
            return count, np.asarray(labels)

        if uses_mask:
            mask_view = mask if mask.dtype == np.uint8 else mask.view(np.uint8)

        with nogil:
            for i in range(n):
                c_key.real = values[i].real
                c_key.imag = values[i].imag
                if (uses_mask and mask_view[i]) or NaNTraits[swiss_complex64_t].AreEqual(c_key, na_val):
                    labels[i] = na_sentinel
                else:
                    if not self.table.get(c_key, &idx):
                        ret = self.table.insert(c_key, <size_t>count)
                        if ret == -1:
                            break
                        if count >= uniques_arr.shape[0]:
                            with gil:
                                uniques.resize(count * 2, refcheck=uniques_view_exits)
                                uniques_arr = uniques
                        uniques_arr[count] = values[i]
                        labels[i] = count
                        count += 1
                    else:
                        labels[i] = <intp_t>idx

        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

        return count, np.asarray(labels)


def ismember_complex64(
    const float complex[:] arr, const float complex[:] values
) -> ndarray:
    """Return boolean array indicating membership of arr elements in values."""
    cdef:
        Py_ssize_t n_arr = len(arr), n_values = len(values)
        SwissTable[swiss_complex64_t, size_t] table
        int ret = 0
        np_uint8_t[::1] result_view

    # Build the lookup table from values
    table.reserve(<size_t>n_values)

    # Pre-allocate result array
    result = np.empty(n_arr, dtype=np.uint8)
    result_view = result

    # Build set from values (key-only, no value storage needed)
    with nogil:
        ret = table.build_set(<swiss_complex64_t*>&values[0], <size_t>n_values)
    if ret == -1:
        raise MemoryError("Failed to insert into Swiss table")

    # Check membership using batch contains
    with nogil:
        table.contains_batch(<swiss_complex64_t*>&arr[0], <size_t>n_arr, &result_view[0])

    return result.view(np.bool_)


def value_count_complex64(const float complex[:] values, bint dropna=True, const uint8_t[:] mask=None):
    """Count occurrences of each unique value (NaN-aware)."""
    cdef:
        Py_ssize_t i, n = len(values), n_unique = 0, na_count = 0, na_add = 0
        SwissTable[swiss_complex64_t, size_t] table
        swiss_complex64_t val
        size_t idx
        int ret = 0
        float complex[::1] keys_arr = np.empty(n, dtype=np.complex64)
        size_t[::1] indices_arr = np.empty(n, dtype=np.uintp)
        int64_t[::1] counts_arr
        bint uses_mask = mask is not None
        bint isna_entry = False

    table.reserve(<size_t>n)

    # Use batch value_count which stores indices for direct access
    if not uses_mask and n > 0 and values.strides[0] == sizeof(swiss_complex64_t):
        with nogil:
            n_unique = table.value_count_batch(
                <swiss_complex64_t*>&values[0], <size_t>n,
                <swiss_complex64_t*>&keys_arr[0], &indices_arr[0], dropna
            )
    else:
        with nogil:
            for i in range(n):
                val.real = values[i].real
                val.imag = values[i].imag
                if uses_mask:
                    if mask[i]:
                        if not dropna:
                            na_count += 1
                        continue
                else:
                    if dropna and NaNTraits[swiss_complex64_t].IsNaN(val):
                        continue

                ret = table.increment(val, &idx)
                if ret == 1:
                    keys_arr[n_unique] = values[i]
                    indices_arr[n_unique] = idx
                    n_unique += 1
                elif ret == -1:
                    break

    if ret == -1:
        raise MemoryError("Failed to insert into Swiss table")

    if na_count > 0:
        na_add = 1

    # Extract counts using stored indices (direct access, no lookup)
    counts_arr = np.empty(n_unique + na_add, dtype=np.int64)
    with nogil:
        for i in range(n_unique):
            counts_arr[i] = <int64_t>table.vals()[indices_arr[i]]

    if na_count > 0:
        keys_arr[n_unique].real = val.real
        keys_arr[n_unique].imag = val.imag
        counts_arr[n_unique] = na_count

    keys = np.asarray(keys_arr)[:n_unique + na_add].copy()
    counts = np.asarray(counts_arr)
    return keys, counts, na_count


def duplicated_complex64(const float complex[:] values, object keep="first",
                         const uint8_t[:] mask=None):
    """
    Return boolean array indicating duplicated values (NaN-aware).

    Parameters
    ----------
    values : ndarray[complex64]
        Array to check for duplicates
    keep : {'first', 'last', False}, default 'first'
        - 'first': Mark duplicates as True except for first occurrence
        - 'last': Mark duplicates as True except for last occurrence
        - False: Mark all duplicates as True
    mask : ndarray[uint8], optional
        If provided, mask[i] == 1 indicates that values[i] is NA/missing.

    Returns
    -------
    ndarray[bool]
        Boolean array indicating duplicates
    """
    cdef:
        Py_ssize_t i, n = len(values), first_na = -1
        SwissTable[swiss_complex64_t, size_t] table
        swiss_complex64_t c_key
        size_t idx
        int ret = 0
        np_uint8_t[::1] result
        bint keep_first = keep == "first"
        bint keep_last = keep == "last"
        bint _keep_none = keep is False
        bint uses_mask = mask is not None
        bint seen_na = False
        bint seen_multiple_na = False

    if keep not in ("first", "last", False):
        raise ValueError('keep must be either "first", "last" or False')

    table.reserve(<size_t>n)

    if not uses_mask and n > 0 and values.strides[0] == sizeof(swiss_complex64_t):
        result = np.empty(n, dtype=np.uint8)
        with nogil:
            if keep_first or keep_last:
                ret = table.duplicated_keep_batch(<swiss_complex64_t*>&values[0], n, keep_first, &result[0])
            else:  # keep == False
                ret = table.duplicated_false_batch(<swiss_complex64_t*>&values[0], n, &result[0])
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

        return np.asarray(result).view(np.bool_)

    result = np.zeros(n, dtype=np.uint8)

    if keep_first:
        with nogil:
            for i in range(n):
                if uses_mask and mask[i]:
                    if seen_na:
                        result[i] = 1
                    else:
                        seen_na = True
                        result[i] = 0
                else:
                    c_key.real = values[i].real
                    c_key.imag = values[i].imag
                    ret = table.insert_key_only(c_key)
                    if ret == 0:
                        result[i] = 1
                    elif ret == -1:
                        break
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    elif keep_last:
        with nogil:
            for i in range(n - 1, -1, -1):
                if uses_mask and mask[i]:
                    if seen_na:
                        result[i] = 1
                    else:
                        seen_na = True
                        result[i] = 0
                else:
                    c_key.real = values[i].real
                    c_key.imag = values[i].imag
                    ret = table.insert_key_only(c_key)
                    if ret == 0:
                        result[i] = 1
                    elif ret == -1:
                        break
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    else:  # keep == False
        with nogil:
            for i in range(n):
                if uses_mask and mask[i]:
                    if not seen_na:
                        first_na = i
                        seen_na = True
                        result[i] = 0
                    elif not seen_multiple_na:
                        result[i] = 1
                        result[first_na] = 1
                        seen_multiple_na = True
                    else:
                        result[i] = 1
                else:
                    c_key.real = values[i].real
                    c_key.imag = values[i].imag
                    ret = table.insert_if_absent(c_key, <size_t>i, &idx)  # noqa: E501
                    if ret == 0:
                        result[idx] = 1
                        result[i] = 1
                    elif ret == -1:
                        break
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    return np.asarray(result).view(np.bool_)


# =============================================================================
# Complex128 Map
# =============================================================================
cdef class SwissComplex128Map(HashTable):
    """
    Swiss Table hash map: complex128 -> size_t

    NaN handling: Complex numbers with NaN components are considered equal.
    """

    def __cinit__(self, size_t size_hint=0):
        self.uses_mask = False
        if size_hint > 0:
            self.table.reserve(size_hint)

    cdef swiss_complex128_t _to_c_complex(self, object key):
        cdef swiss_complex128_t c_key
        c_key.real = key.real
        c_key.imag = key.imag
        return c_key

    def insert(self, key, size_t value):
        cdef swiss_complex128_t c_key = self._to_c_complex(key)
        cdef int ret = self.table.insert(c_key, value)
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")
        return ret == 1

    def get(self, key, default=None):
        cdef swiss_complex128_t c_key = self._to_c_complex(key)
        cdef size_t val
        if self.table.get(c_key, &val):
            return val
        return default

    def __getitem__(self, key):
        cdef swiss_complex128_t c_key = self._to_c_complex(key)
        cdef size_t val
        if self.table.get(c_key, &val):
            return val
        raise KeyError(key)

    def __setitem__(self, key, size_t value):
        cdef swiss_complex128_t c_key = self._to_c_complex(key)
        cdef int ret = self.table.insert(c_key, value)

    cpdef get_item(self, complex128_t key):
        cdef swiss_complex128_t c_key
        c_key.real = key.real
        c_key.imag = key.imag
        cdef size_t val
        if self.table.get(c_key, &val):
            return val
        raise KeyError(key)

    cpdef set_item(self, complex128_t key, Py_ssize_t value):
        cdef swiss_complex128_t c_key
        c_key.real = key.real
        c_key.imag = key.imag
        cdef int ret = self.table.insert(c_key, value)
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    cpdef get_na(self):
        """Extracts the position of na_value from the hashtable.

        Returns
        -------
        The position of the last na value.
        """
        raise NotImplementedError("not implement for get_na")

    cpdef set_na(self, Py_ssize_t val):
        # Caller is responsible for checking for pd.NA
        raise NotImplementedError("not implement for set_na")

    def __contains__(self, key):
        cdef swiss_complex128_t c_key = self._to_c_complex(key)
        return self.table.contains(c_key)

    def __len__(self):
        return self.table.size()

    def __iter__(self):
        cdef size_t index = 0
        cdef swiss_complex128_t c_key
        while index < self.table.capacity():
            index = self.table.iter_next(index)
            if index < self.table.capacity():
                c_key = self.table.key_at(index)
                yield complex(c_key.real, c_key.imag)
                index += 1

    def items(self):
        cdef size_t index = 0
        cdef swiss_complex128_t c_key
        while index < self.table.capacity():
            index = self.table.iter_next(index)
            if index < self.table.capacity():
                c_key = self.table.key_at(index)
                yield (complex(c_key.real, c_key.imag),
                       self.table.val_at(index))
                index += 1

    @property
    def size(self):
        return self.table.size()

    @property
    def capacity(self):
        return self.table.capacity()

    def __repr__(self):
        return f"SwissComplex128Map(size={self.size}, capacity={self.capacity})"

    def map_locations(self, const double complex[:] values, const uint8_t[:] mask = None) -> None:
        """Build table mapping values to their array positions."""
        cdef:
            Py_ssize_t i, n = len(values)
            swiss_complex128_t c_key
            int ret = 0

        if self.uses_mask and mask is not None:
            raise NotImplementedError("mask not implement for swisstable")  # pragma: no cover

        with nogil:
            if n > 0 and values.strides[0] == sizeof(swiss_complex128_t):
                ret = self.table.map_locations(<swiss_complex128_t*>&values[0], <size_t>n)
            else:
                for i in range(n):
                    c_key.real = values[i].real
                    c_key.imag = values[i].imag
                    ret = self.table.insert(c_key, <size_t>i)
                    if ret == -1:
                        break

        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    def lookup(self, const double complex[:] values) -> ndarray:
        """Look up array values in table, returning positions (-1 if not found)."""
        cdef:
            Py_ssize_t i, n = len(values)
            swiss_complex128_t c_key
            size_t loc
            intp_t[::1] locs = np.empty(n, dtype=np.intp)

        with nogil:
            if n > 0 and values.strides[0] == sizeof(swiss_complex128_t):
                self.table.lookup_batch(<swiss_complex128_t*>&values[0], <size_t>n, &locs[0])
            else:
                for i in range(n):
                    c_key.real = values[i].real
                    c_key.imag = values[i].imag
                    if self.table.get(c_key, &loc):
                        locs[i] = <intp_t>loc
                    else:
                        locs[i] = -1

        return np.asarray(locs)

    def unique(self, const double complex[:] values, bint return_inverse=False,
               object mask=None):
        """Calculate unique values and optionally their inverse mapping."""
        cdef:
            Py_ssize_t i, n = len(values), count = 0
            swiss_complex128_t c_key
            size_t idx
            int ret = 0
            double complex[::1] uniques_arr = np.empty(n, dtype=np.complex128)
            intp_t[::1] labels
            bint uses_mask = mask is not None
            bint seen_na = False
            const uint8_t[:] mask_view
            np_uint8_t[::1] result_mask_arr

        if not uses_mask and n > 0 and values.strides[0] == sizeof(swiss_complex128_t):
            if return_inverse:
                labels = np.empty(n, dtype=np.intp)
                with nogil:
                    count = self.table.unique_with_inverse(
                        <swiss_complex128_t*>&values[0], <size_t>n,
                        <swiss_complex128_t*>&uniques_arr[0], &labels[0]
                    )
                if count == -1:
                    raise MemoryError("Failed to insert into Swiss table")
                uniques = np.asarray(uniques_arr)[:count].copy()
                return uniques, np.asarray(labels)
            else:
                with nogil:
                    count = self.table.unique_batch(<swiss_complex128_t*>&values[0], <size_t>n,
                        <swiss_complex128_t*>&uniques_arr[0])
                if count == -1:
                    raise MemoryError("Failed to insert into Swiss table")
                uniques = np.asarray(uniques_arr)[:count].copy()
                return uniques

        if uses_mask:
            mask_view = mask if mask.dtype == np.uint8 else mask.view(np.uint8)
            result_mask_arr = np.empty(n, dtype=np.uint8)

        if return_inverse:
            labels = np.empty(n, dtype=np.intp)
            with nogil:
                for i in range(n):
                    if uses_mask and mask_view[i]:
                        if not seen_na:
                            uniques_arr[count] = values[i]
                            result_mask_arr[count] = 1
                            labels[i] = count
                            seen_na = True
                            count += 1
                        else:
                            for idx in range(<size_t>count):
                                if result_mask_arr[idx] == 1:
                                    labels[i] = <intp_t>idx
                                    break
                    else:
                        c_key.real = values[i].real
                        c_key.imag = values[i].imag
                        ret = self.table.insert_if_absent(c_key, <size_t>count, &idx)  # noqa: E501
                        if ret == 1:
                            uniques_arr[count] = values[i]
                            if uses_mask:
                                result_mask_arr[count] = 0
                            labels[i] = count
                            count += 1
                        elif ret == 0:
                            labels[i] = <intp_t>idx
                        else:
                            break

            if ret == -1:
                raise MemoryError("Failed to insert into Swiss table")

            uniques = np.asarray(uniques_arr)[:count].copy()
            if uses_mask:
                result_mask = np.asarray(result_mask_arr)[:count].copy()
                return uniques, np.asarray(labels), result_mask
            return uniques, np.asarray(labels)
        else:
            with nogil:
                for i in range(n):
                    if uses_mask and mask_view[i]:
                        if not seen_na:
                            uniques_arr[count] = values[i]
                            result_mask_arr[count] = 1
                            seen_na = True
                            count += 1
                    else:
                        c_key.real = values[i].real
                        c_key.imag = values[i].imag
                        ret = self.table.insert_key_only(c_key)
                        if ret == 1:
                            uniques_arr[count] = values[i]
                            if uses_mask:
                                result_mask_arr[count] = 0
                            count += 1
                        elif ret == -1:
                            break

            if ret == -1:
                raise MemoryError("Failed to insert into Swiss table")

            uniques = np.asarray(uniques_arr)[:count].copy()
            if uses_mask:
                return uniques, np.asarray(result_mask_arr)[:count].copy()
            return uniques

    def factorize(self, const double complex[:] values, Py_ssize_t na_sentinel=-1,
                  object na_value=None, object mask=None, bint ignore_na=True):
        """Calculate unique values and labels for categorical encoding (NaN-aware)."""
        cdef:
            Py_ssize_t i, n = len(values), count = 0, na_idx = -1
            swiss_complex128_t c_key
            size_t idx
            int ret = 0
            double complex[::1] uniques_arr = np.empty(n, dtype=np.complex128)
            intp_t[::1] labels = np.empty(n, dtype=np.intp)
            bint uses_mask = mask is not None
            bint seen_na = False
            const uint8_t[:] mask_view
            double real_part, imag_part
            swiss_complex128_t na_val = na_value if na_value is not None else NaNTraits[swiss_complex128_t].NaN()

        # Fast path: no mask - use C batch function for maximum performance
        if not uses_mask and n > 0 and values.strides[0] == sizeof(swiss_complex128_t):
            with nogil:
                if ignore_na:
                    count = self.table.factorize_batch(
                        <swiss_complex128_t*>&values[0], <size_t>n,
                        <swiss_complex128_t*>&uniques_arr[0], &labels[0],
                        na_sentinel, na_val
                    )
                else:
                    count = self.table.unique_with_inverse(
                        <swiss_complex128_t*>&values[0], <size_t>n,
                        <swiss_complex128_t*>&uniques_arr[0], &labels[0]
                    )
            if count == -1:
                raise MemoryError("Failed to insert into Swiss table")
            uniques = np.asarray(uniques_arr)[:count].copy()
            return uniques, np.asarray(labels)

        if uses_mask:
            mask_view = mask if mask.dtype == np.uint8 else mask.view(np.uint8)

        # Pre-reserve capacity to avoid resize checks in hot loop
        self.table.reserve(<size_t>n)

        with nogil:
            for i in range(n):
                c_key.real = values[i].real
                c_key.imag = values[i].imag
                # Check for NA: mask takes precedence, then NaN in components
                if uses_mask and mask_view[i]:
                    if ignore_na:
                        labels[i] = na_sentinel
                    else:
                        if not seen_na:
                            uniques_arr[count] = values[i]
                            na_idx = count
                            labels[i] = count
                            seen_na = True
                            count += 1
                        else:
                            labels[i] = na_idx
                elif ignore_na and NaNTraits[swiss_complex128_t].AreEqual(c_key, na_val):
                    labels[i] = na_sentinel
                else:
                    if not self.table.get(c_key, &idx):
                        ret = self.table.insert(c_key, <size_t>count)
                        if ret == -1:
                            break
                        uniques_arr[count] = values[i]
                        labels[i] = count
                        count += 1
                    else:
                        labels[i] = <intp_t>idx

        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

        uniques = np.asarray(uniques_arr)[:count].copy()
        return uniques, np.asarray(labels)

    def get_labels(self, const double complex[:] values, ndarray uniques, bint uniques_view_exits,
        Py_ssize_t count_prior=0, Py_ssize_t na_sentinel=-1,
        object na_value=None, object mask=None):
        """Calculate unique values and labels for categorical encoding."""
        cdef:
            Py_ssize_t i, n = len(values), count = count_prior
            swiss_complex128_t c_key
            size_t idx
            int ret = 0
            double complex[::1] uniques_arr = uniques
            intp_t[::1] labels = np.empty(n, dtype=np.intp)
            bint uses_na_value = na_value is not None
            bint uses_mask = mask is not None
            bint seen_na = False
            const uint8_t[:] mask_view
            swiss_complex128_t na_val = na_value if uses_na_value else NaNTraits[swiss_complex128_t].NaN()

        # Fast path: no mask - use C batch function for maximum performance
        if not uses_mask and (count_prior + n <= uniques_arr.shape[0]) and (
            n > 0 and values.strides[0] == sizeof(swiss_complex128_t)):
            with nogil:
                count = self.table.factorize_batch(
                    <swiss_complex128_t*>&values[0], <size_t>n,
                    <swiss_complex128_t*>&uniques_arr[0], &labels[0],
                    na_sentinel, na_val, count_prior
                )
            if count == -1:
                raise MemoryError("Failed to insert into Swiss table")
            return count, np.asarray(labels)

        if uses_mask:
            mask_view = mask if mask.dtype == np.uint8 else mask.view(np.uint8)

        with nogil:
            for i in range(n):
                c_key.real = values[i].real
                c_key.imag = values[i].imag
                if (uses_mask and mask_view[i]) or NaNTraits[swiss_complex128_t].AreEqual(c_key, na_val):
                    labels[i] = na_sentinel
                else:
                    if not self.table.get(c_key, &idx):
                        ret = self.table.insert(c_key, <size_t>count)
                        if ret == -1:
                            break
                        if count >= uniques_arr.shape[0]:
                            with gil:
                                uniques.resize(count * 2, refcheck=uniques_view_exits)
                                uniques_arr = uniques
                        uniques_arr[count] = values[i]
                        labels[i] = count
                        count += 1
                    else:
                        labels[i] = <intp_t>idx

        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

        return count, np.asarray(labels)


def ismember_complex128(
    const double complex[:] arr, const double complex[:] values
) -> ndarray:
    """Return boolean array indicating membership of arr elements in values."""
    cdef:
        Py_ssize_t n_arr = len(arr), n_values = len(values)
        SwissTable[swiss_complex128_t, size_t] table
        int ret = 0
        np_uint8_t[::1] result_view

    # Build the lookup table from values
    table.reserve(<size_t>n_values)

    # Pre-allocate result array
    result = np.empty(n_arr, dtype=np.uint8)
    result_view = result

    # Build set from values (key-only, no value storage needed)
    with nogil:
        ret = table.build_set(<swiss_complex128_t*>&values[0], <size_t>n_values)
    if ret == -1:
        raise MemoryError("Failed to insert into Swiss table")

    # Check membership using batch contains
    with nogil:
        table.contains_batch(<swiss_complex128_t*>&arr[0], <size_t>n_arr, &result_view[0])

    return result.view(np.bool_)


def value_count_complex128(const double complex[:] values, bint dropna=True, const uint8_t[:] mask=None):
    """Count occurrences of each unique value (NaN-aware)."""
    cdef:
        Py_ssize_t i, n = len(values), n_unique = 0, na_count = 0, na_add = 0
        SwissTable[swiss_complex128_t, size_t] table
        swiss_complex128_t val
        size_t idx
        int ret = 0
        double complex[::1] keys_arr = np.empty(n, dtype=np.complex128)
        size_t[::1] indices_arr = np.empty(n, dtype=np.uintp)
        int64_t[::1] counts_arr
        bint uses_mask = mask is not None
        bint isna_entry = False

    table.reserve(<size_t>n)

    # Use batch value_count which stores indices for direct access
    if not uses_mask and n > 0 and values.strides[0] == sizeof(swiss_complex128_t):
        with nogil:
            n_unique = table.value_count_batch(
                <swiss_complex128_t*>&values[0], <size_t>n,
                <swiss_complex128_t*>&keys_arr[0], &indices_arr[0], dropna
            )
    else:
        with nogil:
            for i in range(n):
                val.real = values[i].real
                val.imag = values[i].imag
                if uses_mask:
                    if mask[i]:
                        if not dropna:
                            na_count += 1
                        continue
                else:
                    if dropna and NaNTraits[swiss_complex128_t].IsNaN(val):
                        continue

                ret = table.increment(val, &idx)
                if ret == 1:
                    keys_arr[n_unique] = values[i]
                    indices_arr[n_unique] = idx
                    n_unique += 1
                elif ret == -1:
                    break

    if ret == -1:
        raise MemoryError("Failed to insert into Swiss table")

    if na_count > 0:
        na_add = 1

    # Extract counts using stored indices (direct access, no lookup)
    counts_arr = np.empty(n_unique + na_add, dtype=np.int64)
    with nogil:
        for i in range(n_unique):
            counts_arr[i] = <int64_t>table.vals()[indices_arr[i]]

    if na_count > 0:
        keys_arr[n_unique].real = val.real
        keys_arr[n_unique].imag = val.imag
        counts_arr[n_unique] = na_count

    keys = np.asarray(keys_arr)[:n_unique + na_add].copy()
    counts = np.asarray(counts_arr)
    return keys, counts, na_count


def duplicated_complex128(const double complex[:] values, object keep="first",
                          const uint8_t[:] mask=None):
    """
    Return boolean array indicating duplicated values (NaN-aware).

    Parameters
    ----------
    values : ndarray[complex128]
        Array to check for duplicates
    keep : {'first', 'last', False}, default 'first'
        - 'first': Mark duplicates as True except for first occurrence
        - 'last': Mark duplicates as True except for last occurrence
        - False: Mark all duplicates as True
    mask : ndarray[uint8], optional
        If provided, mask[i] == 1 indicates that values[i] is NA/missing.

    Returns
    -------
    ndarray[bool]
        Boolean array indicating duplicates
    """
    cdef:
        Py_ssize_t i, n = len(values), first_na = -1
        SwissTable[swiss_complex128_t, size_t] table
        swiss_complex128_t c_key
        size_t idx
        int ret = 0
        np_uint8_t[::1] result
        bint keep_first = keep == "first"
        bint keep_last = keep == "last"
        bint _keep_none = keep is False
        bint uses_mask = mask is not None
        bint seen_na = False
        bint seen_multiple_na = False

    if keep not in ("first", "last", False):
        raise ValueError('keep must be either "first", "last" or False')

    table.reserve(<size_t>n)

    if not uses_mask and n > 0 and values.strides[0] == sizeof(swiss_complex128_t):
        result = np.empty(n, dtype=np.uint8)
        with nogil:
            if keep_first or keep_last:
                ret = table.duplicated_keep_batch(<swiss_complex128_t*>&values[0], n, keep_first, &result[0])
            else:  # keep == False
                ret = table.duplicated_false_batch(<swiss_complex128_t*>&values[0], n, &result[0])
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

        return np.asarray(result).view(np.bool_)

    result = np.zeros(n, dtype=np.uint8)

    if keep_first:
        with nogil:
            for i in range(n):
                if uses_mask and mask[i]:
                    if seen_na:
                        result[i] = 1
                    else:
                        seen_na = True
                        result[i] = 0
                else:
                    c_key.real = values[i].real
                    c_key.imag = values[i].imag
                    ret = table.insert_key_only(c_key)
                    if ret == 0:
                        result[i] = 1
                    elif ret == -1:
                        break
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    elif keep_last:
        with nogil:
            for i in range(n - 1, -1, -1):
                if uses_mask and mask[i]:
                    if seen_na:
                        result[i] = 1
                    else:
                        seen_na = True
                        result[i] = 0
                else:
                    c_key.real = values[i].real
                    c_key.imag = values[i].imag
                    ret = table.insert_key_only(c_key)
                    if ret == 0:
                        result[i] = 1
                    elif ret == -1:
                        break
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    else:  # keep == False
        with nogil:
            for i in range(n):
                if uses_mask and mask[i]:
                    if not seen_na:
                        first_na = i
                        seen_na = True
                        result[i] = 0
                    elif not seen_multiple_na:
                        result[i] = 1
                        result[first_na] = 1
                        seen_multiple_na = True
                    else:
                        result[i] = 1
                else:
                    c_key.real = values[i].real
                    c_key.imag = values[i].imag
                    ret = table.insert_if_absent(c_key, <size_t>i, &idx)  # noqa: E501
                    if ret == 0:
                        result[idx] = 1
                        result[i] = 1
                    elif ret == -1:
                        break
        if ret == -1:
            raise MemoryError("Failed to insert into Swiss table")

    return np.asarray(result).view(np.bool_)
