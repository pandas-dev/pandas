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
from cpython.object cimport PyObject
from numpy cimport (
    intp_t,
    complex64_t,
    complex128_t,
)
ctypedef intp_t c_ssize_t

# External declarations for all Swiss Table types
cdef extern from "swisstable/swisstable_class.hpp" namespace "pandas::swisstable" nogil:
    cdef cppclass NaNTraits[K]:
        @staticmethod
        K NaN()
        @staticmethod
        bint IsNaN(K key)
        @staticmethod
        bint AreEqual(K a, K b)

    cdef cppclass SwissTable[K, V]:
        SwissTable()

        size_t size() const
        size_t capacity() const
        K *keys()
        V *vals()

        size_t find(K key) const
        int insert(K key, size_t val)
        bint get(K key, size_t *val_out) const
        bint contains(K key) const

        size_t iter_next(size_t index) const
        K key_at(size_t index) const
        size_t val_at(size_t index) const

        int insert_key_only(K key)
        int build_set(const K *keys, size_t n)
        int insert_if_absent(K key, size_t new_val, size_t *val_out)

        int increment(K key, size_t *idx_out)
        bint reserve(size_t want)

        # -------- Batch operations --------
        int map_locations(const K *keys, size_t n)

        void lookup_batch(const K *keys, size_t n, c_ssize_t *locs) const

        void contains_batch(const K *keys, size_t n, uint8_t *result) const

        int64_t unique_batch(const K *keys, size_t n, K *uniques_out)

        int64_t unique_with_inverse(
            const K *keys,
            size_t n,
            K *uniques_out,
            c_ssize_t *labels_out
        )

        int64_t unique_with_inverse(
            const K *keys,
            size_t n,
            K *uniques_out,
            c_ssize_t *labels_out,
            int64_t count_prior
        )

        int64_t factorize_batch(
            const K *keys,
            size_t n,
            K *uniques_out,
            c_ssize_t *labels_out,
            int64_t na_sentinel,
            K nan_key
        )

        int64_t factorize_batch(
            const K *keys,
            size_t n,
            K *uniques_out,
            c_ssize_t *labels_out,
            int64_t na_sentinel,
            K nan_key,
            int64_t count_prior
        )

        int64_t value_count_batch(
            const K *keys,
            size_t n,
            K *keys_out,
            size_t *indices_out,
            bint dropna
        )

        int duplicated_keep_batch(
            const K *keys,
            size_t n,
            bint keep_first,
            uint8_t *result
        )

        int duplicated_false_batch(
            const K *keys,
            size_t n,
            uint8_t *result
        )

    # Complex types
    ctypedef struct swiss_complex64_t:
        float real
        float imag

    ctypedef struct swiss_complex128_t:
        double real
        double imag


# prototypes for sharing
from pandas._libs.hashtable cimport HashTable

cdef class SwissUInt64Map(HashTable):
    cdef SwissTable[uint64_t, size_t] table
    cdef bint uses_mask

    cpdef get_item(self, uint64_t val)
    cpdef set_item(self, uint64_t key, Py_ssize_t val)
    cpdef get_na(self)
    cpdef set_na(self, Py_ssize_t val)

cdef class SwissInt64Map(HashTable):
    cdef SwissTable[int64_t, size_t] table
    cdef bint uses_mask

    cpdef get_item(self, int64_t val)
    cpdef set_item(self, int64_t key, Py_ssize_t val)
    cpdef get_na(self)
    cpdef set_na(self, Py_ssize_t val)

cdef class SwissUInt32Map(HashTable):
    cdef SwissTable[uint32_t, size_t] table
    cdef bint uses_mask

    cpdef get_item(self, uint32_t val)
    cpdef set_item(self, uint32_t key, Py_ssize_t val)
    cpdef get_na(self)
    cpdef set_na(self, Py_ssize_t val)

cdef class SwissInt32Map(HashTable):
    cdef SwissTable[int32_t, size_t] table
    cdef bint uses_mask

    cpdef get_item(self, int32_t val)
    cpdef set_item(self, int32_t key, Py_ssize_t val)
    cpdef get_na(self)
    cpdef set_na(self, Py_ssize_t val)

cdef class SwissUInt16Map(HashTable):
    cdef SwissTable[uint16_t, size_t] table
    cdef bint uses_mask

    cpdef get_item(self, uint16_t val)
    cpdef set_item(self, uint16_t key, Py_ssize_t val)
    cpdef get_na(self)
    cpdef set_na(self, Py_ssize_t val)

cdef class SwissInt16Map(HashTable):
    cdef SwissTable[int16_t, size_t] table
    cdef bint uses_mask

    cpdef get_item(self, int16_t val)
    cpdef set_item(self, int16_t key, Py_ssize_t val)
    cpdef get_na(self)
    cpdef set_na(self, Py_ssize_t val)

cdef class SwissUInt8Map(HashTable):
    cdef SwissTable[uint8_t, size_t] table
    cdef bint uses_mask

    cpdef get_item(self, uint8_t val)
    cpdef set_item(self, uint8_t key, Py_ssize_t val)
    cpdef get_na(self)
    cpdef set_na(self, Py_ssize_t val)

cdef class SwissInt8Map(HashTable):
    cdef SwissTable[int8_t, size_t] table
    cdef bint uses_mask

    cpdef get_item(self, int8_t val)
    cpdef set_item(self, int8_t key, Py_ssize_t val)
    cpdef get_na(self)
    cpdef set_na(self, Py_ssize_t val)

cdef class SwissFloat64Map(HashTable):
    cdef SwissTable[double, size_t] table
    cdef bint uses_mask

    cpdef get_item(self, double val)
    cpdef set_item(self, double key, Py_ssize_t val)
    cpdef get_na(self)
    cpdef set_na(self, Py_ssize_t val)

cdef class SwissFloat32Map(HashTable):
    cdef SwissTable[float, size_t] table
    cdef bint uses_mask

    cpdef get_item(self, float val)
    cpdef set_item(self, float key, Py_ssize_t val)
    cpdef get_na(self)
    cpdef set_na(self, Py_ssize_t val)

cdef class SwissComplex64Map(HashTable):
    cdef SwissTable[swiss_complex64_t, size_t] table
    cdef bint uses_mask

    cpdef get_item(self, complex64_t val)
    cpdef set_item(self, complex64_t key, Py_ssize_t val)
    cpdef get_na(self)
    cpdef set_na(self, Py_ssize_t val)

    cdef swiss_complex64_t _to_c_complex(self, object key)

cdef class SwissComplex128Map(HashTable):
    cdef SwissTable[swiss_complex128_t, size_t] table
    cdef bint uses_mask

    cpdef get_item(self, complex128_t val)
    cpdef set_item(self, complex128_t key, Py_ssize_t val)
    cpdef get_na(self)
    cpdef set_na(self, Py_ssize_t val)

    cdef swiss_complex128_t _to_c_complex(self, object key)
