cimport cython
from cpython.mem cimport (
    PyMem_Free,
    PyMem_Malloc,
)
from cpython.ref cimport (
    Py_INCREF,
    PyObject,
)
from libc.stdlib cimport (
    free,
    malloc,
)

import numpy as np

cimport numpy as cnp
from numpy cimport ndarray

cnp.import_array()


from pandas._libs cimport util
from pandas._libs.dtypes cimport numeric_object_t
from pandas._libs.khash cimport (
    KHASH_TRACE_DOMAIN,
    are_equivalent_float32_t,
    are_equivalent_float64_t,
    are_equivalent_khcomplex64_t,
    are_equivalent_khcomplex128_t,
    kh_needed_n_buckets,
    kh_python_hash_equal,
    kh_python_hash_func,
    khiter_t,
)
from pandas._libs.missing cimport checknull


def get_hashtable_trace_domain():
    return KHASH_TRACE_DOMAIN


def object_hash(obj):
    return kh_python_hash_func(obj)


def objects_are_equal(a, b):
    return kh_python_hash_equal(a, b)


cdef int64_t NPY_NAT = util.get_nat()
SIZE_HINT_LIMIT = (1 << 20) + 7


cdef Py_ssize_t _INIT_VEC_CAP = 128

include "hashtable_class_helper.pxi"
include "hashtable_func_helper.pxi"


# map derived hash-map types onto basic hash-map types:
if np.dtype(np.intp) == np.dtype(np.int64):
    IntpHashTable = Int64HashTable
    unique_label_indices = _unique_label_indices_int64
elif np.dtype(np.intp) == np.dtype(np.int32):
    IntpHashTable = Int32HashTable
    unique_label_indices = _unique_label_indices_int32
else:
    raise ValueError(np.dtype(np.intp))


cdef class Factorizer:
    cdef readonly:
        Py_ssize_t count

    def __cinit__(self, size_hint: int):
        self.count = 0

    def get_count(self) -> int:
        return self.count

    def factorize(self, values, na_sentinel=-1, na_value=None, mask=None) -> np.ndarray:
        raise NotImplementedError


cdef class ObjectFactorizer(Factorizer):
    cdef public:
        PyObjectHashTable table
        ObjectVector uniques

    def __cinit__(self, size_hint: int):
        self.table = PyObjectHashTable(size_hint)
        self.uniques = ObjectVector()

    def factorize(
        self, ndarray[object] values, na_sentinel=-1, na_value=None, mask=None
    ) -> np.ndarray:
        """

        Returns
        -------
        np.ndarray[np.intp]

        Examples
        --------
        Factorize values with nans replaced by na_sentinel

        >>> fac = ObjectFactorizer(3)
        >>> fac.factorize(np.array([1,2,np.nan], dtype='O'), na_sentinel=20)
        array([ 0,  1, 20])
        """
        cdef:
            ndarray[intp_t] labels

        if mask is not None:
            raise NotImplementedError("mask not supported for ObjectFactorizer.")

        if self.uniques.external_view_exists:
            uniques = ObjectVector()
            uniques.extend(self.uniques.to_array())
            self.uniques = uniques
        labels = self.table.get_labels(values, self.uniques,
                                       self.count, na_sentinel, na_value)
        self.count = len(self.uniques)
        return labels

ctypedef fused htfunc_t:
    numeric_object_t
    complex128_t
    complex64_t


cpdef value_count(ndarray[htfunc_t] values, bint dropna, const uint8_t[:] mask=None):
    if htfunc_t is object:
        return value_count_object(values, dropna, mask=mask)

    elif htfunc_t is int8_t:
        return value_count_int8(values, dropna, mask=mask)
    elif htfunc_t is int16_t:
        return value_count_int16(values, dropna, mask=mask)
    elif htfunc_t is int32_t:
        return value_count_int32(values, dropna, mask=mask)
    elif htfunc_t is int64_t:
        return value_count_int64(values, dropna, mask=mask)

    elif htfunc_t is uint8_t:
        return value_count_uint8(values, dropna, mask=mask)
    elif htfunc_t is uint16_t:
        return value_count_uint16(values, dropna, mask=mask)
    elif htfunc_t is uint32_t:
        return value_count_uint32(values, dropna, mask=mask)
    elif htfunc_t is uint64_t:
        return value_count_uint64(values, dropna, mask=mask)

    elif htfunc_t is float64_t:
        return value_count_float64(values, dropna, mask=mask)
    elif htfunc_t is float32_t:
        return value_count_float32(values, dropna, mask=mask)

    elif htfunc_t is complex128_t:
        return value_count_complex128(values, dropna, mask=mask)
    elif htfunc_t is complex64_t:
        return value_count_complex64(values, dropna, mask=mask)

    else:
        raise TypeError(values.dtype)


cpdef duplicated(ndarray[htfunc_t] values,
                 object keep="first",
                 const uint8_t[:] mask=None):
    if htfunc_t is object:
        return duplicated_object(values, keep, mask=mask)

    elif htfunc_t is int8_t:
        return duplicated_int8(values, keep, mask=mask)
    elif htfunc_t is int16_t:
        return duplicated_int16(values, keep, mask=mask)
    elif htfunc_t is int32_t:
        return duplicated_int32(values, keep, mask=mask)
    elif htfunc_t is int64_t:
        return duplicated_int64(values, keep, mask=mask)

    elif htfunc_t is uint8_t:
        return duplicated_uint8(values, keep, mask=mask)
    elif htfunc_t is uint16_t:
        return duplicated_uint16(values, keep, mask=mask)
    elif htfunc_t is uint32_t:
        return duplicated_uint32(values, keep, mask=mask)
    elif htfunc_t is uint64_t:
        return duplicated_uint64(values, keep, mask=mask)

    elif htfunc_t is float64_t:
        return duplicated_float64(values, keep, mask=mask)
    elif htfunc_t is float32_t:
        return duplicated_float32(values, keep, mask=mask)

    elif htfunc_t is complex128_t:
        return duplicated_complex128(values, keep, mask=mask)
    elif htfunc_t is complex64_t:
        return duplicated_complex64(values, keep, mask=mask)

    else:
        raise TypeError(values.dtype)


cpdef ismember(ndarray[htfunc_t] arr, ndarray[htfunc_t] values):
    if htfunc_t is object:
        return ismember_object(arr, values)

    elif htfunc_t is int8_t:
        return ismember_int8(arr, values)
    elif htfunc_t is int16_t:
        return ismember_int16(arr, values)
    elif htfunc_t is int32_t:
        return ismember_int32(arr, values)
    elif htfunc_t is int64_t:
        return ismember_int64(arr, values)

    elif htfunc_t is uint8_t:
        return ismember_uint8(arr, values)
    elif htfunc_t is uint16_t:
        return ismember_uint16(arr, values)
    elif htfunc_t is uint32_t:
        return ismember_uint32(arr, values)
    elif htfunc_t is uint64_t:
        return ismember_uint64(arr, values)

    elif htfunc_t is float64_t:
        return ismember_float64(arr, values)
    elif htfunc_t is float32_t:
        return ismember_float32(arr, values)

    elif htfunc_t is complex128_t:
        return ismember_complex128(arr, values)
    elif htfunc_t is complex64_t:
        return ismember_complex64(arr, values)

    else:
        raise TypeError(values.dtype)


@cython.wraparound(False)
@cython.boundscheck(False)
def mode(ndarray[htfunc_t] values, bint dropna, const uint8_t[:] mask=None):
    # TODO(cython3): use const htfunct_t[:]

    cdef:
        ndarray[htfunc_t] keys
        ndarray[htfunc_t] modes
        ndarray[uint8_t] res_mask = None

        int64_t[::1] counts
        int64_t count, _, max_count = -1
        Py_ssize_t nkeys, k, na_counter, j = 0

    keys, counts, na_counter = value_count(values, dropna, mask=mask)
    nkeys = len(keys)

    modes = np.empty(nkeys, dtype=values.dtype)

    if htfunc_t is not object:
        with nogil:
            for k in range(nkeys):
                count = counts[k]
                if count == max_count:
                    j += 1
                elif count > max_count:
                    max_count = count
                    j = 0
                else:
                    continue

                modes[j] = keys[k]
    else:
        for k in range(nkeys):
            count = counts[k]
            if count == max_count:
                j += 1
            elif count > max_count:
                max_count = count
                j = 0
            else:
                continue

            modes[j] = keys[k]

    if na_counter > 0:
        res_mask = np.zeros(j+1, dtype=np.bool_)
        res_mask[j] = True
    return modes[:j + 1], res_mask
