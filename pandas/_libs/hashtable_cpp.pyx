import cython
import numpy as np

cimport numpy as cnp
from libc.stdint cimport uint32_t
from libc.string cimport memcpy
from libcpp.vector cimport vector

from pandas._libs.khash cimport kh_needed_n_buckets


cdef extern from "<functional>" namespace "std" nogil:
    cdef cppclass hash[T]:
        hash()
        size_t operator()

cdef extern from "pandas/vendored/klib/cpp/khash.hpp" namespace "klib" nogil:
    cdef cppclass KHash[T, Hash, Eq=*, khint_t=*]:
        T *keys
        KHash()
        # TODO: validate we don't need deconstructor
        # ~KHash()
        void exist(khint_t x)
        T &at(khint_t x)
        khint_t get(const T &)
        # TODO: make this khint_t
        # int resize(khint_t)
        int resize(uint32_t)
        khint_t put(const T &, int *)
        # void del(khint_t x)


# TODO: de-duplicate from hashtable.pyx
cdef uint32_t SIZE_HINT_LIMIT = (1 << 20) + 7


@cython.wraparound(False)
@cython.boundscheck(False)
def unique_label_indices(const cnp.npy_intp[:] labels) -> cnp.ndarray:
    """
    Indices of the first occurrences of the unique labels
    *excluding* -1. equivalent to:
        np.unique(labels, return_index=True)[1]
    """
    cdef:
        int ret = 0
        Py_ssize_t i, n = len(labels)
        KHash[cnp.npy_intp, hash[cnp.npy_intp]] *table = (
            new KHash[cnp.npy_intp, hash[cnp.npy_intp]]()
        )
        cnp.ndarray[cnp.npy_intp, ndim=1] arr
        vector[cnp.npy_intp] idx = vector[cnp.npy_intp]()

    table.resize(min(kh_needed_n_buckets(n), SIZE_HINT_LIMIT))

    with nogil:
        for i in range(n):
            table.put(labels[i], &ret)
            if ret != 0:
                # TODO: pandas has a custom resize operation but we
                # rely on C++ stdlib here - how different are they?
                idx.push_back(i)

    # TODO: must be a cleaner way to do this?
    # even arr.data = move(idx.data()) would be better but arr.data is readonly
    arr = np.empty(idx.size(), dtype=np.intp)
    memcpy(arr.data, idx.const_data(), idx.size() * sizeof(cnp.npy_intp))
    arr = arr[np.asarray(labels)[arr].argsort()]

    return arr[1:] if arr.size != 0 and labels[arr[0]] == -1 else arr
