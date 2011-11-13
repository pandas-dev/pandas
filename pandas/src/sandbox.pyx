from numpy cimport *
import numpy as np

import_array()

cdef class ArrayCruncher:

    cdef:
        ndarray arr
        object f
        bint raw
        Py_ssize_t N, K

    def __init__(self, arr, f, axis=0, raw=True):
        self.arr = arr
        self.f = f
        self.raw = raw
        self.N, self.K = arr.shape

    def reduce(self):
        cdef:
            char* dummy_buf
            ndarray arr, result, chunk
            Py_ssize_t i, increment
            flatiter it

        if not self.arr.flags.c_contiguous:
            arr = self.arr.copy('C')
        else:
            arr = self.arr

        increment = self.K * self.arr.dtype.itemsize
        chunk = np.empty(self.K, dtype=arr.dtype)
        result = np.empty(self.N, dtype=arr.dtype)
        it = <flatiter> PyArray_IterNew(result)

        dummy_buf = chunk.data
        chunk.data = arr.data

        for i in range(self.N):
            PyArray_SETITEM(result, PyArray_ITER_DATA(it), self.f(chunk))
            chunk.data = chunk.data + increment
            PyArray_ITER_NEXT(it)

        # so we don't free the wrong memory
        chunk.data = dummy_buf

        return result
