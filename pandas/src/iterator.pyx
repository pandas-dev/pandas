cdef class RowIterator(object):
    cdef:
        ndarray arr, iterbuf
        Py_ssize_t N, K, itemsize
        char* buf

    def __init__(self, ndarray arr):
        self.arr = arr
        self.N, self.K = arr.shape
        self.itemsize = arr.dtype.itemsize
        self.iterbuf = np.empty(self.K, dtype=self.arr.dtype)
        self.buf = self.iterbuf.data

    def __del__(self):
        self.iterbuf.data = self.buf

    def __iter__(self):
        cdef:
            ndarray result = np.empty(self.K, dtype=self.arr.dtype)
            char* buf, arr_buf
            Py_ssize_t i, inc

        buf = result.data
        arr_buf = arr.data

        inc = self.itemsize * self.K

        for i from 0 <= i < self.N:
            result.data = arr_buf
            yield result
            arr_buf = arr_buf + inc

        result.data = buf
