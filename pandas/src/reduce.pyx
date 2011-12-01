from numpy cimport *
import numpy as np

cdef class Reducer:
    '''
    Performs generic reduction operation on a C or Fortran-contiguous ndarray
    while avoiding ndarray construction overhead
    '''
    cdef:
        Py_ssize_t increment, chunksize, nresults
        object arr, dummy, f

    def __init__(self, object arr, object f, axis=1, dummy=None):
        n, k = arr.shape

        if axis == 0:
            if not arr.flags.f_contiguous:
                arr = arr.copy('F')

            self.nresults = k
            self.chunksize = n
            self.increment = n * arr.dtype.itemsize
        else:
            if not arr.flags.c_contiguous:
                arr = arr.copy('C')

            self.nresults = n
            self.chunksize = k
            self.increment = k * arr.dtype.itemsize

        self.f = f
        self.arr = arr
        self.dummy = self._check_dummy(dummy)

    def _check_dummy(self, dummy=None):
        if dummy is None:
            dummy = np.empty(self.chunksize, dtype=self.arr.dtype)
        else:
            if dummy.dtype != self.arr.dtype:
                raise ValueError('Dummy array must be same dtype')
            if len(dummy) != self.chunksize:
                raise ValueError('Dummy array must be length %d' %
                                 self.chunksize)

        return dummy

    def get_result(self):
        cdef:
            char* dummy_buf
            ndarray arr, result, chunk
            Py_ssize_t i
            flatiter it
            object res

        arr = self.arr
        chunk = self.dummy

        dummy_buf = chunk.data
        chunk.data = arr.data
        try:
            for i in range(self.nresults):
                res = self.f(chunk)
                if i == 0:
                    result = self._get_result_array(res)
                    it = <flatiter> PyArray_IterNew(result)

                PyArray_SETITEM(result, PyArray_ITER_DATA(it), res)
                chunk.data = chunk.data + self.increment
                PyArray_ITER_NEXT(it)
        finally:
            # so we don't free the wrong memory
            chunk.data = dummy_buf
        if result.dtype == np.object_:
            result = maybe_convert_objects(result)
        return result

    def _get_result_array(self, object res):
        try:
            assert(not isinstance(res, np.ndarray))
            result = np.empty(self.nresults, dtype='O')
            # if hasattr(res, 'dtype'):
            #     result = np.empty(self.nresults, dtype=res.dtype)
            # else:
            #     result = np.empty(self.nresults, dtype='O')
            result[0] = res
        except Exception:
            raise ValueError('function does not reduce')
        return result

def reduce(arr, f, axis=0, dummy=None):
    reducer = Reducer(arr, f, axis=axis, dummy=dummy)
    return reducer.get_result()
