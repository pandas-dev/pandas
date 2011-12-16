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
            assert(not (isinstance(res, list) and len(res) == len(self.dummy)))

            result = np.empty(self.nresults, dtype='O')
            # if hasattr(res, 'dtype'):
            #     result = np.empty(self.nresults, dtype=res.dtype)
            # else:
            #     result = np.empty(self.nresults, dtype='O')
            result[0] = res
        except Exception:
            raise ValueError('function does not reduce')
        return result

cdef class SeriesGrouper:
    '''
    Performs generic grouping operation while avoiding ndarray construction
    overhead
    '''
    cdef:
        Py_ssize_t nresults, ngroup
        object arr, dummy, f, labels, counts
        bint passed_dummy

    def __init__(self, object arr, object f, object labels, ngroups,
                 dummy=None):
        n = len(arr)

        assert(arr.ndim == 1)

        if not arr.flags.contiguous:
            arr = arr.copy()

        self.labels = labels
        self.f = f
        self.arr = arr

        self.dummy = self._check_dummy(dummy)
        self.passed_dummy = dummy is not None

        self.counts = np.zeros(ngroups, dtype='i4')

        self.ngroups = ngroups

    def _check_dummy(self, dummy=None):
        if dummy is None:
            dummy = np.empty(0, dtype=self.arr.dtype)
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
            ndarray[int32_t] labels, counts
            Py_ssize_t i, group_size, n, lab
            flatiter it
            npy_intp *shape
            object res
            bint initialized = 0
            tuple args
            object kwds

        labels = self.labels
        counts = self.counts

        arr = self.arr
        chunk = self.dummy

        dummy_buf = chunk.data
        chunk.data = arr.data

        shape = chunk.shape
        group_size = 0
        n = len(arr)

        args = cpython.PyTuple_New(1)
        kwds = {}
        cpython.PyTuple_SET_ITEM(args, 0, chunk)
        cpython.Py_INCREF(chunk)

        try:
            for i in range(n):
                group_size += 1

                lab = labels[i]

                if i == n - 1 or lab != labels[i + 1]:
                    chunk.shape[0] = group_size

                    res = cpython.PyObject_Call(self.f, args, kwds)

                    # res = self.f(chunk)
                    if not initialized:
                        result = self._get_result_array(res)
                        it = <flatiter> PyArray_IterNew(result)
                        initialized = 1

                    PyArray_ITER_GOTO1D(it, lab)
                    PyArray_SETITEM(result, PyArray_ITER_DATA(it), res)
                    counts[lab] = group_size

                    chunk.data = chunk.data + group_size
                    group_size = 0
        except:
            raise
        finally:
            # so we don't free the wrong memory
            chunk.shape[0] = 0
            chunk.data = dummy_buf

        if result.dtype == np.object_:
            result = maybe_convert_objects(result)

        return result

    def _get_result_array(self, object res):
        try:
            assert(not isinstance(res, np.ndarray))
            assert(not (isinstance(res, list) and len(res) == len(self.dummy)))

            result = np.empty(self.ngroups, dtype='O')
        except Exception:
            raise ValueError('function does not reduce')
        return result

def reduce(arr, f, axis=0, dummy=None):
    reducer = Reducer(arr, f, axis=axis, dummy=dummy)
    return reducer.get_result()
