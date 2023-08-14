"""
Cython implementations for internal ExtensionArrays.
"""
cimport cython

import numpy as np

cimport numpy as cnp
from cpython cimport PyErr_Clear
from numpy cimport (
    int8_t,
    int64_t,
    ndarray,
    uint8_t,
)

cnp.import_array()

from libc.stdlib cimport (
    free,
    malloc,
)


cdef extern from "pandas/vendored/nanoarrow.h":
    int8_t ArrowBitGet(const uint8_t*, int64_t)
    void ArrowBitSet(uint8_t*, int64_t)
    void ArrowBitClear(uint8_t*, int64_t)


@cython.freelist(16)
cdef class NDArrayBacked:
    """
    Implementing these methods in cython improves performance quite a bit.

    import pandas as pd

    from pandas._libs.arrays import NDArrayBacked as cls

    dti = pd.date_range("2016-01-01", periods=3)
    dta = dti._data
    arr = dta._ndarray

    obj = cls._simple_new(arr, arr.dtype)

    # for foo in [arr, dta, obj]: ...

    %timeit foo.copy()
    299 ns ± 30 ns per loop     # <-- arr underlying ndarray (for reference)
    530 ns ± 9.24 ns per loop   # <-- dta with cython NDArrayBacked
    1.66 µs ± 46.3 ns per loop  # <-- dta without cython NDArrayBacked
    328 ns ± 5.29 ns per loop   # <-- obj with NDArrayBacked.__cinit__
    371 ns ± 6.97 ns per loop   # <-- obj with NDArrayBacked._simple_new

    %timeit foo.T
    125 ns ± 6.27 ns per loop   # <-- arr underlying ndarray (for reference)
    226 ns ± 7.66 ns per loop   # <-- dta with cython NDArrayBacked
    911 ns ± 16.6 ns per loop   # <-- dta without cython NDArrayBacked
    215 ns ± 4.54 ns per loop   # <-- obj with NDArrayBacked._simple_new

    """
    # TODO: implement take in terms of cnp.PyArray_TakeFrom
    # TODO: implement concat_same_type in terms of cnp.PyArray_Concatenate

    # cdef:
    #    readonly ndarray _ndarray
    #    readonly object _dtype

    def __init__(self, ndarray values, object dtype):
        self._ndarray = values
        self._dtype = dtype

    @classmethod
    def _simple_new(cls, ndarray values, object dtype):
        cdef:
            NDArrayBacked obj
        obj = NDArrayBacked.__new__(cls)
        obj._ndarray = values
        obj._dtype = dtype
        return obj

    cpdef NDArrayBacked _from_backing_data(self, ndarray values):
        """
        Construct a new ExtensionArray `new_array` with `arr` as its _ndarray.

        This should round-trip:
            self == self._from_backing_data(self._ndarray)
        """
        # TODO: re-reuse simple_new if/when it can be cpdef
        cdef:
            NDArrayBacked obj
        obj = NDArrayBacked.__new__(type(self))
        obj._ndarray = values
        obj._dtype = self._dtype
        return obj

    cpdef __setstate__(self, state):
        if isinstance(state, dict):
            if "_data" in state:
                data = state.pop("_data")
            elif "_ndarray" in state:
                data = state.pop("_ndarray")
            else:
                raise ValueError  # pragma: no cover
            self._ndarray = data
            self._dtype = state.pop("_dtype")

            for key, val in state.items():
                setattr(self, key, val)
        elif isinstance(state, tuple):
            if len(state) != 3:
                if len(state) == 1 and isinstance(state[0], dict):
                    self.__setstate__(state[0])
                    return
                raise NotImplementedError(state)  # pragma: no cover

            data, dtype = state[:2]
            if isinstance(dtype, np.ndarray):
                dtype, data = data, dtype
            self._ndarray = data
            self._dtype = dtype

            if isinstance(state[2], dict):
                for key, val in state[2].items():
                    setattr(self, key, val)
            else:
                raise NotImplementedError(state)  # pragma: no cover
        else:
            raise NotImplementedError(state)  # pragma: no cover

    def __len__(self) -> int:
        return len(self._ndarray)

    @property
    def shape(self):
        # object cast bc _ndarray.shape is npy_intp*
        return (<object>(self._ndarray)).shape

    @property
    def ndim(self) -> int:
        return self._ndarray.ndim

    @property
    def size(self) -> int:
        # TODO(cython3): use self._ndarray.size
        return cnp.PyArray_SIZE(self._ndarray)

    @property
    def nbytes(self) -> int:
        return cnp.PyArray_NBYTES(self._ndarray)

    def copy(self, order="C"):
        cdef:
            cnp.NPY_ORDER order_code
            int success

        success = cnp.PyArray_OrderConverter(order, &order_code)
        if not success:
            # clear exception so that we don't get a SystemError
            PyErr_Clear()
            # same message used by numpy
            msg = f"order must be one of 'C', 'F', 'A', or 'K' (got '{order}')"
            raise ValueError(msg)

        res_values = cnp.PyArray_NewCopy(self._ndarray, order_code)
        return self._from_backing_data(res_values)

    def delete(self, loc, axis=0):
        res_values = np.delete(self._ndarray, loc, axis=axis)
        return self._from_backing_data(res_values)

    def swapaxes(self, axis1, axis2):
        res_values = cnp.PyArray_SwapAxes(self._ndarray, axis1, axis2)
        return self._from_backing_data(res_values)

    # TODO: pass NPY_MAXDIMS equiv to axis=None?
    def repeat(self, repeats, axis: int | np.integer = 0):
        if axis is None:
            axis = 0
        res_values = cnp.PyArray_Repeat(self._ndarray, repeats, <int>axis)
        return self._from_backing_data(res_values)

    def reshape(self, *args, **kwargs):
        res_values = self._ndarray.reshape(*args, **kwargs)
        return self._from_backing_data(res_values)

    def ravel(self, order="C"):
        # cnp.PyArray_OrderConverter(PyObject* obj, NPY_ORDER* order)
        # res_values = cnp.PyArray_Ravel(self._ndarray, order)
        res_values = self._ndarray.ravel(order)
        return self._from_backing_data(res_values)

    @property
    def T(self):
        res_values = self._ndarray.T
        return self._from_backing_data(res_values)

    def transpose(self, *axes):
        res_values = self._ndarray.transpose(*axes)
        return self._from_backing_data(res_values)

    @classmethod
    def _concat_same_type(cls, to_concat, axis=0):
        # NB: We are assuming at this point that dtypes all match
        new_values = [obj._ndarray for obj in to_concat]
        new_arr = cnp.PyArray_Concatenate(new_values, axis)
        return to_concat[0]._from_backing_data(new_arr)


def _unpickle_bitmaskarray(array):
    bma = BitMaskArray(array)
    return bma


cdef class BitMaskArray:
    cdef:
        Py_ssize_t array_size
        Py_ssize_t array_nbytes
        object array_shape
        uint8_t* validity_buffer
    cdef public:
        object parent

    cdef _setitem_with_integral(self, const int key, const uint8_t value):
        if value:
            ArrowBitSet(self.validity_buffer, key)
        else:
            ArrowBitClear(self.validity_buffer, key)

    cdef void init_from_ndarray(self, const uint8_t[:] arr):
        cdef Py_ssize_t i, arrlen
        self.array_size = arr.size
        self.array_nbytes = self.array_size // 8 + 1
        self.validity_buffer = <uint8_t *>malloc(self.array_nbytes)
        arrlen = len(arr)
        for i in range(arrlen):
            self._setitem_with_integral(i, arr[i])

    def __cinit__(self, data):
        self.parent = None
        if isinstance(data, np.ndarray):
            self.array_shape = data.shape
            self.init_from_ndarray(data.ravel())
        elif isinstance(data, type(self)):
            self.parent = data
            # other attributes are undefined when a parent exists
        else:
            raise TypeError("Unsupported argument to BitMaskArray constructor")

    def __dealloc__(self):
        if not self.parent:
            free(self.validity_buffer)

    def __setitem__(self, key, value):
        cdef int index = 0
        if self.parent is not None:
            self.parent[key] = value
        elif isinstance(key, int) and key >= 0:
            self._setitem_with_integral(key, bool(value))
        else:
            arr = self.to_numpy()
            arr[key] = value
            for val in arr.ravel():
                if val:
                    ArrowBitSet(self.validity_buffer, index)
                else:
                    ArrowBitClear(self.validity_buffer, index)
                index += 1

    def __getitem__(self, key):
        if self.parent is not None:
            return self.parent[key]
        elif isinstance(key, int) and key >= 0:
            return bool(ArrowBitGet(self.validity_buffer, key))
        else:
            return self.to_numpy()[key]

    def __invert__(self):
        if self.parent is not None:
            return ~self.parent
        return ~self.to_numpy()

    def __or__(self, other):
        if self.parent is not None:
            return self.parent | other
        elif isinstance(other, type(self)):
            return self.to_numpy() | other.to_numpy()
        else:
            return self.to_numpy() | other

    def __reduce__(self):
        object_state = (self.to_numpy(),)
        return (_unpickle_bitmaskarray, object_state)

    @property
    def nbytes(self) -> int:
        return self.array_nbytes

    def to_numpy(self) -> ndarray:
        if self.parent is not None:
            return self.parent.to_numpy()

        cdef int i
        cdef ndarray[uint8_t] result = np.empty(self.array_size, dtype=bool)
        for i in range(self.array_size):
            result[i] = self[i]

        return result.reshape(self.array_shape)
