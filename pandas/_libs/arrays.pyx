"""
Cython implementations for internal ExtensionArrays.
"""
cimport cython

import numpy as np

cimport numpy as cnp
from cpython cimport PyErr_Clear
from cpython.slice cimport PySlice_Unpack
from libc.stdlib cimport (
    free,
    malloc,
)
from libc.string cimport memcpy
from numpy cimport (
    int8_t,
    int64_t,
    ndarray,
    uint8_t,
)

from pandas._libs.lib import is_null_slice

cnp.import_array()

cdef extern from "pandas/vendored/nanoarrow.h":
    struct ArrowBuffer:
        uint8_t* data
        int64_t size_bytes
        int64_t capacity_bytes

    struct ArrowBitmap:
        ArrowBuffer buffer
        int64_t size_bits

    void ArrowBitmapInit(ArrowBitmap*)
    void ArrowBitmapReserve(ArrowBitmap*, int64_t)
    void ArrowBitmapAppendUnsafe(ArrowBitmap*, uint8_t, int64_t)
    void ArrowBitmapAppendInt8Unsafe(ArrowBitmap*, const int8_t *, int64_t)
    void ArrowBitmapReset(ArrowBitmap*)
    void ArrowBitsUnpackInt8(const uint8_t*, int64_t, int64_t, int8_t*)
    int8_t ArrowBitGet(const uint8_t*, int64_t)
    void ArrowBitSetTo(uint8_t*, int64_t, uint8_t)
    void ArrowBitsSetTo(uint8_t*, int64_t, int64_t, uint8_t)
    int64_t ArrowBitCountSet(const uint8_t*, int64_t, int64_t)
    void ArrowBitmapReset(ArrowBitmap*)

cdef extern from "pandas/bitmask_algorithms.h":
    void ConcatenateBitmapData(const ArrowBitmap**, size_t, ArrowBitmap*)
    bint BitmapAny(const ArrowBitmap*)
    bint BitmapAll(const ArrowBitmap*)
    bint BitmapOr(const ArrowBitmap*, const ArrowBitmap*, ArrowBitmap*)
    bint BitmapXor(const ArrowBitmap*, const ArrowBitmap*, ArrowBitmap*)
    bint BitmapAnd(const ArrowBitmap*, const ArrowBitmap*, ArrowBitmap*)
    bint BitmapInvert(const ArrowBitmap*, ArrowBitmap*)
    bint BitmapTake(const ArrowBitmap*, const int64_t*, size_t, ArrowBitmap*)
    bint BitmapPutFromBufferMask(ArrowBitmap*, const uint8_t*, size_t, uint8_t)


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


cdef class BitmaskArray:
    cdef:
        ArrowBitmap bitmap
        bint buffer_owner  # set when parent is None, but gives C-level access
        # NumPy compatibility
        cdef Py_ssize_t ndim
        cdef Py_ssize_t[2] shape
        cdef Py_ssize_t[2] strides
        # Buffer protocol support
        int n_consumers
        uint8_t* memview_buf
    cdef public:
        object parent  # assignments gives RC to ensure proper buffer lifecycle

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void init_from_ndarray(self, const uint8_t[::1] arr):
        cdef ArrowBitmap bitmap
        # As long as we have a 1D arr argument we can use .shape[0] to avoid
        # a call to Python via .size
        cdef int64_t nobs = arr.shape[0]
        ArrowBitmapInit(&bitmap)
        ArrowBitmapReserve(&bitmap, nobs)
        ArrowBitmapAppendInt8Unsafe(&bitmap, <const int8_t*>&arr[0], nobs)
        self.bitmap = bitmap
        self.buffer_owner = True

    cdef void init_from_bitmaskarray(self, BitmaskArray bma):
        self.bitmap = bma.bitmap
        self.buffer_owner = False
        self.ndim = bma.ndim
        self.shape[0] = bma.shape[0]
        self.strides[0] = bma.strides[0]
        if self.ndim == 2:
            self.shape[1] = bma.shape[1]
            self.strides[1] = bma.strides[1]

    def __cinit__(self):
        cdef BitmaskArray self_ = self
        self.parent = False
        self_.n_consumers = 0
        self_.memview_buf = NULL

    def __init__(self, data):
        cdef BitmaskArray self_ = self
        if isinstance(data, np.ndarray):
            if not data.flags["C_CONTIGUOUS"]:
                data = np.ascontiguousarray(data)

            self.init_from_ndarray(data.ravel())
            self.parent = None
            self_.ndim = data.ndim
            self_.shape[0] = data.shape[0]
            self_.strides[0] = data.strides[0]
            if (data.ndim == 2):
                self_.shape[1] = data.shape[1]
                self_.strides[1] = data.strides[1]
        elif isinstance(data, type(self)):
            self.init_from_bitmaskarray(data)
            self.parent = data
        else:
            raise TypeError("Unsupported argument to BitmaskArray constructor")

    def __dealloc__(self):
        if self.buffer_owner:
            ArrowBitmapReset(&self.bitmap)

    @staticmethod
    cdef BitmaskArray copy_from_bitmaskarray(BitmaskArray old_bma):
        """
        Constructs a new BitmaskArray from a bitmap pointer. Copies data
        and manages the subsequenty lifecycle of the bitmap.
        """
        # Bypass __init__ calls
        cdef BitmaskArray bma = BitmaskArray.__new__(BitmaskArray)
        cdef uint8_t* buf
        cdef ArrowBitmap bitmap
        # TODO: this leaks a bit into the internals of the nanoarrow bitmap
        # We may want to upstream a BitmapCopy function instead
        ArrowBitmapInit(&bitmap)
        buf = <uint8_t*>malloc(old_bma.bitmap.buffer.size_bytes)
        memcpy(buf, old_bma.bitmap.buffer.data, old_bma.bitmap.buffer.size_bytes)
        bitmap.buffer.size_bytes = old_bma.bitmap.buffer.size_bytes
        bitmap.buffer.capacity_bytes = old_bma.bitmap.buffer.capacity_bytes
        bitmap.size_bits = old_bma.bitmap.size_bits
        bitmap.buffer.data = buf

        bma.bitmap = bitmap
        bma.buffer_owner = True
        bma.ndim = old_bma.ndim
        bma.shape = old_bma.shape
        bma.strides = old_bma.strides
        bma.parent = False

        return bma

    def __len__(self):
        return self.bitmap.size_bits

    def __repr__(self):
        if self.parent:
            par = object.__repr__(self.parent)
        else:
            par = None

        data = self.bytes

        return (
            f"{object.__repr__(self)}\nparent: {par}\ndata: {data}\n"
        )

    @cython.wraparound(False)
    @cython.boundscheck(False)
    @staticmethod
    cdef BitmaskArray c_concatenate(list objs):
        cdef Py_ssize_t i
        cdef int64_t total_bits = 0
        cdef BitmaskArray current_bma
        cdef Py_ssize_t nbitmaps = len(objs)

        cdef BitmaskArray first_bma = objs[0]
        cdef int expected_ndim = first_bma.ndim
        cdef Py_ssize_t expected_stride0 = first_bma.strides[0]
        cdef Py_ssize_t expected_shape1, expected_stride1
        if expected_ndim == 2:
            expected_stride1 = first_bma.strides[1]
            expected_shape1 = first_bma.shape[1]

        cdef Py_ssize_t dim0shape = 0

        cdef ArrowBitmap** bitmaps = <ArrowBitmap**>malloc(
            sizeof(ArrowBitmap*) * nbitmaps
        )

        for i in range(nbitmaps):
            current_bma = <BitmaskArray?>objs[i]
            if (
                    current_bma.ndim != expected_ndim
                    or current_bma.strides[0] != expected_stride0
                    or (
                        expected_ndim == 2 and (
                            current_bma.shape[1] != expected_shape1
                            or current_bma.strides[1] != expected_stride1
                        )
                    )
            ):
                free(bitmaps)
                raise NotImplementedError(
                    "BitmaskArray.concatenate does not support broadcasting"
                )
            total_bits += current_bma.bitmap.size_bits
            bitmaps[i] = &current_bma.bitmap
            dim0shape += current_bma.shape[0]

        # Bypass __init__ calls
        cdef BitmaskArray bma = BitmaskArray.__new__(BitmaskArray)
        cdef ArrowBitmap bitmap

        ArrowBitmapInit(&bitmap)
        ArrowBitmapReserve(&bitmap, total_bits)

        ConcatenateBitmapData(bitmaps, nbitmaps, &bitmap)
        free(bitmaps)

        bma.bitmap = bitmap
        bma.buffer_owner = True

        bma.ndim = expected_ndim
        bma.shape[0] = dim0shape  # only allowed because of axis=0 assumption
        bma.strides[0] = expected_stride0
        if expected_ndim == 2:
            bma.shape[1] = expected_shape1
            bma.strides[1] = expected_stride1

        bma.parent = None

        return bma

    @classmethod
    def concatenate(cls, objs, axis):
        if axis != 0:
            raise NotImplementedError(
                "BitmaskArray.concatenate only implemented for axis=0"
            )

        return BitmaskArray.c_concatenate(objs)

    def __setitem__(self, key, value):
        cdef const uint8_t[:] keymask
        cdef const uint8_t[:] arr1d
        cdef Py_ssize_t i = 0
        cdef Py_ssize_t ckey
        cdef bint cvalue
        cdef BitmaskArray self_ = self

        if isinstance(key, int):
            ckey = key
            cvalue = value
            if ckey >= 0 and ckey < self.bitmap.size_bits:
                ArrowBitSetTo(self.bitmap.buffer.data, ckey, cvalue)
                return

        if is_null_slice(key) and isinstance(value, (int, bool)):
            cvalue = value  # blindly assuming ints are 0 or 1
            ArrowBitsSetTo(
                self.bitmap.buffer.data,
                0,
                self.bitmap.size_bits,
                cvalue
            )
        elif (
                isinstance(key, np.ndarray)
                and key.dtype == bool
                and isinstance(value, (int, bool))
                and len(key) == len(self)
        ):
            keymask = key
            if BitmapPutFromBufferMask(
                    &self_.bitmap,
                    &keymask[0],
                    keymask.shape[0],
                    value
            ) != 0:
                raise ValueError("BitmaskArray.__setitem__ failed!")
        else:
            arr = self.to_numpy()
            arr[key] = value
            arr1d = arr.ravel()
            for i in range(arr1d.shape[0]):
                ArrowBitSetTo(self.bitmap.buffer.data, i, arr1d[i])

    def __getitem__(self, key):
        cdef Py_ssize_t ckey
        cdef Py_ssize_t start, stop, step
        cdef BitmaskArray bma
        cdef ArrowBitmap bitmap
        cdef int64_t nbytes, nbits
        cdef BitmaskArray self_ = self
        # to_numpy can be expensive, so try to avoid for simple cases
        if isinstance(key, int) and self.ndim == 1:
            ckey = key
            if ckey >= 0 and ckey < self.bitmap.size_bits:
                return bool(ArrowBitGet(self.bitmap.buffer.data, ckey))
        elif is_null_slice(key):
            return self
        elif isinstance(key, slice) and self.ndim == 1:
            # fastpath for slices that start at 0 and step 1 at a time
            # towards a positive number.
            # TODO: upstream generic ArrowBitsGet function in nanoarrow
            PySlice_Unpack(key, &start, &stop, &step)
            if start == 0 and stop > 0 and step == 1:
                if stop > self_.bitmap.size_bits:
                    nbits = self_.bitmap.size_bits
                else:
                    nbits = stop

                nbytes = (nbits + 7) // 8

                bma = BitmaskArray.__new__(BitmaskArray)
                ArrowBitmapInit(&bitmap)
                ArrowBitmapReserve(&bitmap, nbits)
                memcpy(bitmap.buffer.data, self_.bitmap.buffer.data, nbytes)
                bitmap.buffer.size_bytes = nbytes
                bitmap.size_bits = nbits

                bma.bitmap = bitmap
                bma.buffer_owner = True
                bma.ndim = self_.ndim
                bma.shape[0] = nbits
                bma.strides = self_.strides
                bma.parent = False

                return bma

        return self.to_numpy()[key]

    def __invert__(self):
        # note that this inverts the entire byte, even if the
        # bitmap only uses a few of the bits within that byte
        # the remaining bits of the byte are of undefined value
        # so be sure to only check bytes we need
        cdef BitmaskArray self_ = self
        cdef BitmaskArray bma = BitmaskArray.__new__(BitmaskArray)
        cdef ArrowBitmap bitmap

        ArrowBitmapInit(&bitmap)
        ArrowBitmapReserve(&bitmap, self_.bitmap.size_bits)

        BitmapInvert(&self_.bitmap, &bitmap)

        bma.bitmap = bitmap
        bma.buffer_owner = True
        bma.ndim = self_.ndim
        bma.shape = self_.shape
        bma.strides = self_.strides
        bma.parent = None

        return bma

    def __and__(self, other):
        cdef ndarray[uint8_t] result
        cdef BitmaskArray other_bma, self_ = self  # self_ required for Cython < 3
        cdef ArrowBitmap bitmap

        if isinstance(other, type(self)):
            # TODO: maybe should return Self here instead of ndarray
            other_bma = other
            if self_.bitmap.size_bits == 0:
                result = np.empty([], dtype=bool)
                if self_.ndim == 2:
                    return result.reshape(self_.shape[0], self_.shape[1])
                return result

            if self_.bitmap.size_bits != other_bma.bitmap.size_bits:
                raise ValueError("bitmaps are not equal size")

            ArrowBitmapInit(&bitmap)
            ArrowBitmapReserve(&bitmap, self_.bitmap.size_bits)
            BitmapAnd(&self_.bitmap, &other_bma.bitmap, &bitmap)

            result = np.empty(self_.bitmap.size_bits, dtype=bool)
            ArrowBitsUnpackInt8(
                bitmap.buffer.data,
                0,
                bitmap.size_bits,
                <int8_t*>&result[0]
            )
            ArrowBitmapReset(&bitmap)

            if self_.ndim == 2:
                return result.reshape(self_.shape[0], self_.shape[1])
            return result

        return self.to_numpy() & other

    def __or__(self, other):
        cdef ndarray[uint8_t] result
        cdef BitmaskArray other_bma, self_ = self  # self_ required for Cython < 3
        cdef ArrowBitmap bitmap

        if isinstance(other, type(self)):
            # TODO: maybe should return Self here instead of ndarray
            other_bma = other
            if self_.bitmap.size_bits == 0:
                result = np.empty([], dtype=bool)
                if self_.ndim == 2:
                    return result.reshape(self_.shape[0], self_.shape[1])
                return result

            if self_.bitmap.size_bits != other_bma.bitmap.size_bits:
                raise ValueError("bitmaps are not equal size")

            ArrowBitmapInit(&bitmap)
            ArrowBitmapReserve(&bitmap, self_.bitmap.size_bits)
            BitmapOr(&self_.bitmap, &other_bma.bitmap, &bitmap)

            result = np.empty(self_.bitmap.size_bits, dtype=bool)
            ArrowBitsUnpackInt8(
                bitmap.buffer.data,
                0,
                bitmap.size_bits,
                <int8_t*>&result[0]
            )
            ArrowBitmapReset(&bitmap)

            if self_.ndim == 2:
                return result.reshape(self_.shape[0], self_.shape[1])
            return result

        return self.to_numpy() | other

    def __xor__(self, other):
        cdef ndarray[uint8_t] result
        cdef BitmaskArray other_bma, self_ = self  # self_ required for Cython < 3
        cdef ArrowBitmap bitmap

        if isinstance(other, type(self)):
            # TODO: maybe should return Self here instead of ndarray
            other_bma = other
            if self_.bitmap.size_bits == 0:
                result = np.empty([], dtype=bool)
                if self_.ndim == 2:
                    return result.reshape(self_.shape[0], self_.shape[1])
                return result

            if self_.bitmap.size_bits != other_bma.bitmap.size_bits:
                raise ValueError("bitmaps are not equal size")

            ArrowBitmapInit(&bitmap)
            ArrowBitmapReserve(&bitmap, self_.bitmap.size_bits)
            BitmapXor(&self_.bitmap, &other_bma.bitmap, &bitmap)

            result = np.empty(self_.bitmap.size_bits, dtype=bool)
            ArrowBitsUnpackInt8(
                bitmap.buffer.data,
                0,
                bitmap.size_bits,
                <int8_t*>&result[0]
            )
            ArrowBitmapReset(&bitmap)
            if self_.ndim == 2:
                return result.reshape(self_.shape[0], self_.shape[1])
            return result

        return self.to_numpy() ^ other

    def __getstate__(self):
        cdef BitmaskArray self_ = self
        state = {
            "parent": self.parent,
            "ndim": self_.ndim,
            "shape0": self_.shape[0],
            "stride0": self_.strides[0],
            "n_consumers": self_.n_consumers,
            "buffer_owner": self_.buffer_owner,
            # Private ArrowBitmap attributes below
            "bitmap.buffer.size_bytes": self_.bitmap.buffer.size_bytes,
            "bitmap.buffer.capacity_bytes": self_.bitmap.buffer.capacity_bytes,
            "bitmap.size_bits": self_.bitmap.size_bits
        }

        if self_.ndim == 2:
            state["shape1"] = self_.shape[1]
            state["stride1"] = self_.strides[1]

        # memview should only exist when n_consumers > 0
        if self_.n_consumers > 0:
            memview_buf_data = bytearray(len(self))
            for i in range(len(self)):
                memview_buf_data[i] = self_.memview_buf[i]

            state["memview_buf_data"] = memview_buf_data

        # Only parents own data
        if self_.buffer_owner:
            bitmap_data = bytearray(self_.bitmap.buffer.size_bytes)
            for i in range(self_.bitmap.buffer.size_bytes):
                bitmap_data[i] = self_.bitmap.buffer.data[i]

            state["bitmap_data"] = bitmap_data

        return state

    def __setstate__(self, state):
        cdef ArrowBitmap bitmap
        cdef BitmaskArray self_ = self, other
        self.parent = state["parent"]
        self_.ndim = state["ndim"]
        self_.shape[0] = state["shape0"]
        self_.strides[0] = state["stride0"]
        self_.n_consumers = state["n_consumers"]
        self_.buffer_owner = state["buffer_owner"]

        nbytes = state["bitmap.buffer.size_bytes"]
        capacity_bytes = state["bitmap.buffer.capacity_bytes"]
        nbits = state["bitmap.size_bits"]

        if self_.ndim == 2:
            self_.shape[1] = state["shape1"]
            self_.strides[1] = state["stride1"]

        if self_.n_consumers > 0:
            self_.memview_buf = <uint8_t *>malloc(nbits)
            memview_buf_data = state["memview_buf_data"]
            for i in range(nbits):
                self_.memview_buf[i] = memview_buf_data[i]

        if not self_.buffer_owner:
            other = self.parent
            self_.bitmap = other.bitmap
            self_.bitmap.size_bits = nbits
            self_.bitmap.buffer.size_bytes = nbytes
            self_.bitmap.buffer.capacity_bytes = capacity_bytes
        else:
            ArrowBitmapInit(&bitmap)

            buf = <uint8_t*>malloc(nbytes)
            data = state["bitmap_data"]
            for i in range(nbytes):
                buf[i] = data[i]

            bitmap.buffer.data = buf
            bitmap.buffer.size_bytes = nbytes
            bitmap.buffer.capacity_bytes = nbytes
            bitmap.size_bits = nbits
            self_.bitmap = bitmap

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __iter__(self):
        cdef Py_ssize_t i
        cdef BitmaskArray self_ = self  # self_ required for Cython < 3
        for i in range(self_.bitmap.size_bits):
            yield bool(ArrowBitGet(self_.bitmap.buffer.data, i))

    def __getbuffer__(self, Py_buffer *buffer, int flags):
        cdef BitmaskArray self_ = self

        if self_.n_consumers == 0:
            self_.memview_buf = <uint8_t*>malloc(self_.bitmap.size_bits)
            ArrowBitsUnpackInt8(
                self_.bitmap.buffer.data,
                0,
                self_.bitmap.size_bits,
                <int8_t*>self_.memview_buf
            )

        buffer.buf = self_.memview_buf
        buffer.format = "?"
        buffer.internal = NULL
        buffer.itemsize = 1
        buffer.len = self_.bitmap.size_bits
        buffer.ndim = self_.ndim
        buffer.obj = self
        buffer.readonly = 1
        buffer.shape = self_.shape
        buffer.strides = self_.strides
        buffer.suboffsets = NULL

        self_.n_consumers += 1

    def __releasebuffer__(self, Py_buffer *buffer):
        cdef BitmaskArray self_ = self
        self_.n_consumers -= 1
        if self_.n_consumers == 0:
            free(self_.memview_buf)

    @property
    def size(self) -> int:
        return self.bitmap.size_bits

    @property
    def nbytes(self) -> int:
        return self.bitmap.buffer.size_bytes

    @property
    def bytes(self):
        cdef Py_ssize_t i, nbytes = self.bitmap.buffer.size_bytes
        arr_bytes = bytearray(nbytes)
        for i in range(nbytes):
            arr_bytes[i] = self.bitmap.buffer.data[i]

        return bytes(arr_bytes)

    @property
    def shape(self):
        """Strictly for NumPy compat in mask_ops"""
        cdef BitmaskArray self_ = self
        if self_.ndim == 1:
            return tuple((self_.shape[0],))
        return tuple((self_.shape[0], self_.shape[1]))

    @property
    def dtype(self):
        """Strictly for NumPy compat in mask_ops"""
        return np.dtype("bool")

    def any(self) -> bool:
        return BitmapAny(&self.bitmap)

    def all(self) -> bool:
        return BitmapAll(&self.bitmap)

    def sum(self) -> int:
        return ArrowBitCountSet(self.bitmap.buffer.data, 0, self.bitmap.size_bits)

    def take_1d(
        self,
        const int64_t[:] indices,
        const int axis=0,
    ):
        cdef BitmaskArray self_ = self
        cdef Py_ssize_t nindices = len(indices)
        if axis != 0:
            raise NotImplementedError(
                "BitmaskArray.take_1d only implemented for axis=0"
            )

        if nindices <= 0:
            raise NotImplementedError(
                "take_1d does not support empty takes"
            )

        cdef ArrowBitmap bitmap
        cdef BitmaskArray bma = BitmaskArray.__new__(BitmaskArray)

        ArrowBitmapInit(&bitmap)
        ArrowBitmapReserve(&bitmap, nindices)

        if BitmapTake(&self_.bitmap, &indices[0], nindices, &bitmap) != 0:
            ArrowBitmapReset(&bitmap)
            raise ValueError("take_1d does not support negative indexing")

        bma.bitmap = bitmap
        bma.buffer_owner = True

        bma.ndim = self_.ndim
        bma.shape[0] = indices.shape[0]
        bma.strides = self_.strides

        bma.parent = None

        return bma

    def copy(self):
        return BitmaskArray.copy_from_bitmaskarray(self)

    def to_numpy(self) -> ndarray:
        cdef BitmaskArray self_ = self
        cdef ndarray[uint8_t] result = np.empty(self_.bitmap.size_bits, dtype=bool)

        ArrowBitsUnpackInt8(
            self_.bitmap.buffer.data,
            0,
            self_.bitmap.size_bits,
            <int8_t*>cnp.PyArray_BYTES(result),
        )

        if self_.ndim == 2:
            return result.reshape(self_.shape[0], self_.shape[1])
        return result
