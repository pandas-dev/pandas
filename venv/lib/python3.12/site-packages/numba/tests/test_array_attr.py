import numpy as np

import unittest
from numba.np.numpy_support import from_dtype
from numba import njit, typeof
from numba.core import types
from numba.tests.support import (TestCase, MemoryLeakMixin,
                                 skip_parfors_unsupported)
from numba.core.errors import TypingError
from numba.experimental import jitclass


def array_dtype(a):
    return a.dtype


def use_dtype(a, b):
    return a.view(b.dtype)


def dtype_eq_int64(a):
    return a.dtype == np.dtype('int64')


def array_itemsize(a):
    return a.itemsize


def array_nbytes(a):
    return a.nbytes


def array_shape(a, i):
    return a.shape[i]


def array_strides(a, i):
    return a.strides[i]


def array_ndim(a):
    return a.ndim


def array_size(a):
    return a.size


def array_flags_contiguous(a):
    return a.flags.contiguous

def array_flags_c_contiguous(a):
    return a.flags.c_contiguous

def array_flags_f_contiguous(a):
    return a.flags.f_contiguous


def nested_array_itemsize(a):
    return a.f.itemsize

def nested_array_nbytes(a):
    return a.f.nbytes

def nested_array_shape(a):
    return a.f.shape


def nested_array_strides(a):
    return a.f.strides


def nested_array_ndim(a):
    return a.f.ndim


def nested_array_size(a):
    return a.f.size


def size_after_slicing_usecase(buf, i):
    sliced = buf[i]
    # Make sure size attribute is not lost
    return sliced.size


def array_ctypes_data(arr):
    return arr.ctypes.data


def array_real(arr):
    return arr.real


def array_imag(arr):
    return arr.imag


class TestArrayAttr(MemoryLeakMixin, TestCase):

    def setUp(self):
        super(TestArrayAttr, self).setUp()
        self.a = np.arange(20, dtype=np.int32).reshape(4, 5)

    def check_unary(self, pyfunc, arr):
        aryty = typeof(arr)
        cfunc = self.get_cfunc(pyfunc, (aryty,))
        expected = pyfunc(arr)
        self.assertPreciseEqual(cfunc(arr), expected)
        # Retry with forced any layout
        cfunc = self.get_cfunc(pyfunc, (aryty.copy(layout='A'),))
        self.assertPreciseEqual(cfunc(arr), expected)

    def check_unary_with_arrays(self, pyfunc,):
        self.check_unary(pyfunc, self.a)
        self.check_unary(pyfunc, self.a.T)
        self.check_unary(pyfunc, self.a[::2])
        # 0-d array
        arr = np.array([42]).reshape(())
        self.check_unary(pyfunc, arr)
        # array with an empty dimension
        arr = np.zeros(0)
        self.check_unary(pyfunc, arr)

        # check with reshape
        self.check_unary(pyfunc, arr.reshape((1, 0, 2)))

    def get_cfunc(self, pyfunc, argspec):
        return njit(argspec)(pyfunc)

    def test_shape(self):
        pyfunc = array_shape
        cfunc = self.get_cfunc(pyfunc, (types.int32[:,:], types.int32))

        for i in range(self.a.ndim):
            self.assertEqual(pyfunc(self.a, i), cfunc(self.a, i))

    def test_strides(self):
        pyfunc = array_strides
        cfunc = self.get_cfunc(pyfunc, (types.int32[:,:], types.int32))

        for i in range(self.a.ndim):
            self.assertEqual(pyfunc(self.a, i), cfunc(self.a, i))

    def test_ndim(self):
        self.check_unary_with_arrays(array_ndim)

    def test_size(self):
        self.check_unary_with_arrays(array_size)

    def test_itemsize(self):
        self.check_unary_with_arrays(array_itemsize)

    def test_nbytes(self):
        self.check_unary_with_arrays(array_nbytes)

    def test_dtype(self):
        pyfunc = array_dtype
        self.check_unary(pyfunc, self.a)
        dtype = np.dtype([('x', np.int8), ('y', np.int8)])
        arr = np.zeros(4, dtype=dtype)
        self.check_unary(pyfunc, arr)

    def test_use_dtype(self):
        # Test using the dtype attribute inside the Numba function itself
        b = np.empty(1, dtype=np.int16)
        pyfunc = use_dtype
        cfunc = self.get_cfunc(pyfunc, (typeof(self.a), typeof(b)))
        expected = pyfunc(self.a, b)
        self.assertPreciseEqual(cfunc(self.a, b), expected)

    def test_dtype_equal(self):
        # Test checking if a dtype is equal to another dtype
        pyfunc = dtype_eq_int64
        self.check_unary(pyfunc, np.empty(1, dtype=np.int16))
        self.check_unary(pyfunc, np.empty(1, dtype=np.int64))

    def test_flags_contiguous(self):
        self.check_unary_with_arrays(array_flags_contiguous)

    def test_flags_c_contiguous(self):
        self.check_unary_with_arrays(array_flags_c_contiguous)

    def test_flags_f_contiguous(self):
        self.check_unary_with_arrays(array_flags_f_contiguous)


class TestNestedArrayAttr(MemoryLeakMixin, unittest.TestCase):
    def setUp(self):
        super(TestNestedArrayAttr, self).setUp()
        dtype = np.dtype([('a', np.int32), ('f', np.int32, (2, 5))])
        self.a = np.recarray(1, dtype)[0]
        self.nbrecord = from_dtype(self.a.dtype)

    def get_cfunc(self, pyfunc):
        return njit((self.nbrecord,))(pyfunc)

    def test_shape(self):
        pyfunc = nested_array_shape
        cfunc = self.get_cfunc(pyfunc)

        self.assertEqual(pyfunc(self.a), cfunc(self.a))

    def test_strides(self):
        pyfunc = nested_array_strides
        cfunc = self.get_cfunc(pyfunc)

        self.assertEqual(pyfunc(self.a), cfunc(self.a))

    def test_ndim(self):
        pyfunc = nested_array_ndim
        cfunc = self.get_cfunc(pyfunc)

        self.assertEqual(pyfunc(self.a), cfunc(self.a))

    def test_nbytes(self):
        pyfunc = nested_array_nbytes
        cfunc = self.get_cfunc(pyfunc)

        self.assertEqual(pyfunc(self.a), cfunc(self.a))

    def test_size(self):
        pyfunc = nested_array_size
        cfunc = self.get_cfunc(pyfunc)

        self.assertEqual(pyfunc(self.a), cfunc(self.a))

    def test_itemsize(self):
        pyfunc = nested_array_itemsize
        cfunc = self.get_cfunc(pyfunc)

        self.assertEqual(pyfunc(self.a), cfunc(self.a))


class TestSlicedArrayAttr(MemoryLeakMixin, unittest.TestCase):
    def test_size_after_slicing(self):
        pyfunc = size_after_slicing_usecase
        cfunc = njit(pyfunc)
        arr = np.arange(2 * 5).reshape(2, 5)
        for i in range(arr.shape[0]):
            self.assertEqual(pyfunc(arr, i), cfunc(arr, i))
        arr = np.arange(2 * 5 * 3).reshape(2, 5, 3)
        for i in range(arr.shape[0]):
            self.assertEqual(pyfunc(arr, i), cfunc(arr, i))


class TestArrayCTypes(MemoryLeakMixin, TestCase):

    _numba_parallel_test_ = False

    def test_array_ctypes_data(self):
        pyfunc = array_ctypes_data
        cfunc = njit(pyfunc)
        arr = np.arange(3)
        self.assertEqual(pyfunc(arr), cfunc(arr))

    @skip_parfors_unsupported
    def test_array_ctypes_ref_error_in_parallel(self):
        # Issue #2887
        from ctypes import CFUNCTYPE, c_void_p, c_int32, c_double, c_bool

        @CFUNCTYPE(c_bool, c_void_p, c_int32, c_void_p)
        def callback(inptr, size, outptr):
            # A ctypes callback that manipulate the incoming pointers.
            try:
                inbuf = (c_double * size).from_address(inptr)
                outbuf = (c_double * 1).from_address(outptr)
                a = np.ndarray(size, buffer=inbuf, dtype=np.float64)
                b = np.ndarray(1, buffer=outbuf, dtype=np.float64)
                b[0] = (a + a.size)[0]
                return True
            except:
                import traceback
                traceback.print_exception()
                return False


        # parallel=True is required to reproduce the error.
        @njit(parallel=True)
        def foo(size):
            arr = np.ones(size)
            out = np.empty(1)
            # Exercise array.ctypes
            inct = arr.ctypes
            outct = out.ctypes
            # The reference to `arr` is dead by now
            status = callback(inct.data, size, outct.data)
            return status, out[0]

        size = 3
        status, got = foo(size)
        self.assertTrue(status)
        self.assertPreciseEqual(got, (np.ones(size) + size)[0])


class TestRealImagAttr(MemoryLeakMixin, TestCase):
    def check_complex(self, pyfunc):
        cfunc = njit(pyfunc)
        # test 1D
        size = 10
        arr = np.arange(size) + np.arange(size) * 10j
        self.assertPreciseEqual(pyfunc(arr), cfunc(arr))
        # test 2D
        arr = arr.reshape(2, 5)
        self.assertPreciseEqual(pyfunc(arr), cfunc(arr))

    def test_complex_real(self):
        self.check_complex(array_real)

    def test_complex_imag(self):
        self.check_complex(array_imag)

    def check_number_real(self, dtype):
        pyfunc = array_real
        cfunc = njit(pyfunc)
        # test 1D
        size = 10
        arr = np.arange(size, dtype=dtype)
        self.assertPreciseEqual(pyfunc(arr), cfunc(arr))
        # test 2D
        arr = arr.reshape(2, 5)
        self.assertPreciseEqual(pyfunc(arr), cfunc(arr))
        # test identity
        self.assertEqual(arr.data, pyfunc(arr).data)
        self.assertEqual(arr.data, cfunc(arr).data)
        # test writable
        real = cfunc(arr)
        self.assertNotEqual(arr[0, 0], 5)
        real[0, 0] = 5
        self.assertEqual(arr[0, 0], 5)

    def test_number_real(self):
        """
        Testing .real of non-complex dtypes
        """
        for dtype in [np.uint8, np.int32, np.float32, np.float64]:
            self.check_number_real(dtype)

    def check_number_imag(self, dtype):
        pyfunc = array_imag
        cfunc = njit(pyfunc)
        # test 1D
        size = 10
        arr = np.arange(size, dtype=dtype)
        self.assertPreciseEqual(pyfunc(arr), cfunc(arr))
        # test 2D
        arr = arr.reshape(2, 5)
        self.assertPreciseEqual(pyfunc(arr), cfunc(arr))
        # test are zeros
        self.assertEqual(cfunc(arr).tolist(), np.zeros_like(arr).tolist())
        # test readonly
        imag = cfunc(arr)
        with self.assertRaises(ValueError) as raises:
            imag[0] = 1
        self.assertEqual('assignment destination is read-only',
                         str(raises.exception))

    def test_number_imag(self):
        """
        Testing .imag of non-complex dtypes
        """
        for dtype in [np.uint8, np.int32, np.float32, np.float64]:
            self.check_number_imag(dtype)

    def test_record_real(self):
        rectyp = np.dtype([('real', np.float32), ('imag', np.complex64)])
        arr = np.zeros(3, dtype=rectyp)
        arr['real'] = np.random.random(arr.size)
        arr['imag'] = np.random.random(arr.size) * 1.3j

        # check numpy behavior
        # .real is identity
        self.assertIs(array_real(arr), arr)
        # .imag is zero_like
        self.assertEqual(array_imag(arr).tolist(), np.zeros_like(arr).tolist())

        # check numba behavior
        # it's most likely a user error, anyway
        jit_array_real = njit(array_real)
        jit_array_imag = njit(array_imag)

        with self.assertRaises(TypingError) as raises:
            jit_array_real(arr)
        self.assertIn("cannot access .real of array of Record",
                      str(raises.exception))

        with self.assertRaises(TypingError) as raises:
            jit_array_imag(arr)
        self.assertIn("cannot access .imag of array of Record",
                      str(raises.exception))

class TestJitclassFlagsSegfault(MemoryLeakMixin, TestCase):
    """Regression test for: https://github.com/numba/numba/issues/4775 """

    def test(self):

        @jitclass(dict())
        class B(object):

            def __init__(self):
                pass

            def foo(self, X):
                X.flags

        Z = B()
        Z.foo(np.ones(4))

if __name__ == '__main__':
    unittest.main()
