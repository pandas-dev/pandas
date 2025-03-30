from itertools import product, cycle
import gc
import sys
import unittest
import warnings

import numpy as np

from numba import jit, njit, typeof
from numba.core import types
from numba.core.errors import TypingError, NumbaValueError
from numba.np.numpy_support import as_dtype, numpy_version
from numba.tests.support import (TestCase, MemoryLeakMixin,
                                 needs_blas, skip_if_numpy_2,
                                 expected_failure_np2)

TIMEDELTA_M = 'timedelta64[M]'
TIMEDELTA_Y = 'timedelta64[Y]'


def np_around_array(arr, decimals, out):
    np.around(arr, decimals, out)

def np_around_binary(val, decimals):
    return np.around(val, decimals)

def np_around_unary(val):
    return np.around(val)

def np_round_array(arr, decimals, out):
    np.round(arr, decimals, out)

def np_round__array(arr, decimals, out):
    np.round_(arr, decimals, out)

def np_round_binary(val, decimals):
    return np.round(val, decimals)

def np_round_unary(val):
    return np.round(val)

def _fixed_np_round(arr, decimals=0, out=None):
    """
    A slightly bugfixed version of np.round().
    """
    if out is not None and arr.dtype.kind == 'c':
        # workaround for https://github.com/numpy/numpy/issues/5779
        _fixed_np_round(arr.real, decimals, out.real)
        _fixed_np_round(arr.imag, decimals, out.imag)
        return out
    else:
        res = np.round(arr, decimals, out)
        if out is None:
            # workaround for https://github.com/numpy/numpy/issues/5780
            def fixup_signed_zero(arg, res):
                if res == 0.0 and arg < 0:
                    return -np.abs(res)
                else:
                    return res
            if isinstance(arr, (complex, np.complexfloating)):
                res = complex(fixup_signed_zero(arr.real, res.real),
                              fixup_signed_zero(arr.imag, res.imag))
            else:
                res = fixup_signed_zero(arr, res)
        return res


def array_T(arr):
    return arr.T

def array_transpose(arr):
    return arr.transpose()

def array_copy(arr):
    return arr.copy()

def np_copy(arr):
    return np.copy(arr)

def np_asfortranarray(arr):
    return np.asfortranarray(arr)

def np_ascontiguousarray(arr):
    return np.ascontiguousarray(arr)

def array_view(arr, newtype):
    return arr.view(newtype)

def array_take(arr, indices):
    return arr.take(indices)

def array_take_kws(arr, indices, axis):
    return arr.take(indices, axis=axis)

def np_arange_1(arg0):
    return np.arange(arg0)

def np_arange_2(arg0, arg1):
    return np.arange(arg0, arg1)

def np_arange_3(arg0, arg1, arg2):
    return np.arange(arg0, arg1, arg2)

def np_arange_4(arg0, arg1, arg2, arg3):
    return np.arange(arg0, arg1, arg2, arg3)

def np_arange_1_stop(arg0, stop):
    return np.arange(arg0, stop=stop)

def np_arange_1_step(arg0, step):
    return np.arange(arg0, step=step)

def np_arange_1_dtype(arg0, dtype):
    return np.arange(arg0, dtype=dtype)

def np_arange_2_step(arg0, arg1, step):
    return np.arange(arg0, arg1, step=step)

def np_arange_2_dtype(arg0, arg1, dtype):
    return np.arange(arg0, arg1, dtype=dtype)

def np_arange_start_stop(start, stop):
    return np.arange(start, stop=stop)

def np_arange_start_stop_step(start, stop, step):
    return np.arange(start, stop=stop, step=step)

def np_arange_start_stop_step_dtype(start, stop, step, dtype):
    return np.arange(start, stop=stop, step=step, dtype=dtype)

def array_fill(arr, val):
    return arr.fill(val)

# XXX Can't pass a dtype as a Dispatcher argument for now
def make_array_view(newtype):
    def array_view(arr):
        return arr.view(newtype)
    return array_view

def array_sliced_view(arr, ):
    return arr[0:4].view(np.float32)[0]

def make_array_astype(newtype):
    def array_astype(arr):
        return arr.astype(newtype)
    return array_astype


def np_frombuffer(b):
    """
    np.frombuffer() on a Python-allocated buffer.
    """
    return np.frombuffer(b)

def np_frombuffer_dtype(b):
    return np.frombuffer(b, dtype=np.complex64)

def np_frombuffer_dtype_str(b):
    return np.frombuffer(b, dtype='complex64')

def np_frombuffer_allocated(shape):
    """
    np.frombuffer() on a Numba-allocated buffer.
    """
    arr = np.ones(shape, dtype=np.int32)
    return np.frombuffer(arr)

def np_frombuffer_allocated_dtype(shape):
    arr = np.ones(shape, dtype=np.int32)
    return np.frombuffer(arr, dtype=np.complex64)

def identity_usecase(a, b):
    return (a is b), (a is not b)

def array_nonzero(a):
    return a.nonzero()

def np_nonzero(a):
    return np.nonzero(a)

def np_where_1(c):
    return np.where(c)

def np_where_3(c, x, y):
    return np.where(c, x, y)

def array_item(a):
    return a.item()

def array_itemset(a, v):
    a.itemset(v)

def array_sum(a, *args):
    return a.sum(*args)

def array_sum_axis_kws(a, axis):
    return a.sum(axis=axis)

def array_sum_dtype_kws(a, dtype):
    return a.sum(dtype=dtype)

def array_sum_axis_dtype_kws(a, dtype, axis):
    return a.sum(axis=axis, dtype=dtype)

def array_sum_axis_dtype_pos(a, a1, a2):
    return a.sum(a1, a2)

def array_sum_const_multi(arr, axis):
    # use np.sum with different constant args multiple times to check
    # for internal compile cache to see if constant-specialization is
    # applied properly.
    a = np.sum(arr, axis=4)
    b = np.sum(arr, 3)
    # the last invocation uses runtime-variable
    c = np.sum(arr, axis)
    # as method
    d = arr.sum(axis=5)
    # negative const axis
    e = np.sum(arr, axis=-1)
    return a, b, c, d, e

def array_sum_const_axis_neg_one(a, axis):
    # use .sum with -1 axis, this is for use with 1D arrays where the above
    # "const_multi" variant would raise errors
    return a.sum(axis=-1)

def array_cumsum(a, *args):
    return a.cumsum(*args)

def array_cumsum_kws(a, axis):
    return a.cumsum(axis=axis)

def array_real(a):
    return np.real(a)

def array_imag(a):
    return np.imag(a)

def np_clip_no_out(a, a_min, a_max):
    return np.clip(a, a_min, a_max)

def np_clip(a, a_min, a_max, out=None):
    return np.clip(a, a_min, a_max, out)

def np_clip_kwargs(a, a_min, a_max, out=None):
    return np.clip(a, a_min, a_max, out=out)

def array_clip(a, a_min=None, a_max=None, out=None):
    return a.clip(a_min, a_max, out)

def array_clip_kwargs(a, a_min=None, a_max=None, out=None):
    return a.clip(a_min, a_max, out=out)

def array_clip_no_out(a, a_min, a_max):
    return a.clip(a_min, a_max)

def array_conj(a):
    return a.conj()

def array_conjugate(a):
    return a.conjugate()

def np_unique(a):
    return np.unique(a)


def array_dot(a, b):
    return a.dot(b)

def array_dot_chain(a, b):
    return a.dot(b).dot(b)

def array_ctor(n, dtype):
    return np.ones(n, dtype=dtype)

class TestArrayMethods(MemoryLeakMixin, TestCase):
    """
    Test various array methods and array-related functions.
    """

    def setUp(self):
        super(TestArrayMethods, self).setUp()

    def check_round_scalar(self, unary_pyfunc, binary_pyfunc):
        base_values = [-3.0, -2.5, -2.25, -1.5, 1.5, 2.25, 2.5, 2.75]
        complex_values = [x * (1 - 1j) for x in base_values]
        int_values = [int(x) for x in base_values]
        argtypes = (types.float64, types.float32, types.int32,
                    types.complex64, types.complex128)
        argvalues = [base_values, base_values, int_values,
                     complex_values, complex_values]

        pyfunc = binary_pyfunc
        for ty, values in zip(argtypes, argvalues):
            cfunc = njit((ty, types.int32))(pyfunc)
            for decimals in (1, 0, -1):
                for v in values:
                    if decimals > 0:
                        v *= 10
                    expected = _fixed_np_round(v, decimals)
                    got = cfunc(v, decimals)
                    self.assertPreciseEqual(got, expected)

        pyfunc = unary_pyfunc
        for ty, values in zip(argtypes, argvalues):
            cfunc = njit((ty,))(pyfunc)
            for v in values:
                expected = _fixed_np_round(v)
                got = cfunc(v)
                self.assertPreciseEqual(got, expected)

    def test_round_scalar(self):
        self.check_round_scalar(np_round_unary, np_round_binary)

    def test_around_scalar(self):
        self.check_round_scalar(np_around_unary, np_around_binary)

    def check_round_array(self, pyfunc):
        def check_round(cfunc, values, inty, outty, decimals):
            # Create input and output arrays of the right type
            arr = values.astype(as_dtype(inty))
            out = np.zeros_like(arr).astype(as_dtype(outty))
            pyout = out.copy()
            _fixed_np_round(arr, decimals, pyout)
            self.memory_leak_setup()
            cfunc(arr, decimals, out)
            self.memory_leak_teardown()
            np.testing.assert_allclose(out, pyout)
            # Output shape mismatch
            with self.assertRaises(ValueError) as raises:
                cfunc(arr, decimals, out[1:])
            self.assertEqual(str(raises.exception),
                             "invalid output shape")

        def check_types(argtypes, outtypes, values):
            for inty, outty in product(argtypes, outtypes):
                argtys = (types.Array(inty, 1, 'A'), types.int32,
                          types.Array(outty, 1, 'A'))
                cfunc = njit(argtys)(pyfunc)
                check_round(cfunc, values, inty, outty, 0)
                check_round(cfunc, values, inty, outty, 1)
                if not isinstance(outty, types.Integer):
                    check_round(cfunc, values * 10, inty, outty, -1)
                else:
                    # Avoid Numpy bug when output is an int:
                    # https://github.com/numpy/numpy/issues/5777
                    pass

        values = np.array([-3.0, -2.5, -2.25, -1.5, 1.5, 2.25, 2.5, 2.75])

        argtypes = (types.float64, types.float32)
        check_types(argtypes, argtypes, values)

        argtypes = (types.complex64, types.complex128)
        check_types(argtypes, argtypes, values * (1 - 1j))

        # Exceptions leak references
        self.disable_leak_check()

    def test_round_array(self):
        self.check_round_array(np_round_array)

    def test_around_array(self):
        self.check_round_array(np_around_array)

    @skip_if_numpy_2
    def test_round__array(self):
        self.check_round_array(np_round__array)

    def test_around_bad_array(self):
        for pyfunc in (np_round_unary, np_around_unary):
            cfunc = jit(nopython=True)(pyfunc)
            msg = '.*The argument "a" must be array-like.*'
            with self.assertRaisesRegex(TypingError, msg):
                cfunc(None)

    def test_around_bad_out(self):
        funcs = [np_round_array, np_around_array]
        if numpy_version < (2, 0):
            funcs.append(np_round__array)
        for py_func in funcs:
            cfunc = jit(nopython=True)(py_func)
            msg = '.*The argument "out" must be an array if it is provided.*'
            with self.assertRaisesRegex(TypingError, msg):
                cfunc(5, 0, out=6)

    def test_array_view(self):

        def run(arr, dtype):
            pyfunc = make_array_view(dtype)
            return njit(pyfunc)(arr)

        def check(arr, dtype):
            expected = arr.view(dtype)
            self.memory_leak_setup()
            got = run(arr, dtype)
            self.assertPreciseEqual(got, expected)
            del got
            self.memory_leak_teardown()

        def check_err(arr, dtype):
            with self.assertRaises(ValueError) as raises:
                run(arr, dtype)
            self.assertEqual(str(raises.exception),
                             "new type not compatible with array")

        def check_err_noncontig_last_axis(arr, dtype):
            # check NumPy interpreted version raises
            msg = ("To change to a dtype of a different size, the last axis "
                   "must be contiguous")
            with self.assertRaises(ValueError) as raises:
                make_array_view(dtype)(arr)
            self.assertEqual(str(raises.exception), msg)
            # check Numba version raises
            with self.assertRaises(ValueError) as raises:
                run(arr, dtype)
            self.assertEqual(str(raises.exception), msg)

        def check_err_0d(arr, dtype):
            # check NumPy interpreted version raises
            msg = ("Changing the dtype of a 0d array is only supported "
                   "if the itemsize is unchanged")
            with self.assertRaises(ValueError) as raises:
                make_array_view(dtype)(arr)
            self.assertEqual(str(raises.exception), msg)
            # check Numba version raises
            with self.assertRaises(ValueError) as raises:
                run(arr, dtype)
            self.assertEqual(str(raises.exception), msg)

        def check_err_smaller_dtype(arr, dtype):
            # check NumPy interpreted version raises
            msg = ("When changing to a smaller dtype, its size must be a "
                   "divisor of the size of original dtype")
            with self.assertRaises(ValueError) as raises:
                make_array_view(dtype)(arr)
            self.assertEqual(str(raises.exception), msg)
            # check Numba version raises
            with self.assertRaises(ValueError) as raises:
                run(arr, dtype)
            self.assertEqual(str(raises.exception), msg)

        def check_err_larger_dtype(arr, dtype):
            # check NumPy interpreted version raises
            msg = ("When changing to a larger dtype, its size must be a "
                   "divisor of the total size in bytes of the last axis "
                   "of the array.")
            with self.assertRaises(ValueError) as raises:
                make_array_view(dtype)(arr)
            self.assertEqual(str(raises.exception), msg)
            # check Numba version raises
            with self.assertRaises(ValueError) as raises:
                run(arr, dtype)
            self.assertEqual(str(raises.exception), msg)

        dt1 = np.dtype([('a', np.int8), ('b', np.int8)])
        dt2 = np.dtype([('u', np.int16), ('v', np.int8)])
        dt3 = np.dtype([('x', np.int16), ('y', np.int16)])

        check_error_larger_dt = check_err_larger_dtype
        check_error_smaller_dt = check_err_smaller_dtype
        check_error_noncontig = check_err_noncontig_last_axis
        check_error_0d = check_err_0d

        # C-contiguous
        arr = np.arange(24, dtype=np.int8)
        check(arr, np.dtype('int16'))
        check(arr, np.int16)
        check(arr, np.int8)
        check(arr, np.float32)
        check(arr, np.complex64)
        check(arr, dt1)
        check(arr, dt2)
        check_error_larger_dt(arr, np.complex128)

        # Last dimension must have a compatible size
        arr = arr.reshape((3, 8))
        check(arr, np.int8)
        check(arr, np.float32)
        check(arr, np.complex64)
        check(arr, dt1)
        check_error_larger_dt(arr, dt2)
        check_error_larger_dt(arr, np.complex128)

        # F-contiguous
        f_arr = np.arange(24, dtype=np.int8).reshape((3, 8)).T
        # neither F or C contiguous
        not_f_or_c_arr = np.zeros((4, 4)).T[::2, ::2]

        check_maybe_error = check_err_noncontig_last_axis

        check(f_arr, np.int8)
        check(not_f_or_c_arr, np.uint64)
        check_maybe_error(f_arr, np.float32)
        check_maybe_error(f_arr, np.complex64)
        check_maybe_error(f_arr, dt1)

        check_error_noncontig(f_arr, dt2)
        check_error_noncontig(f_arr, np.complex128)
        check_error_noncontig(not_f_or_c_arr, np.int8)

        # Non-contiguous: only a type with the same itemsize can be used
        arr = np.arange(16, dtype=np.int32)[::2]
        check(arr, np.uint32)
        check(arr, np.float32)
        check(arr, dt3)
        check_error_noncontig(arr, np.int8)
        check_error_noncontig(arr, np.int16)
        check_error_noncontig(arr, np.int64)
        check_error_noncontig(arr, dt1)
        check_error_noncontig(arr, dt2)

        ## Zero-dim array: only a type with the same itemsize can be used
        arr = np.array([42], dtype=np.int32).reshape(())
        check(arr, np.uint32)
        check(arr, np.float32)
        check(arr, dt3)
        check_error_0d(arr, np.int8)
        check_error_0d(arr, np.int16)
        check_error_0d(arr, np.int64)
        check_error_0d(arr, dt1)
        check_error_0d(arr, dt2)

        # Changing to smaller dtype
        arr = np.array(['abcdef'])
        check_error_smaller_dt(arr, np.complex128)

        # Exceptions leak references
        self.disable_leak_check()

    def test_array_sliced_view(self):
        """
        Test .view() on A layout array but has contiguous innermost dimension.
        """
        pyfunc = array_sliced_view
        cfunc = njit((types.uint8[:],))(pyfunc)

        orig = np.array([1.5, 2], dtype=np.float32)
        byteary = orig.view(np.uint8)

        expect = pyfunc(byteary)
        got = cfunc(byteary)

        self.assertEqual(expect, got)

    def test_array_astype(self):

        def run(arr, dtype):
            pyfunc = make_array_astype(dtype)
            return njit(pyfunc)(arr)

        def check(arr, dtype):
            expected = arr.astype(dtype).copy(order='A')
            got = run(arr, dtype)
            self.assertPreciseEqual(got, expected)

        # C-contiguous
        arr = np.arange(24, dtype=np.int8)
        check(arr, np.dtype('int16'))
        check(arr, np.int32)
        check(arr, np.float32)
        check(arr, np.complex128)
        check(arr, "float32")

        # F-contiguous
        arr = np.arange(24, dtype=np.int8).reshape((3, 8)).T
        check(arr, np.float32)

        # Non-contiguous
        arr = np.arange(16, dtype=np.int32)[::2]
        check(arr, np.uint64)

        # check read only attr does not get copied
        arr = np.arange(16, dtype=np.int32)
        arr.flags.writeable = False
        check(arr, np.int32)

        # Invalid conversion
        dt = np.dtype([('x', np.int8)])
        with self.assertTypingError() as raises:
            check(arr, dt)
        self.assertIn('cannot convert from int32 to Record',
                      str(raises.exception))
        # Check non-Literal string raises
        unicode_val = "float32"
        with self.assertTypingError() as raises:
            @jit(nopython=True)
            def foo(dtype):
                np.array([1]).astype(dtype)
            foo(unicode_val)
        self.assertIn('array.astype if dtype is a string it must be constant',
                      str(raises.exception))

    def check_np_frombuffer(self, pyfunc):

        cfunc = njit(pyfunc)

        def check(buf):
            old_refcnt = sys.getrefcount(buf)
            expected = pyfunc(buf)
            self.memory_leak_setup()
            got = cfunc(buf)
            self.assertPreciseEqual(got, expected)
            del expected
            # Note gc.collect is due to references in `except ... as e` that
            # aren't immediately cleared
            gc.collect()
            self.assertEqual(sys.getrefcount(buf), old_refcnt + 1)
            del got
            gc.collect()
            self.assertEqual(sys.getrefcount(buf), old_refcnt)
            self.memory_leak_teardown()

        b = bytearray(range(16))
        check(b)
        check(bytes(b))
        check(memoryview(b))
        check(np.arange(12))
        b = np.arange(12).reshape((3, 4))
        check(b)

        # Exceptions leak references
        self.disable_leak_check()

        with self.assertRaises(ValueError) as raises:
            cfunc(bytearray(b"xxx"))
        self.assertEqual("buffer size must be a multiple of element size",
                         str(raises.exception))

    def test_np_frombuffer(self):
        self.check_np_frombuffer(np_frombuffer)

    def test_np_frombuffer_dtype(self):
        self.check_np_frombuffer(np_frombuffer_dtype)

    def test_np_frombuffer_dtype_str(self):
        self.check_np_frombuffer(np_frombuffer_dtype_str)

    def test_np_frombuffer_dtype_non_const_str(self):
        @jit(nopython=True)
        def func(buf, dt):
            np.frombuffer(buf, dtype=dt)

        with self.assertRaises(TypingError) as raises:
            func(bytearray(range(16)), 'int32')

        excstr = str(raises.exception)
        msg = ("If np.frombuffer dtype is a string it must be a "
               "string constant.")
        self.assertIn(msg, excstr)

    def test_np_frombuffer_bad_buffer(self):
        @jit(nopython=True)
        def func(buf):
            return np.frombuffer(buf)

        msg = '.*Argument "buffer" must be buffer-like.*'
        with self.assertRaisesRegex(TypingError, msg) as raises:
            func(None)

    def check_layout_dependent_func(self, pyfunc, fac=np.arange):
        def is_same(a, b):
            return a.ctypes.data == b.ctypes.data
        def check_arr(arr):
            cfunc = njit((typeof(arr),))(pyfunc)
            expected = pyfunc(arr)
            got = cfunc(arr)
            self.assertPreciseEqual(expected, got)
            self.assertEqual(is_same(expected, arr), is_same(got, arr))
        arr = fac(24)
        check_arr(arr)
        check_arr(arr.reshape((3, 8)))
        check_arr(arr.reshape((3, 8)).T)
        check_arr(arr.reshape((3, 8))[::2])
        check_arr(arr.reshape((2, 3, 4)))
        check_arr(arr.reshape((2, 3, 4)).T)
        check_arr(arr.reshape((2, 3, 4))[::2])
        arr = np.array([0]).reshape(())
        check_arr(arr)

    def test_array_transpose(self):
        self.check_layout_dependent_func(array_transpose)

    def test_array_T(self):
        self.check_layout_dependent_func(array_T)

    def test_array_copy(self):
        self.check_layout_dependent_func(array_copy)

    def test_np_copy(self):
        self.check_layout_dependent_func(np_copy)

    def check_ascontiguousarray_scalar(self, pyfunc):
        def check_scalar(x):
            cfunc = njit((typeof(x),))(pyfunc)
            expected = pyfunc(x)
            got = cfunc(x)
            self.assertPreciseEqual(expected, got)
        for x in [42, 42.0, 42j, np.float32(42), np.float64(42), True]:
            check_scalar(x)

    def check_bad_array(self, pyfunc):
        msg = '.*The argument "a" must be array-like.*'
        with self.assertRaisesRegex(TypingError, msg) as raises:
            njit((typeof('hello'), ))(pyfunc)

    def test_np_asfortranarray(self):
        self.check_layout_dependent_func(np_asfortranarray)
        self.check_bad_array(np_asfortranarray)
        self.check_ascontiguousarray_scalar(np_asfortranarray)

    def test_np_ascontiguousarray(self):
        self.check_layout_dependent_func(np_ascontiguousarray)
        self.check_bad_array(np_asfortranarray)
        self.check_ascontiguousarray_scalar(np_ascontiguousarray)

    def check_np_frombuffer_allocated(self, pyfunc):

        cfunc = njit(pyfunc)

        def check(shape):
            expected = pyfunc(shape)
            got = cfunc(shape)
            self.assertPreciseEqual(got, expected)

        check((16,))
        check((4, 4))
        check((1, 0, 1))

    def test_np_frombuffer_allocated(self):
        self.check_np_frombuffer_allocated(np_frombuffer_allocated)

    def test_np_frombuffer_allocated2(self):
        self.check_np_frombuffer_allocated(np_frombuffer_allocated_dtype)

    def check_nonzero(self, pyfunc):
        def fac(N):
            np.random.seed(42)
            arr = np.random.random(N)
            arr[arr < 0.3] = 0.0
            arr[arr > 0.7] = float('nan')
            return arr

        def check_arr(arr):
            cfunc = njit((typeof(arr),))(pyfunc)
            expected = pyfunc(arr)
            expected = [a.copy() for a in expected]
            self.assertPreciseEqual(cfunc(arr), expected)

        arr = np.int16([1, 0, -1, 0])
        check_arr(arr)
        arr = np.bool_([1, 0, 1])
        check_arr(arr)

        arr = fac(24)
        check_arr(arr)
        check_arr(arr.reshape((3, 8)))
        check_arr(arr.reshape((3, 8)).T)
        check_arr(arr.reshape((3, 8))[::2])
        check_arr(arr.reshape((2, 3, 4)))
        check_arr(arr.reshape((2, 3, 4)).T)
        check_arr(arr.reshape((2, 3, 4))[::2])

        arr = np.array(["Hello", "", "world"])
        check_arr(arr)

        for v in (0.0, 1.5, float('nan')):
            arr = np.array([v]).reshape(())
            if numpy_version < (2, 1):
                check_arr(arr)
            else:
                with self.assertRaises((ValueError, TypingError)) as raises:
                    njit((typeof(arr),))(pyfunc)
                msg = "Calling nonzero on 0d arrays is not allowed. Use " \
                      "np.atleast_1d(scalar).nonzero() instead."
                self.assertIn(msg, str(raises.exception))

    def test_array_nonzero(self):
        self.check_nonzero(array_nonzero)

    def test_np_nonzero(self):
        self.check_nonzero(np_nonzero)

    def test_np_where_1(self):
        self.check_nonzero(np_where_1)

    def test_np_where_3(self):
        pyfunc = np_where_3
        def fac(N):
            np.random.seed(42)
            arr = np.random.random(N)
            arr[arr < 0.3] = 0.0
            arr[arr > 0.7] = float('nan')
            return arr

        layouts = cycle(['C', 'F', 'A'])
        _types = [np.int32, np.int64, np.float32, np.float64, np.complex64,
                  np.complex128]

        np.random.seed(42)

        def check_arr(arr, layout=False):
            np.random.shuffle(_types)
            if layout != False:
                x = np.zeros_like(arr, dtype=_types[0], order=layout)
                y = np.zeros_like(arr, dtype=_types[1], order=layout)
                arr = arr.copy(order=layout)
            else:
                x = np.zeros_like(arr, dtype=_types[0], order=next(layouts))
                y = np.zeros_like(arr, dtype=_types[1], order=next(layouts))
            x.fill(4)
            y.fill(9)
            cfunc = njit((typeof(arr), typeof(x), typeof(y)))(pyfunc)
            expected = pyfunc(arr, x, y)
            got = cfunc(arr, x, y)
            self.assertPreciseEqual(got, expected)

        def check_scal(scal):
            x = 4
            y = 5
            np.random.shuffle(_types)
            x = _types[0](4)
            y = _types[1](5)
            cfunc = njit((typeof(scal), typeof(x), typeof(y)))(pyfunc)
            expected = pyfunc(scal, x, y)
            got = cfunc(scal, x, y)
            self.assertPreciseEqual(got, expected)

        arr = np.int16([1, 0, -1, 0])
        check_arr(arr)
        arr = np.bool_([1, 0, 1])
        check_arr(arr)

        arr = fac(24)
        check_arr(arr)
        check_arr(arr.reshape((3, 8)))
        check_arr(arr.reshape((3, 8)).T)
        check_arr(arr.reshape((3, 8))[::2])
        check_arr(arr.reshape((2, 3, 4)))
        check_arr(arr.reshape((2, 3, 4)).T)
        check_arr(arr.reshape((2, 3, 4))[::2])

        check_arr(arr.reshape((2, 3, 4)), layout='F')
        check_arr(arr.reshape((2, 3, 4)).T, layout='F')
        check_arr(arr.reshape((2, 3, 4))[::2], layout='F')

        for v in (0.0, 1.5, float('nan')):
            arr = np.array([v]).reshape(())
            check_arr(arr)

        for x in (0, 1, True, False, 2.5, 0j):
            check_scal(x)

    def test_np_where_3_broadcast_x_y_scalar(self):
        pyfunc = np_where_3
        cfunc = jit(nopython=True)(pyfunc)

        def check_ok(args):
            expected = pyfunc(*args)
            got = cfunc(*args)
            self.assertPreciseEqual(got, expected)

        def a_variations():
            a = np.linspace(-2, 4, 20)
            self.random.shuffle(a)
            yield a
            yield a.reshape(2, 5, 2)
            yield a.reshape(4, 5, order='F')
            yield a.reshape(2, 5, 2)[::-1]

        for a in a_variations():
            params = (a > 0, 0, 1)
            check_ok(params)

            params = (a < 0, np.nan, 1 + 4j)
            check_ok(params)

            params = (a > 1, True, False)
            check_ok(params)

    def test_np_where_3_broadcast_x_or_y_scalar(self):
        pyfunc = np_where_3
        cfunc = jit(nopython=True)(pyfunc)

        def check_ok(args):
            condition, x, y = args

            expected = pyfunc(condition, x, y)
            got = cfunc(condition, x, y)
            self.assertPreciseEqual(got, expected)

            # swap x and y
            expected = pyfunc(condition, y, x)
            got = cfunc(condition, y, x)
            self.assertPreciseEqual(got, expected)

        def array_permutations():
            x = np.arange(9).reshape(3, 3)
            yield x
            yield x * 1.1
            yield np.asfortranarray(x)
            yield x[::-1]
            yield np.linspace(-10, 10, 60).reshape(3, 4, 5) * 1j

        def scalar_permutations():
            yield 0
            yield 4.3
            yield np.nan
            yield True
            yield 8 + 4j

        for x in array_permutations():
            for y in scalar_permutations():
                x_mean = np.mean(x)
                condition = x > x_mean
                params = (condition, x, y)
                check_ok(params)

    def test_np_where_numpy_basic(self):
        # https://github.com/numpy/numpy/blob/fe2bb380fd9a084b622ff3f00cb6f245e8c1a10e/numpy/core/tests/test_multiarray.py#L8670-L8694
        pyfunc = np_where_3
        cfunc = jit(nopython=True)(pyfunc)

        # skipping unsupported dtypes:
        # np.longdouble, np.clongdouble
        dts = [bool, np.int16, np.int32, np.int64, np.double, np.complex128]
        for dt in dts:
            c = np.ones(53, dtype=bool)
            np.testing.assert_equal(cfunc( c, dt(0), dt(1)), dt(0))
            np.testing.assert_equal(cfunc(~c, dt(0), dt(1)), dt(1))
            np.testing.assert_equal(cfunc(True, dt(0), dt(1)), dt(0))
            np.testing.assert_equal(cfunc(False, dt(0), dt(1)), dt(1))
            d = np.ones_like(c).astype(dt)
            e = np.zeros_like(d)
            r = d.astype(dt)
            c[7] = False
            r[7] = e[7]
            np.testing.assert_equal(cfunc(c, e, e), e)
            np.testing.assert_equal(cfunc(c, d, e), r)
            np.testing.assert_equal(cfunc(c, d, e[0]), r)
            np.testing.assert_equal(cfunc(c, d[0], e), r)
            np.testing.assert_equal(cfunc(c[::2], d[::2], e[::2]), r[::2])
            np.testing.assert_equal(cfunc(c[1::2], d[1::2], e[1::2]), r[1::2])
            np.testing.assert_equal(cfunc(c[::3], d[::3], e[::3]), r[::3])
            np.testing.assert_equal(cfunc(c[1::3], d[1::3], e[1::3]), r[1::3])
            np.testing.assert_equal(cfunc(c[::-2], d[::-2], e[::-2]), r[::-2])
            np.testing.assert_equal(cfunc(c[::-3], d[::-3], e[::-3]), r[::-3])
            np.testing.assert_equal(cfunc(c[1::-3], d[1::-3], e[1::-3]), r[1::-3])

    def test_np_where_numpy_ndim(self):
        # https://github.com/numpy/numpy/blob/fe2bb380fd9a084b622ff3f00cb6f245e8c1a10e/numpy/core/tests/test_multiarray.py#L8737-L8749
        pyfunc = np_where_3
        cfunc = jit(nopython=True)(pyfunc)

        c = [True, False]
        a = np.zeros((2, 25))
        b = np.ones((2, 25))
        r = cfunc(np.array(c)[:,np.newaxis], a, b)
        np.testing.assert_array_equal(r[0], a[0])
        np.testing.assert_array_equal(r[1], b[0])

        a = a.T
        b = b.T
        r = cfunc(c, a, b)
        np.testing.assert_array_equal(r[:,0], a[:,0])
        np.testing.assert_array_equal(r[:,1], b[:,0])

    def test_np_where_numpy_dtype_mix(self):
        # https://github.com/numpy/numpy/blob/fe2bb380fd9a084b622ff3f00cb6f245e8c1a10e/numpy/core/tests/test_multiarray.py#L8751-L8773
        pyfunc = np_where_3
        cfunc = jit(nopython=True)(pyfunc)

        c = np.array([False, True, False, False, False, False, True, False,
                     False, False, True, False])
        a = np.uint32(1)
        b = np.array([5., 0., 3., 2., -1., -4., 0., -10., 10., 1., 0., 3.],
                      dtype=np.float64)
        r = np.array([5., 1., 3., 2., -1., -4., 1., -10., 10., 1., 1., 3.],
                     dtype=np.float64)
        np.testing.assert_equal(cfunc(c, a, b), r)

        a = a.astype(np.float32)
        b = b.astype(np.int64)
        np.testing.assert_equal(cfunc(c, a, b), r)

        # non bool mask
        c = c.astype(int)
        c[c != 0] = 34242324
        np.testing.assert_equal(cfunc(c, a, b), r)
        # invert
        tmpmask = c != 0
        c[c == 0] = 41247212
        c[tmpmask] = 0
        np.testing.assert_equal(cfunc(c, b, a), r)

    def test_np_where_numpy_test_error(self):
        # https://github.com/numpy/numpy/blob/fe2bb380fd9a084b622ff3f00cb6f245e8c1a10e/numpy/core/tests/test_multiarray.py#L8794-L8799
        pyfunc = np_where_3
        cfunc = jit(nopython=True)(pyfunc)

        c = [True, True]
        a = np.ones((4, 5))
        b = np.ones((5, 5))

        self.disable_leak_check()
        with self.assertRaisesRegex(ValueError, "objects cannot be broadcast"):
            cfunc(c, a, b)

        with self.assertRaisesRegex(ValueError, "objects cannot be broadcast"):
            cfunc(c[0], a, b)

    def test_np_where_invalid_inputs(self):
        pyfunc = np_where_3
        cfunc = jit(nopython=True)(pyfunc)

        msg = 'The argument "condition" must be array-like'
        with self.assertRaisesRegex(TypingError, msg):
            cfunc(None, 2, 3)

        msg = 'The argument "x" must be array-like if provided'
        with self.assertRaisesRegex(TypingError, msg):
            cfunc(1, 'hello', 3)

        msg = 'The argument "y" must be array-like if provided'
        with self.assertRaisesRegex(TypingError, msg):
            cfunc(1, 2, 'world')

        # None values are not yet supported in np.where
        msg = 'Argument "x" or "y" cannot be None'
        with self.assertRaisesRegex(TypingError, msg):
            cfunc(1, None, None)

    def test_arange_1_arg(self):

        all_pyfuncs = (
            np_arange_1,
            lambda x: np.arange(x, 10),
            lambda x: np.arange(7, step=max(1, abs(x)))
        )

        for pyfunc in all_pyfuncs:
            cfunc = jit(nopython=True)(pyfunc)

            def check_ok(arg0):
                expected = pyfunc(arg0)
                got = cfunc(arg0)
                np.testing.assert_allclose(expected, got)

            check_ok(0)
            check_ok(1)
            check_ok(4)
            check_ok(5.5)
            check_ok(-3)
            check_ok(complex(4, 4))
            check_ok(np.int8(0))

    def test_arange_2_arg(self):
        def check_ok(arg0, arg1, pyfunc, cfunc):
            expected = pyfunc(arg0, arg1)
            got = cfunc(arg0, arg1)
            np.testing.assert_allclose(expected, got)

        all_pyfuncs = (
            np_arange_2,
            np_arange_start_stop,
            np_arange_1_stop,
            np_arange_1_step,
            lambda x, y: np.arange(x, y, 5),
            lambda x, y: np.arange(2, y, step=x),
        )

        for pyfunc in all_pyfuncs:
            cfunc = jit(nopython=True)(pyfunc)

            check_ok(-1, 5, pyfunc, cfunc)
            check_ok(-8, -1, pyfunc, cfunc)
            check_ok(4, 0.5, pyfunc, cfunc)
            check_ok(0.5, 4, pyfunc, cfunc)
            check_ok(3, None, pyfunc, cfunc)
            if numpy_version < (2, 0):
                check_ok(complex(1, 1), complex(4, 4), pyfunc, cfunc)
                check_ok(complex(4, 4), complex(1, 1), pyfunc, cfunc)

        pyfunc = np_arange_1_dtype
        cfunc = jit(nopython=True)(pyfunc)

        check_ok(5, np.float32, pyfunc, cfunc)
        check_ok(2.0, np.int32, pyfunc, cfunc)
        check_ok(7, None, pyfunc, cfunc)
        check_ok(np.int8(0), None, pyfunc, cfunc)

        if numpy_version < (2, 0):
            check_ok(10, np.complex128, pyfunc, cfunc)
            check_ok(np.complex64(10), np.complex128, pyfunc, cfunc)

    def test_arange_3_arg(self):
        windows64 = sys.platform.startswith('win32') and sys.maxsize > 2 ** 32

        def check_ok(arg0, arg1, arg2, pyfunc, cfunc, check_dtype=False):
            expected = pyfunc(arg0, arg1, arg2)
            got = cfunc(arg0, arg1, arg2)
            np.testing.assert_allclose(expected, got)
            # windows 64 cannot differentiate between a python int and a
            # np.int64 which means the result from numba is int64 more often
            # than in NumPy.
            if not windows64:
                self.assertEqual(expected.dtype, got.dtype)

        for pyfunc in (np_arange_3, np_arange_2_step, np_arange_start_stop_step):
            cfunc = jit(nopython=True)(pyfunc)

            check_ok(0, 5, 1, pyfunc, cfunc)
            check_ok(-8, -1, 3, pyfunc, cfunc)
            check_ok(0, -10, -2, pyfunc, cfunc)
            check_ok(0.5, 4, 2, pyfunc, cfunc)
            check_ok(0, 1, 0.1, pyfunc, cfunc)
            check_ok(3, 6, None, pyfunc, cfunc)
            check_ok(3, None, None, pyfunc, cfunc)
            check_ok(np.int8(0), np.int8(5), np.int8(1), pyfunc, cfunc)
            check_ok(np.int8(0), np.int16(5), np.int32(1), pyfunc, cfunc)
            # check upcasting logic, this matters most on windows
            i8 = np.int8
            check_ok(i8(0), i8(5), i8(1), pyfunc, cfunc, True) # C int
            check_ok(np.int64(0), i8(5), i8(1), pyfunc, cfunc, True) # int64
            if numpy_version < (2, 0):
                check_ok(0, complex(4, 4), complex(1, 1), pyfunc, cfunc)

        pyfunc = np_arange_2_dtype
        cfunc = jit(nopython=True)(pyfunc)

        check_ok(1, 5, np.float32, pyfunc, cfunc)
        check_ok(2.0, 8, np.int32, pyfunc, cfunc)
        check_ok(1, 7, None, pyfunc, cfunc)
        check_ok(np.int8(0), np.int32(5), None, pyfunc, cfunc, True)
        if numpy_version < (2, 0):
            check_ok(-2, 10, np.complex128, pyfunc, cfunc)
            check_ok(3, np.complex64(10), np.complex128, pyfunc, cfunc)

    def test_arange_4_arg(self):
        for pyfunc in (np_arange_4, np_arange_start_stop_step_dtype):
            cfunc = jit(nopython=True)(pyfunc)

            def check_ok(arg0, arg1, arg2, arg3):
                expected = pyfunc(arg0, arg1, arg2, arg3)
                got = cfunc(arg0, arg1, arg2, arg3)
                np.testing.assert_allclose(expected, got)

            check_ok(0, 5, 1, np.float64)
            check_ok(-8, -1, 3, np.int32)
            check_ok(0, -10, -2, np.float32)
            check_ok(0.5, 4, 2, None)
            check_ok(3, 6, None, None)
            check_ok(3, None, None, None)
            if numpy_version < (2, 0):
                check_ok(0, 1, 0.1, np.complex128)
                check_ok(0, complex(4, 4), complex(1, 1), np.complex128)

    def test_arange_throws(self):
        # Exceptions leak references
        self.disable_leak_check()

        bad_funcs_1 = [
            lambda x: np.arange(stop=x),
            lambda x: np.arange(step=x),
            lambda x: np.arange(dtype=x),
        ]
        bad_funcs_2 = [
            lambda x, y: np.arange(stop=x, step=y),
            lambda x, y: np.arange(stop=x, dtype=y),
        ]

        for pyfunc in bad_funcs_1:
            with self.assertRaises(TypingError) as raises:
                cfunc = jit(nopython=True)(pyfunc)
                cfunc(2)
        for pyfunc in bad_funcs_2:
            with self.assertRaises(TypingError) as raises:
                cfunc = jit(nopython=True)(pyfunc)
                cfunc(2, 6)

        # check step size = 0, this is nonsense
        pyfunc = np_arange_3
        cfunc = jit(nopython=True)(pyfunc)
        for f in (pyfunc, cfunc,):
            for inputs in [(1, np.int16(2), 0), (1, 2, 0)]:
                # there's a different error depending on whether any of the
                # input values are np scalars
                permitted_errors = (ZeroDivisionError, ValueError)
                with self.assertRaises(permitted_errors) as raises:
                    # this will raise RuntimeWarning's about zero division
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        f(*inputs)
                    self.assertIn("Maximum allowed size exceeded",
                                str(raises.exception))

    def test_arange_accuracy(self):
        # Checking arange reasonably replicates NumPy's algorithm
        # see https://github.com/numba/numba/issues/6768
        @jit(nopython=True)
        def foo(step):
            return np.arange(0, 1 + step, step)

        x = 0.010101010101010102
        self.assertPreciseEqual(foo(x), foo.py_func(x))

    def test_item(self):
        pyfunc = array_item
        cfunc = jit(nopython=True)(pyfunc)

        def check_ok(arg):
            expected = pyfunc(arg)
            got = cfunc(arg)
            self.assertPreciseEqual(got, expected)

        def check_err(arg):
            with self.assertRaises(ValueError) as raises:
                cfunc(arg)
            self.assertIn("item(): can only convert an array of size 1 to a Python scalar",
                          str(raises.exception))

        # Exceptions leak references
        self.disable_leak_check()

        # Test on different kinds of scalars and 1-item arrays
        check_ok(np.float32([1.5]))
        check_ok(np.complex128([[1.5j]]))
        check_ok(np.array(1.5))
        check_ok(np.bool_(True))
        check_ok(np.float32(1.5))

        check_err(np.array([1, 2]))
        check_err(np.array([]))

    @skip_if_numpy_2
    def test_itemset(self):
        pyfunc = array_itemset
        cfunc = jit(nopython=True)(pyfunc)

        def check_ok(a, v):
            expected = a.copy()
            got = a.copy()
            pyfunc(expected, v)
            cfunc(got, v)
            self.assertPreciseEqual(got, expected)

        def check_err(a):
            with self.assertRaises(ValueError) as raises:
                cfunc(a, 42)
            self.assertIn("itemset(): can only write to an array of size 1",
                          str(raises.exception))

        # Exceptions leak references
        self.disable_leak_check()

        # Test on different kinds of 1-item arrays
        check_ok(np.float32([1.5]), 42)
        check_ok(np.complex128([[1.5j]]), 42)
        check_ok(np.array(1.5), 42)

        check_err(np.array([1, 2]))
        check_err(np.array([]))

    def test_sum(self):
        """ test sum over a whole range of dtypes, no axis or dtype parameter
        """
        pyfunc = array_sum
        cfunc = jit(nopython=True)(pyfunc)
        all_dtypes = [np.float64, np.float32, np.int64, np.int32,
                      np.complex64, np.complex128, np.timedelta64]
        all_test_arrays = [
            [np.ones((7, 6, 5, 4, 3), arr_dtype),
             np.ones(1, arr_dtype),
             np.ones((7, 3), arr_dtype) * -5]
            for arr_dtype in all_dtypes]

        unsigned_dtypes = [np.uint32, np.uint64, np.bool_]
        all_test_arrays = [
            [np.ones((7, 6, 5, 4, 3), arr_dtype),
             np.ones(1, arr_dtype)]
            for arr_dtype in unsigned_dtypes]

        for arr_list in all_test_arrays:
            for arr in arr_list:
                with self.subTest("Test np.sum with {} input ".format(arr.dtype)):
                    self.assertPreciseEqual(pyfunc(arr), cfunc(arr))

    def test_sum_axis_kws1(self):
        """ test sum with axis parameter over a whole range of dtypes  """
        pyfunc = array_sum_axis_kws
        cfunc = jit(nopython=True)(pyfunc)
        all_dtypes = [np.float64, np.float32, np.int64, np.complex64,
                      np.complex128, TIMEDELTA_M]
        all_test_arrays = [
            [np.ones((7, 6, 5, 4, 3), arr_dtype),
             np.ones(1, arr_dtype),
             np.ones((7, 3), arr_dtype) * -5]
            for arr_dtype in all_dtypes]

        unsigned_dtypes = [np.uint64, np.bool_]
        all_test_arrays += [
            [np.ones((7, 6, 5, 4, 3), arr_dtype),
             np.ones(1, arr_dtype)]
            for arr_dtype in unsigned_dtypes]

        for arr_list in all_test_arrays:
            for arr in arr_list:
                for axis in (0, 1, 2):
                    if axis > len(arr.shape)-1:
                        continue
                    with self.subTest("Testing np.sum(axis) with {} "
                                      "input ".format(arr.dtype)):
                        self.assertPreciseEqual(pyfunc(arr, axis=axis),
                                                cfunc(arr, axis=axis))

    def test_sum_axis_kws2(self):
        """  testing uint32 and int32 separately

        uint32 and int32 must be tested separately because Numpy's current
        behaviour is different in 64bits Windows (accumulates as int32)
        and 64bits Linux (accumulates as int64), while Numba has decided to always
        accumulate as int64, when the OS is 64bits. No testing has been done
        for behaviours in 32 bits platforms.
        """
        pyfunc = array_sum_axis_kws
        cfunc = jit(nopython=True)(pyfunc)
        all_dtypes = [np.int32]
        # expected return dtypes in Numba
        out_dtypes = {np.dtype('int32'): np.int64, np.dtype('uint32'): np.uint64,
                      np.dtype('int64'): np.int64,
                      np.dtype(TIMEDELTA_M): np.dtype(TIMEDELTA_M)}
        all_test_arrays = [
            [np.ones((7, 6, 5, 4, 3), arr_dtype),
             np.ones(1, arr_dtype),
             np.ones((7, 3), arr_dtype) * -5]
            for arr_dtype in all_dtypes]

        unsigned_dtypes = [np.uint32]
        all_test_arrays += [
            [np.ones((7, 6, 5, 4, 3), arr_dtype),
             np.ones(1, arr_dtype)]
            for arr_dtype in unsigned_dtypes]

        for arr_list in all_test_arrays:
            for arr in arr_list:
                for axis in (0, 1, 2):
                    if axis > len(arr.shape)-1:
                        continue
                    with self.subTest("Testing np.sum(axis) with {} "
                                      "input ".format(arr.dtype)):
                        npy_res = pyfunc(arr, axis=axis)
                        numba_res = cfunc(arr, axis=axis)
                        if isinstance(numba_res, np.ndarray):
                            self.assertPreciseEqual(
                                npy_res.astype(out_dtypes[arr.dtype]),
                                numba_res.astype(out_dtypes[arr.dtype]))
                        else:
                            # the results are scalars
                            self.assertEqual(npy_res, numba_res)

    def test_sum_dtype_kws(self):
        """ test sum with dtype parameter over a whole range of dtypes """
        pyfunc = array_sum_dtype_kws
        cfunc = jit(nopython=True)(pyfunc)
        all_dtypes = [np.float64, np.float32, np.int64, np.int32,
                      np.complex64, np.complex128]
        all_test_arrays = [
            [np.ones((7, 6, 5, 4, 3), arr_dtype),
             np.ones(1, arr_dtype),
             np.ones((7, 3), arr_dtype) * -5]
            for arr_dtype in all_dtypes]

        unsigned_dtypes = [np.uint32, np.uint64, np.bool_]
        all_test_arrays = [
            [np.ones((7, 6, 5, 4, 3), arr_dtype),
             np.ones(1, arr_dtype)]
            for arr_dtype in unsigned_dtypes]

        out_dtypes = {np.dtype('float64'): [np.float64],
                      np.dtype('float32'): [np.float64, np.float32],
                      np.dtype('int64'): [np.float64, np.int64, np.float32],
                      np.dtype('int32'): [np.float64, np.int64, np.float32, np.int32],
                      np.dtype('uint32'): [np.float64, np.int64, np.float32],
                      np.dtype('uint64'): [np.float64, np.int64],
                      np.dtype('bool'): [np.float64, np.int64, np.float32, np.int32, np.bool_],
                      np.dtype('complex64'): [np.complex64, np.complex128],
                      np.dtype('complex128'): [np.complex128]}

        for arr_list in all_test_arrays:
            for arr in arr_list:
                for out_dtype in out_dtypes[arr.dtype]:
                    subtest_str = ("Testing np.sum with {} input and {} output"
                                   .format(arr.dtype, out_dtype))
                    with self.subTest(subtest_str):
                            self.assertPreciseEqual(pyfunc(arr, dtype=out_dtype),
                                                    cfunc(arr, dtype=out_dtype))

    def test_sum_axis_dtype_kws(self):
        """ test sum with axis and dtype parameters over a whole range of dtypes """
        pyfunc = array_sum_axis_dtype_kws
        cfunc = jit(nopython=True)(pyfunc)
        all_dtypes = [np.float64, np.float32, np.int64, np.int32,
                      np.complex64, np.complex128]
        all_test_arrays = [
            [np.ones((7, 6, 5, 4, 3), arr_dtype),
             np.ones(1, arr_dtype),
             np.ones((7, 3), arr_dtype) * -5]
            for arr_dtype in all_dtypes]

        unsigned_dtypes = [np.uint32, np.uint64, np.bool_]
        all_test_arrays = [
            [np.ones((7, 6, 5, 4, 3), arr_dtype),
             np.ones(1, arr_dtype)]
            for arr_dtype in unsigned_dtypes]

        out_dtypes = {np.dtype('float64'): [np.float64],
                      np.dtype('float32'): [np.float64, np.float32],
                      np.dtype('int64'): [np.float64, np.int64, np.float32],
                      np.dtype('int32'): [np.float64, np.int64, np.float32, np.int32],
                      np.dtype('uint32'): [np.float64, np.int64, np.float32],
                      np.dtype('uint64'): [np.float64, np.uint64],
                      np.dtype('bool'): [np.float64, np.int64, np.float32, np.int32, np.bool_],
                      np.dtype('complex64'): [np.complex64, np.complex128],
                      np.dtype('complex128'): [np.complex128]}

        for arr_list in all_test_arrays:
            for arr in arr_list:
                for out_dtype in out_dtypes[arr.dtype]:
                    for axis in (0, 1, 2):
                        if axis > len(arr.shape) - 1:
                            continue
                        subtest_str = ("Testing np.sum with {} input and {} output "
                                       .format(arr.dtype, out_dtype))
                        with self.subTest(subtest_str):
                            py_res = pyfunc(arr, axis=axis, dtype=out_dtype)
                            nb_res = cfunc(arr, axis=axis, dtype=out_dtype)
                            self.assertPreciseEqual(py_res, nb_res)

    def test_sum_axis_dtype_pos_arg(self):
        """ testing that axis and dtype inputs work when passed as positional """
        pyfunc = array_sum_axis_dtype_pos
        cfunc = jit(nopython=True)(pyfunc)
        dtype = np.float64
        # OK
        a = np.ones((7, 6, 5, 4, 3))
        self.assertPreciseEqual(pyfunc(a, 1, dtype),
                                cfunc(a,  1, dtype))

        self.assertPreciseEqual(pyfunc(a, 2, dtype),
                                cfunc(a, 2, dtype))

    def test_sum_1d_kws(self):
        # check 1d reduces to scalar
        pyfunc = array_sum_axis_kws
        cfunc = jit(nopython=True)(pyfunc)
        a = np.arange(10.)
        self.assertPreciseEqual(pyfunc(a, axis=0), cfunc(a, axis=0))
        pyfunc = array_sum_const_axis_neg_one
        cfunc = jit(nopython=True)(pyfunc)
        a = np.arange(10.)
        self.assertPreciseEqual(pyfunc(a, axis=-1), cfunc(a, axis=-1))

    def test_sum_const(self):
        pyfunc = array_sum_const_multi
        cfunc = jit(nopython=True)(pyfunc)

        arr = np.ones((3, 4, 5, 6, 7, 8))
        axis = 1
        self.assertPreciseEqual(pyfunc(arr, axis), cfunc(arr, axis))
        axis = 2
        self.assertPreciseEqual(pyfunc(arr, axis), cfunc(arr, axis))

    def test_sum_exceptions(self):
        # Exceptions leak references
        self.disable_leak_check()
        pyfunc = array_sum
        cfunc = jit(nopython=True)(pyfunc)

        a = np.ones((7, 6, 5, 4, 3))
        b = np.ones((4, 3))
        # BAD: axis > dimensions
        with self.assertRaises(ValueError):
            cfunc(b, 2)
        # BAD: negative axis
        with self.assertRaises(ValueError):
            cfunc(a, -1)
        # BAD: axis greater than 3
        with self.assertRaises(ValueError):
            cfunc(a, 4)

    def test_sum_const_negative(self):
        # Exceptions leak references
        self.disable_leak_check()

        @jit(nopython=True)
        def foo(arr):
            return arr.sum(axis=-3)

        # ndim == 4, axis == -3, OK
        a = np.ones((1, 2, 3, 4))
        self.assertPreciseEqual(foo(a), foo.py_func(a))
        # ndim == 3, axis == -3, OK
        a = np.ones((1, 2, 3))
        self.assertPreciseEqual(foo(a), foo.py_func(a))
        # ndim == 2, axis == -3, BAD
        a = np.ones((1, 2))
        with self.assertRaises(NumbaValueError) as raises:
            foo(a)
        errmsg = "'axis' entry (-1) is out of bounds"
        self.assertIn(errmsg, str(raises.exception))
        with self.assertRaises(ValueError) as raises:
            foo.py_func(a)
        self.assertIn("out of bounds", str(raises.exception))

    def test_cumsum(self):
        pyfunc = array_cumsum
        cfunc = jit(nopython=True)(pyfunc)
        # OK
        a = np.ones((2, 3))
        self.assertPreciseEqual(pyfunc(a), cfunc(a))
        # BAD: with axis
        with self.assertRaises(TypingError):
            cfunc(a, 1)
        # BAD: with kw axis
        pyfunc = array_cumsum_kws
        cfunc = jit(nopython=True)(pyfunc)
        with self.assertRaises(TypingError):
            cfunc(a, axis=1)

    def test_take(self):
        pyfunc = array_take
        cfunc = jit(nopython=True)(pyfunc)

        def check(arr, ind):
            expected = pyfunc(arr, ind)
            got = cfunc(arr, ind)
            self.assertPreciseEqual(expected, got)
            if hasattr(expected, 'order'):
                self.assertEqual(expected.order == got.order)

        # need to check:
        # 1. scalar index
        # 2. 1d array index
        # 3. nd array index, >2d and F order
        # 4. reflected list
        # 5. tuples

        test_indices = []
        test_indices.append(1)
        test_indices.append(5)
        test_indices.append(11)
        test_indices.append(-2)
        test_indices.append(np.array([1, 5, 1, 11, 3]))
        test_indices.append(np.array([[1, 5, 1], [11, 3, 0]], order='F'))
        test_indices.append(np.array([[[1, 5, 1], [11, 3, 0]]]))
        test_indices.append(np.array([[[[1, 5]], [[11, 0]], [[1, 2]]]]))
        test_indices.append([1, 5, 1, 11, 3])
        test_indices.append((1, 5, 1))
        test_indices.append(((1, 5, 1), (11, 3, 2)))
        test_indices.append((((1,), (5,), (1,)), ((11,), (3,), (2,))))

        layouts = cycle(['C', 'F', 'A'])

        for dt in [np.float64, np.int64, np.complex128]:
            A = np.arange(12, dtype=dt).reshape((4, 3), order=next(layouts))
            for ind in test_indices:
                check(A, ind)

        #check illegal access raises
        A = np.arange(12, dtype=dt).reshape((4, 3), order=next(layouts))
        szA = A.size
        illegal_indices = [szA, -szA - 1, np.array(szA), np.array(-szA - 1),
                           [szA], [-szA - 1]]
        for x in illegal_indices:
            with self.assertRaises(IndexError):
                cfunc(A, x) #  oob raises

        # check float indexing raises
        with self.assertRaises(TypingError):
            cfunc(A, [1.7])

        #exceptions leak refs
        self.disable_leak_check()

    def test_fill(self):
        pyfunc = array_fill
        cfunc = jit(nopython=True)(pyfunc)
        def check(arr, val):
            expected = np.copy(arr)
            erv = pyfunc(expected, val)
            self.assertTrue(erv is None)
            got = np.copy(arr)
            grv = cfunc(got, val)
            self.assertTrue(grv is None)
            # check mutation is the same
            self.assertPreciseEqual(expected, got)

        # scalar
        A = np.arange(1)
        for x in [np.float64, np.bool_]:
            check(A, x(10))

        # 2d
        A = np.arange(12).reshape(3, 4)
        for x in [np.float64, np.bool_]:
            check(A, x(10))

        # 4d
        A = np.arange(48, dtype=np.complex64).reshape(2, 3, 4, 2)
        for x in [np.float64, np.complex128, np.bool_]:
            check(A, x(10))

    def test_real(self):
        pyfunc = array_real
        cfunc = jit(nopython=True)(pyfunc)

        x = np.linspace(-10, 10)
        np.testing.assert_equal(pyfunc(x), cfunc(x))

        x, y = np.meshgrid(x, x)
        z = x + 1j*y
        np.testing.assert_equal(pyfunc(z), cfunc(z))

    def test_imag(self):
        pyfunc = array_imag
        cfunc = jit(nopython=True)(pyfunc)

        x = np.linspace(-10, 10)
        np.testing.assert_equal(pyfunc(x), cfunc(x))

        x, y = np.meshgrid(x, x)
        z = x + 1j*y
        np.testing.assert_equal(pyfunc(z), cfunc(z))

    def _lower_clip_result_test_util(self, func, a, a_min, a_max):
        # verifies that type-inference is working on the return value
        # this used to trigger issue #3489
        def lower_clip_result(a):
            return np.expm1(func(a, a_min, a_max))

        np.testing.assert_almost_equal(
            lower_clip_result(a),
            jit(nopython=True)(lower_clip_result)(a))

    def test_clip(self):
        has_out = (np_clip, np_clip_kwargs, array_clip, array_clip_kwargs)
        has_no_out = (np_clip_no_out, array_clip_no_out)
        # TODO: scalars are not tested (issue #3469)
        for a in (np.linspace(-10, 10, 101),
                  np.linspace(-10, 10, 40).reshape(5, 2, 4)):
            for pyfunc in has_out + has_no_out:
                cfunc = jit(nopython=True)(pyfunc)

                msg = "array_clip: must set either max or min"
                with self.assertRaisesRegex(ValueError, msg):
                    cfunc(a, None, None)

                np.testing.assert_equal(pyfunc(a, 0, None), cfunc(a, 0, None))
                np.testing.assert_equal(pyfunc(a, None, 0), cfunc(a, None, 0))

                np.testing.assert_equal(pyfunc(a, -5, 5), cfunc(a, -5, 5))

                if pyfunc in has_out:
                    pyout = np.empty_like(a)
                    cout = np.empty_like(a)
                    np.testing.assert_equal(pyfunc(a, -5, 5, pyout),
                                            cfunc(a, -5, 5, cout))
                    np.testing.assert_equal(pyout, cout)

                self._lower_clip_result_test_util(cfunc, a, -5, 5)

    def test_clip_array_min_max(self):
        has_out = (np_clip, np_clip_kwargs, array_clip, array_clip_kwargs)
        has_no_out = (np_clip_no_out, array_clip_no_out)
        # TODO: scalars are not tested (issue #3469)
        a = np.linspace(-10, 10, 40).reshape(5, 2, 4)
        a_min_arr = np.arange(-8, 0).astype(a.dtype).reshape(2, 4)
        a_max_arr = np.arange(0, 8).astype(a.dtype).reshape(2, 4)
        mins = [0, -5, a_min_arr, None]
        maxs = [0, 5, a_max_arr, None]
        for pyfunc in has_out + has_no_out:
            cfunc = jit(nopython=True)(pyfunc)

            for a_min in mins:
                for a_max in maxs:

                    if a_min is None and a_max is None:
                        msg = "array_clip: must set either max or min"
                        with self.assertRaisesRegex(ValueError, msg):
                            cfunc(a, None, None)
                        continue

                    np.testing.assert_equal(pyfunc(a, a_min, a_max), cfunc(a, a_min, a_max))

                    if pyfunc in has_out:
                        pyout = np.empty_like(a)
                        cout = np.empty_like(a)
                        np.testing.assert_equal(pyfunc(a, a_min, a_max, pyout),
                                                cfunc(a, a_min, a_max, cout))
                        np.testing.assert_equal(pyout, cout)

                    self._lower_clip_result_test_util(cfunc, a, a_min, a_max)

    def test_clip_bad_array(self):
        cfunc = jit(nopython=True)(np_clip)
        msg = '.*The argument "a" must be array-like.*'
        with self.assertRaisesRegex(TypingError, msg):
            cfunc(None, 0, 10)

    def test_clip_bad_min(self):
        cfunc = jit(nopython=True)(np_clip)
        msg = '.*The argument "a_min" must be a number.*'
        with self.assertRaisesRegex(TypingError, msg):
            cfunc(1, 'a', 10)

    def test_clip_bad_max(self):
        cfunc = jit(nopython=True)(np_clip)
        msg = '.*The argument "a_max" must be a number.*'
        with self.assertRaisesRegex(TypingError, msg):
            cfunc(1, 1, 'b')

    def test_clip_bad_out(self):
        cfunc = jit(nopython=True)(np_clip)
        msg = '.*The argument "out" must be an array if it is provided.*'
        with self.assertRaisesRegex(TypingError, msg):
            cfunc(5, 1, 10, out=6)

    def test_clip_no_broadcast(self):
        self.disable_leak_check()
        cfunc = jit(nopython=True)(np_clip)
        msg = ".*shape mismatch: objects cannot be broadcast to a single shape.*"
        a = np.linspace(-10, 10, 40).reshape(5, 2, 4)
        a_min_arr = np.arange(-5, 0).astype(a.dtype).reshape(5, 1)
        a_max_arr = np.arange(0, 5).astype(a.dtype).reshape(5, 1)
        min_max = [(0, a_max_arr), (-5, a_max_arr),
                   (a_min_arr, a_max_arr),
                   (a_min_arr, 0), (a_min_arr, 5)]
        for a_min, a_max in min_max:
            with self.assertRaisesRegex(ValueError, msg):
                cfunc(a, a_min, a_max)

    def test_conj(self):
        for pyfunc in [array_conj, array_conjugate]:
            cfunc = jit(nopython=True)(pyfunc)

            x = np.linspace(-10, 10)
            np.testing.assert_equal(pyfunc(x), cfunc(x))

            x, y = np.meshgrid(x, x)
            z = x + 1j*y
            np.testing.assert_equal(pyfunc(z), cfunc(z))

    def test_unique(self):
        pyfunc = np_unique
        cfunc = jit(nopython=True)(pyfunc)

        def check(a):
            np.testing.assert_equal(pyfunc(a), cfunc(a))

        check(np.array([[1, 1, 3], [3, 4, 5]]))
        check(np.array(np.zeros(5)))
        check(np.array([[3.1, 3.1], [1.7, 2.29], [3.3, 1.7]]))
        check(np.array([]))

    @needs_blas
    def test_array_dot(self):
        # just ensure that the dot impl dispatches correctly, do
        # not test dot itself, this is done in test_linalg.
        pyfunc = array_dot
        cfunc = jit(nopython=True)(pyfunc)
        a = np.arange(20.).reshape(4, 5)
        b = np.arange(5.)
        np.testing.assert_equal(pyfunc(a, b), cfunc(a, b))

        # check that chaining works
        pyfunc = array_dot_chain
        cfunc = jit(nopython=True)(pyfunc)
        a = np.arange(16.).reshape(4, 4)
        np.testing.assert_equal(pyfunc(a, a), cfunc(a, a))

    def test_array_ctor_with_dtype_arg(self):
        # Test using np.dtype and np.generic (i.e. np.dtype.type) has args
        pyfunc = array_ctor
        cfunc = jit(nopython=True)(pyfunc)
        n = 2
        args = n, np.int32
        np.testing.assert_array_equal(pyfunc(*args), cfunc(*args))
        args = n, np.dtype('int32')
        np.testing.assert_array_equal(pyfunc(*args), cfunc(*args))
        args = n, np.float32
        np.testing.assert_array_equal(pyfunc(*args), cfunc(*args))
        args = n, np.dtype('f4')
        np.testing.assert_array_equal(pyfunc(*args), cfunc(*args))

class TestArrayComparisons(TestCase):

    def test_identity(self):
        def check(a, b, expected):
            cfunc = njit((typeof(a), typeof(b)))(pyfunc)
            self.assertPreciseEqual(cfunc(a, b),
                                    (expected, not expected))

        pyfunc = identity_usecase

        arr = np.zeros(10, dtype=np.int32).reshape((2, 5))
        check(arr, arr, True)
        check(arr, arr[:], True)
        check(arr, arr.copy(), False)
        check(arr, arr.view('uint32'), False)
        check(arr, arr.T, False)
        check(arr, arr[:-1], False)

    # Other comparison operators ('==', etc.) are tested in test_ufuncs


if __name__ == '__main__':
    unittest.main()
