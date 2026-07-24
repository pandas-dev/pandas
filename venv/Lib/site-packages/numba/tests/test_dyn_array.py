import contextlib
import sys
import numpy as np
import random
import re
import threading
import gc

from numba.core.errors import TypingError
from numba import njit
from numba.core import types, utils, config
from numba.tests.support import MemoryLeakMixin, TestCase, tag, skip_if_32bit
from numba.core.utils import PYVERSION
import unittest


nrtjit = njit(_nrt=True, nogil=True)


def np_concatenate1(a, b, c):
    return np.concatenate((a, b, c))

def np_concatenate2(a, b, c, axis):
    return np.concatenate((a, b, c), axis=axis)

def np_stack1(a, b, c):
    return np.stack((a, b, c))

def np_stack2(a, b, c, axis):
    return np.stack((a, b, c), axis=axis)

def np_hstack(a, b, c):
    return np.hstack((a, b, c))

def np_vstack(a, b, c):
    return np.vstack((a, b, c))

def np_row_stack(a, b, c):
    return np.row_stack((a, b, c))

def np_dstack(a, b, c):
    return np.dstack((a, b, c))

def np_column_stack(a, b, c):
    return np.column_stack((a, b, c))


class BaseTest(TestCase):

    def check_outputs(self, pyfunc, argslist, exact=True):
        cfunc = nrtjit(pyfunc)
        for args in argslist:
            expected = pyfunc(*args)
            ret = cfunc(*args)
            self.assertEqual(ret.size, expected.size)
            self.assertEqual(ret.dtype, expected.dtype)
            self.assertStridesEqual(ret, expected)
            if exact:
                np.testing.assert_equal(expected, ret)
            else:
                np.testing.assert_allclose(expected, ret)


class NrtRefCtTest(MemoryLeakMixin):
    def assert_array_nrt_refct(self, arr, expect):
        self.assertEqual(arr.base.refcount, expect)


class TestDynArray(NrtRefCtTest, TestCase):

    def test_empty_0d(self):
        @nrtjit
        def foo():
            arr = np.empty(())
            arr[()] = 42
            return arr

        arr = foo()
        self.assert_array_nrt_refct(arr, 1)
        np.testing.assert_equal(42, arr)
        self.assertEqual(arr.size, 1)
        self.assertEqual(arr.shape, ())
        self.assertEqual(arr.dtype, np.dtype(np.float64))
        self.assertEqual(arr.strides, ())
        arr.fill(123)  # test writability
        np.testing.assert_equal(123, arr)
        del arr

    def test_empty_1d(self):
        @nrtjit
        def foo(n):
            arr = np.empty(n)
            for i in range(n):
                arr[i] = i

            return arr

        n = 3
        arr = foo(n)
        self.assert_array_nrt_refct(arr, 1)
        np.testing.assert_equal(np.arange(n), arr)
        self.assertEqual(arr.size, n)
        self.assertEqual(arr.shape, (n,))
        self.assertEqual(arr.dtype, np.dtype(np.float64))
        self.assertEqual(arr.strides, (np.dtype(np.float64).itemsize,))
        arr.fill(123)  # test writability
        np.testing.assert_equal(123, arr)
        del arr

    def test_empty_2d(self):
        def pyfunc(m, n):
            arr = np.empty((m, n), np.int32)
            for i in range(m):
                for j in range(n):
                    arr[i, j] = i + j

            return arr

        cfunc = nrtjit(pyfunc)
        m = 4
        n = 3
        expected_arr = pyfunc(m, n)
        got_arr = cfunc(m, n)
        self.assert_array_nrt_refct(got_arr, 1)
        np.testing.assert_equal(expected_arr, got_arr)

        self.assertEqual(expected_arr.size, got_arr.size)
        self.assertEqual(expected_arr.shape, got_arr.shape)
        self.assertEqual(expected_arr.strides, got_arr.strides)

        del got_arr

    def test_empty_3d(self):
        def pyfunc(m, n, p):
            arr = np.empty((m, n, p), np.int32)
            for i in range(m):
                for j in range(n):
                    for k in range(p):
                        arr[i, j, k] = i + j + k

            return arr

        cfunc = nrtjit(pyfunc)
        m = 4
        n = 3
        p = 2
        expected_arr = pyfunc(m, n, p)
        got_arr = cfunc(m, n, p)
        self.assert_array_nrt_refct(got_arr, 1)
        np.testing.assert_equal(expected_arr, got_arr)

        self.assertEqual(expected_arr.size, got_arr.size)
        self.assertEqual(expected_arr.shape, got_arr.shape)
        self.assertEqual(expected_arr.strides, got_arr.strides)

        del got_arr

    def test_empty_2d_sliced(self):
        def pyfunc(m, n, p):
            arr = np.empty((m, n), np.int32)
            for i in range(m):
                for j in range(n):
                    arr[i, j] = i + j

            return arr[p]

        cfunc = nrtjit(pyfunc)
        m = 4
        n = 3
        p = 2
        expected_arr = pyfunc(m, n, p)
        got_arr = cfunc(m, n, p)
        self.assert_array_nrt_refct(got_arr, 1)
        np.testing.assert_equal(expected_arr, got_arr)

        self.assertEqual(expected_arr.size, got_arr.size)
        self.assertEqual(expected_arr.shape, got_arr.shape)
        self.assertEqual(expected_arr.strides, got_arr.strides)

        del got_arr

    def test_return_global_array(self):
        y = np.ones(4, dtype=np.float32)
        initrefct = sys.getrefcount(y)

        def return_external_array():
            return y

        cfunc = nrtjit(return_external_array)
        out = cfunc()

        # out reference by cfunc
        self.assertEqual(initrefct + 1, sys.getrefcount(y))

        np.testing.assert_equal(y, out)
        np.testing.assert_equal(y, np.ones(4, dtype=np.float32))
        np.testing.assert_equal(out, np.ones(4, dtype=np.float32))

        del out
        gc.collect()
        # out is only referenced by cfunc
        self.assertEqual(initrefct + 1, sys.getrefcount(y))

        del cfunc
        gc.collect()
        # y is no longer referenced by cfunc
        self.assertEqual(initrefct, sys.getrefcount(y))

    def test_return_global_array_sliced(self):
        y = np.ones(4, dtype=np.float32)

        def return_external_array():
            return y[2:]

        cfunc = nrtjit(return_external_array)
        out = cfunc()
        self.assertIsNone(out.base)

        yy = y[2:]
        np.testing.assert_equal(yy, out)
        np.testing.assert_equal(yy, np.ones(2, dtype=np.float32))
        np.testing.assert_equal(out, np.ones(2, dtype=np.float32))

    def test_array_pass_through(self):
        def pyfunc(y):
            return y

        arr = np.ones(4, dtype=np.float32)

        cfunc = nrtjit(pyfunc)
        expected = cfunc(arr)
        got = pyfunc(arr)

        np.testing.assert_equal(expected, arr)
        np.testing.assert_equal(expected, got)
        self.assertIs(expected, arr)
        self.assertIs(expected, got)

    def test_array_pass_through_sliced(self):
        def pyfunc(y):
            return y[y.size // 2:]

        arr = np.ones(4, dtype=np.float32)

        initrefct = sys.getrefcount(arr)

        cfunc = nrtjit(pyfunc)
        got = cfunc(arr)
        self.assertEqual(initrefct + 1, sys.getrefcount(arr))
        expected = pyfunc(arr)
        self.assertEqual(initrefct + 2, sys.getrefcount(arr))

        np.testing.assert_equal(expected, arr[arr.size // 2])
        np.testing.assert_equal(expected, got)

        del expected
        self.assertEqual(initrefct + 1, sys.getrefcount(arr))
        del got
        self.assertEqual(initrefct, sys.getrefcount(arr))

    def test_ufunc_with_allocated_output(self):

        def pyfunc(a, b):
            out = np.empty(a.shape)
            np.add(a, b, out)
            return out

        cfunc = nrtjit(pyfunc)

        # 1D case
        arr_a = np.random.random(10)
        arr_b = np.random.random(10)

        np.testing.assert_equal(pyfunc(arr_a, arr_b),
                                cfunc(arr_a, arr_b))

        self.assert_array_nrt_refct(cfunc(arr_a, arr_b), 1)

        # 2D case
        arr_a = np.random.random(10).reshape(2, 5)
        arr_b = np.random.random(10).reshape(2, 5)

        np.testing.assert_equal(pyfunc(arr_a, arr_b),
                                cfunc(arr_a, arr_b))

        self.assert_array_nrt_refct(cfunc(arr_a, arr_b), 1)

        # 3D case
        arr_a = np.random.random(70).reshape(2, 5, 7)
        arr_b = np.random.random(70).reshape(2, 5, 7)

        np.testing.assert_equal(pyfunc(arr_a, arr_b),
                                cfunc(arr_a, arr_b))

        self.assert_array_nrt_refct(cfunc(arr_a, arr_b), 1)

    def test_allocation_mt(self):
        """
        This test exercises the array allocation in multithreaded usecase.
        This stress the freelist inside NRT.
        """

        def pyfunc(inp):
            out = np.empty(inp.size)

            # Zero fill
            for i in range(out.size):
                out[i] = 0

            for i in range(inp[0]):
                # Allocate inside a loop
                tmp = np.empty(inp.size)
                # Write to tmp
                for j in range(tmp.size):
                    tmp[j] = inp[j]
                # out = tmp + i
                for j in range(tmp.size):
                    out[j] += tmp[j] + i

            return out

        cfunc = nrtjit(pyfunc)
        size = 10  # small array size so that the computation is short
        arr = np.random.randint(1, 10, size)
        frozen_arr = arr.copy()

        np.testing.assert_equal(pyfunc(arr), cfunc(arr))
        # Ensure we did not modify the input
        np.testing.assert_equal(frozen_arr, arr)

        workers = []
        inputs = []
        outputs = []

        # Make wrapper to store the output
        def wrapped(inp, out):
            out[:] = cfunc(inp)

        # Create a lot of worker threads to create contention
        for i in range(100):
            arr = np.random.randint(1, 10, size)
            out = np.empty_like(arr)
            thread = threading.Thread(target=wrapped,
                                      args=(arr, out),
                                      name="worker{0}".format(i))
            workers.append(thread)
            inputs.append(arr)
            outputs.append(out)

        # Launch worker threads
        for thread in workers:
            thread.start()

        # Join worker threads
        for thread in workers:
            thread.join()

        # Check result
        for inp, out in zip(inputs, outputs):
            np.testing.assert_equal(pyfunc(inp), out)

    def test_refct_mt(self):
        """
        This test exercises the refct in multithreaded code
        """

        def pyfunc(n, inp):
            out = np.empty(inp.size)
            for i in range(out.size):
                out[i] = inp[i] + 1
            # Use swap to trigger many refct ops
            for i in range(n):
                out, inp = inp, out
            return out

        cfunc = nrtjit(pyfunc)
        size = 10
        input = np.arange(size, dtype=float)
        expected_refct = sys.getrefcount(input)
        swapct = random.randrange(1000)
        expected = pyfunc(swapct, input)
        np.testing.assert_equal(expected, cfunc(swapct, input))
        # The following checks can discover a reference count error
        del expected
        self.assertEqual(expected_refct, sys.getrefcount(input))

        workers = []
        outputs = []
        swapcts = []

        # Make wrapper to store the output
        def wrapped(n, input, out):
            out[:] = cfunc(n, input)

        # Create worker threads
        for i in range(100):
            out = np.empty(size)
            # All thread shares the same input
            swapct = random.randrange(1000)
            thread = threading.Thread(target=wrapped,
                                      args=(swapct, input, out),
                                      name="worker{0}".format(i))
            workers.append(thread)
            outputs.append(out)
            swapcts.append(swapct)

        # Launch worker threads
        for thread in workers:
            thread.start()

        # Join worker threads
        for thread in workers:
            thread.join()

        # Check result
        for swapct, out in zip(swapcts, outputs):
            np.testing.assert_equal(pyfunc(swapct, input), out)

        del outputs, workers
        # The following checks can discover a reference count error
        self.assertEqual(expected_refct, sys.getrefcount(input))

    @skip_if_32bit
    def test_invalid_size_array(self):

        @njit
        def foo(x):
            np.empty(x)

        # Exceptions leak references
        self.disable_leak_check()

        with self.assertRaises(MemoryError) as raises:
            foo(types.size_t.maxval // 8 // 2)

        self.assertIn("Allocation failed", str(raises.exception))

    def test_swap(self):

        def pyfunc(x, y, t):
            """Swap array x and y for t number of times
            """
            for i in range(t):
                x, y = y, x

            return x, y


        cfunc = nrtjit(pyfunc)

        x = np.random.random(100)
        y = np.random.random(100)

        t = 100

        initrefct = sys.getrefcount(x), sys.getrefcount(y)
        expect, got = pyfunc(x, y, t), cfunc(x, y, t)
        self.assertIsNone(got[0].base)
        self.assertIsNone(got[1].base)
        np.testing.assert_equal(expect, got)
        del expect, got
        self.assertEqual(initrefct, (sys.getrefcount(x), sys.getrefcount(y)))

    def test_return_tuple_of_array(self):

        def pyfunc(x):
            y = np.empty(x.size)
            for i in range(y.size):
                y[i] = x[i] + 1
            return x, y

        cfunc = nrtjit(pyfunc)

        x = np.random.random(5)
        initrefct = sys.getrefcount(x)
        expected_x, expected_y = pyfunc(x)
        got_x, got_y = cfunc(x)
        self.assertIs(x, expected_x)
        self.assertIs(x, got_x)
        np.testing.assert_equal(expected_x, got_x)
        np.testing.assert_equal(expected_y, got_y)
        del expected_x, got_x
        self.assertEqual(initrefct, sys.getrefcount(x))

        self.assertEqual(sys.getrefcount(expected_y), sys.getrefcount(got_y))

    def test_return_tuple_of_array_created(self):

        def pyfunc(x):
            y = np.empty(x.size)
            for i in range(y.size):
                y[i] = x[i] + 1
            out = y, y
            return out

        cfunc = nrtjit(pyfunc)

        x = np.random.random(5)
        expected_x, expected_y = pyfunc(x)
        got_x, got_y = cfunc(x)
        np.testing.assert_equal(expected_x, got_x)
        np.testing.assert_equal(expected_y, got_y)

        if PYVERSION in ((3, 14), ):
            expected_refcount = 1
        elif PYVERSION in ((3, 10), (3, 11), (3, 12), (3, 13)):
            expected_refcount = 2
        else:
            raise NotImplementedError(PYVERSION)

        self.assertEqual(expected_refcount, sys.getrefcount(got_y))
        self.assertEqual(expected_refcount, sys.getrefcount(got_x))

    def test_issue_with_return_leak(self):
        """
        Dispatcher returns a new reference.
        It need to workaround it for now.
        """
        @nrtjit
        def inner(out):
            return out

        def pyfunc(x):
            return inner(x)

        cfunc = nrtjit(pyfunc)

        arr = np.arange(10)
        old_refct = sys.getrefcount(arr)

        if PYVERSION in ((3, 14), ):
            self.assertEqual(2, sys.getrefcount(pyfunc(arr)))
            self.assertEqual(2, sys.getrefcount(cfunc(arr)))
        elif PYVERSION in ((3, 10), (3, 11), (3, 12), (3, 13)):
            self.assertEqual(old_refct, sys.getrefcount(pyfunc(arr)))
            self.assertEqual(old_refct, sys.getrefcount(cfunc(arr)))
        else:
            raise NotImplementedError(PYVERSION)

        self.assertEqual(old_refct, sys.getrefcount(arr))


class ConstructorBaseTest(NrtRefCtTest):

    def check_0d(self, pyfunc):
        cfunc = nrtjit(pyfunc)
        expected = pyfunc()
        ret = cfunc()
        self.assert_array_nrt_refct(ret, 1)
        self.assertEqual(ret.size, expected.size)
        self.assertEqual(ret.shape, expected.shape)
        self.assertEqual(ret.dtype, expected.dtype)
        self.assertEqual(ret.strides, expected.strides)
        self.check_result_value(ret, expected)
        # test writability
        expected = np.empty_like(ret) # np.full_like was not added until Numpy 1.8
        expected.fill(123)
        ret.fill(123)
        np.testing.assert_equal(ret, expected)

    def check_1d(self, pyfunc):
        cfunc = nrtjit(pyfunc)
        n = 3
        expected = pyfunc(n)
        ret = cfunc(n)
        self.assert_array_nrt_refct(ret, 1)
        self.assertEqual(ret.size, expected.size)
        self.assertEqual(ret.shape, expected.shape)
        self.assertEqual(ret.dtype, expected.dtype)
        self.assertEqual(ret.strides, expected.strides)
        self.check_result_value(ret, expected)
        # test writability
        expected = np.empty_like(ret) # np.full_like was not added until Numpy 1.8
        expected.fill(123)
        ret.fill(123)
        np.testing.assert_equal(ret, expected)
        # errors
        with self.assertRaises(ValueError) as cm:
            cfunc(-1)
        self.assertEqual(str(cm.exception), "negative dimensions not allowed")

    def check_2d(self, pyfunc):
        cfunc = nrtjit(pyfunc)
        m, n = 2, 3
        expected = pyfunc(m, n)
        ret = cfunc(m, n)
        self.assert_array_nrt_refct(ret, 1)
        self.assertEqual(ret.size, expected.size)
        self.assertEqual(ret.shape, expected.shape)
        self.assertEqual(ret.dtype, expected.dtype)
        self.assertEqual(ret.strides, expected.strides)
        self.check_result_value(ret, expected)
        # test writability
        expected = np.empty_like(ret)  # np.full_like was not added until Numpy 1.8
        expected.fill(123)
        ret.fill(123)
        np.testing.assert_equal(ret, expected)
        # errors
        with self.assertRaises(ValueError) as cm:
            cfunc(2, -1)
        self.assertEqual(str(cm.exception), "negative dimensions not allowed")

    def check_alloc_size(self, pyfunc):
        """Checks that pyfunc will error, not segfaulting due to array size."""
        cfunc = nrtjit(pyfunc)
        with self.assertRaises(ValueError) as e:
            cfunc()
        self.assertIn(
            "array is too big",
            str(e.exception)
        )


class TestNdZeros(ConstructorBaseTest, TestCase):

    def setUp(self):
        super(TestNdZeros, self).setUp()
        self.pyfunc = np.zeros

    def check_result_value(self, ret, expected):
        np.testing.assert_equal(ret, expected)

    def test_0d(self):
        pyfunc = self.pyfunc
        def func():
            return pyfunc(())
        self.check_0d(func)

    def test_1d(self):
        pyfunc = self.pyfunc
        def func(n):
            return pyfunc(n)
        self.check_1d(func)

    def test_1d_dtype(self):
        pyfunc = self.pyfunc
        def func(n):
            return pyfunc(n, np.int32)
        self.check_1d(func)

    def test_1d_dtype_instance(self):
        # dtype as numpy dtype, not as scalar class
        pyfunc = self.pyfunc
        _dtype = np.dtype('int32')
        def func(n):
            return pyfunc(n, _dtype)
        self.check_1d(func)

    def test_1d_dtype_str(self):
        pyfunc = self.pyfunc
        _dtype = 'int32'
        def func(n):
            return pyfunc(n, _dtype)
        self.check_1d(func)

        def func(n):
            return pyfunc(n, 'complex128')
        self.check_1d(func)

    def test_1d_dtype_str_alternative_spelling(self):
        # like test_1d_dtype_str but using the shorthand type spellings
        pyfunc = self.pyfunc
        _dtype = 'i4'
        def func(n):
            return pyfunc(n, _dtype)
        self.check_1d(func)

        def func(n):
            return pyfunc(n, 'c8')
        self.check_1d(func)

    def test_1d_dtype_str_structured_dtype(self):
        # test_1d_dtype_str but using a structured dtype
        pyfunc = self.pyfunc
        _dtype = "i4, (2,3)f8"
        def func(n):
            return pyfunc(n, _dtype)
        self.check_1d(func)

    def test_1d_dtype_non_const_str(self):
        pyfunc = self.pyfunc

        @njit
        def func(n, dt):
            return pyfunc(n, dt)

        with self.assertRaises(TypingError) as raises:
            func(5, 'int32')

        excstr = str(raises.exception)
        msg = (f"If np.{self.pyfunc.__name__} dtype is a string it must be a "
               "string constant.")
        self.assertIn(msg, excstr)

    def test_1d_dtype_invalid_str(self):
        pyfunc = self.pyfunc

        @njit
        def func(n):
            return pyfunc(n, 'ABCDEF')

        with self.assertRaises(TypingError) as raises:
            func(5)

        excstr = str(raises.exception)
        self.assertIn("Invalid NumPy dtype specified: 'ABCDEF'", excstr)

    def test_2d(self):
        pyfunc = self.pyfunc
        def func(m, n):
            return pyfunc((m, n))
        self.check_2d(func)

    def test_2d_shape_dtypes(self):
        # Test for issue #4575
        pyfunc = self.pyfunc
        def func1(m, n):
            return pyfunc((np.int16(m), np.int32(n)))
        self.check_2d(func1)
        # Using a 64-bit value checks that 32 bit systems will downcast to intp
        def func2(m, n):
            return pyfunc((np.int64(m), np.int8(n)))
        self.check_2d(func2)
        # Make sure an error is thrown if we can't downcast safely
        if config.IS_32BITS:
            cfunc = nrtjit(lambda m, n: pyfunc((m, n)))
            with self.assertRaises(ValueError):
                cfunc(np.int64(1 << (32 - 1)), 1)

    def test_2d_dtype_kwarg(self):
        pyfunc = self.pyfunc
        def func(m, n):
            return pyfunc((m, n), dtype=np.complex64)
        self.check_2d(func)

    def test_2d_dtype_str_kwarg(self):
        pyfunc = self.pyfunc
        def func(m, n):
            return pyfunc((m, n), dtype='complex64')
        self.check_2d(func)

    def test_2d_dtype_str_kwarg_alternative_spelling(self):
        # as test_2d_dtype_str_kwarg but with the numpy shorthand type spelling
        pyfunc = self.pyfunc
        def func(m, n):
            return pyfunc((m, n), dtype='c8')
        self.check_2d(func)

    def test_alloc_size(self):
        pyfunc = self.pyfunc
        width = types.intp.bitwidth
        def gen_func(shape, dtype):
            return lambda : pyfunc(shape, dtype)
        # Under these values numba will segfault, but that's another issue
        self.check_alloc_size(gen_func(1 << width - 2, np.intp))
        self.check_alloc_size(gen_func((1 << width - 8, 64), np.intp))


class TestNdOnes(TestNdZeros):

    def setUp(self):
        super(TestNdOnes, self).setUp()
        self.pyfunc = np.ones

    @unittest.expectedFailure
    def test_1d_dtype_str_structured_dtype(self):
        super().test_1d_dtype_str_structured_dtype()


class TestNdFull(ConstructorBaseTest, TestCase):

    def check_result_value(self, ret, expected):
        np.testing.assert_equal(ret, expected)

    def test_0d(self):
        def func():
            return np.full((), 4.5)
        self.check_0d(func)

    def test_1d(self):
        def func(n):
            return np.full(n, 4.5)
        self.check_1d(func)

    def test_1d_dtype(self):
        def func(n):
            return np.full(n, 4.5, np.bool_)
        self.check_1d(func)

    def test_1d_dtype_instance(self):
        dtype = np.dtype('bool')
        def func(n):
            return np.full(n, 4.5, dtype)
        self.check_1d(func)

    def test_1d_dtype_str(self):
        def func(n):
            return np.full(n, 4.5, 'bool_')
        self.check_1d(func)

    def test_1d_dtype_str_alternative_spelling(self):
        # like test_1d_dtype_str but using the shorthand type spelling
        def func(n):
            return np.full(n, 4.5, '?')
        self.check_1d(func)

    def test_1d_dtype_non_const_str(self):

        @njit
        def func(n, fv, dt):
            return np.full(n, fv, dt)

        with self.assertRaises(TypingError) as raises:
            func((5,), 4.5, 'int32')

        excstr = str(raises.exception)
        msg = ("If np.full dtype is a string it must be a "
               "string constant.")
        self.assertIn(msg, excstr)

    def test_1d_dtype_invalid_str(self):

        @njit
        def func(n, fv):
            return np.full(n, fv, 'ABCDEF')

        with self.assertRaises(TypingError) as raises:
            func((5,), 4.5)

        excstr = str(raises.exception)
        self.assertIn("Invalid NumPy dtype specified: 'ABCDEF'", excstr)

    def test_2d(self):
        def func(m, n):
            return np.full((m, n), 4.5)
        self.check_2d(func)

    def test_2d_dtype_kwarg(self):
        def func(m, n):
            return np.full((m, n), 1 + 4.5j, dtype=np.complex64)
        self.check_2d(func)

    def test_2d_dtype_from_type(self):
        # tests issue #2862
        def func(m, n):
            return np.full((m, n), np.int32(1))
        self.check_2d(func)

        # Complex uses `.real`, imaginary part dropped
        def func(m, n):
            return np.full((m, n), np.complex128(1))
        self.check_2d(func)

        # and that if a dtype is specified, this influences the return type
        def func(m, n):
            return np.full((m, n), 1, dtype=np.int8)
        self.check_2d(func)

    def test_2d_shape_dtypes(self):
        # Test for issue #4575
        def func1(m, n):
            return np.full((np.int16(m), np.int32(n)), 4.5)
        self.check_2d(func1)
        # Using a 64-bit value checks that 32 bit systems will downcast to intp
        def func2(m, n):
            return np.full((np.int64(m), np.int8(n)), 4.5)
        self.check_2d(func2)
        # Make sure an error is thrown if we can't downcast safely
        if config.IS_32BITS:
            cfunc = nrtjit(lambda m, n: np.full((m, n), 4.5))
            with self.assertRaises(ValueError):
                cfunc(np.int64(1 << (32 - 1)), 1)

    def test_alloc_size(self):
        width = types.intp.bitwidth
        def gen_func(shape, value):
            return lambda : np.full(shape, value)
        # Under these values numba will segfault, but that's another issue
        self.check_alloc_size(gen_func(1 << width - 2, 1))
        self.check_alloc_size(gen_func((1 << width - 8, 64), 1))


class ConstructorLikeBaseTest(object):

    def mutate_array(self, arr):
        try:
            arr.fill(42)
        except (TypeError, ValueError):
            # Try something else (e.g. Numpy 1.6 with structured dtypes)
            fill_value = b'x' * arr.dtype.itemsize
            arr.fill(fill_value)

    def check_like(self, pyfunc, dtype):
        def check_arr(arr):
            expected = pyfunc(arr)
            ret = cfunc(arr)
            self.assertEqual(ret.size, expected.size)
            self.assertEqual(ret.dtype, expected.dtype)
            self.assertStridesEqual(ret, expected)
            self.check_result_value(ret, expected)
            # test writability
            self.mutate_array(ret)
            self.mutate_array(expected)
            np.testing.assert_equal(ret, expected)

        orig = np.linspace(0, 5, 6).astype(dtype)
        cfunc = nrtjit(pyfunc)

        for shape in (6, (2, 3), (1, 2, 3), (3, 1, 2), ()):
            if shape == ():
                arr = orig[-1:].reshape(())
            else:
                arr = orig.reshape(shape)
            check_arr(arr)
            # Non-contiguous array
            if arr.ndim > 0:
                check_arr(arr[::2])
            # Check new array doesn't inherit readonly flag
            arr.flags['WRITEABLE'] = False
            # verify read-only
            with self.assertRaises(ValueError):
                arr[0] = 1
            check_arr(arr)

        # Scalar argument => should produce a 0-d array
        check_arr(orig[0])


class TestNdEmptyLike(ConstructorLikeBaseTest, TestCase):

    def setUp(self):
        super(TestNdEmptyLike, self).setUp()
        self.pyfunc = np.empty_like

    def check_result_value(self, ret, expected):
        pass

    def test_like(self):
        pyfunc = self.pyfunc
        def func(arr):
            return pyfunc(arr)
        self.check_like(func, np.float64)

    def test_like_structured(self):
        dtype = np.dtype([('a', np.int16), ('b', np.float32)])
        pyfunc = self.pyfunc
        def func(arr):
            return pyfunc(arr)
        self.check_like(func, dtype)

    def test_like_dtype(self):
        pyfunc = self.pyfunc
        def func(arr):
            return pyfunc(arr, np.int32)
        self.check_like(func, np.float64)

    def test_like_dtype_instance(self):
        dtype = np.dtype('int32')
        pyfunc = self.pyfunc
        def func(arr):
            return pyfunc(arr, dtype)
        self.check_like(func, np.float64)

    def test_like_dtype_structured(self):
        dtype = np.dtype([('a', np.int16), ('b', np.float32)])
        pyfunc = self.pyfunc
        def func(arr):
            return pyfunc(arr, dtype)
        self.check_like(func, np.float64)

    def test_like_dtype_kwarg(self):
        pyfunc = self.pyfunc
        def func(arr):
            return pyfunc(arr, dtype=np.int32)
        self.check_like(func, np.float64)

    def test_like_dtype_str_kwarg(self):
        pyfunc = self.pyfunc
        def func(arr):
            return pyfunc(arr, dtype='int32')
        self.check_like(func, np.float64)

    def test_like_dtype_str_kwarg_alternative_spelling(self):
        pyfunc = self.pyfunc
        def func(arr):
            return pyfunc(arr, dtype='i4')
        self.check_like(func, np.float64)

    def test_like_dtype_non_const_str(self):
        pyfunc = self.pyfunc

        @njit
        def func(n, dt):
            return pyfunc(n, dt)

        with self.assertRaises(TypingError) as raises:
            func(np.ones(4), 'int32')

        excstr = str(raises.exception)
        msg = (f"If np.{self.pyfunc.__name__} dtype is a string it must be a "
               "string constant.")
        self.assertIn(msg, excstr)
        self.assertIn(
            '{}(array(float64, 1d, C), unicode_type)'.format(pyfunc.__name__),
            excstr)

    def test_like_dtype_invalid_str(self):
        pyfunc = self.pyfunc

        @njit
        def func(n):
            return pyfunc(n, 'ABCDEF')

        with self.assertRaises(TypingError) as raises:
            func(np.ones(4))

        excstr = str(raises.exception)
        self.assertIn("Invalid NumPy dtype specified: 'ABCDEF'", excstr)


class TestNdZerosLike(TestNdEmptyLike):

    def setUp(self):
        super(TestNdZerosLike, self).setUp()
        self.pyfunc = np.zeros_like

    def check_result_value(self, ret, expected):
        np.testing.assert_equal(ret, expected)

    def test_like_structured(self):
        super(TestNdZerosLike, self).test_like_structured()

    def test_like_dtype_structured(self):
        super(TestNdZerosLike, self).test_like_dtype_structured()


class TestNdOnesLike(TestNdZerosLike):

    def setUp(self):
        super(TestNdOnesLike, self).setUp()
        self.pyfunc = np.ones_like
        self.expected_value = 1

    # Not supported yet.

    @unittest.expectedFailure
    def test_like_structured(self):
        super(TestNdOnesLike, self).test_like_structured()

    @unittest.expectedFailure
    def test_like_dtype_structured(self):
        super(TestNdOnesLike, self).test_like_dtype_structured()


class TestNdFullLike(ConstructorLikeBaseTest, TestCase):

    def check_result_value(self, ret, expected):
        np.testing.assert_equal(ret, expected)

    def test_like(self):
        def func(arr):
            return np.full_like(arr, 3.5)
        self.check_like(func, np.float64)

    # Not supported yet.
    @unittest.expectedFailure
    def test_like_structured(self):
        dtype = np.dtype([('a', np.int16), ('b', np.float32)])
        def func(arr):
            return np.full_like(arr, 4.5)
        self.check_like(func, dtype)

    def test_like_dtype(self):
        def func(arr):
            return np.full_like(arr, 4.5, np.bool_)
        self.check_like(func, np.float64)

    def test_like_dtype_instance(self):
        dtype = np.dtype('bool')
        def func(arr):
            return np.full_like(arr, 4.5, dtype)
        self.check_like(func, np.float64)

    def test_like_dtype_kwarg(self):
        def func(arr):
            return np.full_like(arr, 4.5, dtype=np.bool_)
        self.check_like(func, np.float64)

    def test_like_dtype_str_kwarg(self):
        def func(arr):
            return np.full_like(arr, 4.5, 'bool_')
        self.check_like(func, np.float64)

    def test_like_dtype_str_kwarg_alternative_spelling(self):
        def func(arr):
            return np.full_like(arr, 4.5, dtype='?')
        self.check_like(func, np.float64)

    def test_like_dtype_non_const_str_kwarg(self):

        @njit
        def func(arr, fv, dt):
            return np.full_like(arr, fv, dt)

        with self.assertRaises(TypingError) as raises:
            func(np.ones(3,), 4.5, 'int32')

        excstr = str(raises.exception)
        msg = ("If np.full_like dtype is a string it must be a "
               "string constant.")
        self.assertIn(msg, excstr)

    def test_like_dtype_invalid_str(self):

        @njit
        def func(arr, fv):
            return np.full_like(arr, fv, "ABCDEF")

        with self.assertRaises(TypingError) as raises:
            func(np.ones(4), 3.4)

        excstr = str(raises.exception)
        self.assertIn("Invalid NumPy dtype specified: 'ABCDEF'", excstr)


class TestNdIdentity(BaseTest):

    def check_identity(self, pyfunc):
        self.check_outputs(pyfunc, [(3,)])

    def test_identity(self):
        def func(n):
            return np.identity(n)
        self.check_identity(func)

    def test_identity_dtype(self):
        for dtype in (np.complex64, np.int16, np.bool_, np.dtype('bool'),
                      'bool_'):
            def func(n):
                return np.identity(n, dtype)
            self.check_identity(func)

    def test_like_dtype_non_const_str_kwarg(self):

        @njit
        def func(n, dt):
            return np.identity(n, dt)

        with self.assertRaises(TypingError) as raises:
            func(4, 'int32')

        excstr = str(raises.exception)
        msg = ("If np.identity dtype is a string it must be a "
               "string constant.")
        self.assertIn(msg, excstr)



class TestNdEye(BaseTest):

    def test_eye_n(self):
        def func(n):
            return np.eye(n)
        self.check_outputs(func, [(1,), (3,)])

    def test_eye_n_dtype(self):
        # check None option, dtype class, instance of dtype class
        for dt in (None, np.complex128, np.complex64(1)):
            def func(n, dtype=dt):
                return np.eye(n, dtype=dtype)
            self.check_outputs(func, [(1,), (3,)])

    def test_eye_n_m(self):
        def func(n, m):
            return np.eye(n, m)
        self.check_outputs(func, [(1, 2), (3, 2), (0, 3)])

    def check_eye_n_m_k(self, func):
        self.check_outputs(func, [(1, 2, 0),
                                  (3, 4, 1),
                                  (3, 4, -1),
                                  (4, 3, -2),
                                  (4, 3, -5),
                                  (4, 3, 5)])

    def test_eye_n_m_k(self):
        def func(n, m, k):
            return np.eye(n, m, k)
        self.check_eye_n_m_k(func)

    def test_eye_n_m_k_dtype(self):
        def func(n, m, k):
            return np.eye(N=n, M=m, k=k, dtype=np.int16)
        self.check_eye_n_m_k(func)

    def test_eye_n_m_k_dtype_instance(self):
        dtype = np.dtype('int16')
        def func(n, m, k):
            return np.eye(N=n, M=m, k=k, dtype=dtype)
        self.check_eye_n_m_k(func)


class TestNdDiag(TestCase):

    def setUp(self):
        v = np.array([1, 2, 3])
        hv = np.array([[1, 2, 3]])
        vv = np.transpose(hv)
        self.vectors = [v, hv, vv]
        a3x4 = np.arange(12).reshape(3, 4)
        a4x3 = np.arange(12).reshape(4, 3)
        self.matricies = [a3x4, a4x3]
        def func(q):
            return np.diag(q)
        self.py = func
        self.jit = nrtjit(func)

        def func_kwarg(q, k=0):
            return np.diag(q, k=k)
        self.py_kw = func_kwarg
        self.jit_kw = nrtjit(func_kwarg)

    def check_diag(self, pyfunc, nrtfunc, *args, **kwargs):
        expected = pyfunc(*args, **kwargs)
        computed = nrtfunc(*args, **kwargs)
        self.assertEqual(computed.size, expected.size)
        self.assertEqual(computed.dtype, expected.dtype)
        # NOTE: stride not tested as np returns a RO view, nb returns new data
        np.testing.assert_equal(expected, computed)

    # create a diag matrix from a vector
    def test_diag_vect_create(self):
        for d in self.vectors:
            self.check_diag(self.py, self.jit, d)

    # create a diag matrix from a vector at a given offset
    def test_diag_vect_create_kwarg(self):
        for k in range(-10, 10):
            for d in self.vectors:
                self.check_diag(self.py_kw, self.jit_kw, d, k=k)

    # extract the diagonal
    def test_diag_extract(self):
        for d in self.matricies:
            self.check_diag(self.py, self.jit, d)

    # extract a diagonal at a given offset
    def test_diag_extract_kwarg(self):
        for k in range(-4, 4):
            for d in self.matricies:
                self.check_diag(self.py_kw, self.jit_kw, d, k=k)

    # check error handling
    def test_error_handling(self):
        d = np.array([[[1.]]])
        cfunc = nrtjit(self.py)

        # missing arg
        with self.assertRaises(TypeError):
            cfunc()

        # > 2d
        with self.assertRaises(TypingError):
            cfunc(d)
        with self.assertRaises(TypingError):
            dfunc = nrtjit(self.py_kw)
            dfunc(d, k=3)

    def test_bad_shape(self):
        cfunc = nrtjit(self.py)
        msg = '.*The argument "v" must be array-like.*'
        with self.assertRaisesRegex(TypingError, msg) as raises:
            cfunc(None)

class TestLinspace(BaseTest):

    def test_linspace_2(self):
        def pyfunc(n, m):
            return np.linspace(n, m)
        self.check_outputs(pyfunc,
                           [(0, 4), (1, 100), (-3.5, 2.5), (-3j, 2+3j),
                            (2, 1), (1+0.5j, 1.5j)])

    def test_linspace_3(self):
        def pyfunc(n, m, p):
            return np.linspace(n, m, p)
        self.check_outputs(pyfunc,
                           [(0, 4, 9), (1, 4, 3), (-3.5, 2.5, 8),
                            (-3j, 2+3j, 7), (2, 1, 0),
                            (1+0.5j, 1.5j, 5), (1, 1e100, 1)])

    def test_linspace_accuracy(self):
        # Checking linspace reasonably replicates NumPy's algorithm
        # see https://github.com/numba/numba/issues/6768
        @nrtjit
        def foo(n, m, p):
            return np.linspace(n, m, p)

        n, m, p = 0.0, 1.0, 100
        self.assertPreciseEqual(foo(n, m, p), foo.py_func(n, m, p))


class TestNpyEmptyKeyword(TestCase):
    def _test_with_dtype_kw(self, dtype):
        def pyfunc(shape):
            return np.empty(shape, dtype=dtype)

        shapes = [1, 5, 9]

        cfunc = nrtjit(pyfunc)
        for s in shapes:
            expected = pyfunc(s)
            got = cfunc(s)
            self.assertEqual(expected.dtype, got.dtype)
            self.assertEqual(expected.shape, got.shape)

    def test_with_dtype_kws(self):
        for dtype in [np.int32, np.float32, np.complex64, np.dtype('complex64')]:
            self._test_with_dtype_kw(dtype)

    def _test_with_shape_and_dtype_kw(self, dtype):
        def pyfunc(shape):
            return np.empty(shape=shape, dtype=dtype)

        shapes = [1, 5, 9]

        cfunc = nrtjit(pyfunc)
        for s in shapes:
            expected = pyfunc(s)
            got = cfunc(s)
            self.assertEqual(expected.dtype, got.dtype)
            self.assertEqual(expected.shape, got.shape)

    def test_with_shape_and_dtype_kws(self):
        for dtype in [np.int32, np.float32, np.complex64, np.dtype('complex64')]:
            self._test_with_shape_and_dtype_kw(dtype)

    def test_empty_no_args(self):

        def pyfunc():
            return np.empty()

        cfunc = nrtjit(pyfunc)

        # Trigger the compilation
        # That will cause a TypingError due to missing shape argument
        with self.assertRaises(TypingError):
            cfunc()


class TestNpArray(MemoryLeakMixin, BaseTest):

    def test_0d(self):
        def pyfunc(arg):
            return np.array(arg)

        cfunc = nrtjit(pyfunc)
        got = cfunc(42)
        self.assertPreciseEqual(got, np.array(42, dtype=np.intp))
        got = cfunc(2.5)
        self.assertPreciseEqual(got, np.array(2.5))

    def test_0d_with_dtype(self):
        def pyfunc(arg):
            return np.array(arg, dtype=np.int16)

        self.check_outputs(pyfunc, [(42,), (3.5,)])

    def test_1d(self):
        def pyfunc(arg):
            return np.array(arg)

        cfunc = nrtjit(pyfunc)
        # A list
        got = cfunc([2, 3, 42])
        self.assertPreciseEqual(got, np.intp([2, 3, 42]))
        # A heterogeneous tuple
        got = cfunc((1.0, 2.5j, 42))
        self.assertPreciseEqual(got, np.array([1.0, 2.5j, 42]))
        # An empty tuple
        got = cfunc(())
        self.assertPreciseEqual(got, np.float64(()))

    def test_1d_with_dtype(self):
        def pyfunc(arg):
            return np.array(arg, dtype=np.float32)

        self.check_outputs(pyfunc,
                           [([2, 42],),
                            ([3.5, 1.0],),
                            ((1, 3.5, 42),),
                            ((),),
                            ])

    def test_1d_with_str_dtype(self):
        def pyfunc(arg):
            return np.array(arg, dtype='float32')

        self.check_outputs(pyfunc,
                           [([2, 42],),
                            ([3.5, 1.0],),
                            ((1, 3.5, 42),),
                            ((),),
                            ])

    def test_1d_with_non_const_str_dtype(self):

        @njit
        def func(arg, dt):
            return np.array(arg, dtype=dt)

        with self.assertRaises(TypingError) as raises:
            func((5, 3), 'int32')

        excstr = str(raises.exception)
        msg = (f"If np.array dtype is a string it must be a "
               "string constant.")
        self.assertIn(msg, excstr)

    def test_2d(self):
        def pyfunc(arg):
            return np.array(arg)

        cfunc = nrtjit(pyfunc)
        # A list of tuples
        got = cfunc([(1, 2), (3, 4)])
        self.assertPreciseEqual(got, np.intp([[1, 2], [3, 4]]))
        got = cfunc([(1, 2.5), (3, 4.5)])
        self.assertPreciseEqual(got, np.float64([[1, 2.5], [3, 4.5]]))
        # A tuple of lists
        got = cfunc(([1, 2], [3, 4]))
        self.assertPreciseEqual(got, np.intp([[1, 2], [3, 4]]))
        got = cfunc(([1, 2], [3.5, 4.5]))
        self.assertPreciseEqual(got, np.float64([[1, 2], [3.5, 4.5]]))
        # A tuple of tuples
        got = cfunc(((1.5, 2), (3.5, 4.5)))
        self.assertPreciseEqual(got, np.float64([[1.5, 2], [3.5, 4.5]]))
        got = cfunc(((), ()))
        self.assertPreciseEqual(got, np.float64(((), ())))

    def test_2d_with_dtype(self):
        def pyfunc(arg):
            return np.array(arg, dtype=np.int32)

        cfunc = nrtjit(pyfunc)
        got = cfunc([(1, 2.5), (3, 4.5)])
        self.assertPreciseEqual(got, np.int32([[1, 2], [3, 4]]))

    def test_raises(self):

        def pyfunc(arg):
            return np.array(arg)

        cfunc = nrtjit(pyfunc)

        @contextlib.contextmanager
        def check_raises(msg):
            with self.assertRaises(TypingError) as raises:
                yield
            self.assertIn(msg, str(raises.exception))

        with check_raises(('array(float64, 1d, C) not allowed in a '
                           'homogeneous sequence')):
            cfunc(np.array([1.]))

        with check_raises(('type Tuple(int64, reflected list(int64)<iv=None>) '
                          'does not have a regular shape')):
            cfunc((np.int64(1), [np.int64(2)]))

        with check_raises(
                "cannot convert Tuple(int64, Record(a[type=int32;offset=0],"
                "b[type=float32;offset=4];8;False)) to a homogeneous type",
            ):
            st = np.dtype([('a', 'i4'), ('b', 'f4')])
            val = np.zeros(1, dtype=st)[0]
            cfunc(((1, 2), (np.int64(1), val)))

    def test_bad_array(self):
        @njit
        def func(obj):
            return np.array(obj)

        msg = '.*The argument "object" must be array-like.*'
        with self.assertRaisesRegex(TypingError, msg) as raises:
            func(None)

    def test_bad_dtype(self):
        @njit
        def func(obj, dt):
            return np.array(obj, dt)

        msg = '.*The argument "dtype" must be a data-type if it is provided.*'
        with self.assertRaisesRegex(TypingError, msg) as raises:
            func(5, 4)


class TestNpConcatenate(MemoryLeakMixin, TestCase):
    """
    Tests for np.concatenate().
    """

    def _3d_arrays(self):
        a = np.arange(24).reshape((4, 3, 2))
        b = a + 10
        c = (b + 10).copy(order='F')
        d = (c + 10)[::-1]
        e = (d + 10)[...,::-1]
        return a, b, c, d, e

    @contextlib.contextmanager
    def assert_invalid_sizes_over_dim(self, axis):
        with self.assertRaises(ValueError) as raises:
            yield
        self.assertIn("input sizes over dimension %d do not match" % axis,
                      str(raises.exception))

    def test_3d(self):
        pyfunc = np_concatenate2
        cfunc = nrtjit(pyfunc)

        def check(a, b, c, axis):
            for ax in (axis, -3 + axis):
                expected = pyfunc(a, b, c, axis=ax)
                got = cfunc(a, b, c, axis=ax)
                self.assertPreciseEqual(got, expected)

        def check_all_axes(a, b, c):
            for axis in range(3):
                check(a, b, c, axis)

        a, b, c, d, e = self._3d_arrays()

        # Inputs with equal sizes
        # C, C, C
        check_all_axes(a, b, b)
        # C, C, F
        check_all_axes(a, b, c)
        # F, F, F
        check_all_axes(a.T, b.T, a.T)
        # F, F, C
        check_all_axes(a.T, b.T, c.T)
        # F, F, A
        check_all_axes(a.T, b.T, d.T)
        # A, A, A
        # (note Numpy may select the layout differently for other inputs)
        check_all_axes(d.T, e.T, d.T)

        # Inputs with compatible sizes
        check(a[1:], b, c[::-1], axis=0)
        check(a, b[:,1:], c, axis=1)
        check(a, b, c[:,:,1:], axis=2)

        # Different but compatible dtypes
        check_all_axes(a, b.astype(np.float64), b)

        # Exceptions leak references
        self.disable_leak_check()

        # Incompatible sizes
        for axis in (1, 2, -2, -1):
            with self.assert_invalid_sizes_over_dim(0):
                cfunc(a[1:], b, b, axis)
        for axis in (0, 2, -3, -1):
            with self.assert_invalid_sizes_over_dim(1):
                cfunc(a, b[:,1:], b, axis)

    def test_3d_no_axis(self):
        pyfunc = np_concatenate1
        cfunc = nrtjit(pyfunc)

        def check(a, b, c):
            expected = pyfunc(a, b, c)
            got = cfunc(a, b, c)
            self.assertPreciseEqual(got, expected)

        a, b, c, d, e = self._3d_arrays()

        # Inputs with equal sizes
        # C, C, C
        check(a, b, b)
        # C, C, F
        check(a, b, c)
        # F, F, F
        check(a.T, b.T, a.T)
        # F, F, C
        check(a.T, b.T, c.T)
        # F, F, A
        check(a.T, b.T, d.T)
        # A, A, A
        # (note Numpy may select the layout differently for other inputs)
        check(d.T, e.T, d.T)

        # Inputs with compatible sizes
        check(a[1:], b, c[::-1])

        # Exceptions leak references
        self.disable_leak_check()

        # Incompatible sizes
        with self.assert_invalid_sizes_over_dim(1):
            cfunc(a, b[:,1:], b)

    def test_typing_errors(self):
        pyfunc = np_concatenate1
        cfunc = nrtjit(pyfunc)

        a = np.arange(15)
        b = a.reshape((3, 5))
        c = a.astype(np.dtype([('x', np.int8)]))
        d = np.array(42)

        # Different dimensionalities
        with self.assertTypingError() as raises:
            cfunc(a, b, b)
        self.assertIn("all the input arrays must have same number of dimensions",
                      str(raises.exception))

        # Incompatible dtypes
        with self.assertTypingError() as raises:
            cfunc(a, c, c)
        self.assertIn("input arrays must have compatible dtypes",
                      str(raises.exception))

        # 0-d arrays
        with self.assertTypingError() as raises:
            cfunc(d, d, d)
        self.assertIn("zero-dimensional arrays cannot be concatenated",
                      str(raises.exception))

        # non-tuple input
        with self.assertTypingError() as raises:
            cfunc(c, 1, c)
        self.assertIn('expecting a non-empty tuple of arrays', str(raises.exception))


@unittest.skipUnless(hasattr(np, "stack"), "this Numpy doesn't have np.stack()")
class TestNpStack(MemoryLeakMixin, TestCase):
    """
    Tests for np.stack().
    """

    def _3d_arrays(self):
        a = np.arange(24).reshape((4, 3, 2))
        b = a + 10
        c = (b + 10).copy(order='F')
        d = (c + 10)[::-1]
        e = (d + 10)[...,::-1]
        return a, b, c, d, e

    @contextlib.contextmanager
    def assert_invalid_sizes(self):
        with self.assertRaises(ValueError) as raises:
            yield
        self.assertIn("all input arrays must have the same shape",
                      str(raises.exception))

    def check_stack(self, pyfunc, cfunc, args):
        expected = pyfunc(*args)
        got = cfunc(*args)
        # Numba doesn't choose the same layout as Numpy.
        # We would like to check the result is contiguous, but we can't
        # rely on the "flags" attribute when there are 1-sized
        # dimensions.
        self.assertEqual(got.shape, expected.shape)
        self.assertPreciseEqual(got.flatten(), expected.flatten())

    def check_3d(self, pyfunc, cfunc, generate_starargs):
        def check(a, b, c, args):
            self.check_stack(pyfunc, cfunc, (a, b, c) + args)

        def check_all_axes(a, b, c):
            for args in generate_starargs():
                check(a, b, c, args)

        a, b, c, d, e = self._3d_arrays()

        # C, C, C
        check_all_axes(a, b, b)
        # C, C, F
        check_all_axes(a, b, c)
        # F, F, F
        check_all_axes(a.T, b.T, a.T)
        # F, F, C
        check_all_axes(a.T, b.T, c.T)
        # F, F, A
        check_all_axes(a.T, b.T, d.T)
        # A, A, A
        check_all_axes(d.T, e.T, d.T)

        # Different but compatible dtypes
        check_all_axes(a, b.astype(np.float64), b)

    def check_runtime_errors(self, cfunc, generate_starargs):
        # Exceptions leak references
        self.assert_no_memory_leak()
        self.disable_leak_check()

        # Inputs have different shapes
        a, b, c, d, e = self._3d_arrays()
        with self.assert_invalid_sizes():
            args = next(generate_starargs())
            cfunc(a[:-1], b, c, *args)

    def test_3d(self):
        """
        stack(3d arrays, axis)
        """
        pyfunc = np_stack2
        cfunc = nrtjit(pyfunc)

        def generate_starargs():
            for axis in range(3):
                yield (axis,)
                yield (-3 + axis,)

        self.check_3d(pyfunc, cfunc, generate_starargs)
        self.check_runtime_errors(cfunc, generate_starargs)

    def test_3d_no_axis(self):
        """
        stack(3d arrays)
        """
        pyfunc = np_stack1
        cfunc = nrtjit(pyfunc)

        def generate_starargs():
            yield()

        self.check_3d(pyfunc, cfunc, generate_starargs)
        self.check_runtime_errors(cfunc, generate_starargs)

    def test_0d(self):
        """
        stack(0d arrays)
        """
        pyfunc = np_stack1
        cfunc = nrtjit(pyfunc)

        a = np.array(42)
        b = np.array(-5j)
        c = np.array(True)

        self.check_stack(pyfunc, cfunc, (a, b, c))

    def check_xxstack(self, pyfunc, cfunc):
        """
        3d and 0d tests for hstack(), vstack(), dstack().
        """
        def generate_starargs():
            yield()

        self.check_3d(pyfunc, cfunc, generate_starargs)
        # 0d
        a = np.array(42)
        b = np.array(-5j)
        c = np.array(True)
        self.check_stack(pyfunc, cfunc, (a, b, a))

    def test_hstack(self):
        pyfunc = np_hstack
        cfunc = nrtjit(pyfunc)

        self.check_xxstack(pyfunc, cfunc)
        # 1d
        a = np.arange(5)
        b = np.arange(6) + 10
        self.check_stack(pyfunc, cfunc, (a, b, b))
        # 2d
        a = np.arange(6).reshape((2, 3))
        b = np.arange(8).reshape((2, 4)) + 100
        self.check_stack(pyfunc, cfunc, (a, b, a))

    def test_vstack(self):
        # Since np.row_stack is an alias for np.vstack, it does not need a
        # separate Numba implementation. For every test for np.vstack, the same
        # test for np.row_stack has been added.
        functions = [np_vstack, np_row_stack]
        for pyfunc in functions:
            cfunc = nrtjit(pyfunc)

            self.check_xxstack(pyfunc, cfunc)
            # 1d
            a = np.arange(5)
            b = a + 10
            self.check_stack(pyfunc, cfunc, (a, b, b))
            # 2d
            a = np.arange(6).reshape((3, 2))
            b = np.arange(8).reshape((4, 2)) + 100
            self.check_stack(pyfunc, cfunc, (a, b, b))

    def test_dstack(self):
        pyfunc = np_dstack
        cfunc = nrtjit(pyfunc)

        self.check_xxstack(pyfunc, cfunc)
        # 1d
        a = np.arange(5)
        b = a + 10
        self.check_stack(pyfunc, cfunc, (a, b, b))
        # 2d
        a = np.arange(12).reshape((3, 4))
        b = a + 100
        self.check_stack(pyfunc, cfunc, (a, b, b))

    def test_column_stack(self):
        pyfunc = np_column_stack
        cfunc = nrtjit(pyfunc)

        a = np.arange(4)
        b = a + 10
        c = np.arange(12).reshape((4, 3))
        self.check_stack(pyfunc, cfunc, (a, b, c))

        # Exceptions leak references
        self.assert_no_memory_leak()
        self.disable_leak_check()

        # Invalid dims
        a = np.array(42)
        with self.assertTypingError():
            cfunc((a, a, a))
        a = a.reshape((1, 1, 1))
        with self.assertTypingError():
            cfunc((a, a, a))

    def test_bad_arrays(self):
        for pyfunc in (np_stack1, np_hstack, np_vstack, np_dstack, np_column_stack):
            cfunc = nrtjit(pyfunc)
            c = np.arange(12).reshape((4, 3))

            # non-tuple input
            with self.assertTypingError() as raises:
                cfunc(c, 1, c)
            self.assertIn('expecting a non-empty tuple of arrays', str(raises.exception))


def benchmark_refct_speed():
    def pyfunc(x, y, t):
        """Swap array x and y for t number of times
        """
        for i in range(t):
            x, y = y, x
        return x, y

    cfunc = nrtjit(pyfunc)

    x = np.random.random(100)
    y = np.random.random(100)
    t = 10000

    def bench_pyfunc():
        pyfunc(x, y, t)

    def bench_cfunc():
        cfunc(x, y, t)

    python_time = utils.benchmark(bench_pyfunc)
    numba_time = utils.benchmark(bench_cfunc)
    print(python_time)
    print(numba_time)


if __name__ == "__main__":
    unittest.main()
