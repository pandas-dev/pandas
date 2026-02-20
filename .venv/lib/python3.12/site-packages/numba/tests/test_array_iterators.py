import itertools

import numpy as np

from numba import jit, njit, typeof
from numba.core import types
from numba.tests.support import TestCase, MemoryLeakMixin
import unittest


def array_iter(arr):
    total = 0
    for i, v in enumerate(arr):
        total += i * v
    return total

def array_iter_items(arr):
    return list(iter(arr))

def array_view_iter(arr, idx):
    total = 0
    for i, v in enumerate(arr[idx]):
        total += i * v
    return total

def array_flat(arr, out):
    for i, v in enumerate(arr.flat):
        out[i] = v

def array_flat_getitem(arr, ind):
    return arr.flat[ind]

def array_flat_setitem(arr, ind, val):
    arr.flat[ind] = val

def array_flat_sum(arr):
    s = 0
    for i, v in enumerate(arr.flat):
        s = s + (i + 1) * v
    return s

def array_flat_len(arr):
    return len(arr.flat)

def array_ndenumerate_zero_dim(arr):
    for i, v in np.ndenumerate(arr):
        return (i, v)
    return None

def array_ndenumerate_sum(arr):
    s = 0
    for (i, j), v in np.ndenumerate(arr):
        s = s + (i + 1) * (j + 1) * v
    return s

def np_ndindex_empty():
    s = 0
    for ind in np.ndindex(()):
        s += s + len(ind) + 1
    return s

def np_ndindex(x, y):
    s = 0
    n = 0
    for i, j in np.ndindex(x, y):
        s = s + (i + 1) * (j + 1)
    return s

def np_ndindex_array(arr):
    s = 0
    n = 0
    for indices in np.ndindex(arr.shape):
        for i, j in enumerate(indices):
            s = s + (i + 1) * (j + 1)
    return s

def np_nditer1(a):
    res = []
    for u in np.nditer(a):
        res.append(u.item())
    return res

def np_nditer2(a, b):
    res = []
    for u, v in np.nditer((a, b)):
        res.append((u.item(), v.item()))
    return res

def np_nditer3(a, b, c):
    res = []
    for u, v, w in np.nditer((a, b, c)):
        res.append((u.item(), v.item(), w.item()))
    return res

def iter_next(arr):
    it = iter(arr)
    it2 = iter(arr)
    return next(it), next(it), next(it2)


#
# Test premature free (see issue #2112).
# The following test allocates an array ``x`` inside the body.
# The compiler will put a ``del x`` right after the last use of ``x``,
# which is right after the creation of the array iterator and
# before the loop is entered.  If the iterator does not incref the array,
# the iterator will be reading garbage data of free'ed memory.
#

def array_flat_premature_free(size):
    x = np.arange(size)
    res = np.zeros_like(x, dtype=np.intp)
    for i, v in enumerate(x.flat):
        res[i] = v
    return res

def array_ndenumerate_premature_free(size):
    x = np.arange(size)
    res = np.zeros_like(x, dtype=np.intp)
    for i, v in np.ndenumerate(x):
        res[i] = v
    return res


class TestArrayIterators(MemoryLeakMixin, TestCase):
    """
    Test array.flat, np.ndenumerate(), etc.
    """

    def setUp(self):
        super(TestArrayIterators, self).setUp()

    def check_array_iter_1d(self, arr):
        pyfunc = array_iter
        cfunc = njit((typeof(arr),))(pyfunc)
        expected = pyfunc(arr)
        self.assertPreciseEqual(cfunc(arr), expected)

    def check_array_iter_items(self, arr):
        pyfunc = array_iter_items
        cfunc = njit((typeof(arr),))(pyfunc)
        expected = pyfunc(arr)
        self.assertPreciseEqual(cfunc(arr), expected)

    def check_array_view_iter(self, arr, index):
        pyfunc = array_view_iter
        cfunc = njit((typeof(arr), typeof(index),))(pyfunc)
        expected = pyfunc(arr, index)
        self.assertPreciseEqual(cfunc(arr, index), expected)

    def check_array_flat(self, arr, arrty=None):
        out = np.zeros(arr.size, dtype=arr.dtype)
        nb_out = out.copy()
        if arrty is None:
            arrty = typeof(arr)

        cfunc = njit((arrty, typeof(out),))(array_flat)

        array_flat(arr, out)
        cfunc(arr, nb_out)

        self.assertPreciseEqual(out, nb_out)

    def check_array_unary(self, arr, arrty, func):
        cfunc = njit((arrty,))(func)
        self.assertPreciseEqual(cfunc(arr), func(arr))

    def check_array_ndenumerate_sum(self, arr, arrty):
        self.check_array_unary(arr, arrty, array_ndenumerate_sum)

    def test_array_iter(self):
        # Test iterating over arrays
        arr = np.arange(6)
        self.check_array_iter_1d(arr)
        self.check_array_iter_items(arr)
        arr = arr[::2]
        self.assertFalse(arr.flags.c_contiguous)
        self.assertFalse(arr.flags.f_contiguous)
        self.check_array_iter_1d(arr)
        self.check_array_iter_items(arr)
        arr = np.bool_([1, 0, 0, 1])
        self.check_array_iter_1d(arr)
        self.check_array_iter_items(arr)
        arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        self.check_array_iter_items(arr)
        self.check_array_iter_items(arr.T)

    def test_array_iter_yielded_order(self):
        # See issue #5692
        @jit(nopython=True)
        def foo(arr):
            t = []
            for y1 in arr:
                for y2 in y1:
                    t.append(y2.ravel())
            return t

        # 'F' ordered
        arr = np.arange(24).reshape((2, 3, 4), order='F')
        expected = foo.py_func(arr)
        got = foo(arr)
        self.assertPreciseEqual(expected, got)

        # 'A' ordered, outer strided
        arr = np.arange(64).reshape((4, 8, 2), order='F')[::2, :, :]
        expected = foo.py_func(arr)
        got = foo(arr)
        self.assertPreciseEqual(expected, got)

        # 'A' ordered, middle strided
        arr = np.arange(64).reshape((4, 8, 2), order='F')[:, ::2, :]
        expected = foo.py_func(arr)
        got = foo(arr)
        self.assertPreciseEqual(expected, got)

        # 'A' ordered, inner strided
        arr = np.arange(64).reshape((4, 8, 2), order='F')[:, :, ::2]
        expected = foo.py_func(arr)
        got = foo(arr)
        self.assertPreciseEqual(expected, got)

        @jit(nopython=True)
        def flag_check(arr):
            out = []
            for sub in arr:
                out.append((sub, sub.flags.c_contiguous,
                            sub.flags.f_contiguous))
            return out

        arr = np.arange(10).reshape((2, 5), order='F')
        expected = flag_check.py_func(arr)
        got = flag_check(arr)

        self.assertEqual(len(expected), len(got))
        ex_arr, e_flag_c, e_flag_f = expected[0]
        go_arr, g_flag_c, g_flag_f = got[0]
        np.testing.assert_allclose(ex_arr, go_arr)
        self.assertEqual(e_flag_c, g_flag_c)
        self.assertEqual(e_flag_f, g_flag_f)

    def test_array_view_iter(self):
        # Test iterating over a 1d view over a 2d array
        arr = np.arange(12).reshape((3, 4))
        self.check_array_view_iter(arr, 1)
        self.check_array_view_iter(arr.T, 1)
        arr = arr[::2]
        self.check_array_view_iter(arr, 1)
        arr = np.bool_([1, 0, 0, 1]).reshape((2, 2))
        self.check_array_view_iter(arr, 1)

    def test_array_flat_3d(self):
        arr = np.arange(24).reshape(4, 2, 3)

        arrty = typeof(arr)
        self.assertEqual(arrty.ndim, 3)
        self.assertEqual(arrty.layout, 'C')
        self.assertTrue(arr.flags.c_contiguous)
        # Test with C-contiguous array
        self.check_array_flat(arr)
        # Test with Fortran-contiguous array
        arr = arr.transpose()
        self.assertFalse(arr.flags.c_contiguous)
        self.assertTrue(arr.flags.f_contiguous)
        self.assertEqual(typeof(arr).layout, 'F')
        self.check_array_flat(arr)
        # Test with non-contiguous array
        arr = arr[::2]
        self.assertFalse(arr.flags.c_contiguous)
        self.assertFalse(arr.flags.f_contiguous)
        self.assertEqual(typeof(arr).layout, 'A')
        self.check_array_flat(arr)
        # Boolean array
        arr = np.bool_([1, 0, 0, 1] * 2).reshape((2, 2, 2))
        self.check_array_flat(arr)

    def test_array_flat_empty(self):
        # Test .flat with various shapes of empty arrays, contiguous
        # and non-contiguous (see issue #846).

        # Define a local checking function, Numba's `typeof` ends up aliasing
        # 0d C and F ordered arrays, so the check needs to go via the compile
        # result entry point to bypass type checking.
        def check(arr, arrty):
            cfunc = njit((arrty,))(array_flat_sum)
            cres = cfunc.overloads[(arrty,)]
            got = cres.entry_point(arr)
            expected = cfunc.py_func(arr)
            self.assertPreciseEqual(expected, got)

        arr = np.zeros(0, dtype=np.int32)
        arr = arr.reshape(0, 2)
        arrty = types.Array(types.int32, 2, layout='C')
        check(arr, arrty)
        arrty = types.Array(types.int32, 2, layout='F')
        check(arr, arrty)
        arrty = types.Array(types.int32, 2, layout='A')
        check(arr, arrty)
        arr = arr.reshape(2, 0)
        arrty = types.Array(types.int32, 2, layout='C')
        check(arr, arrty)
        arrty = types.Array(types.int32, 2, layout='F')
        check(arr, arrty)
        arrty = types.Array(types.int32, 2, layout='A')
        check(arr, arrty)

    def test_array_flat_getitem(self):
        # Test indexing of array.flat object
        pyfunc = array_flat_getitem
        cfunc = njit(pyfunc)
        def check(arr, ind):
            expected = pyfunc(arr, ind)
            self.assertEqual(cfunc(arr, ind), expected)

        arr = np.arange(24).reshape(4, 2, 3)
        for i in range(arr.size):
            check(arr, i)
        arr = arr.T
        for i in range(arr.size):
            check(arr, i)
        arr = arr[::2]
        for i in range(arr.size):
            check(arr, i)
        arr = np.array([42]).reshape(())
        for i in range(arr.size):
            check(arr, i)
        # Boolean array
        arr = np.bool_([1, 0, 0, 1])
        for i in range(arr.size):
            check(arr, i)
        arr = arr[::2]
        for i in range(arr.size):
            check(arr, i)

    def test_array_flat_setitem(self):
        # Test indexing of array.flat object
        pyfunc = array_flat_setitem
        cfunc = njit(pyfunc)
        def check(arr, ind):
            # Use np.copy() to keep the layout
            expected = np.copy(arr)
            got = np.copy(arr)
            pyfunc(expected, ind, 123)
            cfunc(got, ind, 123)
            self.assertPreciseEqual(got, expected)

        arr = np.arange(24).reshape(4, 2, 3)
        for i in range(arr.size):
            check(arr, i)
        arr = arr.T
        for i in range(arr.size):
            check(arr, i)
        arr = arr[::2]
        for i in range(arr.size):
            check(arr, i)
        arr = np.array([42]).reshape(())
        for i in range(arr.size):
            check(arr, i)
        # Boolean array
        arr = np.bool_([1, 0, 0, 1])
        for i in range(arr.size):
            check(arr, i)
        arr = arr[::2]
        for i in range(arr.size):
            check(arr, i)

    def test_array_flat_len(self):
        # Test len(array.flat)
        pyfunc = array_flat_len
        cfunc = njit(array_flat_len)
        def check(arr):
            expected = pyfunc(arr)
            self.assertPreciseEqual(cfunc(arr), expected)

        arr = np.arange(24).reshape(4, 2, 3)
        check(arr)
        arr = arr.T
        check(arr)
        arr = arr[::2]
        check(arr)
        arr = np.array([42]).reshape(())
        check(arr)

    def test_array_flat_premature_free(self):
        cfunc = njit((types.intp,))(array_flat_premature_free)
        expect = array_flat_premature_free(6)
        got = cfunc(6)
        self.assertTrue(got.sum())
        self.assertPreciseEqual(expect, got)

    def test_array_ndenumerate_2d(self):
        arr = np.arange(12).reshape(4, 3)
        arrty = typeof(arr)
        self.assertEqual(arrty.ndim, 2)
        self.assertEqual(arrty.layout, 'C')
        self.assertTrue(arr.flags.c_contiguous)
        # Test with C-contiguous array
        self.check_array_ndenumerate_sum(arr, arrty)
        # Test with Fortran-contiguous array
        arr = arr.transpose()
        self.assertFalse(arr.flags.c_contiguous)
        self.assertTrue(arr.flags.f_contiguous)
        arrty = typeof(arr)
        self.assertEqual(arrty.layout, 'F')
        self.check_array_ndenumerate_sum(arr, arrty)
        # Test with non-contiguous array
        arr = arr[::2]
        self.assertFalse(arr.flags.c_contiguous)
        self.assertFalse(arr.flags.f_contiguous)
        arrty = typeof(arr)
        self.assertEqual(arrty.layout, 'A')
        self.check_array_ndenumerate_sum(arr, arrty)
        # Boolean array
        arr = np.bool_([1, 0, 0, 1]).reshape((2, 2))
        self.check_array_ndenumerate_sum(arr, typeof(arr))

    def test_array_ndenumerate_empty(self):
        # Define a local checking function, Numba's `typeof` ends up aliasing
        # 0d C and F ordered arrays, so the check needs to go via the compile
        # result entry point to bypass type checking.
        def check(arr, arrty):
            cfunc = njit((arrty,))(array_ndenumerate_sum)
            cres = cfunc.overloads[(arrty,)]
            got = cres.entry_point(arr)
            expected = cfunc.py_func(arr)
            np.testing.assert_allclose(expected, got)

        arr = np.zeros(0, dtype=np.int32)
        arr = arr.reshape(0, 2)
        arrty = types.Array(types.int32, 2, layout='C')
        check(arr, arrty)
        arrty = types.Array(types.int32, 2, layout='F')
        check(arr, arrty)
        arrty = types.Array(types.int32, 2, layout='A')
        check(arr, arrty)
        arr = arr.reshape(2, 0)
        arrty = types.Array(types.int32, 2, layout='C')
        check(arr, arrty)
        arrty = types.Array(types.int32, 2, layout='F')
        check(arr, arrty)
        arrty = types.Array(types.int32, 2, layout='A')
        check(arr, arrty)

    def test_array_ndenumerate_premature_free(self):
        cfunc = njit((types.intp,))(array_ndenumerate_premature_free)
        expect = array_ndenumerate_premature_free(6)
        got = cfunc(6)
        self.assertTrue(got.sum())
        self.assertPreciseEqual(expect, got)

    def test_array_ndenumerate_zero_dim(self):
        func = array_ndenumerate_zero_dim
        cfunc = njit(func)
        arr = np.array(97)
        expected = func(arr)
        got = cfunc(arr)
        self.assertPreciseEqual(got, expected)
        self.assertEqual(expected[0], ())
        self.assertEqual(expected[1], 97)

    def test_np_ndindex(self):
        func = np_ndindex
        cfunc = njit((types.int32, types.int32,))(func)
        self.assertPreciseEqual(cfunc(3, 4), func(3, 4))
        self.assertPreciseEqual(cfunc(3, 0), func(3, 0))
        self.assertPreciseEqual(cfunc(0, 3), func(0, 3))
        self.assertPreciseEqual(cfunc(0, 0), func(0, 0))

    def test_np_ndindex_array(self):
        func = np_ndindex_array
        arr = np.arange(12, dtype=np.int32) + 10
        self.check_array_unary(arr, typeof(arr), func)
        arr = arr.reshape((4, 3))
        self.check_array_unary(arr, typeof(arr), func)
        arr = arr.reshape((2, 2, 3))
        self.check_array_unary(arr, typeof(arr), func)

    def test_np_ndindex_empty(self):
        func = np_ndindex_empty
        cfunc = njit((),)(func)
        self.assertPreciseEqual(cfunc(), func())

    def test_iter_next(self):
        # This also checks memory management with iter() and next()
        func = iter_next
        arr = np.arange(12, dtype=np.int32) + 10
        self.check_array_unary(arr, typeof(arr), func)


class TestNdIter(MemoryLeakMixin, TestCase):
    """
    Test np.nditer()
    """

    def inputs(self):
        # All those inputs are compatible with a (3, 4) main shape

        # scalars
        yield np.float32(100)

        # 0-d arrays
        yield np.array(102, dtype=np.int16)

        # 1-d arrays
        yield np.arange(4).astype(np.complex64)
        yield np.arange(8)[::2]

        # 2-d arrays
        a = np.arange(12).reshape((3, 4))
        yield a
        yield a.copy(order='F')
        a = np.arange(24).reshape((6, 4))[::2]
        yield a

    def basic_inputs(self):
        yield np.arange(4).astype(np.complex64)
        yield np.arange(8)[::2]
        a = np.arange(12).reshape((3, 4))
        yield a
        yield a.copy(order='F')

    def check_result(self, got, expected):
        self.assertEqual(set(got), set(expected), (got, expected))

    def test_nditer1(self):
        pyfunc = np_nditer1
        cfunc = jit(nopython=True)(pyfunc)
        for a in self.inputs():
            expected = pyfunc(a)
            got = cfunc(a)
            self.check_result(got, expected)

    def test_nditer2(self):
        pyfunc = np_nditer2
        cfunc = jit(nopython=True)(pyfunc)
        for a, b in itertools.product(self.inputs(), self.inputs()):
            expected = pyfunc(a, b)
            got = cfunc(a, b)
            self.check_result(got, expected)

    def test_nditer3(self):
        pyfunc = np_nditer3
        cfunc = jit(nopython=True)(pyfunc)
        # Use a restricted set of inputs, to shorten test time
        inputs = self.basic_inputs
        for a, b, c in itertools.product(inputs(), inputs(), inputs()):
            expected = pyfunc(a, b, c)
            got = cfunc(a, b, c)
            self.check_result(got, expected)

    def test_errors(self):
        # Incompatible shapes
        pyfunc = np_nditer2
        cfunc = jit(nopython=True)(pyfunc)

        self.disable_leak_check()

        def check_incompatible(a, b):
            with self.assertRaises(ValueError) as raises:
                cfunc(a, b)
            self.assertIn("operands could not be broadcast together",
                          str(raises.exception))

        check_incompatible(np.arange(2), np.arange(3))
        a = np.arange(12).reshape((3, 4))
        b = np.arange(3)
        check_incompatible(a, b)


if __name__ == '__main__':
    unittest.main()
