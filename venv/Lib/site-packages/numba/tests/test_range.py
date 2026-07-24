import unittest

import numpy as np

from numba import jit, njit
from numba.core import types
from numba.tests.support import TestCase

from numba.core.inline_closurecall import length_of_iterator

def loop1(n):
    s = 0
    for i in range(n):
        s += i
    return s


def loop2(a, b):
    s = 0
    for i in range(a, b):
        s += i
    return s


def loop3(a, b, c):
    s = 0
    for i in range(a, b, c):
        s += i
    return s


def range_len1(n):
    return len(range(n))

def range_len2(a, b):
    return len(range(a, b))

def range_len3(a, b, c):
    return len(range(a, b, c))
def range_iter_len1(a):
    return length_of_iterator(iter(range(a)))

def range_iter_len2(a):
    return length_of_iterator(iter(a))

def range_attrs(start, stop, step):
    r1 = range(start)
    r2 = range(start, stop)
    r3 = range(start, stop, step)
    tmp = []
    for r in (r1, r2, r3):
        tmp.append((r.start, r.stop, r.step))
    return tmp

def range_contains(val, start, stop, step):
    r1 = range(start)
    r2 = range(start, stop)
    r3 = range(start, stop, step)
    return [val in r for r in (r1, r2, r3)]


class TestRange(unittest.TestCase):

    def test_loop1_int16(self):
        pyfunc = loop1
        cfunc = njit((types.int16,))(pyfunc)
        self.assertTrue(cfunc(5), pyfunc(5))

    def test_loop2_int16(self):
        pyfunc = loop2
        cfunc = njit((types.int16, types.int16))(pyfunc)
        self.assertTrue(cfunc(1, 6), pyfunc(1, 6))

    def test_loop3_int32(self):
        pyfunc = loop3
        cfunc = njit((types.int32, types.int32, types.int32))(pyfunc)
        arglist = [
            (1, 2, 1),
            (2, 8, 3),
            (-10, -11, -10),
            (-10, -10, -2),
        ]
        for args in arglist:
            self.assertEqual(cfunc(*args), pyfunc(*args))

    def test_range_len1(self):
        pyfunc = range_len1
        typelist = [types.int16, types.int32, types.int64]
        arglist = [5, 0, -5]
        for typ in typelist:
            cfunc = njit((typ,))(pyfunc)
            for arg in arglist:
                self.assertEqual(cfunc(typ(arg)), pyfunc(typ(arg)))

    def test_range_len2(self):
        pyfunc = range_len2
        typelist = [types.int16, types.int32, types.int64]
        arglist = [(1,6), (6,1), (-5, -1)]
        for typ in typelist:
            cfunc = njit((typ, typ))(pyfunc)
            for args in arglist:
                args_ = tuple(typ(x) for x in args)
                self.assertEqual(cfunc(*args_), pyfunc(*args_))

    def test_range_len3(self):
        pyfunc = range_len3
        typelist = [types.int16, types.int32, types.int64]
        arglist = [
            (1, 2, 1),
            (2, 8, 3),
            (-10, -11, -10),
            (-10, -10, -2),
        ]
        for typ in typelist:
            cfunc = njit((typ, typ, typ))(pyfunc)
            for args in arglist:
                args_ = tuple(typ(x) for x in args)
                self.assertEqual(cfunc(*args_), pyfunc(*args_))

    def test_range_iter_len1(self):
        range_func = range_len1
        range_iter_func = range_iter_len1
        typelist = [types.int16, types.int32, types.int64]
        arglist = [5, 0, -5]
        for typ in typelist:
            cfunc = njit((typ,))(range_iter_func)
            for arg in arglist:
                self.assertEqual(cfunc(typ(arg)), range_func(typ(arg)))

    def test_range_iter_list(self):
        range_iter_func = range_iter_len2
        cfunc = njit((types.List(types.intp, reflected=True),))(range_iter_func)
        arglist = [1, 2, 3, 4, 5]
        self.assertEqual(cfunc(arglist), len(arglist))

    def test_range_attrs(self):
        pyfunc = range_attrs
        arglist = [(0, 0, 1),
                   (0, -1, 1),
                   (-1, 1, 1),
                   (-1, 4, 1),
                   (-1, 4, 10),
                   (5, -5, -2),]

        cfunc = njit((types.int64, types.int64, types.int64),)(pyfunc)
        for arg in arglist:
            self.assertEqual(cfunc(*arg), pyfunc(*arg))

    def test_range_contains(self):
        pyfunc = range_contains
        arglist = [(0, 0, 1),
                   (-1, 0, 1),
                   (1, 0, -1),
                   (0, -1, 1),
                   (0, 1, -1),
                   (-1, 1, 1),
                   (-1, 4, 1),
                   (-1, 4, 10),
                   (5, -5, -2),]

        bool_vals = [True, False]
        int_vals = [-10, -6, -5, -4, -2, -1, 0,
                     1, 2, 4, 5, 6, 10]
        float_vals = [-1.1, -1.0, 0.0, 1.0, 1.1]
        complex_vals = [1 + 0j, 1 + 1j, 1.1 + 0j, 1.0 + 1.1j]

        vallist = (bool_vals + int_vals + float_vals
                + complex_vals)

        cfunc = njit(pyfunc)
        for arg in arglist:
            for val in vallist:
                self.assertEqual(cfunc(val, *arg), pyfunc(val, *arg))

        non_numeric_vals = [{'a': 1}, [1, ], 'abc', (1,)]

        cfunc_obj = jit(pyfunc, forceobj=True)
        for arg in arglist:
            for val in non_numeric_vals:
                self.assertEqual(cfunc_obj(val, *arg), pyfunc(val, *arg))


@njit
def my_arange(start, stop, step):
    x = np.zeros(len(range(start, stop, step)), dtype=np.uint64)
    i = 0
    for v in range(start, stop, step):
        x[i] = v
        i += 1
    return x


class TestRangeNumpy(TestCase):
    def test_range_safe_cast_mixed(self):
        """Test that mixing `uint64` and `int64` works."""
        a = my_arange(np.uint64(6), np.uint64(0), np.int64(-1))
        self.assertPreciseEqual(a, np.arange(6, 0, -1, dtype=np.uint64))


if __name__ == '__main__':
    unittest.main()
