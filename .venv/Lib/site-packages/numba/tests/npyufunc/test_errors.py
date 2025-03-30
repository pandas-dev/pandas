import contextlib
import sys

import numpy as np

from numba import vectorize, guvectorize

from numba.tests.support import (TestCase, CheckWarningsMixin,
                                 skip_macos_fenv_errors)
import unittest


def sqrt(val):
    if val < 0.0:
        raise ValueError('Value must be positive')
    return val ** 0.5


def gufunc_foo(inp, n, out):
    for i in range(inp.shape[0]):
        if inp[i] < 0:
            raise ValueError('Value must be positive')
        out[i] = inp[i] * n[0]

def truediv(a, b):
    return a / b

def floordiv(a, b):
    return a // b

def remainder(a, b):
    return a % b

def power(a, b):
    return a ** b


class TestExceptions(TestCase):
    """
    Test raising exceptions inside ufuncs.
    """

    def check_ufunc_raise(self, **vectorize_args):
        f = vectorize(['float64(float64)'], **vectorize_args)(sqrt)
        arr = np.array([1, 4, -2, 9, -1, 16], dtype=np.float64)
        out = np.zeros_like(arr)
        with self.assertRaises(ValueError) as cm:
            f(arr, out)
        self.assertIn('Value must be positive', str(cm.exception))
        # All values were computed except for the ones giving an error
        self.assertEqual(list(out), [1, 2, 0, 3, 0, 4])

    def test_ufunc_raise(self):
        self.check_ufunc_raise(nopython=True)

    def test_ufunc_raise_objmode(self):
        self.check_ufunc_raise(forceobj=True)

    def check_gufunc_raise(self, **vectorize_args):
        f = guvectorize(['int32[:], int32[:], int32[:]'], '(n),()->(n)',
                        **vectorize_args)(gufunc_foo)
        arr = np.array([1, 2, -3, 4], dtype=np.int32)
        out = np.zeros_like(arr)
        with self.assertRaises(ValueError) as cm:
            f(arr, 2, out)
        # The gufunc bailed out after the error
        self.assertEqual(list(out), [2, 4, 0, 0])

    def test_gufunc_raise(self):
        self.check_gufunc_raise(nopython=True)

    def test_gufunc_raise_objmode(self):
        self.check_gufunc_raise(forceobj=True)

class TestFloatingPointExceptions(TestCase, CheckWarningsMixin):
    """
    Test floating-point exceptions inside ufuncs.

    Note the warnings emitted by Numpy reflect IEEE-754 semantics.
    """

    def check_truediv_real(self, dtype):
        """
        Test 1 / 0 and 0 / 0.
        """
        f = vectorize(nopython=True)(truediv)
        a = np.array([5., 6., 0., 8.], dtype=dtype)
        b = np.array([1., 0., 0., 4.], dtype=dtype)
        expected = np.array([5., float('inf'), float('nan'), 2.])
        with self.check_warnings(["divide by zero encountered",
                                  "invalid value encountered"]):
            res = f(a, b)
            self.assertPreciseEqual(res, expected)

    def test_truediv_float(self):
        self.check_truediv_real(np.float64)

    def test_truediv_integer(self):
        self.check_truediv_real(np.int32)

    def check_divmod_float(self, pyfunc, values, messages):
        """
        Test 1 // 0 and 0 // 0.
        """
        f = vectorize(nopython=True)(pyfunc)
        a = np.array([5., 6., 0., 9.])
        b = np.array([1., 0., 0., 4.])
        expected = np.array(values)
        with self.check_warnings(messages):
            res = f(a, b)
            self.assertPreciseEqual(res, expected)

    def test_floordiv_float(self):
        self.check_divmod_float(floordiv,
                                [5.0, float('inf'), float('nan'), 2.0],
                                ["divide by zero encountered",
                                 "invalid value encountered"])

    @skip_macos_fenv_errors
    def test_remainder_float(self):
        self.check_divmod_float(remainder,
                                [0.0, float('nan'), float('nan'), 1.0],
                                ["invalid value encountered"])

    def check_divmod_int(self, pyfunc, values):
        """
        Test 1 % 0 and 0 % 0.
        """
        f = vectorize(nopython=True)(pyfunc)
        a = np.array([5, 6, 0, 9])
        b = np.array([1, 0, 0, 4])
        expected = np.array(values)
        # No warnings raised because LLVM makes it difficult
        with self.check_warnings([]):
            res = f(a, b)
            self.assertPreciseEqual(res, expected)

    def test_floordiv_int(self):
        self.check_divmod_int(floordiv, [5, 0, 0, 2])

    def test_remainder_int(self):
        self.check_divmod_int(remainder, [0, 0, 0, 1])

    def test_power_float(self):
        """
        Test 0 ** -1 and 2 ** <big number>.
        """
        f = vectorize(nopython=True)(power)
        a = np.array([5., 0., 2., 8.])
        b = np.array([1., -1., 1e20, 4.])
        expected = np.array([5., float('inf'), float('inf'), 4096.])
        with self.check_warnings(["divide by zero encountered",
                                  "overflow encountered"]):
            res = f(a, b)
            self.assertPreciseEqual(res, expected)

    def test_power_integer(self):
        """
        Test 0 ** -1.
        Note 2 ** <big number> returns an undefined value (depending
        on the algorithm).
        """
        dtype = np.int64
        f = vectorize(["int64(int64, int64)"], nopython=True)(power)
        a = np.array([5, 0, 6], dtype=dtype)
        b = np.array([1, -1, 2], dtype=dtype)
        expected = np.array([5, -2**63, 36], dtype=dtype)
        with self.check_warnings([]):
            res = f(a, b)
            self.assertPreciseEqual(res, expected)


if __name__ == "__main__":
    unittest.main()
