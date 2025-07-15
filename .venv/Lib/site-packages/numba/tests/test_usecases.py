import itertools
import unittest
import numpy as np

from numba import jit, njit
from numba.core import types
from numba.tests import usecases
from numba.tests.support import TestCase


class TestUsecases(TestCase):
    # NOTE: All these test cases are run in subprocesses to achieve total
    # isolation.

    @TestCase.run_test_in_subprocess
    def test_andor(self):
        pyfunc = usecases.andor
        cfunc = njit((types.int32, types.int32))(pyfunc)

        # Argument boundaries
        xs = -1, 0, 1, 9, 10, 11
        ys = -1, 0, 1, 9, 10, 11

        for args in itertools.product(xs, ys):
            self.assertEqual(pyfunc(*args), cfunc(*args), "args %s" % (args,))

    @TestCase.run_test_in_subprocess
    def test_sum1d(self):
        pyfunc = usecases.sum1d
        cfunc = njit((types.int32, types.int32))(pyfunc)

        ss = -1, 0, 1, 100, 200
        es = -1, 0, 1, 100, 200

        for args in itertools.product(ss, es):
            self.assertEqual(pyfunc(*args), cfunc(*args), args)

    @TestCase.run_test_in_subprocess
    def test_sum1d_pyobj(self):
        pyfunc = usecases.sum1d
        cfunc = jit((types.int32, types.int32), forceobj=True)(pyfunc)

        ss = -1, 0, 1, 100, 200
        es = -1, 0, 1, 100, 200

        for args in itertools.product(ss, es):
            self.assertEqual(pyfunc(*args), cfunc(*args), args)

    @TestCase.run_test_in_subprocess
    def test_sum2d(self):
        pyfunc = usecases.sum2d
        cfunc = njit((types.int32, types.int32))(pyfunc)

        ss = -1, 0, 1, 100, 200
        es = -1, 0, 1, 100, 200

        for args in itertools.product(ss, es):
            self.assertEqual(pyfunc(*args), cfunc(*args), args)

    @TestCase.run_test_in_subprocess
    def test_while_count(self):
        pyfunc = usecases.while_count
        cfunc = njit((types.int32, types.int32))(pyfunc)

        ss = -1, 0, 1, 100, 200
        es = -1, 0, 1, 100, 200

        for args in itertools.product(ss, es):
            self.assertEqual(pyfunc(*args), cfunc(*args), args)

    @TestCase.run_test_in_subprocess
    def test_copy_arrays(self):
        pyfunc = usecases.copy_arrays
        arraytype = types.Array(types.int32, 1, 'A')
        cfunc = njit((arraytype, arraytype))(pyfunc)

        nda = 0, 1, 10, 100

        for nd in nda:
            a = np.arange(nd, dtype='int32')
            b = np.empty_like(a)
            args = a, b

            cfunc(*args)
            self.assertPreciseEqual(a, b, msg=str(args))

    @TestCase.run_test_in_subprocess
    def test_copy_arrays2d(self):
        pyfunc = usecases.copy_arrays2d
        arraytype = types.Array(types.int32, 2, 'A')
        cfunc = njit((arraytype, arraytype))(pyfunc)

        nda = (0, 0), (1, 1), (2, 5), (4, 25)

        for nd in nda:
            d1, d2 = nd
            a = np.arange(d1 * d2, dtype='int32').reshape(d1, d2)
            b = np.empty_like(a)
            args = a, b

            cfunc(*args)
            self.assertPreciseEqual(a, b, msg=str(args))

    @TestCase.run_test_in_subprocess
    def test_string_concat(self):
        pyfunc = usecases.string_concat
        cfunc = jit((types.int32, types.int32), forceobj=True)(pyfunc)

        xs = -1, 0, 1
        ys = -1, 0, 1

        for x, y in itertools.product(xs, ys):
            args = x, y
            self.assertEqual(pyfunc(*args), cfunc(*args), args)

    @TestCase.run_test_in_subprocess
    def test_string_len(self):
        pyfunc = usecases.string_len
        cfunc = jit((types.pyobject,), forceobj=True)(pyfunc)

        test_str = '123456'
        self.assertEqual(pyfunc(test_str), cfunc(test_str))
        test_str = '1'
        self.assertEqual(pyfunc(test_str), cfunc(test_str))
        test_str = ''
        self.assertEqual(pyfunc(test_str), cfunc(test_str))

    @TestCase.run_test_in_subprocess
    def test_string_slicing(self):
        pyfunc = usecases.string_slicing
        cfunc = jit((types.pyobject,) * 3, forceobj=True)(pyfunc)

        test_str = '123456'
        self.assertEqual(pyfunc(test_str, 0, 3), cfunc(test_str, 0, 3))
        self.assertEqual(pyfunc(test_str, 1, 5), cfunc(test_str, 1, 5))
        self.assertEqual(pyfunc(test_str, 2, 3), cfunc(test_str, 2, 3))

    @TestCase.run_test_in_subprocess
    def test_string_conversion(self):
        pyfunc = usecases.string_conversion

        cfunc = jit((types.int32,), forceobj=True)(pyfunc)
        self.assertEqual(pyfunc(1), cfunc(1))

        cfunc = jit((types.float32,), forceobj=True)(pyfunc)
        self.assertEqual(pyfunc(1.1), cfunc(1.1))

    @TestCase.run_test_in_subprocess
    def test_string_comparisons(self):
        import operator
        pyfunc = usecases.string_comparison
        cfunc = jit((types.pyobject, types.pyobject, types.pyobject),
                    forceobj=True)(pyfunc)

        test_str1 = '123'
        test_str2 = '123'
        op = operator.eq
        self.assertEqual(pyfunc(test_str1, test_str2, op),
            cfunc(test_str1, test_str2, op))

        test_str1 = '123'
        test_str2 = '456'
        op = operator.eq
        self.assertEqual(pyfunc(test_str1, test_str2, op),
            cfunc(test_str1, test_str2, op))

        test_str1 = '123'
        test_str2 = '123'
        op = operator.ne
        self.assertEqual(pyfunc(test_str1, test_str2, op),
            cfunc(test_str1, test_str2, op))

        test_str1 = '123'
        test_str2 = '456'
        op = operator.ne
        self.assertEqual(pyfunc(test_str1, test_str2, op),
            cfunc(test_str1, test_str2, op))

    @TestCase.run_test_in_subprocess
    def test_blackscholes_cnd(self):
        pyfunc = usecases.blackscholes_cnd
        cfunc = njit((types.float32,))(pyfunc)

        ds = -0.5, 0, 0.5

        for d in ds:
            args = (d,)
            self.assertEqual(pyfunc(*args), cfunc(*args), args)


if __name__ == '__main__':
    unittest.main()
