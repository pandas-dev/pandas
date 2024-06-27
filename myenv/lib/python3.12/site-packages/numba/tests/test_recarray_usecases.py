import sys
import unittest

import numpy as np

from numba import njit
from numba.core import types
from numba.tests.support import captured_stdout, TestCase
from numba.np import numpy_support


def usecase1(arr1, arr2):
    """Base on https://github.com/numba/numba/issues/370

    Modified to add test-able side effect.
    """
    n1 = arr1.size
    n2 = arr2.size

    for i1 in range(n1):
        st1 = arr1[i1]
        for i2 in range(n2):
            st2 = arr2[i2]
            st2.row += st1.p * st2.p + st1.row - st1.col

        st1.p += st2.p
        st1.col -= st2.col


def usecase2(x, N):
    """
    Base on test1 of https://github.com/numba/numba/issues/381
    """
    for k in range(N):
        y = x[k]
        print(y.f1, y.s1, y.f2)


def usecase3(x, N):
    """
    Base on test2 of https://github.com/numba/numba/issues/381
    """
    for k in range(N):
        print(x.f1[k], x.s1[k], x.f2[k])


def usecase4(x, N):
    """
    Base on test3 of https://github.com/numba/numba/issues/381
    """
    for k in range(N):
        y = x[k]
        print(y.f1, x.s1[k], y.f2)


def usecase5(x, N):
    """
    Base on test4 of https://github.com/numba/numba/issues/381
    """
    for k in range(N):
        print(x[k].f1, x.s1[k], x[k].f2)


class TestRecordUsecase(TestCase):

    def setUp(self):
        fields = [('f1', '<f8'), ('s1', '|S3'), ('f2', '<f8')]
        self.unaligned_dtype = np.dtype(fields)
        self.aligned_dtype = np.dtype(fields, align=True)

    def test_usecase1(self):
        pyfunc = usecase1

        # This is an unaligned dtype
        mystruct_dt = np.dtype([('p', np.float64),
                           ('row', np.float64),
                           ('col', np.float64)])
        mystruct = numpy_support.from_dtype(mystruct_dt)

        cfunc = njit((mystruct[:], mystruct[:]))(pyfunc)

        st1 = np.recarray(3, dtype=mystruct_dt)
        st2 = np.recarray(3, dtype=mystruct_dt)

        st1.p = np.arange(st1.size) + 1
        st1.row = np.arange(st1.size) + 1
        st1.col = np.arange(st1.size) + 1

        st2.p = np.arange(st2.size) + 1
        st2.row = np.arange(st2.size) + 1
        st2.col = np.arange(st2.size) + 1

        expect1 = st1.copy()
        expect2 = st2.copy()

        got1 = expect1.copy()
        got2 = expect2.copy()

        pyfunc(expect1, expect2)
        cfunc(got1, got2)

        np.testing.assert_equal(expect1, got1)
        np.testing.assert_equal(expect2, got2)

    def _setup_usecase2to5(self, dtype):
        N = 5
        a = np.recarray(N, dtype=dtype)
        a.f1 = np.arange(N)
        a.f2 = np.arange(2, N + 2)
        a.s1 = np.array(['abc'] * a.shape[0], dtype='|S3')
        return a

    def _test_usecase2to5(self, pyfunc, dtype):
        array = self._setup_usecase2to5(dtype)
        record_type = numpy_support.from_dtype(dtype)
        cfunc = njit((record_type[:], types.intp))(pyfunc)

        with captured_stdout():
            pyfunc(array, len(array))
            expect = sys.stdout.getvalue()

        with captured_stdout():
            cfunc(array, len(array))
            got = sys.stdout.getvalue()

        self.assertEqual(expect, got)

    def test_usecase2(self):
        self._test_usecase2to5(usecase2, self.unaligned_dtype)
        self._test_usecase2to5(usecase2, self.aligned_dtype)

    def test_usecase3(self):
        self._test_usecase2to5(usecase3, self.unaligned_dtype)
        self._test_usecase2to5(usecase3, self.aligned_dtype)

    def test_usecase4(self):
        self._test_usecase2to5(usecase4, self.unaligned_dtype)
        self._test_usecase2to5(usecase4, self.aligned_dtype)

    def test_usecase5(self):
        self._test_usecase2to5(usecase5, self.unaligned_dtype)
        self._test_usecase2to5(usecase5, self.aligned_dtype)


if __name__ == '__main__':
    unittest.main()
