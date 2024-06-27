import numpy as np

import unittest
from numba import njit
from numba.core import types, errors
from numba.tests.support import TestCase


def setitem_slice(a, start, stop, step, scalar):
    a[start:stop:step] = scalar


def usecase(obs, nPoints):
    center = nPoints // 2
    obs[0:center] = np.arange(center)
    obs[center] = 321
    obs[(center + 1):] = np.arange(nPoints - center - 1)


class TestStoreSlice(TestCase):

    def test_usecase(self):
        n = 10
        obs_got = np.zeros(n)
        obs_expected = obs_got.copy()

        cfunc = njit((types.float64[:], types.intp))(usecase)
        cfunc(obs_got, n)
        usecase(obs_expected, n)

        self.assertPreciseEqual(obs_got, obs_expected)

    def test_array_slice_setitem(self):
        n = 10
        argtys = (types.int64[:], types.int64, types.int64, types.int64,
                  types.int64)
        cfunc = njit(argtys)(setitem_slice)
        a = np.arange(n, dtype=np.int64)
        # tuple is (start, stop, step, scalar)
        tests = ((2, 6, 1, 7),
                 (2, 6, -1, 7),
                 (-2, len(a), 2, 77),
                 (-2, 2 * len(a), 2, 77),
                 (-2, -6, 3, 88),
                 (-2, -6, -3, 9999),
                 (-6, -2, 4, 88),
                 (-6, -2, -4, 88),
                 (16, 20, 2, 88),
                 (16, 20, -2, 88),
                 )

        for start, stop, step, scalar in tests:
            a = np.arange(n, dtype=np.int64)
            b = np.arange(n, dtype=np.int64)
            cfunc(a, start, stop, step, scalar)
            setitem_slice(b, start, stop, step, scalar)
            self.assertPreciseEqual(a, b)

        # test if step = 0
        a = np.arange(n, dtype=np.int64)
        with self.assertRaises(ValueError) as cm:
            cfunc(a, 3, 6, 0, 88)
        self.assertEqual(str(cm.exception), "slice step cannot be zero")


if __name__ == '__main__':
    unittest.main()

