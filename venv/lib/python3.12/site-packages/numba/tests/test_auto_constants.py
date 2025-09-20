import math
import sys

import numpy as np

from numba import njit
import numba.tests.usecases as uc
import unittest


class TestAutoConstants(unittest.TestCase):
    def test_numpy_nan(self):

        @njit
        def f():
            return np.nan

        self.assertTrue(math.isnan(f()))
        self.assertTrue(math.isnan(f.py_func()))

    def test_sys_constant(self):

        @njit
        def f():
            return sys.hexversion

        self.assertEqual(f(), f.py_func())

    def test_module_string_constant(self):

        @njit
        def f():
            return uc._GLOBAL_STR
        self.assertEqual(f(), f.py_func())


if __name__ == '__main__':
    unittest.main()
