"""
Tests issues or edge cases for producing invalid NRT refct
"""


import gc

import numpy as np

import unittest
from numba import njit
from numba.core.runtime import rtsys
from numba.tests.support import TestCase, EnableNRTStatsMixin


class TestNrtRefCt(EnableNRTStatsMixin, TestCase):

    def setUp(self):
        # Clean up any NRT-backed objects hanging in a dead reference cycle
        gc.collect()
        super(TestNrtRefCt, self).setUp()

    def test_no_return(self):
        """
        Test issue #1291
        """

        @njit
        def foo(n):
            for i in range(n):
                temp = np.zeros(2)
            return 0

        n = 10
        init_stats = rtsys.get_allocation_stats()
        foo(n)
        cur_stats = rtsys.get_allocation_stats()
        self.assertEqual(cur_stats.alloc - init_stats.alloc, n)
        self.assertEqual(cur_stats.free - init_stats.free, n)

    def test_escaping_var_init_in_loop(self):
        """
        Test issue #1297
        """

        @njit
        def g(n):

            x = np.zeros((n, 2))

            for i in range(n):
                y = x[i]

            for i in range(n):
                y = x[i]

            return 0

        init_stats = rtsys.get_allocation_stats()
        g(10)
        cur_stats = rtsys.get_allocation_stats()
        self.assertEqual(cur_stats.alloc - init_stats.alloc, 1)
        self.assertEqual(cur_stats.free - init_stats.free, 1)

    def test_invalid_computation_of_lifetime(self):
        """
        Test issue #1573
        """
        @njit
        def if_with_allocation_and_initialization(arr1, test1):
            tmp_arr = np.zeros_like(arr1)

            for i in range(tmp_arr.shape[0]):
                pass

            if test1:
                np.zeros_like(arr1)

            return tmp_arr

        arr = np.random.random((5, 5))  # the values are not consumed

        init_stats = rtsys.get_allocation_stats()
        if_with_allocation_and_initialization(arr, False)
        cur_stats = rtsys.get_allocation_stats()
        self.assertEqual(cur_stats.alloc - init_stats.alloc,
                         cur_stats.free - init_stats.free)

    def test_del_at_beginning_of_loop(self):
        """
        Test issue #1734
        """
        @njit
        def f(arr):
            res = 0

            for i in (0, 1):
                # `del t` is issued here before defining t.  It must be
                # correctly handled by the lowering phase.
                t = arr[i]
                if t[i] > 1:
                    res += t[i]

            return res

        arr = np.ones((2, 2))
        init_stats = rtsys.get_allocation_stats()
        f(arr)
        cur_stats = rtsys.get_allocation_stats()
        self.assertEqual(cur_stats.alloc - init_stats.alloc,
                         cur_stats.free - init_stats.free)


if __name__ == '__main__':
    unittest.main()
