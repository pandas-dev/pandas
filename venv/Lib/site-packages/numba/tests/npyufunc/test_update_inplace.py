# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import, division

import unittest

import numpy as np
from numba import guvectorize
from numba.tests.support import TestCase


def py_replace_2nd(x_t, y_1):
    for t in range(0, x_t.shape[0], 2):
        x_t[t] = y_1[0]


def py_update_3(x0_t, x1_t, x2_t, y_1):
    for t in range(0, x0_t.shape[0]):
        x0_t[t] = y_1[0]
        x1_t[t] = 2 * y_1[0]
        x2_t[t] = 3 * y_1[0]


class TestUpdateInplace(TestCase):

    def _run_test_for_gufunc(self, gufunc, py_func, expect_f4_to_pass=True,
                             z=2):
        for dtype, expect_to_pass in [('f8', True), ('f4', expect_f4_to_pass)]:
            inputs = [np.zeros(10, dtype) for _ in range(gufunc.nin - 1)]
            ex_inputs = [x_t.copy() for x_t in inputs]

            gufunc(*inputs, z)
            py_func(*ex_inputs, np.array([z]))

            for i, (x_t, ex_x_t) in enumerate(zip(inputs, ex_inputs)):
                if expect_to_pass:
                    np.testing.assert_equal(x_t, ex_x_t, err_msg='input %s' % i)
                else:
                    self.assertFalse((x_t == ex_x_t).all(), msg='input %s' % i)

    def test_update_inplace(self):
        # test without writable_args
        gufunc = guvectorize(['void(f8[:], f8[:])'], '(t),()',
                             nopython=True)(py_replace_2nd)
        self._run_test_for_gufunc(gufunc, py_replace_2nd,
                                  expect_f4_to_pass=False)

        # test with writable_args
        gufunc = guvectorize(['void(f8[:], f8[:])'], '(t),()',
                             nopython=True, writable_args=(0,))(py_replace_2nd)
        self._run_test_for_gufunc(gufunc, py_replace_2nd)

        # test with writable_args as strings
        gufunc = guvectorize(['void(f8[:], f8[:])'], '(t),()', nopython=True,
                             writable_args=('x_t',))(py_replace_2nd)
        self._run_test_for_gufunc(gufunc, py_replace_2nd)

    def test_update_inplace_with_cache(self):
        # test with writable_args
        gufunc = guvectorize(['void(f8[:], f8[:])'], '(t),()',
                             nopython=True, writable_args=(0,),
                             cache=True)(py_replace_2nd)
        # 2nd time it is loaded from cache
        gufunc = guvectorize(['void(f8[:], f8[:])'], '(t),()',
                             nopython=True, writable_args=(0,),
                             cache=True)(py_replace_2nd)
        self._run_test_for_gufunc(gufunc, py_replace_2nd)

    def test_update_inplace_parallel(self):
        # test with writable_args
        gufunc = guvectorize(['void(f8[:], f8[:])'], '(t),()',
                             nopython=True, writable_args=(0,),
                             target='parallel')(py_replace_2nd)
        self._run_test_for_gufunc(gufunc, py_replace_2nd)

    def test_update_inplace_3(self):
        # test without writable_args
        gufunc = guvectorize(['void(f8[:], f8[:], f8[:], f8[:])'],
                             '(t),(t),(t),()',
                             nopython=True)(py_update_3)
        self._run_test_for_gufunc(gufunc, py_update_3, expect_f4_to_pass=False)

        # test with writable_args
        gufunc = guvectorize(['void(f8[:], f8[:], f8[:], f8[:])'],
                             '(t),(t),(t),()', nopython=True,
                             writable_args=(0, 1, 2))(py_update_3)
        self._run_test_for_gufunc(gufunc, py_update_3)

        # test with writable_args as mix of strings and ints
        gufunc = guvectorize(['void(f8[:], f8[:], f8[:], f8[:])'],
                             '(t),(t),(t),()', nopython=True,
                             writable_args=('x0_t', 'x1_t', 2))(py_update_3)
        self._run_test_for_gufunc(gufunc, py_update_3)

    def test_exceptions(self):
        # check that len(writable_args) <= nin
        with self.assertRaises(ValueError):
            guvectorize(['void(f8[:], f8[:])'], '(t),()', nopython=True,
                        writable_args=(0, 1, 2, 5))(py_replace_2nd)

        # check that all values in writable_args are between 0 and nin
        with self.assertRaises(ValueError):
            guvectorize(['void(f8[:], f8[:])'], '(t),()',
                        nopython=True, writable_args=(5,))(py_replace_2nd)

        with self.assertRaises(ValueError):
            guvectorize(['void(f8[:], f8[:])'], '(t),()',
                        nopython=True, writable_args=(-1,))(py_replace_2nd)

        # check that exception is raised when passing non-existing argument name
        with self.assertRaises(RuntimeError):
            guvectorize(['void(f8[:], f8[:])'], '(t),()',
                        nopython=True, writable_args=('z_t',))(py_replace_2nd)

        # writable_args are not supported for target='cuda'
        with self.assertRaises(TypeError):
            guvectorize(['void(f8[:], f8[:])'], '(t),()',
                        nopython=True, writable_args=(0,),
                        target='cuda')(py_replace_2nd)


if __name__ == '__main__':
    unittest.main()
