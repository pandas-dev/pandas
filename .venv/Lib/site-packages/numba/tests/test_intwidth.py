import unittest

import math
import sys

from numba import jit
from numba.core import utils
from numba.tests.support import TestCase, tag


max_uint64 = 18446744073709551615

def usecase_uint64_global():
    return max_uint64

def usecase_uint64_constant():
    return 18446744073709551615

def usecase_uint64_func():
    return max(18446744073709551614, 18446744073709551615)

def usecase_int64_pos():
    return 9223372036854775807

def usecase_int64_neg():
    return -9223372036854775808

def usecase_int64_func():
    return (max(9223372036854775807, -9223372036854775808)
            + min(9223372036854775807, -9223372036854775808))


class IntWidthTest(TestCase):

    def check_nullary_func(self, pyfunc, **kwargs):
        cfunc = jit(**kwargs)(pyfunc)
        self.assertPreciseEqual(cfunc(), pyfunc())

    def test_global_uint64(self, nopython=False):
        pyfunc = usecase_uint64_global
        self.check_nullary_func(pyfunc, nopython=nopython)

    def test_global_uint64_npm(self):
        self.test_global_uint64(nopython=True)

    def test_constant_uint64(self, nopython=False):
        pyfunc = usecase_uint64_constant
        self.check_nullary_func(pyfunc, nopython=nopython)

    def test_constant_uint64_npm(self):
        self.test_constant_uint64(nopython=True)

    def test_constant_uint64_function_call(self, nopython=False):
        pyfunc = usecase_uint64_func
        self.check_nullary_func(pyfunc, nopython=nopython)

    def test_constant_uint64_function_call_npm(self):
        self.test_constant_uint64_function_call(nopython=True)

    def test_bit_length(self):
        f = utils.bit_length
        self.assertEqual(f(0x7f), 7)
        self.assertEqual(f(-0x7f), 7)
        self.assertEqual(f(0x80), 8)
        self.assertEqual(f(-0x80), 7)
        self.assertEqual(f(0xff), 8)
        self.assertEqual(f(-0xff), 8)
        self.assertEqual(f(0x100), 9)
        self.assertEqual(f(-0x100), 8)
        self.assertEqual(f(-0x101), 9)
        self.assertEqual(f(0x7fffffff), 31)
        self.assertEqual(f(-0x7fffffff), 31)
        self.assertEqual(f(-0x80000000), 31)
        self.assertEqual(f(0x80000000), 32)
        self.assertEqual(f(0xffffffff), 32)
        self.assertEqual(f(0xffffffffffffffff), 64)
        self.assertEqual(f(0x10000000000000000), 65)

    def test_constant_int64(self, nopython=False):
        self.check_nullary_func(usecase_int64_pos, nopython=nopython)
        self.check_nullary_func(usecase_int64_neg, nopython=nopython)
        self.check_nullary_func(usecase_int64_func, nopython=nopython)

    def test_constant_int64_npm(self):
        self.test_constant_int64(nopython=True)


if __name__ == '__main__':
    unittest.main()

