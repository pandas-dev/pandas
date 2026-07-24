"""
This file will be copied to a temporary directory in order to
exercise caching compiled C callbacks.

See test_cfunc.py.
"""

import sys

from numba import cfunc, jit
from numba.tests.support import TestCase, captured_stderr


Z = 1

add_sig = "float64(float64, float64)"

div_sig = "float64(int64, int64)"


@cfunc(add_sig, cache=True, nopython=True)
def add_usecase(x, y):
    return x + y + Z

@cfunc(add_sig, nopython=True)
def add_nocache_usecase(x, y):
    return x + y + Z

@cfunc(div_sig, cache=True, nopython=True)
def div_usecase(a, b):
    return a / b


@jit(nopython=True)
def inner(x, y):
    return x + y + Z

@cfunc(add_sig, cache=True, nopython=True)
def outer(x, y):
    return inner(-y, x)


class _TestModule(TestCase):
    """
    Tests for functionality of this module's cfuncs.
    Note this does not define any "test_*" method, instead check_module()
    should be called by hand.
    """

    def check_module(self, mod):
        f = mod.add_usecase
        self.assertPreciseEqual(f.ctypes(2.0, 3.0), 6.0)
        f = mod.add_nocache_usecase
        self.assertPreciseEqual(f.ctypes(2.0, 3.0), 6.0)
        f = mod.outer
        self.assertPreciseEqual(f.ctypes(5.0, 2.0), 4.0)

        f = mod.div_usecase
        with captured_stderr() as err:
            self.assertPreciseEqual(f.ctypes(7, 2), 3.5)
        self.assertEqual(err.getvalue(), "")
        with captured_stderr() as err:
            f.ctypes(7, 0)
        err = err.getvalue()
        self.assertIn("ZeroDivisionError", err)


def self_test():
    mod = sys.modules[__name__]
    _TestModule().check_module(mod)
