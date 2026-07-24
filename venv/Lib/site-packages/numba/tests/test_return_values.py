"""
Test return values
"""


import math

import unittest
from numba import jit
from numba.core import types
from numba.core.errors import TypingError, NumbaTypeError


enable_pyobj_flags = {'forceobj': True}
no_pyobj_flags = {'nopython': True}


def get_nopython_func():
    return abs

def get_pyobj_func():
    return open

def get_module_func():
    return math.floor


class TestReturnValues(unittest.TestCase):

    def test_nopython_func(self, flags=enable_pyobj_flags):
        # Test returning func that is supported in nopython mode
        pyfunc = get_nopython_func
        cfunc = jit((), **flags)(pyfunc)
        if flags == enable_pyobj_flags:
            result = cfunc()
            self.assertEqual(result, abs)
        else:
            self.fail("Unexpected successful compilation.")

    def test_nopython_func_npm(self):
        with self.assertRaises(NumbaTypeError):
            self.test_nopython_func(flags=no_pyobj_flags)

    def test_pyobj_func(self, flags=enable_pyobj_flags):
        # Test returning func that is only supported in object mode
        pyfunc = get_pyobj_func
        cfunc = jit((), **flags)(pyfunc)
        if flags == enable_pyobj_flags:
            result = cfunc()
            self.assertEqual(result, open)
        else:
            self.fail("Unexpected successful compilation.")

    def test_pyobj_func_npm(self):
        with self.assertRaises(TypingError):
            self.test_pyobj_func(flags=no_pyobj_flags)

    def test_module_func(self, flags=enable_pyobj_flags):
        # Test returning imported func that is only supported in object mode
        pyfunc = get_module_func
        cfunc = jit((), **flags)(pyfunc)
        if flags == enable_pyobj_flags:
            result = cfunc()
            self.assertEqual(result, math.floor)
        else:
            self.fail("Unexpected successful compilation.")

    def test_module_func_npm(self):
        with self.assertRaises(NumbaTypeError):
            self.test_module_func(flags=no_pyobj_flags)


if __name__ == '__main__':
    unittest.main()
