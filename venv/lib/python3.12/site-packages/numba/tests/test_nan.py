import unittest
from numba import jit
from numba.core import types

enable_pyobj_flags = {'forceobj': True}
no_pyobj_flags = {'nopython': True}


def isnan(x):
    return x != x


def isequal(x):
    return x == x


class TestNaN(unittest.TestCase):

    def test_nans(self, flags=enable_pyobj_flags):
        pyfunc = isnan
        cfunc = jit((types.float64,), **flags)(pyfunc)

        self.assertTrue(cfunc(float('nan')))
        self.assertFalse(cfunc(1.0))

        pyfunc = isequal
        cfunc = jit((types.float64,), **flags)(pyfunc)

        self.assertFalse(cfunc(float('nan')))
        self.assertTrue(cfunc(1.0))

    def test_nans_npm(self):
        self.test_nans(flags=no_pyobj_flags)


if __name__ == '__main__':
    unittest.main()
