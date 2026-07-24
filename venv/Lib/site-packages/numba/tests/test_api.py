import warnings

import numba
from numba import jit, njit

from numba.tests.support import TestCase, always_test
import unittest


class TestNumbaModule(TestCase):
    """
    Test the APIs exposed by the top-level `numba` module.
    """

    def check_member(self, name):
        self.assertTrue(hasattr(numba, name), name)
        self.assertIn(name, numba.__all__)

    @always_test
    def test_numba_module(self):
        # jit
        self.check_member("jit")
        self.check_member("vectorize")
        self.check_member("guvectorize")
        self.check_member("njit")
        # errors
        self.check_member("NumbaError")
        self.check_member("TypingError")
        # types
        self.check_member("int32")
        # misc
        numba.__version__  # not in __all__


class TestJitDecorator(TestCase):
    """
    Test the jit and njit decorators
    """
    def test_jit_nopython_forceobj(self):
        with self.assertRaises(ValueError) as cm:
            jit(nopython=True, forceobj=True)
        self.assertIn(
            "Only one of 'nopython' or 'forceobj' can be True.",
            str(cm.exception)
        )

        def py_func(x):
            return x

        jit_func = jit(nopython=True)(py_func)
        jit_func(1)
        # Check length of nopython_signatures to check
        # which mode the function was compiled in
        self.assertEqual(len(jit_func.nopython_signatures), 1)

        jit_func = jit(forceobj=True)(py_func)
        jit_func(1)
        self.assertEqual(len(jit_func.nopython_signatures), 0)

    def test_njit_nopython_forceobj(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', RuntimeWarning)
            njit(forceobj=True)
        self.assertEqual(len(w), 1)
        self.assertIn(
            'forceobj is set for njit and is ignored', str(w[0].message)
        )

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', RuntimeWarning)
            njit(nopython=True)
        self.assertEqual(len(w), 1)
        self.assertIn(
            'nopython is set for njit and is ignored', str(w[0].message)
        )

        def py_func(x):
            return x

        jit_func = njit(nopython=True)(py_func)
        jit_func(1)
        self.assertEqual(len(jit_func.nopython_signatures), 1)

        jit_func = njit(forceobj=True)(py_func)
        jit_func(1)
        # Since forceobj is ignored this has to compile in nopython mode
        self.assertEqual(len(jit_func.nopython_signatures), 1)


if __name__ == '__main__':
    unittest.main()
