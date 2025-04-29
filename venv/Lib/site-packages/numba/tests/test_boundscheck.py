import numpy as np

from numba.cuda.testing import SerialMixin
from numba import typeof, cuda, njit
from numba.core.types import float64
from numba.tests.support import TestCase, MemoryLeakMixin
from numba.core import config
import unittest


def basic_array_access(a):
    return a[10]


def slice_array_access(a):
    # The first index (slice) is not bounds checked
    return a[10:, 10]


def fancy_array_access(x):
    a = np.array([1, 2, 3])
    return x[a]


def fancy_array_modify(x):
    a = np.array([1, 2, 3])
    x[a] = 0
    return x


class TestBoundsCheckNoError(MemoryLeakMixin, TestCase):

    @TestCase.run_test_in_subprocess(envvars={'NUMBA_BOUNDSCHECK': ''})
    def test_basic_array_boundscheck(self):
        self.assertIsNone(config.BOUNDSCHECK)

        a = np.arange(5)
        # Check the numpy behavior to make sure the test is correct
        with self.assertRaises(IndexError):
            # TODO: When we raise the same error message as numpy, test that
            # they are the same
            basic_array_access(a)

        at = typeof(a)
        noboundscheck = njit((at,))(basic_array_access)
        # Check that the default flag doesn't raise
        noboundscheck(a)
        # boundscheck(a) is tested in TestBoundsCheckError below

    @TestCase.run_test_in_subprocess(envvars={'NUMBA_BOUNDSCHECK': ''})
    def test_slice_array_boundscheck(self):
        self.assertIsNone(config.BOUNDSCHECK)

        a = np.ones((5, 5))
        b = np.ones((5, 20))
        with self.assertRaises(IndexError):
            # TODO: When we raise the same error message as numpy, test that
            # they are the same
            slice_array_access(a)
        # Out of bounds on a slice doesn't raise
        slice_array_access(b)

        at = typeof(a)
        rt = float64[:]
        noboundscheck = njit(rt(at))(slice_array_access)
        boundscheck = njit(rt(at), boundscheck=True)(slice_array_access)
        # Check that the default flag doesn't raise
        noboundscheck(a)
        noboundscheck(b)
        # boundscheck(a) is tested in TestBoundsCheckError below

        # Doesn't raise
        boundscheck(b)

    @TestCase.run_test_in_subprocess(envvars={'NUMBA_BOUNDSCHECK': ''})
    def test_fancy_indexing_boundscheck(self):
        self.assertIsNone(config.BOUNDSCHECK)

        a = np.arange(3)
        b = np.arange(4)

        # Check the numpy behavior to ensure the test is correct.
        with self.assertRaises(IndexError):
            # TODO: When we raise the same error message as numpy, test that
            # they are the same
            fancy_array_access(a)
        fancy_array_access(b)

        at = typeof(a)
        rt = at.dtype[:]
        noboundscheck = njit(rt(at))(fancy_array_access)
        boundscheck = njit(rt(at), boundscheck=True)(fancy_array_access)
        # Check that the default flag doesn't raise
        noboundscheck(a)
        noboundscheck(b)
        # boundscheck(a) is tested in TestBoundsCheckError below

        # Doesn't raise
        boundscheck(b)


class TestNoCudaBoundsCheck(SerialMixin, TestCase):
    @unittest.skipIf(not cuda.is_available(), "NO CUDA")
    @TestCase.run_test_in_subprocess(envvars={'NUMBA_BOUNDSCHECK': '1'})
    def test_no_cuda_boundscheck(self):
        self.assertTrue(config.BOUNDSCHECK)
        with self.assertRaises(NotImplementedError):
            @cuda.jit(boundscheck=True)
            def func():
                pass

        # Make sure we aren't raising "not supported" error if we aren't
        # requesting bounds checking anyway. Related pull request: #5257
        @cuda.jit(boundscheck=False)
        def func3():
            pass

        @cuda.jit
        def func2(x, a):
            a[1] = x[1]

        a = np.ones((1,))
        x = np.zeros((1,))
        # Out of bounds but doesn't raise (it does raise in the simulator,
        # so skip there)
        if not config.ENABLE_CUDASIM:
            func2[1, 1](x, a)


# This is a separate test because the jitted functions that raise exceptions
# have memory leaks.
class TestBoundsCheckError(TestCase):
    @TestCase.run_test_in_subprocess(envvars={'NUMBA_BOUNDSCHECK': ''})
    def test_basic_array_boundscheck(self):
        self.assertIsNone(config.BOUNDSCHECK)

        a = np.arange(5)
        # Check the numpy behavior to make sure the test is correct
        with self.assertRaises(IndexError):
            # TODO: When we raise the same error message as numpy, test that
            # they are the same
            basic_array_access(a)

        at = typeof(a)
        boundscheck = njit((at,), boundscheck=True)(basic_array_access)

        with self.assertRaises(IndexError):
            boundscheck(a)

    @TestCase.run_test_in_subprocess(envvars={'NUMBA_BOUNDSCHECK': ''})
    def test_slice_array_boundscheck(self):
        self.assertIsNone(config.BOUNDSCHECK)

        a = np.ones((5, 5))
        b = np.ones((5, 20))
        with self.assertRaises(IndexError):
            # TODO: When we raise the same error message as numpy, test that
            # they are the same
            slice_array_access(a)
        # Out of bounds on a slice doesn't raise
        slice_array_access(b)

        at = typeof(a)
        rt = float64[:]
        boundscheck = njit(rt(at), boundscheck=True)(slice_array_access)

        with self.assertRaises(IndexError):
            boundscheck(a)

    @TestCase.run_test_in_subprocess(envvars={'NUMBA_BOUNDSCHECK': ''})
    def test_fancy_indexing_boundscheck(self):
        self.assertIsNone(config.BOUNDSCHECK)

        a = np.arange(3)
        b = np.arange(4)

        # Check the numpy behavior to ensure the test is correct.
        with self.assertRaises(IndexError):
            # TODO: When we raise the same error message as numpy, test that
            # they are the same
            fancy_array_access(a)
        fancy_array_access(b)

        at = typeof(a)
        rt = at.dtype[:]
        boundscheck = njit(rt(at), boundscheck=True)(fancy_array_access)

        with self.assertRaises(IndexError):
            boundscheck(a)

    @TestCase.run_test_in_subprocess(envvars={'NUMBA_BOUNDSCHECK': ''})
    def test_fancy_indexing_with_modification_boundscheck(self):
        self.assertIsNone(config.BOUNDSCHECK)

        a = np.arange(3)
        b = np.arange(4)

        # Check the numpy behavior to ensure the test is correct.
        with self.assertRaises(IndexError):
            # TODO: When we raise the same error message as numpy, test that
            # they are the same
            fancy_array_modify(a)
        fancy_array_modify(b)

        at = typeof(a)
        rt = at.dtype[:]
        boundscheck = njit(rt(at), boundscheck=True)(fancy_array_modify)

        with self.assertRaises(IndexError):
            boundscheck(a)


class TestBoundsEnvironmentVariable(TestCase):
    def setUp(self):
        @njit
        def default(x):
            return x[1]

        @njit(boundscheck=False)
        def off(x):
            return x[1]

        @njit(boundscheck=True)
        def on(x):
            return x[1]

        self.default = default
        self.off = off
        self.on = on

    @TestCase.run_test_in_subprocess(envvars={'NUMBA_BOUNDSCHECK': ''})
    def test_boundscheck_unset(self):
        self.assertIsNone(config.BOUNDSCHECK)

        a = np.array([1])

        # Doesn't raise
        self.default(a)
        self.off(a)

        with self.assertRaises(IndexError):
            self.on(a)

    @TestCase.run_test_in_subprocess(envvars={'NUMBA_BOUNDSCHECK': '1'})
    def test_boundscheck_enabled(self):
        self.assertTrue(config.BOUNDSCHECK)

        a = np.array([1])

        with self.assertRaises(IndexError):
            self.default(a)
            self.off(a)
            self.on(a)

    @TestCase.run_test_in_subprocess(envvars={'NUMBA_BOUNDSCHECK': '0'})
    def test_boundscheck_disabled(self):
        self.assertFalse(config.BOUNDSCHECK)

        a = np.array([1])

        # Doesn't raise
        self.default(a)
        self.off(a)
        self.on(a)


if __name__ == '__main__':
    unittest.main()
