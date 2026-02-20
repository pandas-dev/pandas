import numpy as np

from numba import cuda
import unittest
from numba.cuda.testing import CUDATestCase


@cuda.jit
def foo(x):
    i = cuda.grid(1)
    if i < x.size:
        x[i] += 1


class TestForAll(CUDATestCase):
    def test_forall_1(self):
        arr = np.arange(11)
        orig = arr.copy()
        foo.forall(arr.size)(arr)
        np.testing.assert_array_almost_equal(arr, orig + 1)

    def test_forall_2(self):
        @cuda.jit("void(float32, float32[:], float32[:])")
        def bar(a, x, y):
            i = cuda.grid(1)
            if i < x.size:
                y[i] = a * x[i] + y[i]

        x = np.arange(13, dtype=np.float32)
        y = np.arange(13, dtype=np.float32)
        oldy = y.copy()
        a = 1.234
        bar.forall(y.size)(a, x, y)
        np.testing.assert_array_almost_equal(y, a * x + oldy, decimal=3)

    def test_forall_no_work(self):
        # Ensure that forall doesn't launch a kernel with no blocks when called
        # with 0 elements. See Issue #5017.
        arr = np.arange(11)
        foo.forall(0)(arr)

    def test_forall_negative_work(self):
        # Ensure that forall doesn't allow the creation of a forall with a
        # negative element count.
        with self.assertRaises(ValueError) as raises:
            foo.forall(-1)
        self.assertIn("Can't create ForAll with negative task count",
                      str(raises.exception))


if __name__ == '__main__':
    unittest.main()
