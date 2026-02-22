import numpy as np
from numba import cuda
from numba.core.config import ENABLE_CUDASIM
from numba.cuda.testing import CUDATestCase
import unittest

# Avoid recompilation of the sum_reduce function by keeping it at global scope
sum_reduce = cuda.Reduce(lambda a, b: a + b)


class TestReduction(CUDATestCase):
    def _sum_reduce(self, n):
        A = (np.arange(n, dtype=np.float64) + 1)
        expect = A.sum()
        got = sum_reduce(A)
        self.assertEqual(expect, got)

    def test_sum_reduce(self):
        if ENABLE_CUDASIM:
            # Minimal test set for the simulator (which only wraps
            # functools.reduce)
            test_sizes = [ 1, 16 ]
        else:
            # Tests around the points where blocksize changes, and around larger
            # powers of two, sums of powers of two, and some "random" sizes
            test_sizes = [ 1, 15, 16, 17, 127, 128, 129, 1023, 1024,
                           1025, 1536, 1048576, 1049600, 1049728, 34567 ]
        # Avoid recompilation by keeping sum_reduce here
        for n in test_sizes:
            self._sum_reduce(n)

    def test_empty_array_host(self):
        A = (np.arange(0, dtype=np.float64) + 1)
        expect = A.sum()
        got = sum_reduce(A)
        self.assertEqual(expect, got)

    def test_empty_array_device(self):
        A = (np.arange(0, dtype=np.float64) + 1)
        dA = cuda.to_device(A)
        expect = A.sum()
        got = sum_reduce(dA)
        self.assertEqual(expect, got)

    def test_prod_reduce(self):
        prod_reduce = cuda.reduce(lambda a, b: a * b)
        A = (np.arange(64, dtype=np.float64) + 1)
        expect = A.prod()
        got = prod_reduce(A, init=1)
        np.testing.assert_allclose(expect, got)

    def test_max_reduce(self):
        max_reduce = cuda.Reduce(lambda a, b: max(a, b))
        A = (np.arange(3717, dtype=np.float64) + 1)
        expect = A.max()
        got = max_reduce(A, init=0)
        self.assertEqual(expect, got)

    def test_non_identity_init(self):
        init = 3
        A = (np.arange(10, dtype=np.float64) + 1)
        expect = A.sum() + init
        got = sum_reduce(A, init=init)
        self.assertEqual(expect, got)

    def test_result_on_device(self):
        A = (np.arange(10, dtype=np.float64) + 1)
        got = cuda.to_device(np.zeros(1, dtype=np.float64))
        expect = A.sum()
        res = sum_reduce(A, res=got)
        self.assertIsNone(res)
        self.assertEqual(expect, got[0])


if __name__ == '__main__':
    unittest.main()
