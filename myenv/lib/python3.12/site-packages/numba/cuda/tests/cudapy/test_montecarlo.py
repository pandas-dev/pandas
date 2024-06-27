import math
from numba import cuda
from numba.cuda.testing import unittest, CUDATestCase


class TestCudaMonteCarlo(CUDATestCase):
    def test_montecarlo(self):
        """Just make sure we can compile this
        """

        @cuda.jit(
            'void(double[:], double[:], double, double, double, double[:])')
        def step(last, paths, dt, c0, c1, normdist):
            i = cuda.grid(1)
            if i >= paths.shape[0]:
                return
            noise = normdist[i]
            paths[i] = last[i] * math.exp(c0 * dt + c1 * noise)


if __name__ == '__main__':
    unittest.main()
