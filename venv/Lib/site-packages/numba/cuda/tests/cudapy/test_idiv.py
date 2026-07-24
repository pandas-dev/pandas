import numpy as np
from numba import cuda, float32, float64, int32, void
from numba.cuda.testing import unittest, CUDATestCase


class TestCudaIDiv(CUDATestCase):
    def test_inplace_div(self):

        @cuda.jit(void(float32[:, :], int32, int32))
        def div(grid, l_x, l_y):
            for x in range(l_x):
                for y in range(l_y):
                    grid[x, y] /= 2.0

        x = np.ones((2, 2), dtype=np.float32)
        grid = cuda.to_device(x)
        div[1, 1](grid, 2, 2)
        y = grid.copy_to_host()
        self.assertTrue(np.all(y == 0.5))

    def test_inplace_div_double(self):

        @cuda.jit(void(float64[:, :], int32, int32))
        def div_double(grid, l_x, l_y):
            for x in range(l_x):
                for y in range(l_y):
                    grid[x, y] /= 2.0

        x = np.ones((2, 2), dtype=np.float64)
        grid = cuda.to_device(x)
        div_double[1, 1](grid, 2, 2)
        y = grid.copy_to_host()
        self.assertTrue(np.all(y == 0.5))


if __name__ == '__main__':
    unittest.main()
