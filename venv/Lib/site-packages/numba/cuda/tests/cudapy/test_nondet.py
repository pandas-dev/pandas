import numpy as np
from numba import cuda, float32, void
from numba.cuda.testing import unittest, CUDATestCase


def generate_input(n):
    A = np.array(np.arange(n * n).reshape(n, n), dtype=np.float32)
    B = np.array(np.arange(n) + 0, dtype=A.dtype)
    return A, B


class TestCudaNonDet(CUDATestCase):
    def test_for_pre(self):
        """Test issue with loop not running due to bad sign-extension at the for
        loop precondition.
        """

        @cuda.jit(void(float32[:, :], float32[:, :], float32[:]))
        def diagproduct(c, a, b):
            startX, startY = cuda.grid(2)
            gridX = cuda.gridDim.x * cuda.blockDim.x
            gridY = cuda.gridDim.y * cuda.blockDim.y
            height = c.shape[0]
            width = c.shape[1]

            for x in range(startX, width, (gridX)):
                for y in range(startY, height, (gridY)):
                    c[y, x] = a[y, x] * b[x]

        N = 8

        A, B = generate_input(N)

        F = np.empty(A.shape, dtype=A.dtype)

        blockdim = (32, 8)
        griddim = (1, 1)

        dA = cuda.to_device(A)
        dB = cuda.to_device(B)
        dF = cuda.to_device(F, copy=False)
        diagproduct[griddim, blockdim](dF, dA, dB)

        E = np.dot(A, np.diag(B))
        np.testing.assert_array_almost_equal(dF.copy_to_host(), E)


if __name__ == '__main__':
    unittest.main()
