import numpy as np
from numba import cuda, float32, int32, void
from numba.cuda.testing import unittest, CUDATestCase


class TestCudaPy2Div(CUDATestCase):
    def test_py2_div_issue(self):
        @cuda.jit(void(float32[:], float32[:], float32[:], int32))
        def preCalc(y, yA, yB, numDataPoints):
            i = cuda.grid(1)
            # k is unused, but may be part of the trigger for the bug this
            # tests for.
            k = i % numDataPoints  # noqa: F841

            ans = float32(1.001 * float32(i))

            y[i] = ans
            yA[i] = ans * 1.0
            yB[i] = ans / 1.0

        numDataPoints = 15

        y = np.zeros(numDataPoints, dtype=np.float32)
        yA = np.zeros(numDataPoints, dtype=np.float32)
        yB = np.zeros(numDataPoints, dtype=np.float32)
        preCalc[1, 15](y, yA, yB, numDataPoints)

        self.assertTrue(np.all(y == yA))
        self.assertTrue(np.all(y == yB))


if __name__ == '__main__':
    unittest.main()
