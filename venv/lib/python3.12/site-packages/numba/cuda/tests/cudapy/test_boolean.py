import numpy as np
from numba.cuda.testing import unittest, CUDATestCase
from numba import cuda


def boolean_func(A, vertial):
    if vertial:
        A[0] = 123
    else:
        A[0] = 321


class TestCudaBoolean(CUDATestCase):
    def test_boolean(self):
        func = cuda.jit('void(float64[:], bool_)')(boolean_func)
        A = np.array([0], dtype='float64')
        func[1, 1](A, True)
        self.assertTrue(A[0] == 123)
        func[1, 1](A, False)
        self.assertTrue(A[0] == 321)


if __name__ == '__main__':
    unittest.main()
