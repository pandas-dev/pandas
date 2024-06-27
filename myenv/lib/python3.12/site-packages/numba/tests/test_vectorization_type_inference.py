from numba import vectorize, jit, bool_, double, int_, float32, typeof, int8
import unittest
import numpy as np


def add(a, b):
    return a + b


def func(dtypeA, dtypeB):
    A = np.arange(10, dtype=dtypeA)
    B = np.arange(10, dtype=dtypeB)
    return typeof(vector_add(A, B))


class TestVectTypeInfer(unittest.TestCase):

    def test_type_inference(self):
        """This is testing numpy ufunc dispatch machinery
        """
        global vector_add
        vector_add = vectorize([
            bool_(double, int_),
            double(double, double),
            float32(double, float32),
        ])(add)

        def numba_type_equal(a, b):
            self.assertEqual(a.dtype, b.dtype)
            self.assertEqual(a.ndim, b.ndim)

        numba_type_equal(func(np.dtype(np.float64), np.dtype('i')), bool_[:])
        numba_type_equal(func(np.dtype(np.float64), np.dtype(np.float64)),
                         double[:])
        # This is because the double(double, double) matches first
        numba_type_equal(func(np.dtype(np.float64), np.dtype(np.float32)),
                         double[:])


if __name__ == '__main__':
    unittest.main()
