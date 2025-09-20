import numpy as np
from numba import vectorize
from numba.cuda.testing import skip_on_cudasim, CUDATestCase
import unittest


@skip_on_cudasim('ufunc API unsupported in the simulator')
class TestVectorizeComplex(CUDATestCase):
    def test_vectorize_complex(self):
        @vectorize(['complex128(complex128)'], target='cuda')
        def vcomp(a):
            return a * a + 1.

        A = np.arange(5, dtype=np.complex128)
        B = vcomp(A)
        self.assertTrue(np.allclose(A * A + 1., B))


if __name__ == '__main__':
    unittest.main()
