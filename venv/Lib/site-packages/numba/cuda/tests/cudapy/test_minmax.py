import numpy as np

from numba import cuda, float64
from numba.cuda.testing import unittest, CUDATestCase, skip_on_cudasim


def builtin_max(A, B, C):
    i = cuda.grid(1)

    if i >= len(C):
        return

    C[i] = float64(max(A[i], B[i]))


def builtin_min(A, B, C):
    i = cuda.grid(1)

    if i >= len(C):
        return

    C[i] = float64(min(A[i], B[i]))


@skip_on_cudasim('Tests PTX emission')
class TestCudaMinMax(CUDATestCase):
    def _run(
            self,
            kernel,
            numpy_equivalent,
            ptx_instruction,
            dtype_left,
            dtype_right,
            n=5):
        kernel = cuda.jit(kernel)

        c = np.zeros(n, dtype=np.float64)
        a = np.arange(n, dtype=dtype_left) + .5
        b = np.full(n, fill_value=2, dtype=dtype_right)

        kernel[1, c.shape](a, b, c)
        np.testing.assert_allclose(c, numpy_equivalent(a, b))

        ptx = next(p for p in kernel.inspect_asm().values())
        self.assertIn(ptx_instruction, ptx)

    def test_max_f8f8(self):
        self._run(
            builtin_max,
            np.maximum,
            'max.f64',
            np.float64,
            np.float64)

    def test_max_f4f8(self):
        self._run(
            builtin_max,
            np.maximum,
            'max.f64',
            np.float32,
            np.float64)

    def test_max_f8f4(self):
        self._run(
            builtin_max,
            np.maximum,
            'max.f64',
            np.float64,
            np.float32)

    def test_max_f4f4(self):
        self._run(
            builtin_max,
            np.maximum,
            'max.f32',
            np.float32,
            np.float32)

    def test_min_f8f8(self):
        self._run(
            builtin_min,
            np.minimum,
            'min.f64',
            np.float64,
            np.float64)

    def test_min_f4f8(self):
        self._run(
            builtin_min,
            np.minimum,
            'min.f64',
            np.float32,
            np.float64)

    def test_min_f8f4(self):
        self._run(
            builtin_min,
            np.minimum,
            'min.f64',
            np.float64,
            np.float32)

    def test_min_f4f4(self):
        self._run(
            builtin_min,
            np.minimum,
            'min.f32',
            np.float32,
            np.float32)


if __name__ == '__main__':
    unittest.main()
