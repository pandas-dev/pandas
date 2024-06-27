import math
import numpy as np
from numba import cuda, float64, int8, int32, void
from numba.cuda.testing import unittest, CUDATestCase


def cu_mat_power(A, power, power_A):
    y, x = cuda.grid(2)

    m, n = power_A.shape
    if x >= n or y >= m:
        return

    power_A[y, x] = math.pow(A[y, x], int32(power))


def cu_mat_power_binop(A, power, power_A):
    y, x = cuda.grid(2)

    m, n = power_A.shape
    if x >= n or y >= m:
        return

    power_A[y, x] = A[y, x] ** power


def vec_pow(r, x, y):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = pow(x[i], y[i])


def vec_pow_binop(r, x, y):
    i = cuda.grid(1)

    if i < len(r):
        r[i] = x[i] ** y[i]


def vec_pow_inplace_binop(r, x):
    i = cuda.grid(1)

    if i < len(r):
        r[i] **= x[i]


def random_complex(N):
    np.random.seed(123)
    return (np.random.random(1) + np.random.random(1) * 1j)


class TestCudaPowi(CUDATestCase):
    def test_powi(self):
        dec = cuda.jit(void(float64[:, :], int8, float64[:, :]))
        kernel = dec(cu_mat_power)

        power = 2
        A = np.arange(10, dtype=np.float64).reshape(2, 5)
        Aout = np.empty_like(A)
        kernel[1, A.shape](A, power, Aout)
        self.assertTrue(np.allclose(Aout, A ** power))

    def test_powi_binop(self):
        dec = cuda.jit(void(float64[:, :], int8, float64[:, :]))
        kernel = dec(cu_mat_power_binop)

        power = 2
        A = np.arange(10, dtype=np.float64).reshape(2, 5)
        Aout = np.empty_like(A)
        kernel[1, A.shape](A, power, Aout)
        self.assertTrue(np.allclose(Aout, A ** power))

    # Relative tolerance kwarg is provided because 1.0e-7 (the default for
    # assert_allclose) is a bit tight for single precision.
    def _test_cpow(self, dtype, func, rtol=1.0e-7):
        N = 32
        x = random_complex(N).astype(dtype)
        y = random_complex(N).astype(dtype)
        r = np.zeros_like(x)

        cfunc = cuda.jit(func)
        cfunc[1, N](r, x, y)
        np.testing.assert_allclose(r, x ** y, rtol=rtol)

        # Checks special cases
        x = np.asarray([0.0j, 1.0j], dtype=dtype)
        y = np.asarray([0.0j, 1.0], dtype=dtype)
        r = np.zeros_like(x)

        cfunc[1, 2](r, x, y)
        np.testing.assert_allclose(r, x ** y, rtol=rtol)

    def test_cpow_complex64_pow(self):
        self._test_cpow(np.complex64, vec_pow, rtol=3.0e-7)

    def test_cpow_complex64_binop(self):
        self._test_cpow(np.complex64, vec_pow_binop, rtol=3.0e-7)

    def test_cpow_complex128_pow(self):
        self._test_cpow(np.complex128, vec_pow)

    def test_cpow_complex128_binop(self):
        self._test_cpow(np.complex128, vec_pow_binop)

    def _test_cpow_inplace_binop(self, dtype, rtol=1.0e-7):
        N = 32
        x = random_complex(N).astype(dtype)
        y = random_complex(N).astype(dtype)
        r = x ** y

        cfunc = cuda.jit(vec_pow_inplace_binop)
        cfunc[1, N](x, y)
        np.testing.assert_allclose(x, r, rtol=rtol)

    def test_cpow_complex64_inplace_binop(self):
        self._test_cpow_inplace_binop(np.complex64, rtol=3.0e-7)

    def test_cpow_complex128_inplace_binop(self):
        self._test_cpow_inplace_binop(np.complex128, rtol=3.0e-7)


if __name__ == '__main__':
    unittest.main()
