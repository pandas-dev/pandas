import pytest
import numpy as np
from numpy.testing import assert_allclose
from pytest import raises as assert_raises

from numpy import array, transpose, zeros_like, empty
from scipy.linalg import (cholesky, cholesky_banded, cho_solve_banded,
     cho_factor, cho_solve, LinAlgError)

from scipy.linalg._testutils import assert_no_overwrite
from .test_basic import parametrize_overwrite_arg


class TestCholesky:

    def test_simple(self):
        a = [[8, 2, 3], [2, 9, 3], [3, 3, 6]]
        c = cholesky(a)
        assert_allclose(c.T @ c, a, atol=1e-14)
        c = transpose(c)
        a = c @ c.T
        assert_allclose(cholesky(a, lower=1), c, atol=1e-14)

    def test_check_finite(self):
        a = [[8, 2, 3], [2, 9, 3], [3, 3, 6]]
        c = cholesky(a, check_finite=False)
        assert_allclose(c.T @ c, a, atol=1e-14)
        c = transpose(c)
        a = c @ c.T
        assert_allclose(cholesky(a, lower=1, check_finite=False), c, atol=1e-14)

    def test_simple_complex(self):
        m = array([[3+1j, 3+4j, 5], [0, 2+2j, 2+7j], [0, 0, 7+4j]])
        a = np.conj(m.T) @ m
        c = cholesky(a)
        a1 = np.conj(c.T) @ c
        assert_allclose(a, a1)
        c = transpose(c)
        a = c @ np.conj(c.T)
        assert_allclose(cholesky(a, lower=1), c, atol=1e-14)

    @pytest.mark.parametrize("n", [5, 10, 20])
    @pytest.mark.parametrize("dtype", [np.float32, np.float64,
                                        np.complex64, np.complex128])
    @pytest.mark.parametrize("order", ["F", "C"])
    def test_random(self, n, dtype, order):
        atol = 5e-4 if dtype in [np.float32, np.complex64] else 1e-12
        rng = np.random.default_rng(seed=12345)

        for k in range(2):
            m = rng.normal(size=(n, n)).astype(dtype)
            if np.issubdtype(dtype, np.complexfloating):
                m += 1j * rng.normal(size=(n, n)).astype(dtype)
            for i in range(n):
                m[i, i] = 20*(.1+m[i, i])

            a = np.asarray(np.conj(m.T) @ m, order=order)

            c = cholesky(a, lower=0)
            a1 = np.conj(c.T) @ c
            assert_allclose(a, a1, atol=atol)

            c = transpose(c)
            a = c @ np.conj(c.T)
            assert_allclose(cholesky(a, lower=1), c, atol=atol)

    @pytest.mark.xslow
    def test_int_overflow(self):
       # regression test for
       # https://github.com/scipy/scipy/issues/17436
       # the problem was an int overflow in zeroing out
       # the unused triangular part
       n = 47_000
       x = np.eye(n, dtype=np.float64, order='F')
       x[:4, :4] = np.array([[4, -2, 3, -1],
                             [-2, 4, -3, 1],
                             [3, -3, 5, 0],
                             [-1, 1, 0, 5]])

       cholesky(x, check_finite=False, overwrite_a=True)  # should not segfault

    @pytest.mark.parametrize('dt', [int, float, np.float32, complex, np.complex64])
    @pytest.mark.parametrize('dt_b', [int, float, np.float32, complex, np.complex64])
    def test_empty(self, dt, dt_b):
        a = empty((0, 0), dtype=dt)

        c = cholesky(a)
        assert c.shape == (0, 0)
        assert c.dtype == cholesky(np.eye(2, dtype=dt)).dtype

        c_and_lower = (c, True)
        b = np.asarray([], dtype=dt_b)
        x = cho_solve(c_and_lower, b)
        assert x.shape == (0,)
        assert x.dtype == cho_solve((np.eye(2, dtype=dt), True),
                                     np.ones(2, dtype=dt_b)).dtype

        b = empty((0, 0), dtype=dt_b)
        x = cho_solve(c_and_lower, b)
        assert x.shape == (0, 0)
        assert x.dtype == cho_solve((np.eye(2, dtype=dt), True),
                                     np.ones(2, dtype=dt_b)).dtype

        a1 = array([])
        a2 = array([[]])
        a3 = []
        a4 = [[]]
        for x in ([a1, a2, a3, a4]):
            assert_raises(ValueError, cholesky, x)

    @pytest.mark.parametrize("shape", [(3,), (2, 3), (3, 2)])
    def test_invalid_shape(self, shape):
        rng = np.random.default_rng(seed=12345)
        a = rng.normal(size=shape)

        with pytest.raises(ValueError, match="Expected a square matrix"):
            cholesky(a)

    @pytest.mark.parametrize("dtype", [np.float32, np.float64,
                             np.complex64, np.complex128])
    @pytest.mark.parametrize("order", ["C", "F"])
    def test_non_posdef(self, dtype, order):
        n = 10
        rng = np.random.default_rng(seed=12345)

        lmbd = np.diag(-np.abs(rng.normal(size=(n,))).astype(dtype))
        x = rng.normal(size=(n, n)).astype(dtype)
        if np.issubdtype(dtype, np.complexfloating):
            x += 1j * rng.normal(size=(n, n)).astype(dtype)

        a = np.asarray(x @ lmbd @ np.conj(x.T), order=order)

        with pytest.raises(LinAlgError, match="Internal potrf return info"):
            cholesky(a)

    @parametrize_overwrite_arg
    @pytest.mark.parametrize("dtype", [int, float])
    @pytest.mark.parametrize("order", ["C", "F"])
    def test_overwrite_args(self, overwrite_kw, dtype, order):
        rng = np.random.default_rng(seed=12345)
        n = 5

        x = rng.normal(size=(n, n))
        x = x.astype(dtype)
        a = x @ x.T
        a = a.astype(dtype, order=order)
        a_ref = np.copy(a)

        c = cholesky(a, **overwrite_kw)
        overwrite_a = overwrite_kw.get("overwrite_a", False)
        a_inplace = (overwrite_a and (dtype is not int)
                     and (a.flags["F_CONTIGUOUS"] or a.flags["C_CONTIGUOUS"]))

        assert np.shares_memory(a, c) == a_inplace
        assert np.all(a == a_ref) != a_inplace

    @pytest.mark.parametrize("lower", [True, False])
    def test_overwrite_cleaning(self, lower):
        # Since only half the array is copied over, this checks if the other
        # half is correctly put to 0 when the input array is passed in instead.
        rng = np.random.default_rng(seed=12345)
        n = 5

        x = rng.normal(size=(n, n))
        a = np.asfortranarray(x @ x.T)
        a_ref = np.copy(a)

        c = cholesky(a, lower=lower, overwrite_a=True)

        if lower:
            assert_allclose(np.triu(a, k=1), np.zeros_like(a), atol=1e-14)
            assert_allclose(c @ c.T, a_ref, atol=1e-14)
        else:
            assert_allclose(np.tril(a, k=-1), np.zeros_like(a), atol=1e-14)
            assert_allclose(c.T @ c, a_ref, atol=1e-14)

    @pytest.mark.parametrize("lower", [True, False])
    @pytest.mark.parametrize("order", ['C', 'F', 'noncontig_C', 'noncontig_F'])
    @pytest.mark.parametrize("overwrite_a", [True, False])
    @pytest.mark.parametrize("dtype", [np.float64, np.complex128])
    def test_nan_in_unused_triangle(self, lower, order, overwrite_a, dtype):
        # The unused triangle should never be read, so NaN there must not
        # affect the result.
        rng = np.random.default_rng(seed=12345)
        n = 5
        x = rng.normal(size=(n, n))
        if np.issubdtype(dtype, np.complexfloating):
            x = x + 1j * rng.normal(size=(n, n))
        a = (x @ x.conj().T + n * np.eye(n)).astype(dtype)
        c_ref = cholesky(a, lower=lower)

        # Prepare the array in the desired layout
        if order.startswith('noncontig'):
            # Unequal strides to probe strided path thoroughly
            big_order = 'F' if order.endswith('_F') else 'C'
            a_big = np.zeros((2*n, 3*n), dtype=dtype, order=big_order)
            a_big[::2, ::3] = a
            a_nan = a_big[::2, ::3]  # non-contiguous view
            assert not a_nan.flags['C_CONTIGUOUS']
            assert not a_nan.flags['F_CONTIGUOUS']
        else:
            a_nan = np.array(a, order=order)

        # Put NaN in the unused triangle
        if lower:
            a_nan[np.triu_indices(n, 1)] = np.nan
        else:
            a_nan[np.tril_indices(n, -1)] = np.nan

        c = cholesky(a_nan, lower=lower, overwrite_a=overwrite_a,
                     check_finite=False)
        assert_allclose(c, c_ref, atol=1e-14)
        # Verify reconstruction
        if lower:
            assert_allclose(c @ c.conj().T, a, atol=1e-13)
        else:
            assert_allclose(c.conj().T @ c, a, atol=1e-13)


class TestCholeskyBanded:
    """Tests for cholesky_banded() and cho_solve_banded."""

    def test_check_finite(self):
        # Symmetric positive definite banded matrix `a`
        a = array([[4.0, 1.0, 0.0, 0.0],
                   [1.0, 4.0, 0.5, 0.0],
                   [0.0, 0.5, 4.0, 0.2],
                   [0.0, 0.0, 0.2, 4.0]])
        # Banded storage form of `a`.
        ab = array([[-1.0, 1.0, 0.5, 0.2],
                    [4.0, 4.0, 4.0, 4.0]])
        c = cholesky_banded(ab, lower=False, check_finite=False)
        ufac = zeros_like(a)
        ufac[list(range(4)), list(range(4))] = c[-1]
        ufac[(0, 1, 2), (1, 2, 3)] = c[0, 1:]
        assert_allclose(a, ufac.T @ ufac, atol=1e-14)

        b = array([0.0, 0.5, 4.2, 4.2])
        x = cho_solve_banded((c, False), b, check_finite=False)
        assert_allclose(x, [0.0, 0.0, 1.0, 1.0], atol=1e-14)

    def test_upper_real(self):
        # Symmetric positive definite banded matrix `a`
        a = array([[4.0, 1.0, 0.0, 0.0],
                   [1.0, 4.0, 0.5, 0.0],
                   [0.0, 0.5, 4.0, 0.2],
                   [0.0, 0.0, 0.2, 4.0]])
        # Banded storage form of `a`.
        ab = array([[-1.0, 1.0, 0.5, 0.2],
                    [4.0, 4.0, 4.0, 4.0]])
        c = cholesky_banded(ab, lower=False)
        ufac = zeros_like(a)
        ufac[list(range(4)), list(range(4))] = c[-1]
        ufac[(0, 1, 2), (1, 2, 3)] = c[0, 1:]
        assert_allclose(a, ufac.T @ ufac, atol=1e-14)

        b = array([0.0, 0.5, 4.2, 4.2])
        x = cho_solve_banded((c, False), b)
        assert_allclose(x, [0.0, 0.0, 1.0, 1.0], atol=1e-14)

    def test_upper_complex(self):
        # Hermitian positive definite banded matrix `a`
        a = array([[4.0, 1.0, 0.0, 0.0],
                   [1.0, 4.0, 0.5, 0.0],
                   [0.0, 0.5, 4.0, -0.2j],
                   [0.0, 0.0, 0.2j, 4.0]])
        # Banded storage form of `a`.
        ab = array([[-1.0, 1.0, 0.5, -0.2j],
                    [4.0, 4.0, 4.0, 4.0]])
        c = cholesky_banded(ab, lower=False)
        ufac = zeros_like(a)
        ufac[list(range(4)), list(range(4))] = c[-1]
        ufac[(0, 1, 2), (1, 2, 3)] = c[0, 1:]
        assert_allclose(a, np.conj(ufac.T) @ ufac, atol=1e-14)

        b = array([0.0, 0.5, 4.0-0.2j, 0.2j + 4.0])
        x = cho_solve_banded((c, False), b)
        assert_allclose(x, [0.0, 0.0, 1.0, 1.0], atol=1e-14)

    def test_lower_real(self):
        # Symmetric positive definite banded matrix `a`
        a = array([[4.0, 1.0, 0.0, 0.0],
                   [1.0, 4.0, 0.5, 0.0],
                   [0.0, 0.5, 4.0, 0.2],
                   [0.0, 0.0, 0.2, 4.0]])
        # Banded storage form of `a`.
        ab = array([[4.0, 4.0, 4.0, 4.0],
                    [1.0, 0.5, 0.2, -1.0]])
        c = cholesky_banded(ab, lower=True)
        lfac = zeros_like(a)
        lfac[list(range(4)), list(range(4))] = c[0]
        lfac[(1, 2, 3), (0, 1, 2)] = c[1, :3]
        assert_allclose(a, lfac @ lfac.T, atol=1e-14)

        b = array([0.0, 0.5, 4.2, 4.2])
        x = cho_solve_banded((c, True), b)
        assert_allclose(x, [0.0, 0.0, 1.0, 1.0], atol=1e-14)

    def test_lower_complex(self):
        # Hermitian positive definite banded matrix `a`
        a = array([[4.0, 1.0, 0.0, 0.0],
                   [1.0, 4.0, 0.5, 0.0],
                   [0.0, 0.5, 4.0, -0.2j],
                   [0.0, 0.0, 0.2j, 4.0]])
        # Banded storage form of `a`.
        ab = array([[4.0, 4.0, 4.0, 4.0],
                    [1.0, 0.5, 0.2j, -1.0]])
        c = cholesky_banded(ab, lower=True)
        lfac = zeros_like(a)
        lfac[list(range(4)), list(range(4))] = c[0]
        lfac[(1, 2, 3), (0, 1, 2)] = c[1, :3]
        assert_allclose(a, lfac @ np.conj(lfac.T), atol=1e-14)

        b = array([0.0, 0.5j, 3.8j, 3.8])
        x = cho_solve_banded((c, True), b)
        assert_allclose(x, [0.0, 0.0, 1.0j, 1.0], atol=1e-14)

    @pytest.mark.parametrize('dt', [int, float, np.float32, complex, np.complex64])
    @pytest.mark.parametrize('dt_b', [int, float, np.float32, complex, np.complex64])
    def test_empty(self, dt, dt_b):
        ab = empty((0, 0), dtype=dt)

        cb = cholesky_banded(ab)
        assert cb.shape == (0, 0)

        m = cholesky_banded(np.array([[0, 0], [1, 1]], dtype=dt))
        assert cb.dtype == m.dtype

        cb_and_lower = (cb, True)
        b = np.asarray([], dtype=dt_b)
        x = cho_solve_banded(cb_and_lower, b)
        assert x.shape == (0,)

        dtype_nonempty = cho_solve_banded((m, True), np.ones(2, dtype=dt_b)).dtype
        assert x.dtype == dtype_nonempty

        b = empty((0, 0), dtype=dt_b)
        x = cho_solve_banded(cb_and_lower, b)
        assert x.shape == (0, 0)
        assert x.dtype == dtype_nonempty


class TestOverwrite:
    def test_cholesky(self):
        assert_no_overwrite(cholesky, [(3, 3)])

    def test_cho_factor(self):
        assert_no_overwrite(cho_factor, [(3, 3)])

    def test_cho_solve(self):
        x = array([[2, -1, 0], [-1, 2, -1], [0, -1, 2]])
        xcho = cho_factor(x)
        assert_no_overwrite(lambda b: cho_solve(xcho, b), [(3,)])

    def test_cholesky_banded(self):
        assert_no_overwrite(cholesky_banded, [(2, 3)])

    def test_cho_solve_banded(self):
        x = array([[0, -1, -1], [2, 2, 2]])
        xcho = cholesky_banded(x)
        assert_no_overwrite(lambda b: cho_solve_banded((xcho, False), b),
                            [(3,)])


class TestChoFactor:
    @pytest.mark.parametrize('dt', [int, float, np.float32, complex, np.complex64])
    def test_empty(self, dt):
        a = np.empty((0, 0), dtype=dt)
        x, lower = cho_factor(a)

        assert x.shape == (0, 0)

        xx, lower = cho_factor(np.eye(2, dtype=dt))
        assert x.dtype == xx.dtype

    @pytest.mark.parametrize("n", [5, 10, 20])
    @pytest.mark.parametrize("dtype", [np.float32, np.float64,
                                np.complex64, np.complex128])
    @pytest.mark.parametrize("order", ["C", "F"])
    @pytest.mark.parametrize("lower", [True, False])
    def test_random(self, n, dtype, order, lower):
        atol = 1e-5 if dtype in [np.float32, np.complex64] else 1e-12
        rng = np.random.default_rng(seed=12345)

        lmbd = np.diag(np.abs(rng.normal(size=(n,)).astype(dtype)))
        x = rng.normal(size=(n, n)).astype(dtype)
        if np.issubdtype(dtype, np.complexfloating):
            x += 1j * rng.normal(size=(n, n)).astype(dtype)

        a = np.asarray(x @ lmbd @ np.conj(x.T), order=order)
        c = cholesky(a, lower=lower)
        c_f, lower_out = cho_factor(a, lower=lower)

        assert lower_out == lower

        # Theoretically exactly equal
        # For `cho_factor` the other half of the matrix is irrelevant.
        if lower:
            l = np.tril(c_f)
            assert_allclose(c, l, atol=atol)
            assert_allclose(a, l @ np.conj(l.T), atol=atol)
        else:
            u = np.triu(c_f)
            assert_allclose(c, u, atol=atol)
            assert_allclose(a, np.conj(u.T) @ u, atol=atol)

    @pytest.mark.parametrize("shape", [(3,), (3, 2), (2, 3)])
    def test_invalid_shape(self, shape):
        rng = np.random.default_rng(seed=12345)
        a = rng.normal(size=shape)

        with pytest.raises(ValueError, match="Expected a square matrix or batch"):
            cho_factor(a)

    @pytest.mark.parametrize("dtype", [np.float32, np.float64,
                             np.complex64, np.complex128])
    @pytest.mark.parametrize("order", ["C", "F"])
    def test_non_posdef(self, dtype, order):
        n = 10
        rng = np.random.default_rng(seed=12345)

        lmbd = np.diag(-np.abs(rng.normal(size=(n,))).astype(dtype))
        x = rng.normal(size=(n, n)).astype(dtype)
        if np.issubdtype(dtype, np.complexfloating):
            x += 1j * rng.normal(size=(n, n)).astype(dtype)

        a = np.asarray(x @ lmbd @ np.conj(x.T), order=order)

        with pytest.raises(LinAlgError, match="Internal potrf return info"):
            cho_factor(a)

    @parametrize_overwrite_arg
    @pytest.mark.parametrize("dtype", [int, float])
    @pytest.mark.parametrize("order", ["C", "F"])
    @pytest.mark.parametrize("lower", [True, False])
    def test_overwrite_args(self, overwrite_kw, dtype, order, lower):
        rng = np.random.default_rng(seed=12345)
        n = 5

        x = rng.normal(size=(n, n))
        x = x.astype(dtype)
        a = x @ x.T
        a = a.astype(dtype, order=order)
        a_ref = np.copy(a)

        c, lower = cho_factor(a, **overwrite_kw, lower=lower)
        overwrite_a = overwrite_kw.get("overwrite_a", False)
        a_inplace = (overwrite_a and (dtype is not int)
                     and (a.flags["F_CONTIGUOUS"] or a.flags["C_CONTIGUOUS"]))

        assert np.shares_memory(a, c) == a_inplace
        assert np.all(a == a_ref) != a_inplace
