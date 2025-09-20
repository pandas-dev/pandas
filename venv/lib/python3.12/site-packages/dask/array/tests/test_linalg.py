from __future__ import annotations

import sys

import pytest

pytest.importorskip("numpy")
pytest.importorskip("scipy")

import numpy as np
import scipy.linalg
from packaging.version import Version

import dask.array as da
from dask.array.linalg import qr, sfqr, svd, svd_compressed, tsqr
from dask.array.numpy_compat import _np_version
from dask.array.utils import assert_eq, same_keys, svd_flip


@pytest.mark.parametrize(
    "m,n,chunks,error_type",
    [
        (20, 10, 10, None),  # tall-skinny regular blocks
        (20, 10, (3, 10), None),  # tall-skinny regular fat layers
        (20, 10, ((8, 4, 8), 10), None),  # tall-skinny irregular fat layers
        (40, 10, ((15, 5, 5, 8, 7), 10), None),  # tall-skinny non-uniform chunks (why?)
        (128, 2, (16, 2), None),  # tall-skinny regular thin layers; recursion_depth=1
        (
            129,
            2,
            (16, 2),
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=2 --> 17x2
        (
            130,
            2,
            (16, 2),
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=2 --> 18x2 next
        (
            131,
            2,
            (16, 2),
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=2 --> 18x2 next
        (300, 10, (40, 10), None),  # tall-skinny regular thin layers; recursion_depth=2
        (300, 10, (30, 10), None),  # tall-skinny regular thin layers; recursion_depth=3
        (300, 10, (20, 10), None),  # tall-skinny regular thin layers; recursion_depth=4
        (10, 5, 10, None),  # single block tall
        (5, 10, 10, None),  # single block short
        (10, 10, 10, None),  # single block square
        (10, 40, (10, 10), ValueError),  # short-fat regular blocks
        (10, 40, (10, 15), ValueError),  # short-fat irregular blocks
        (
            10,
            40,
            (10, (15, 5, 5, 8, 7)),
            ValueError,
        ),  # short-fat non-uniform chunks (why?)
        (20, 20, 10, ValueError),  # 2x2 regular blocks
    ],
)
def test_tsqr(m, n, chunks, error_type):
    mat = np.random.default_rng().random((m, n))
    data = da.from_array(mat, chunks=chunks, name="A")

    # qr
    m_q = m
    n_q = min(m, n)
    m_r = n_q
    n_r = n

    # svd
    m_u = m
    n_u = min(m, n)
    n_s = n_q
    m_vh = n_q
    n_vh = n
    d_vh = max(m_vh, n_vh)  # full matrix returned

    if error_type is None:
        # test QR
        q, r = tsqr(data)
        assert_eq((m_q, n_q), q.shape)  # shape check
        assert_eq((m_r, n_r), r.shape)  # shape check
        assert_eq(mat, da.dot(q, r))  # accuracy check
        assert_eq(np.eye(n_q, n_q), da.dot(q.T, q))  # q must be orthonormal
        assert_eq(r, da.triu(r.rechunk(r.shape[0])))  # r must be upper triangular

        # test SVD
        u, s, vh = tsqr(data, compute_svd=True)
        s_exact = np.linalg.svd(mat)[1]
        assert_eq(s, s_exact)  # s must contain the singular values
        assert_eq((m_u, n_u), u.shape)  # shape check
        assert_eq((n_s,), s.shape)  # shape check
        assert_eq((d_vh, d_vh), vh.shape)  # shape check
        assert_eq(np.eye(n_u, n_u), da.dot(u.T, u))  # u must be orthonormal
        assert_eq(np.eye(d_vh, d_vh), da.dot(vh, vh.T))  # vh must be orthonormal
        assert_eq(mat, da.dot(da.dot(u, da.diag(s)), vh[:n_q]))  # accuracy check
    else:
        with pytest.raises(error_type):
            q, r = tsqr(data)
        with pytest.raises(error_type):
            u, s, vh = tsqr(data, compute_svd=True)


@pytest.mark.parametrize(
    "m_min,n_max,chunks,vary_rows,vary_cols,error_type",
    [
        (10, 5, (10, 5), True, False, None),  # single block tall
        (10, 5, (10, 5), False, True, None),  # single block tall
        (10, 5, (10, 5), True, True, None),  # single block tall
        (40, 5, (10, 5), True, False, None),  # multiple blocks tall
        (40, 5, (10, 5), False, True, None),  # multiple blocks tall
        (40, 5, (10, 5), True, True, None),  # multiple blocks tall
        (
            300,
            10,
            (40, 10),
            True,
            False,
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=2
        (
            300,
            10,
            (30, 10),
            True,
            False,
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=3
        (
            300,
            10,
            (20, 10),
            True,
            False,
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=4
        (
            300,
            10,
            (40, 10),
            False,
            True,
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=2
        (
            300,
            10,
            (30, 10),
            False,
            True,
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=3
        (
            300,
            10,
            (20, 10),
            False,
            True,
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=4
        (
            300,
            10,
            (40, 10),
            True,
            True,
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=2
        (
            300,
            10,
            (30, 10),
            True,
            True,
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=3
        (
            300,
            10,
            (20, 10),
            True,
            True,
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=4
    ],
)
def test_tsqr_uncertain(m_min, n_max, chunks, vary_rows, vary_cols, error_type):
    mat = np.random.default_rng().random((m_min * 2, n_max))
    m, n = m_min * 2, n_max
    mat[0:m_min, 0] += 1
    _c0 = mat[:, 0]
    _r0 = mat[0, :]
    c0 = da.from_array(_c0, chunks=m_min, name="c")
    r0 = da.from_array(_r0, chunks=n_max, name="r")
    data = da.from_array(mat, chunks=chunks, name="A")
    if vary_rows:
        data = data[c0 > 0.5, :]
        mat = mat[_c0 > 0.5, :]
        m = mat.shape[0]
    if vary_cols:
        data = data[:, r0 > 0.5]
        mat = mat[:, _r0 > 0.5]
        n = mat.shape[1]

    # qr
    m_q = m
    n_q = min(m, n)
    m_r = n_q
    n_r = n

    # svd
    m_u = m
    n_u = min(m, n)
    n_s = n_q
    m_vh = n_q
    n_vh = n
    d_vh = max(m_vh, n_vh)  # full matrix returned

    if error_type is None:
        # test QR
        q, r = tsqr(data)
        q = q.compute()  # because uncertainty
        r = r.compute()
        assert_eq((m_q, n_q), q.shape)  # shape check
        assert_eq((m_r, n_r), r.shape)  # shape check
        assert_eq(mat, np.dot(q, r))  # accuracy check
        assert_eq(np.eye(n_q, n_q), np.dot(q.T, q))  # q must be orthonormal
        assert_eq(r, np.triu(r))  # r must be upper triangular

        # test SVD
        u, s, vh = tsqr(data, compute_svd=True)
        u = u.compute()  # because uncertainty
        s = s.compute()
        vh = vh.compute()
        s_exact = np.linalg.svd(mat)[1]
        assert_eq(s, s_exact)  # s must contain the singular values
        assert_eq((m_u, n_u), u.shape)  # shape check
        assert_eq((n_s,), s.shape)  # shape check
        assert_eq((d_vh, d_vh), vh.shape)  # shape check
        assert_eq(np.eye(n_u, n_u), np.dot(u.T, u))  # u must be orthonormal
        assert_eq(np.eye(d_vh, d_vh), np.dot(vh, vh.T))  # vh must be orthonormal
        assert_eq(mat, np.dot(np.dot(u, np.diag(s)), vh[:n_q]))  # accuracy check
    else:
        with pytest.raises(error_type):
            q, r = tsqr(data)
        with pytest.raises(error_type):
            u, s, vh = tsqr(data, compute_svd=True)


def test_tsqr_zero_height_chunks():
    m_q = 10
    n_q = 5
    m_r = 5
    n_r = 5

    # certainty
    mat = np.random.default_rng().random((10, 5))
    x = da.from_array(mat, chunks=((4, 0, 1, 0, 5), (5,)))
    q, r = da.linalg.qr(x)
    assert_eq((m_q, n_q), q.shape)  # shape check
    assert_eq((m_r, n_r), r.shape)  # shape check
    assert_eq(mat, da.dot(q, r))  # accuracy check
    assert_eq(np.eye(n_q, n_q), da.dot(q.T, q))  # q must be orthonormal
    assert_eq(r, da.triu(r.rechunk(r.shape[0])))  # r must be upper triangular

    # uncertainty
    mat2 = np.vstack([mat, -(np.ones((10, 5)))])
    v2 = mat2[:, 0]
    x2 = da.from_array(mat2, chunks=5)
    c = da.from_array(v2, chunks=5)
    x = x2[c >= 0, :]  # remove the ones added above to yield mat
    q, r = da.linalg.qr(x)
    q = q.compute()  # because uncertainty
    r = r.compute()
    assert_eq((m_q, n_q), q.shape)  # shape check
    assert_eq((m_r, n_r), r.shape)  # shape check
    assert_eq(mat, np.dot(q, r))  # accuracy check
    assert_eq(np.eye(n_q, n_q), np.dot(q.T, q))  # q must be orthonormal
    assert_eq(r, np.triu(r))  # r must be upper triangular


@pytest.mark.parametrize(
    "m,n,chunks,error_type",
    [
        (20, 10, 10, ValueError),  # tall-skinny regular blocks
        (20, 10, (3, 10), ValueError),  # tall-skinny regular fat layers
        (20, 10, ((8, 4, 8), 10), ValueError),  # tall-skinny irregular fat layers
        (
            40,
            10,
            ((15, 5, 5, 8, 7), 10),
            ValueError,
        ),  # tall-skinny non-uniform chunks (why?)
        (
            128,
            2,
            (16, 2),
            ValueError,
        ),  # tall-skinny regular thin layers; recursion_depth=1
        (
            129,
            2,
            (16, 2),
            ValueError,
        ),  # tall-skinny regular thin layers; recursion_depth=2 --> 17x2
        (
            130,
            2,
            (16, 2),
            ValueError,
        ),  # tall-skinny regular thin layers; recursion_depth=2 --> 18x2 next
        (
            131,
            2,
            (16, 2),
            ValueError,
        ),  # tall-skinny regular thin layers; recursion_depth=2 --> 18x2 next
        (
            300,
            10,
            (40, 10),
            ValueError,
        ),  # tall-skinny regular thin layers; recursion_depth=2
        (
            300,
            10,
            (30, 10),
            ValueError,
        ),  # tall-skinny regular thin layers; recursion_depth=3
        (
            300,
            10,
            (20, 10),
            ValueError,
        ),  # tall-skinny regular thin layers; recursion_depth=4
        (10, 5, 10, None),  # single block tall
        (5, 10, 10, None),  # single block short
        (10, 10, 10, None),  # single block square
        (10, 40, (10, 10), None),  # short-fat regular blocks
        (10, 40, (10, 15), None),  # short-fat irregular blocks
        (10, 40, (10, (15, 5, 5, 8, 7)), None),  # short-fat non-uniform chunks (why?)
        (20, 20, 10, ValueError),  # 2x2 regular blocks
    ],
)
def test_sfqr(m, n, chunks, error_type):
    mat = np.random.default_rng().random((m, n))
    data = da.from_array(mat, chunks=chunks, name="A")
    m_q = m
    n_q = min(m, n)
    m_r = n_q
    n_r = n
    m_qtq = n_q

    if error_type is None:
        q, r = sfqr(data)
        assert_eq((m_q, n_q), q.shape)  # shape check
        assert_eq((m_r, n_r), r.shape)  # shape check
        assert_eq(mat, da.dot(q, r))  # accuracy check
        assert_eq(np.eye(m_qtq, m_qtq), da.dot(q.T, q))  # q must be orthonormal
        assert_eq(r, da.triu(r.rechunk(r.shape[0])))  # r must be upper triangular
    else:
        with pytest.raises(error_type):
            q, r = sfqr(data)


@pytest.mark.parametrize(
    "m,n,chunks,error_type",
    [
        (20, 10, 10, None),  # tall-skinny regular blocks
        (20, 10, (3, 10), None),  # tall-skinny regular fat layers
        (20, 10, ((8, 4, 8), 10), None),  # tall-skinny irregular fat layers
        (40, 10, ((15, 5, 5, 8, 7), 10), None),  # tall-skinny non-uniform chunks (why?)
        (128, 2, (16, 2), None),  # tall-skinny regular thin layers; recursion_depth=1
        (
            129,
            2,
            (16, 2),
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=2 --> 17x2
        (
            130,
            2,
            (16, 2),
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=2 --> 18x2 next
        (
            131,
            2,
            (16, 2),
            None,
        ),  # tall-skinny regular thin layers; recursion_depth=2 --> 18x2 next
        (300, 10, (40, 10), None),  # tall-skinny regular thin layers; recursion_depth=2
        (300, 10, (30, 10), None),  # tall-skinny regular thin layers; recursion_depth=3
        (300, 10, (20, 10), None),  # tall-skinny regular thin layers; recursion_depth=4
        (10, 5, 10, None),  # single block tall
        (5, 10, 10, None),  # single block short
        (10, 10, 10, None),  # single block square
        (10, 40, (10, 10), None),  # short-fat regular blocks
        (10, 40, (10, 15), None),  # short-fat irregular blocks
        (10, 40, (10, (15, 5, 5, 8, 7)), None),  # short-fat non-uniform chunks (why?)
        (20, 20, 10, NotImplementedError),  # 2x2 regular blocks
    ],
)
def test_qr(m, n, chunks, error_type):
    mat = np.random.default_rng().random((m, n))
    data = da.from_array(mat, chunks=chunks, name="A")
    m_q = m
    n_q = min(m, n)
    m_r = n_q
    n_r = n
    m_qtq = n_q

    if error_type is None:
        q, r = qr(data)
        assert_eq((m_q, n_q), q.shape)  # shape check
        assert_eq((m_r, n_r), r.shape)  # shape check
        assert_eq(mat, da.dot(q, r))  # accuracy check
        assert_eq(np.eye(m_qtq, m_qtq), da.dot(q.T, q))  # q must be orthonormal
        assert_eq(r, da.triu(r.rechunk(r.shape[0])))  # r must be upper triangular
    else:
        with pytest.raises(error_type):
            q, r = qr(data)


def test_linalg_consistent_names():
    m, n = 20, 10
    mat = np.random.default_rng().random((m, n))
    data = da.from_array(mat, chunks=(10, n), name="A")

    q1, r1 = qr(data)
    q2, r2 = qr(data)
    assert same_keys(q1, q2)
    assert same_keys(r1, r2)

    u1, s1, v1 = svd(data)
    u2, s2, v2 = svd(data)
    assert same_keys(u1, u2)
    assert same_keys(s1, s2)
    assert same_keys(v1, v2)


@pytest.mark.parametrize("m,n", [(10, 20), (15, 15), (20, 10)])
def test_dask_svd_self_consistent(m, n):
    a = np.random.default_rng().random((m, n))
    d_a = da.from_array(a, chunks=(3, n), name="A")

    d_u, d_s, d_vt = da.linalg.svd(d_a)
    u, s, vt = da.compute(d_u, d_s, d_vt)

    for d_e, e in zip([d_u, d_s, d_vt], [u, s, vt]):
        assert d_e.shape == e.shape
        assert d_e.dtype == e.dtype


@pytest.mark.parametrize("iterator", ["power", "QR"])
def test_svd_compressed_compute(iterator):
    x = da.ones((100, 100), chunks=(10, 10))
    u, s, v = da.linalg.svd_compressed(
        x, k=2, iterator=iterator, n_power_iter=1, compute=True, seed=123
    )
    uu, ss, vv = da.linalg.svd_compressed(
        x, k=2, iterator=iterator, n_power_iter=1, seed=123
    )

    assert len(v.dask) < len(vv.dask)
    assert_eq(v, vv)


@pytest.mark.parametrize("iterator", [("power", 2), ("QR", 2)])
def test_svd_compressed(iterator):
    m, n = 100, 50
    r = 5
    a = da.random.default_rng().random((m, n), chunks=(m, n))

    # calculate approximation and true singular values
    u, s, vt = svd_compressed(
        a, 2 * r, iterator=iterator[0], n_power_iter=iterator[1], seed=4321
    )  # worst case
    s_true = scipy.linalg.svd(a.compute(), compute_uv=False)

    # compute the difference with original matrix
    norm = scipy.linalg.norm((a - (u[:, :r] * s[:r]) @ vt[:r, :]).compute(), 2)

    # ||a-a_hat||_2 <= (1+tol)s_{k+1}: based on eq. 1.10/1.11:
    # Halko, Nathan, Per-Gunnar Martinsson, and Joel A. Tropp.
    # "Finding structure with randomness: Probabilistic algorithms for constructing
    # approximate matrix decompositions." SIAM review 53.2 (2011): 217-288.
    frac = norm / s_true[r + 1] - 1
    # Tolerance determined via simulation to be slightly above max norm of difference matrix in 10k samples.
    # See https://github.com/dask/dask/pull/6799#issuecomment-726631175 for more details.
    tol = 0.4
    assert frac < tol

    assert_eq(np.eye(r, r), da.dot(u[:, :r].T, u[:, :r]))  # u must be orthonormal
    assert_eq(np.eye(r, r), da.dot(vt[:r, :], vt[:r, :].T))  # v must be orthonormal


@pytest.mark.parametrize(
    "input_dtype, output_dtype", [(np.float32, np.float32), (np.float64, np.float64)]
)
def test_svd_compressed_dtype_preservation(input_dtype, output_dtype):
    x = da.random.default_rng().random((50, 50), chunks=(50, 50)).astype(input_dtype)
    u, s, vt = svd_compressed(x, 1, seed=4321)
    assert u.dtype == s.dtype == vt.dtype == output_dtype


@pytest.mark.parametrize("chunks", [(10, 50), (50, 10), (-1, -1)])
@pytest.mark.parametrize("dtype", [np.float32, np.float64])
def test_svd_dtype_preservation(chunks, dtype):
    x = da.random.default_rng().random((50, 50), chunks=chunks).astype(dtype)
    u, s, v = svd(x)
    assert u.dtype == s.dtype == v.dtype == dtype


def test_svd_compressed_deterministic():
    m, n = 30, 25
    x = da.random.default_rng(1234).random(size=(m, n), chunks=(5, 5))
    u, s, vt = svd_compressed(x, 3, seed=1234)
    u2, s2, vt2 = svd_compressed(x, 3, seed=1234)

    assert all(da.compute((u == u2).all(), (s == s2).all(), (vt == vt2).all()))


@pytest.mark.parametrize("m", [5, 10, 15, 20])
@pytest.mark.parametrize("n", [5, 10, 15, 20])
@pytest.mark.parametrize("k", [5])
@pytest.mark.parametrize("chunks", [(5, 10), (10, 5)])
def test_svd_compressed_shapes(m, n, k, chunks):
    x = da.random.default_rng().random(size=(m, n), chunks=chunks)
    u, s, v = svd_compressed(x, k, n_power_iter=1, compute=True, seed=1)
    u, s, v = da.compute(u, s, v)
    r = min(m, n, k)
    assert u.shape == (m, r)
    assert s.shape == (r,)
    assert v.shape == (r, n)


def _check_lu_result(p, l, u, A):
    assert np.allclose(p.dot(l).dot(u), A)

    # check triangulars
    assert_eq(l, da.tril(l), check_graph=False)
    assert_eq(u, da.triu(u), check_graph=False)


def test_lu_1():
    A1 = np.array([[7, 3, -1, 2], [3, 8, 1, -4], [-1, 1, 4, -1], [2, -4, -1, 6]])

    A2 = np.array(
        [
            [7, 0, 0, 0, 0, 0],
            [0, 8, 0, 0, 0, 0],
            [0, 0, 4, 0, 0, 0],
            [0, 0, 0, 6, 0, 0],
            [0, 0, 0, 0, 3, 0],
            [0, 0, 0, 0, 0, 5],
        ]
    )
    # without shuffle
    for A, chunk in zip([A1, A2], [2, 2]):
        dA = da.from_array(A, chunks=(chunk, chunk))
        p, l, u = scipy.linalg.lu(A)
        dp, dl, du = da.linalg.lu(dA)
        assert_eq(p, dp, check_graph=False)
        assert_eq(l, dl, check_graph=False)
        assert_eq(u, du, check_graph=False)
        _check_lu_result(dp, dl, du, A)

    A3 = np.array(
        [
            [7, 3, 2, 1, 4, 1],
            [7, 11, 5, 2, 5, 2],
            [21, 25, 16, 10, 16, 5],
            [21, 41, 18, 13, 16, 11],
            [14, 46, 23, 24, 21, 22],
            [0, 56, 29, 17, 14, 8],
        ]
    )

    # with shuffle
    for A, chunk in zip([A3], [2]):
        dA = da.from_array(A, chunks=(chunk, chunk))
        p, l, u = scipy.linalg.lu(A)
        dp, dl, du = da.linalg.lu(dA)
        _check_lu_result(dp, dl, du, A)


@pytest.mark.slow
@pytest.mark.parametrize("size", [10, 20, 30, 50])
@pytest.mark.filterwarnings("ignore:Increasing:dask.array.core.PerformanceWarning")
def test_lu_2(size):
    rng = np.random.default_rng(10)
    A = rng.integers(0, 10, (size, size))

    dA = da.from_array(A, chunks=(5, 5))
    dp, dl, du = da.linalg.lu(dA)
    _check_lu_result(dp, dl, du, A)


@pytest.mark.slow
@pytest.mark.parametrize("size", [50, 100, 200])
def test_lu_3(size):
    rng = np.random.default_rng(10)
    A = rng.integers(0, 10, (size, size))

    dA = da.from_array(A, chunks=(25, 25))
    dp, dl, du = da.linalg.lu(dA)
    _check_lu_result(dp, dl, du, A)


def test_lu_errors():
    rng = np.random.default_rng()

    A = rng.integers(0, 11, (10, 10, 10))
    dA = da.from_array(A, chunks=(5, 5, 5))
    pytest.raises(ValueError, lambda: da.linalg.lu(dA))

    A = rng.integers(0, 11, (10, 8))
    dA = da.from_array(A, chunks=(5, 4))
    pytest.raises(ValueError, lambda: da.linalg.lu(dA))

    A = rng.integers(0, 11, (20, 20))
    dA = da.from_array(A, chunks=(5, 4))
    pytest.raises(ValueError, lambda: da.linalg.lu(dA))


@pytest.mark.parametrize(("shape", "chunk"), [(20, 10), (50, 10), (70, 20)])
def test_solve_triangular_vector(shape, chunk):
    rng = np.random.default_rng(1)

    A = rng.integers(1, 11, (shape, shape))
    b = rng.integers(1, 11, shape)

    # upper
    Au = np.triu(A)
    dAu = da.from_array(Au, (chunk, chunk))
    db = da.from_array(b, chunk)
    res = da.linalg.solve_triangular(dAu, db)
    assert_eq(res, scipy.linalg.solve_triangular(Au, b))
    assert_eq(dAu.dot(res), b.astype(float), rtol=1e-4)

    # lower
    Al = np.tril(A)
    dAl = da.from_array(Al, (chunk, chunk))
    db = da.from_array(b, chunk)
    res = da.linalg.solve_triangular(dAl, db, lower=True)
    assert_eq(res, scipy.linalg.solve_triangular(Al, b, lower=True))
    assert_eq(dAl.dot(res), b.astype(float))


@pytest.mark.parametrize(("shape", "chunk"), [(20, 10), (50, 10), (50, 20)])
def test_solve_triangular_matrix(shape, chunk):
    rng = np.random.default_rng(1)

    A = rng.integers(1, 10, (shape, shape))
    b = rng.integers(1, 10, (shape, 5))

    # upper
    Au = np.triu(A)
    dAu = da.from_array(Au, (chunk, chunk))
    db = da.from_array(b, (chunk, 5))
    res = da.linalg.solve_triangular(dAu, db)
    assert_eq(res, scipy.linalg.solve_triangular(Au, b))
    assert_eq(dAu.dot(res), b.astype(float))

    # lower
    Al = np.tril(A)
    dAl = da.from_array(Al, (chunk, chunk))
    db = da.from_array(b, (chunk, 5))
    res = da.linalg.solve_triangular(dAl, db, lower=True)
    assert_eq(res, scipy.linalg.solve_triangular(Al, b, lower=True))
    assert_eq(dAl.dot(res), b.astype(float))


@pytest.mark.parametrize(("shape", "chunk"), [(20, 10), (50, 10), (50, 20)])
def test_solve_triangular_matrix2(shape, chunk):
    rng = np.random.default_rng(1)

    A = rng.integers(1, 10, (shape, shape))
    b = rng.integers(1, 10, (shape, shape))

    # upper
    Au = np.triu(A)
    dAu = da.from_array(Au, (chunk, chunk))
    db = da.from_array(b, (chunk, chunk))
    res = da.linalg.solve_triangular(dAu, db)
    assert_eq(res, scipy.linalg.solve_triangular(Au, b))
    assert_eq(dAu.dot(res), b.astype(float))

    # lower
    Al = np.tril(A)
    dAl = da.from_array(Al, (chunk, chunk))
    db = da.from_array(b, (chunk, chunk))
    res = da.linalg.solve_triangular(dAl, db, lower=True)
    assert_eq(res, scipy.linalg.solve_triangular(Al, b, lower=True))
    assert_eq(dAl.dot(res), b.astype(float))


def test_solve_triangular_errors():
    A = np.random.default_rng().integers(0, 10, (10, 10, 10))
    b = np.random.default_rng().integers(1, 10, 10)
    dA = da.from_array(A, chunks=(5, 5, 5))
    db = da.from_array(b, chunks=5)
    pytest.raises(ValueError, lambda: da.linalg.solve_triangular(dA, db))

    A = np.random.default_rng().integers(0, 10, (10, 10))
    b = np.random.default_rng().integers(1, 10, 10)
    dA = da.from_array(A, chunks=(3, 3))
    db = da.from_array(b, chunks=5)
    pytest.raises(ValueError, lambda: da.linalg.solve_triangular(dA, db))


@pytest.mark.parametrize(("shape", "chunk"), [(20, 10), (50, 10)])
def test_solve(shape, chunk):
    rng = np.random.default_rng(1)

    A = rng.integers(1, 10, (shape, shape))
    dA = da.from_array(A, (chunk, chunk))

    # vector
    b = rng.integers(1, 10, shape)
    db = da.from_array(b, chunk)

    res = da.linalg.solve(dA, db)
    assert_eq(res, scipy.linalg.solve(A, b), check_graph=False)
    assert_eq(dA.dot(res), b.astype(float), check_graph=False)

    # tall-and-skinny matrix
    b = rng.integers(1, 10, (shape, 5))
    db = da.from_array(b, (chunk, 5))

    res = da.linalg.solve(dA, db)
    assert_eq(res, scipy.linalg.solve(A, b), check_graph=False)
    assert_eq(dA.dot(res), b.astype(float), check_graph=False)

    # matrix
    b = rng.integers(1, 10, (shape, shape))
    db = da.from_array(b, (chunk, chunk))

    res = da.linalg.solve(dA, db)
    assert_eq(res, scipy.linalg.solve(A, b), check_graph=False)
    assert_eq(dA.dot(res), b.astype(float), check_graph=False)


@pytest.mark.parametrize(("shape", "chunk"), [(20, 10), (50, 10)])
def test_inv(shape, chunk):
    rng = np.random.default_rng(1)

    A = rng.integers(1, 10, (shape, shape))
    dA = da.from_array(A, (chunk, chunk))

    res = da.linalg.inv(dA)
    assert_eq(res, scipy.linalg.inv(A), check_graph=False)
    assert_eq(dA.dot(res), np.eye(shape, dtype=float), check_graph=False)


def _get_symmat(size):
    rng = np.random.default_rng(1)
    A = rng.integers(1, 21, (size, size))
    lA = np.tril(A)
    return lA.dot(lA.T)


# `sym_pos` kwarg was deprecated in scipy 1.9.0
# ref: https://github.com/dask/dask/issues/9335
def _scipy_linalg_solve(a, b, assume_a):
    if Version(scipy.__version__) >= Version("1.9.0"):
        return scipy.linalg.solve(a=a, b=b, assume_a=assume_a)
    elif assume_a == "pos":
        return scipy.linalg.solve(a=a, b=b, sym_pos=True)
    else:
        return scipy.linalg.solve(a=a, b=b)


@pytest.mark.parametrize(("shape", "chunk"), [(20, 10), (30, 6)])
def test_solve_assume_a(shape, chunk):
    rng = np.random.default_rng(1)

    A = _get_symmat(shape)
    dA = da.from_array(A, (chunk, chunk))

    # vector
    b = rng.integers(1, 10, shape)
    db = da.from_array(b, chunk)

    res = da.linalg.solve(dA, db, assume_a="pos")
    assert_eq(res, _scipy_linalg_solve(A, b, assume_a="pos"), check_graph=False)
    assert_eq(dA.dot(res), b.astype(float), check_graph=False)

    # tall-and-skinny matrix
    b = rng.integers(1, 10, (shape, 5))
    db = da.from_array(b, (chunk, 5))

    res = da.linalg.solve(dA, db, assume_a="pos")
    assert_eq(res, _scipy_linalg_solve(A, b, assume_a="pos"), check_graph=False)
    assert_eq(dA.dot(res), b.astype(float), check_graph=False)

    # matrix
    b = rng.integers(1, 10, (shape, shape))
    db = da.from_array(b, (chunk, chunk))

    res = da.linalg.solve(dA, db, assume_a="pos")
    assert_eq(res, _scipy_linalg_solve(A, b, assume_a="pos"), check_graph=False)
    assert_eq(dA.dot(res), b.astype(float), check_graph=False)

    with pytest.warns(FutureWarning, match="sym_pos keyword is deprecated"):
        res = da.linalg.solve(dA, db, sym_pos=True)
        assert_eq(res, _scipy_linalg_solve(A, b, assume_a="pos"), check_graph=False)
        assert_eq(dA.dot(res), b.astype(float), check_graph=False)

    with pytest.warns(FutureWarning, match="sym_pos keyword is deprecated"):
        res = da.linalg.solve(dA, db, sym_pos=False)
        assert_eq(res, _scipy_linalg_solve(A, b, assume_a="gen"), check_graph=False)
        assert_eq(dA.dot(res), b.astype(float), check_graph=False)


@pytest.mark.parametrize(("shape", "chunk"), [(20, 10), (12, 3), (30, 3), (30, 6)])
def test_cholesky(shape, chunk):
    A = _get_symmat(shape)
    dA = da.from_array(A, (chunk, chunk))
    assert_eq(
        da.linalg.cholesky(dA).compute(),
        scipy.linalg.cholesky(A),
        check_graph=False,
        check_chunks=False,
    )
    assert_eq(
        da.linalg.cholesky(dA, lower=True),
        scipy.linalg.cholesky(A, lower=True),
        check_graph=False,
        check_chunks=False,
    )


@pytest.mark.parametrize("iscomplex", [False, True])
@pytest.mark.parametrize(("nrow", "ncol", "chunk"), [(20, 10, 5), (100, 10, 10)])
def test_lstsq(nrow, ncol, chunk, iscomplex):
    rng = np.random.default_rng(1)
    A = rng.integers(1, 20, (nrow, ncol))
    b = rng.integers(1, 20, nrow)
    if iscomplex:
        A = A + 1.0j * rng.integers(1, 20, A.shape)
        b = b + 1.0j * rng.integers(1, 20, b.shape)

    dA = da.from_array(A, (chunk, ncol))
    db = da.from_array(b, chunk)

    x, r, rank, s = np.linalg.lstsq(A, b, rcond=-1)
    dx, dr, drank, ds = da.linalg.lstsq(dA, db)

    assert_eq(dx, x)
    assert_eq(dr, r)
    assert drank.compute() == rank
    assert_eq(ds, s)

    # reduce rank causes multicollinearity, only compare rank
    A[:, 1] = A[:, 2]
    dA = da.from_array(A, (chunk, ncol))
    db = da.from_array(b, chunk)
    x, r, rank, s = np.linalg.lstsq(
        A, b, rcond=np.finfo(np.double).eps * max(nrow, ncol)
    )
    assert rank == ncol - 1
    dx, dr, drank, ds = da.linalg.lstsq(dA, db)
    assert drank.compute() == rank

    # 2D case
    A = rng.integers(1, 20, (nrow, ncol))
    b2D = rng.integers(1, 20, (nrow, ncol // 2))
    if iscomplex:
        A = A + 1.0j * rng.integers(1, 20, A.shape)
        b2D = b2D + 1.0j * rng.integers(1, 20, b2D.shape)
    dA = da.from_array(A, (chunk, ncol))
    db2D = da.from_array(b2D, (chunk, ncol // 2))
    x, r, rank, s = np.linalg.lstsq(A, b2D, rcond=-1)
    dx, dr, drank, ds = da.linalg.lstsq(dA, db2D)

    assert_eq(dx, x)
    assert_eq(dr, r)
    assert drank.compute() == rank
    assert_eq(ds, s)


def test_no_chunks_svd():
    x = np.random.default_rng().random((100, 10))
    u, s, v = np.linalg.svd(x, full_matrices=False)

    for chunks in [((np.nan,) * 10, (10,)), ((np.nan,) * 10, (np.nan,))]:
        dx = da.from_array(x, chunks=(10, 10))
        dx._chunks = chunks

        du, ds, dv = da.linalg.svd(dx)

        assert_eq(s, ds)
        assert_eq(u.dot(np.diag(s)).dot(v), du.dot(da.diag(ds)).dot(dv))
        assert_eq(du.T.dot(du), np.eye(10))
        assert_eq(dv.T.dot(dv), np.eye(10))

        dx = da.from_array(x, chunks=(10, 10))
        dx._chunks = ((np.nan,) * 10, (np.nan,))
        assert_eq(abs(v), abs(dv))
        assert_eq(abs(u), abs(du))


@pytest.mark.parametrize("shape", [(10, 20), (10, 10), (20, 10)])
@pytest.mark.parametrize("chunks", [(-1, -1), (10, -1), (-1, 10)])
@pytest.mark.parametrize("dtype", ["f4", "f8"])
def test_svd_flip_correction(shape, chunks, dtype):
    # Verify that sign-corrected SVD results can still
    # be used to reconstruct inputs
    x = da.random.default_rng().random(size=shape, chunks=chunks).astype(dtype)
    u, s, v = da.linalg.svd(x)

    # Choose precision in evaluation based on float precision
    decimal = 9 if np.dtype(dtype).itemsize > 4 else 6

    # Validate w/ dask inputs
    uf, vf = svd_flip(u, v)
    assert uf.dtype == u.dtype
    assert vf.dtype == v.dtype
    np.testing.assert_almost_equal(np.asarray(np.dot(uf * s, vf)), x, decimal=decimal)

    # Validate w/ numpy inputs
    uc, vc = svd_flip(*da.compute(u, v))
    assert uc.dtype == u.dtype
    assert vc.dtype == v.dtype
    np.testing.assert_almost_equal(np.asarray(np.dot(uc * s, vc)), x, decimal=decimal)


@pytest.mark.parametrize("dtype", ["f2", "f4", "f8", "f16", "c8", "c16", "c32"])
@pytest.mark.parametrize("u_based", [True, False])
def test_svd_flip_sign(dtype, u_based):
    try:
        x = np.array(
            [[1, -1, 1, -1], [1, -1, 1, -1], [-1, 1, 1, -1], [-1, 1, 1, -1]],
            dtype=dtype,
        )
    except TypeError:
        pytest.skip("128-bit floats not supported by NumPy")
    u, v = svd_flip(x, x.T, u_based_decision=u_based)
    assert u.dtype == x.dtype
    assert v.dtype == x.dtype
    # Verify that all singular vectors have same
    # sign except for the last one (i.e. last column)
    y = x.copy()
    y[:, -1] *= y.dtype.type(-1)
    assert_eq(u, y)
    assert_eq(v, y.T)


@pytest.mark.parametrize("chunks", [(10, -1), (-1, 10), (9, -1), (-1, 9)])
@pytest.mark.parametrize("shape", [(10, 100), (100, 10), (10, 10)])
def test_svd_supported_array_shapes(chunks, shape):
    # Test the following cases for tall-skinny, short-fat and square arrays:
    # - no chunking
    # - chunking that contradicts shape (e.g. a 10x100 array with 9x100 chunks)
    # - chunking that aligns with shape (e.g. a 10x100 array with 10x9 chunks)
    x = np.random.default_rng().random(shape)
    dx = da.from_array(x, chunks=chunks)

    du, ds, dv = da.linalg.svd(dx)
    du, dv = da.compute(du, dv)

    nu, ns, nv = np.linalg.svd(x, full_matrices=False)

    # Correct signs before comparison
    du, dv = svd_flip(du, dv)
    nu, nv = svd_flip(nu, nv)

    assert_eq(du, nu)
    assert_eq(ds, ns)
    assert_eq(dv, nv)


def test_svd_incompatible_chunking():
    with pytest.raises(
        NotImplementedError, match="Array must be chunked in one dimension only"
    ):
        x = da.random.default_rng().random((10, 10), chunks=(5, 5))
        da.linalg.svd(x)


@pytest.mark.parametrize("ndim", [0, 1, 3])
def test_svd_incompatible_dimensions(ndim):
    with pytest.raises(ValueError, match="Array must be 2D"):
        x = da.random.default_rng().random((10,) * ndim, chunks=(-1,) * ndim)
        da.linalg.svd(x)


@pytest.mark.xfail(
    sys.platform == "darwin" and _np_version < Version("1.22"),
    reason="https://github.com/dask/dask/issues/7189",
    strict=False,
)
@pytest.mark.parametrize(
    "shape, chunks, axis",
    [[(5,), (2,), None], [(5,), (2,), 0], [(5,), (2,), (0,)], [(5, 6), (2, 2), None]],
)
@pytest.mark.parametrize("norm", [None, 1, -1, np.inf, -np.inf])
@pytest.mark.parametrize("keepdims", [False, True])
def test_norm_any_ndim(shape, chunks, axis, norm, keepdims):
    a = np.random.default_rng().random(shape)
    d = da.from_array(a, chunks=chunks)

    a_r = np.linalg.norm(a, ord=norm, axis=axis, keepdims=keepdims)
    d_r = da.linalg.norm(d, ord=norm, axis=axis, keepdims=keepdims)

    assert_eq(a_r, d_r)


@pytest.mark.xfail(
    _np_version < Version("1.23"),
    reason="https://github.com/numpy/numpy/pull/17709",
    strict=False,
)
@pytest.mark.parametrize("precision", ["single", "double"])
@pytest.mark.parametrize("isreal", [True, False])
@pytest.mark.parametrize("keepdims", [False, True])
@pytest.mark.parametrize("norm", [None, 1, -1, np.inf, -np.inf])
def test_norm_any_prec(norm, keepdims, precision, isreal):
    shape, chunks, axis = (5,), (2,), None

    precs_r = {"single": "float32", "double": "float64"}
    precs_c = {"single": "complex64", "double": "complex128"}

    dtype = precs_r[precision] if isreal else precs_c[precision]

    a = np.random.default_rng().random(shape).astype(dtype)
    d = da.from_array(a, chunks=chunks)
    d_a = np.linalg.norm(a, ord=norm, axis=axis, keepdims=keepdims)
    d_r = da.linalg.norm(d, ord=norm, axis=axis, keepdims=keepdims)

    assert d_r.dtype == precs_r[precision]
    assert d_r.dtype == d_a.dtype


@pytest.mark.slow
@pytest.mark.xfail(
    sys.platform == "darwin" and _np_version < Version("1.22"),
    reason="https://github.com/dask/dask/issues/7189",
    strict=False,
)
@pytest.mark.parametrize(
    "shape, chunks",
    [
        [(5,), (2,)],
        [(5, 3), (2, 2)],
        [(4, 5, 3), (2, 2, 2)],
        [(4, 5, 2, 3), (2, 2, 2, 2)],
        [(2, 5, 2, 4, 3), (2, 2, 2, 2, 2)],
    ],
)
@pytest.mark.parametrize("norm", [None, 1, -1, np.inf, -np.inf])
@pytest.mark.parametrize("keepdims", [False, True])
def test_norm_any_slice(shape, chunks, norm, keepdims):
    a = np.random.default_rng().random(shape)
    d = da.from_array(a, chunks=chunks)

    for firstaxis in range(len(shape)):
        for secondaxis in range(len(shape)):
            if firstaxis != secondaxis:
                axis = (firstaxis, secondaxis)
            else:
                axis = firstaxis
            a_r = np.linalg.norm(a, ord=norm, axis=axis, keepdims=keepdims)
            d_r = da.linalg.norm(d, ord=norm, axis=axis, keepdims=keepdims)
            assert_eq(a_r, d_r)


@pytest.mark.parametrize(
    "shape, chunks, axis", [[(5,), (2,), None], [(5,), (2,), 0], [(5,), (2,), (0,)]]
)
@pytest.mark.parametrize("norm", [0, 2, -2, 0.5])
@pytest.mark.parametrize("keepdims", [False, True])
def test_norm_1dim(shape, chunks, axis, norm, keepdims):
    a = np.random.default_rng().random(shape)
    d = da.from_array(a, chunks=chunks)

    a_r = np.linalg.norm(a, ord=norm, axis=axis, keepdims=keepdims)
    d_r = da.linalg.norm(d, ord=norm, axis=axis, keepdims=keepdims)
    assert_eq(a_r, d_r)


@pytest.mark.parametrize(
    "shape, chunks, axis",
    [[(5, 6), (2, 2), None], [(5, 6), (2, 2), (0, 1)], [(5, 6), (2, 2), (1, 0)]],
)
@pytest.mark.parametrize("norm", ["fro", "nuc", 2, -2])
@pytest.mark.parametrize("keepdims", [False, True])
def test_norm_2dim(shape, chunks, axis, norm, keepdims):
    a = np.random.default_rng().random(shape)
    d = da.from_array(a, chunks=chunks)

    # Need one chunk on last dimension for svd.
    if norm == "nuc" or norm == 2 or norm == -2:
        d = d.rechunk({-1: -1})

    a_r = np.linalg.norm(a, ord=norm, axis=axis, keepdims=keepdims)
    d_r = da.linalg.norm(d, ord=norm, axis=axis, keepdims=keepdims)

    assert_eq(a_r, d_r)


@pytest.mark.parametrize(
    "shape, chunks, axis",
    [[(3, 2, 4), (2, 2, 2), (1, 2)], [(2, 3, 4, 5), (2, 2, 2, 2), (-1, -2)]],
)
@pytest.mark.parametrize("norm", ["nuc", 2, -2])
@pytest.mark.parametrize("keepdims", [False, True])
def test_norm_implemented_errors(shape, chunks, axis, norm, keepdims):
    a = np.random.default_rng().random(shape)
    d = da.from_array(a, chunks=chunks)
    if len(shape) > 2 and len(axis) == 2:
        with pytest.raises(NotImplementedError):
            da.linalg.norm(d, ord=norm, axis=axis, keepdims=keepdims)
