from __future__ import annotations

import numpy as np
import pytest
from packaging.version import Version

pytestmark = pytest.mark.gpu

import dask.array as da
from dask.array.utils import assert_eq

cupy = pytest.importorskip("cupy")
cupy_version = Version(cupy.__version__)


@pytest.mark.skipif(
    cupy_version < Version("6.1.0"),
    reason="Requires CuPy 6.1.0+ (with https://github.com/cupy/cupy/pull/2209)",
)
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
    mat = cupy.random.default_rng().random((m, n))
    data = da.from_array(mat, chunks=chunks, name="A", asarray=False)

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
        q, r = da.linalg.tsqr(data)
        assert_eq((m_q, n_q), q.shape)  # shape check
        assert_eq((m_r, n_r), r.shape)  # shape check
        assert_eq(mat, da.dot(q, r))  # accuracy check
        assert_eq(cupy.eye(n_q, n_q), da.dot(q.T, q))  # q must be orthonormal
        assert_eq(r, np.triu(r.rechunk(r.shape[0])))  # r must be upper triangular

        # test SVD
        u, s, vh = da.linalg.tsqr(data, compute_svd=True)
        s_exact = np.linalg.svd(mat)[1]
        assert_eq(s, s_exact)  # s must contain the singular values
        assert_eq((m_u, n_u), u.shape)  # shape check
        assert_eq((n_s,), s.shape)  # shape check
        assert_eq((d_vh, d_vh), vh.shape)  # shape check
        assert_eq(
            np.eye(n_u, n_u), da.dot(u.T, u), check_type=False
        )  # u must be orthonormal
        assert_eq(
            np.eye(d_vh, d_vh), da.dot(vh, vh.T), check_type=False
        )  # vh must be orthonormal
        assert_eq(mat, da.dot(da.dot(u, da.diag(s)), vh[:n_q]))  # accuracy check
    else:
        with pytest.raises(error_type):
            q, r = da.linalg.tsqr(data)
        with pytest.raises(error_type):
            u, s, vh = da.linalg.tsqr(data, compute_svd=True)


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
    mat = cupy.random.default_rng().random((m_min * 2, n_max))
    m, n = m_min * 2, n_max
    mat[0:m_min, 0] += 1
    _c0 = mat[:, 0]
    _r0 = mat[0, :]
    c0 = da.from_array(_c0, chunks=m_min, name="c", asarray=False)
    r0 = da.from_array(_r0, chunks=n_max, name="r", asarray=False)
    data = da.from_array(mat, chunks=chunks, name="A", asarray=False)
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
        q, r = da.linalg.tsqr(data)
        q = q.compute()  # because uncertainty
        r = r.compute()
        assert_eq((m_q, n_q), q.shape)  # shape check
        assert_eq((m_r, n_r), r.shape)  # shape check
        assert_eq(mat, np.dot(q, r))  # accuracy check
        assert_eq(
            np.eye(n_q, n_q), np.dot(q.T, q), check_type=False
        )  # q must be orthonormal
        assert_eq(r, np.triu(r))  # r must be upper triangular

        # test SVD
        u, s, vh = da.linalg.tsqr(data, compute_svd=True)
        u = u.compute()  # because uncertainty
        s = s.compute()
        vh = vh.compute()
        s_exact = np.linalg.svd(mat)[1]
        assert_eq(s, s_exact)  # s must contain the singular values
        assert_eq((m_u, n_u), u.shape)  # shape check
        assert_eq((n_s,), s.shape)  # shape check
        assert_eq((d_vh, d_vh), vh.shape)  # shape check
        assert_eq(
            np.eye(n_u, n_u), np.dot(u.T, u), check_type=False
        )  # u must be orthonormal
        assert_eq(
            np.eye(d_vh, d_vh), np.dot(vh, vh.T), check_type=False
        )  # vh must be orthonormal
        assert_eq(
            mat, np.dot(np.dot(u, np.diag(s)), vh[:n_q]), check_type=False
        )  # accuracy check
    else:
        with pytest.raises(error_type):
            q, r = da.linalg.tsqr(data)
        with pytest.raises(error_type):
            u, s, vh = da.linalg.tsqr(data, compute_svd=True)


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
        q, r = da.linalg.sfqr(data)
        assert_eq((m_q, n_q), q.shape)  # shape check
        assert_eq((m_r, n_r), r.shape)  # shape check
        assert_eq(mat, da.dot(q, r))  # accuracy check
        assert_eq(np.eye(m_qtq, m_qtq), da.dot(q.T, q))  # q must be orthonormal
        assert_eq(r, da.triu(r.rechunk(r.shape[0])))  # r must be upper triangular
    else:
        with pytest.raises(error_type):
            q, r = da.linalg.sfqr(data)


@pytest.mark.parametrize("iscomplex", [False, True])
@pytest.mark.parametrize(("nrow", "ncol", "chunk"), [(20, 10, 5), (100, 10, 10)])
def test_lstsq(nrow, ncol, chunk, iscomplex):
    rng = cupy.random.default_rng(1)
    A = rng.integers(1, 20, (nrow, ncol))
    b = rng.integers(1, 20, nrow)
    if iscomplex:
        A = A + 1.0j * rng.integers(1, 20, A.shape)
        b = b + 1.0j * rng.integers(1, 20, b.shape)

    dA = da.from_array(A, (chunk, ncol))
    db = da.from_array(b, chunk)

    x, r, rank, s = cupy.linalg.lstsq(A, b, rcond=-1)
    dx, dr, drank, ds = da.linalg.lstsq(dA, db)

    assert_eq(dx, x)
    assert_eq(dr, r)
    assert drank.compute() == rank
    assert_eq(ds, s)

    # reduce rank causes multicollinearity, only compare rank
    A[:, 1] = A[:, 2]
    dA = da.from_array(A, (chunk, ncol))
    db = da.from_array(b, chunk)
    x, r, rank, s = cupy.linalg.lstsq(
        A, b, rcond=cupy.finfo(cupy.double).eps * max(nrow, ncol)
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
    x, r, rank, s = cupy.linalg.lstsq(A, b2D, rcond=-1)
    dx, dr, drank, ds = da.linalg.lstsq(dA, db2D)

    assert_eq(dx, x)
    assert_eq(dr, r)
    assert drank.compute() == rank
    assert_eq(ds, s)


def _get_symmat(size):
    rng = cupy.random.default_rng(1)
    A = rng.integers(1, 21, (size, size))
    lA = cupy.tril(A)
    return lA.dot(lA.T)


@pytest.mark.parametrize(("shape", "chunk"), [(20, 10), (12, 3), (30, 3), (30, 6)])
def test_cholesky(shape, chunk):
    scipy_linalg = pytest.importorskip("scipy.linalg")

    A = _get_symmat(shape)
    dA = da.from_array(A, (chunk, chunk))

    # Need to take the transpose because default in `cupy.linalg.cholesky` is
    # to return lower triangle
    assert_eq(
        da.linalg.cholesky(dA),
        cupy.linalg.cholesky(A).T,
        check_graph=False,
        check_chunks=False,
    )
    assert_eq(
        da.linalg.cholesky(dA, lower=True).map_blocks(cupy.asnumpy),
        scipy_linalg.cholesky(cupy.asnumpy(A), lower=True),
        check_graph=False,
        check_chunks=False,
    )
