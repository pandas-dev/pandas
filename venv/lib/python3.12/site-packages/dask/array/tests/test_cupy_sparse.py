from __future__ import annotations

import pytest

pytestmark = pytest.mark.gpu

import dask.array as da
from dask.array.utils import assert_eq

cupy = pytest.importorskip("cupy")
cupyx = pytest.importorskip("cupyx")


def test_sparse_hstack_vstack_csr():
    pytest.importorskip("cupyx")
    x = cupy.arange(24, dtype=cupy.float32).reshape(4, 6)

    sp = da.from_array(x, chunks=(2, 3), asarray=False, fancy=False)
    sp = sp.map_blocks(cupyx.scipy.sparse.csr_matrix, dtype=cupy.float32)

    y = sp.compute()

    assert cupyx.scipy.sparse.isspmatrix(y)
    assert_eq(x, y.todense())


@pytest.mark.parametrize("axis", [0, 1])
def test_sparse_concatenate(axis):
    pytest.importorskip("cupyx")

    rng = da.random.default_rng(cupy.random.default_rng())
    meta = cupyx.scipy.sparse.csr_matrix((0, 0))

    xs = []
    ys = []
    for _ in range(2):
        x = rng.random((1000, 10), chunks=(100, 10))
        x[x < 0.9] = 0
        xs.append(x)
        ys.append(x.map_blocks(cupyx.scipy.sparse.csr_matrix, meta=meta))

    z = da.concatenate(ys, axis=axis)
    z = z.compute()

    if axis == 0:
        sp_concatenate = cupyx.scipy.sparse.vstack
    elif axis == 1:
        sp_concatenate = cupyx.scipy.sparse.hstack
    z_expected = sp_concatenate(
        [cupyx.scipy.sparse.csr_matrix(e.compute()) for e in xs]
    )

    assert (z.toarray() == z_expected.toarray()).all()


@pytest.mark.parametrize("sp_format", ["csr", "csc"])
@pytest.mark.parametrize(
    "input_sizes",
    [
        {
            "x_shape": (4, 8),
            "y_shape": (8, 6),
            "x_chunks": (2, 4),
            "y_chunks": (4, 3),
        },
        {
            "x_shape": (4, 4),
            "y_shape": (4, 4),
            "x_chunks": (2, 2),
            "y_chunks": (2, 2),
        },
        {
            "x_shape": (4, 4),
            "y_shape": (4, 4),
            "x_chunks": (4, 2),
            "y_chunks": (2, 4),
        },
    ],
)
def test_sparse_dot(sp_format, input_sizes):
    pytest.importorskip("cupyx")

    if sp_format == "csr":
        sp_matrix = cupyx.scipy.sparse.csr_matrix
    elif sp_format == "csc":
        sp_matrix = cupyx.scipy.sparse.csc_matrix
    dtype = "f"
    density = 0.3
    rng = cupy.random.default_rng()
    x_shape, x_chunks = input_sizes["x_shape"], input_sizes["x_chunks"]
    y_shape, y_chunks = input_sizes["y_shape"], input_sizes["y_chunks"]
    x = rng.random(x_shape, dtype=dtype)
    y = rng.random(y_shape, dtype=dtype)
    x[x < 1 - density] = 0
    y[y < 1 - density] = 0
    z = x.dot(y)

    da_x = da.from_array(x, chunks=x_chunks, asarray=False, fancy=False)
    da_y = da.from_array(y, chunks=y_chunks, asarray=False, fancy=False)
    da_x = da_x.map_blocks(sp_matrix, meta=sp_matrix(cupy.array([0], dtype=dtype)))
    da_y = da_y.map_blocks(sp_matrix, meta=sp_matrix(cupy.array([0], dtype=dtype)))
    da_z = da.dot(da_x, da_y).compute()

    assert cupyx.scipy.sparse.isspmatrix(da_z)
    assert_eq(z, da_z.todense())
