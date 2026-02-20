from __future__ import annotations

import numpy as np
import pytest
from packaging.version import Version

import dask
import dask.array as da
from dask.array.reductions import nannumel, numel
from dask.array.utils import assert_eq

sparse = pytest.importorskip("sparse")
SPARSE_VERSION = Version(sparse.__version__)
if sparse:
    # Test failures on older versions of Numba.
    # Conda-Forge provides 0.35.0 on windows right now, causing failures like
    # searchsorted() got an unexpected keyword argument 'side'
    pytest.importorskip("numba", minversion="0.40.0")


functions = [
    lambda x: x,
    lambda x: da.expm1(x),
    lambda x: 2 * x,
    lambda x: x / 2,
    lambda x: x**2,
    lambda x: x + x,
    lambda x: x * x,
    lambda x: x[0],
    lambda x: x[:, 1],
    lambda x: x[:1, None, 1:3],
    lambda x: x.T,
    lambda x: da.transpose(x, (1, 2, 0)),
    lambda x: da.nanmean(x),
    lambda x: da.nanmean(x, axis=1),
    lambda x: da.nanmax(x),
    lambda x: da.nanmin(x),
    lambda x: da.nanprod(x),
    lambda x: da.nanstd(x),
    lambda x: da.nanvar(x),
    lambda x: da.nansum(x),
    # These nan* variants are are not implemented by sparse.COO
    # lambda x: da.median(x, axis=0),
    # lambda x: da.nanargmax(x),
    # lambda x: da.nanargmin(x),
    # lambda x: da.nancumprod(x, axis=0),
    # lambda x: da.nancumsum(x, axis=0),
    lambda x: x.sum(),
    lambda x: x.moment(order=0),
    lambda x: x.mean(),
    lambda x: x.mean(axis=1),
    lambda x: x.std(),
    lambda x: x.var(),
    lambda x: x.dot(np.arange(x.shape[-1])),
    lambda x: x.dot(np.eye(x.shape[-1])),
    lambda x: da.tensordot(x, np.ones(x.shape[:2]), axes=[(0, 1), (0, 1)]),
    lambda x: x.sum(axis=0),
    lambda x: x.max(axis=0),
    lambda x: x.sum(axis=(1, 2)),
    lambda x: x.astype(np.complex128),
    lambda x: x.map_blocks(lambda x: x * 2),
    lambda x: x.map_overlap(lambda x: x * 2, depth=0, trim=True, boundary="none"),
    lambda x: x.map_overlap(lambda x: x * 2, depth=0, trim=False, boundary="none"),
    lambda x: x.round(1),
    lambda x: x.reshape((x.shape[0] * x.shape[1], x.shape[2])),
    lambda x: abs(x),
    lambda x: x > 0.5,
    lambda x: x.rechunk((4, 4, 4)),
    lambda x: x.rechunk((2, 2, 1)),
    lambda x: np.isneginf(x),
    lambda x: np.isposinf(x),
    lambda x: np.zeros_like(x),
    lambda x: np.ones_like(x),
    lambda x: np.full_like(x, fill_value=2),
]


@pytest.mark.parametrize("func", functions)
def test_basic(func):
    x = da.random.random((2, 3, 4), chunks=(1, 2, 2))
    x[x < 0.8] = 0

    y = x.map_blocks(sparse.COO.from_numpy)

    xx = func(x)
    yy = func(y)

    assert_eq(xx, yy, check_meta=False)

    if yy.shape:
        zz = yy.compute()
        if not isinstance(zz, sparse.COO):
            assert (zz != 1).sum() > np.prod(zz.shape) / 2  # mostly dense


def test_tensordot():
    x = da.random.random((2, 3, 4), chunks=(1, 2, 2))
    x[x < 0.8] = 0
    y = da.random.random((4, 3, 2), chunks=(2, 2, 1))
    y[y < 0.8] = 0

    xx = x.map_blocks(sparse.COO.from_numpy)
    yy = y.map_blocks(sparse.COO.from_numpy)

    assert_eq(da.tensordot(x, y, axes=(2, 0)), da.tensordot(xx, yy, axes=(2, 0)))
    assert_eq(da.tensordot(x, y, axes=(1, 1)), da.tensordot(xx, yy, axes=(1, 1)))
    assert_eq(
        da.tensordot(x, y, axes=((1, 2), (1, 0))),
        da.tensordot(xx, yy, axes=((1, 2), (1, 0))),
    )


def test_metadata():
    y = da.random.random((10, 10), chunks=(5, 5))
    y[y < 0.8] = 0
    z = sparse.COO.from_numpy(y.compute())
    y = y.map_blocks(sparse.COO.from_numpy)

    assert isinstance(y._meta, sparse.COO)
    assert isinstance((y + 1)._meta, sparse.COO)
    assert isinstance(y.sum(axis=0)._meta, sparse.COO)
    assert isinstance(y.var(axis=0)._meta, sparse.COO)
    assert isinstance(y[:5, ::2]._meta, sparse.COO)
    assert isinstance(y.rechunk((2, 2))._meta, sparse.COO)
    assert isinstance((y - z)._meta, sparse.COO)
    assert isinstance(y.persist()._meta, sparse.COO)
    assert isinstance(np.concatenate([y, y])._meta, sparse.COO)
    assert isinstance(np.concatenate([y, y[:0], y])._meta, sparse.COO)
    assert isinstance(np.stack([y, y])._meta, sparse.COO)
    assert isinstance(np.stack([y[:0], y[:0]])._meta, sparse.COO)
    assert isinstance(np.concatenate([y[:0], y[:0]])._meta, sparse.COO)


def test_html_repr():
    pytest.importorskip("jinja2")
    y = da.random.random((10, 10), chunks=(5, 5))
    y[y < 0.8] = 0
    y = y.map_blocks(sparse.COO.from_numpy)

    text = y._repr_html_()

    assert "COO" in text
    assert "sparse" in text
    assert "Bytes" not in text


def test_from_delayed_meta():
    def f():
        return sparse.COO.from_numpy(np.eye(3))

    d = dask.delayed(f)()
    x = da.from_delayed(d, shape=(3, 3), meta=sparse.COO.from_numpy(np.eye(1)))
    assert isinstance(x._meta, sparse.COO)
    assert_eq(x, x)


def test_from_array():
    x = sparse.COO.from_numpy(np.eye(10))
    d = da.from_array(x, chunks=(5, 5))

    assert isinstance(d._meta, sparse.COO)
    assert_eq(d, d)
    assert isinstance(d.compute(), sparse.COO)


def test_map_blocks():
    x = da.eye(10, chunks=5)
    y = x.map_blocks(sparse.COO.from_numpy, meta=sparse.COO.from_numpy(np.eye(1)))
    assert isinstance(y._meta, sparse.COO)
    assert_eq(y, y)


def test_meta_from_array():
    x = sparse.COO.from_numpy(np.eye(1))
    y = da.utils.meta_from_array(x, ndim=2)
    assert isinstance(y, sparse.COO)


@pytest.mark.parametrize("numel", [numel, nannumel])
@pytest.mark.parametrize("axis", [0, (0, 1), None])
@pytest.mark.parametrize("keepdims", [True, False])
def test_numel(numel, axis, keepdims):
    x = np.random.random((2, 3, 4))
    x[x < 0.8] = 0
    x[x > 0.9] = np.nan

    xs = sparse.COO.from_numpy(x, fill_value=0.0)

    assert_eq(
        numel(x, axis=axis, keepdims=keepdims), numel(xs, axis=axis, keepdims=keepdims)
    )
