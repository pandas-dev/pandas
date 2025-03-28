from __future__ import annotations

import numpy as np
import pytest

import dask.array as da
from dask.array.utils import assert_eq

xr = pytest.importorskip("xarray")


def test_mean():
    y = da.mean(xr.DataArray([1, 2, 3.0]))
    assert isinstance(y, da.Array)
    assert_eq(y, y)


def test_asarray():
    y = da.asarray(xr.DataArray([1, 2, 3.0]))
    assert isinstance(y, da.Array)
    assert_eq(y, y)


def test_asanyarray():
    y = da.asanyarray(xr.DataArray([1, 2, 3.0]))
    assert isinstance(y, da.Array)
    assert_eq(y, y)


def test_asarray_xarray_intersphinx_workaround():
    # test that the intersphinx workaround in https://github.com/pydata/xarray/issues/4279 works
    module = xr.DataArray.__module__
    try:
        xr.DataArray.__module__ = "xarray"
        y = da.asarray(xr.DataArray([1, 2, 3.0]))
        assert isinstance(y, da.Array)
        assert type(y._meta).__name__ == "ndarray"
        assert_eq(y, y)
    finally:
        xr.DataArray.__module__ = module


def test_fft():
    # Regression test for https://github.com/dask/dask/issues/9679
    coord = da.arange(8, chunks=-1)
    data = da.random.random((8, 8), chunks=-1) + 1
    x = xr.DataArray(data, coords={"x": coord, "y": coord}, dims=["x", "y"])
    result = da.fft.fft(x)
    expected = da.fft.fft(x.data)
    assert_eq(result, expected)


def test_polyfit_reshaping():
    # Regression test for https://github.com/pydata/xarray/issues/4554
    arr = xr.DataArray(da.ones((10, 20, 30), chunks=(1, 5, 30)), dims=["z", "y", "x"])
    result = arr.polyfit("x", 1)
    assert result.polyfit_coefficients.chunks == ((2,), (1,) * 10, (5,) * 4)


def test_positional_indexer_multiple_variables():
    n = 200
    ds = xr.Dataset(
        data_vars=dict(
            a=(
                ["x", "y", "time"],
                da.random.randint(1, 100, (10, 10, n), chunks=(-1, -1, n // 2)),
            ),
            b=(
                ["x", "y", "time"],
                da.random.randint(1, 100, (10, 10, n), chunks=(-1, -1, n // 2)),
            ),
        ),
        coords=dict(
            x=list(range(10)),
            y=list(range(10)),
            time=np.arange(n),
        ),
    )
    indexer = np.arange(n)
    np.random.shuffle(indexer)
    result = ds.isel(time=indexer)
    graph = result.__dask_graph__()
    assert len({k for k in graph if "shuffle-taker" in k}) == 4
    assert len({k for k in graph if "shuffle-sorter" in k}) == 2
