from __future__ import annotations

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
