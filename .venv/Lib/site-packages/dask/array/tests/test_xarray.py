from __future__ import annotations

from functools import partial

import numpy as np
import pytest

import dask
import dask.array as da
from dask.array.utils import assert_eq
from dask.utils import tmpdir

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


@pytest.mark.parametrize("compute", [True, False])
def test_xarray_blockwise_fusion_store(compute):
    def custom_scheduler_get(dsk, keys, expected, **kwargs):
        dsk = dsk.__dask_graph__()
        assert (
            len(dsk) == expected
        ), f"False number of tasks got {len(dsk)} but expected {expected}"
        return [42 for _ in keys]

    # First test that this mocking stuff works as expecged
    with pytest.raises(AssertionError, match="False number of tasks"):
        scheduler = partial(custom_scheduler_get, expected=42)
        dask.compute(da.ones(10), scheduler=scheduler)

    coord = da.arange(8, chunks=-1)
    data = da.random.random((8, 8), chunks=-1) + 1
    x = xr.DataArray(data, coords={"x": coord, "y": coord}, dims=["x", "y"])

    y = ((x + 1) * 2) / 2 - 1

    # Everything fused into one compute task
    # one finalize Alias
    expected = 2
    scheduler = partial(custom_scheduler_get, expected=expected)
    dask.compute(y, scheduler=scheduler)

    with tmpdir() as dirname:
        if compute:
            with dask.config.set(scheduler=scheduler):
                y.to_zarr(dirname, compute=True)
        else:
            # There's a delayed finalize store smashed on top which won't be fused by
            # default
            expected += 1
            scheduler = partial(custom_scheduler_get, expected=expected)
            stored = y.to_zarr(dirname, compute=False)
            dask.compute(stored, scheduler=scheduler)


@pytest.mark.parametrize("wrap_xarray", [False, True])
def test_shared_tasks(wrap_xarray):
    total_calls = 0

    def my_func(a: np.ndarray) -> np.ndarray:
        nonlocal total_calls
        total_calls += 1
        return a

    in1 = da.zeros((5, 5), chunks=2)
    res = da.map_blocks(
        my_func, in1, meta=np.array((), dtype=in1.dtype), dtype=in1.dtype
    )
    out1 = res + 1
    out2 = res + 2
    if wrap_xarray:
        out1 = xr.DataArray(out1)
        out2 = xr.DataArray(out2)

    comp_res = dask.compute(out1, out2)

    if wrap_xarray:
        assert isinstance(comp_res[0], xr.DataArray)
        assert isinstance(comp_res[0].data, np.ndarray)
    else:
        assert isinstance(comp_res[0], np.ndarray)
    assert total_calls == in1.blocks.size


def test_slicing():
    ds = xr.Dataset(
        {"z": (("t", "p", "y", "x"), np.ones((1, 1, 31, 40)))},
    )
    ds = ds.chunk()
    subset = ds.isel(t=[0], p=0).z[:, ::10, ::10][:, ::-1, :]
    subset2 = ds.isel(t=[0]).isel(p=0).z[:, ::10, ::10][:, ::-1, :]
    assert_eq(subset.data, subset2.data)
