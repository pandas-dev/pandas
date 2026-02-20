"""isort:skip_file"""

from __future__ import annotations

import pickle
from typing import TYPE_CHECKING, Any

import numpy as np
import pytest

if TYPE_CHECKING:
    import dask
    import dask.array as da
    import distributed
else:
    dask = pytest.importorskip("dask")
    da = pytest.importorskip("dask.array")
    distributed = pytest.importorskip("distributed")

import contextlib

from dask.distributed import Client, Lock
from distributed.client import futures_of
from distributed.utils_test import (
    cleanup,  # noqa: F401
    client,  # noqa: F401
    cluster,
    cluster_fixture,  # noqa: F401
    gen_cluster,
    loop,  # noqa: F401
    loop_in_thread,  # noqa: F401
)

import xarray as xr
from xarray.backends.locks import HDF5_LOCK, CombinedLock, SerializableLock
from xarray.tests import (
    assert_allclose,
    assert_identical,
    has_h5netcdf,
    has_netCDF4,
    has_scipy,
    requires_cftime,
    requires_netCDF4,
    requires_zarr,
)
from xarray.tests.test_backends import (
    ON_WINDOWS,
    create_tmp_file,
)
from xarray.tests.test_dataset import create_test_data


@pytest.fixture
def tmp_netcdf_filename(tmpdir):
    return str(tmpdir.join("testfile.nc"))


ENGINES = []
if has_scipy:
    ENGINES.append("scipy")
if has_netCDF4:
    ENGINES.append("netcdf4")
if has_h5netcdf:
    ENGINES.append("h5netcdf")

NC_FORMATS = {
    "netcdf4": [
        "NETCDF3_CLASSIC",
        "NETCDF3_64BIT_OFFSET",
        "NETCDF3_64BIT_DATA",
        "NETCDF4_CLASSIC",
        "NETCDF4",
    ],
    "scipy": ["NETCDF3_CLASSIC", "NETCDF3_64BIT"],
    "h5netcdf": ["NETCDF4"],
}

ENGINES_AND_FORMATS = [
    ("netcdf4", "NETCDF3_CLASSIC"),
    ("netcdf4", "NETCDF4_CLASSIC"),
    ("netcdf4", "NETCDF4"),
    ("h5netcdf", "NETCDF4"),
    ("scipy", "NETCDF3_64BIT"),
]


@pytest.mark.parametrize("engine,nc_format", ENGINES_AND_FORMATS)
@pytest.mark.parametrize("compute", [True, False])
def test_dask_distributed_netcdf_roundtrip(
    loop,  # noqa: F811
    tmp_netcdf_filename,
    engine,
    nc_format,
    compute,
):
    if engine not in ENGINES:
        pytest.skip("engine not available")

    chunks = {"dim1": 4, "dim2": 3, "dim3": 6}

    with cluster() as (s, [_a, _b]):
        with Client(s["address"], loop=loop):
            original = create_test_data().chunk(chunks)

            if engine == "scipy":
                with pytest.raises(NotImplementedError):
                    original.to_netcdf(
                        tmp_netcdf_filename, engine=engine, format=nc_format
                    )
                return

            result = original.to_netcdf(
                tmp_netcdf_filename, engine=engine, format=nc_format, compute=compute
            )
            if not compute:
                result.compute()

            with xr.open_dataset(
                tmp_netcdf_filename, chunks=chunks, engine=engine
            ) as restored:
                assert isinstance(restored.var1.data, da.Array)
                computed = restored.compute()
                assert_allclose(original, computed)


@requires_netCDF4
def test_dask_distributed_write_netcdf_with_dimensionless_variables(
    loop,  # noqa: F811
    tmp_netcdf_filename,
):
    with cluster() as (s, [_a, _b]):
        with Client(s["address"], loop=loop):
            original = xr.Dataset({"x": da.zeros(())})
            original.to_netcdf(tmp_netcdf_filename)

            with xr.open_dataset(tmp_netcdf_filename) as actual:
                assert actual.x.shape == ()


@requires_cftime
@requires_netCDF4
@pytest.mark.parametrize("parallel", (True, False))
def test_open_mfdataset_can_open_files_with_cftime_index(parallel, tmp_path):
    T = xr.date_range("20010101", "20010501", calendar="360_day", use_cftime=True)
    Lon = np.arange(100)
    data = np.random.random((T.size, Lon.size))
    da = xr.DataArray(data, coords={"time": T, "Lon": Lon}, name="test")
    file_path = tmp_path / "test.nc"
    da.to_netcdf(file_path)
    with cluster() as (s, [_a, _b]):
        with Client(s["address"]):
            with xr.open_mfdataset(file_path, parallel=parallel) as tf:
                assert_identical(tf["test"], da)


@requires_cftime
@requires_netCDF4
@pytest.mark.parametrize("parallel", (True, False))
def test_open_mfdataset_multiple_files_parallel_distributed(parallel, tmp_path):
    lon = np.arange(100)
    time = xr.date_range("20010101", periods=100, calendar="360_day", use_cftime=True)
    data = np.random.random((time.size, lon.size))
    da = xr.DataArray(data, coords={"time": time, "lon": lon}, name="test")

    fnames = []
    for i in range(0, 100, 10):
        fname = tmp_path / f"test_{i}.nc"
        da.isel(time=slice(i, i + 10)).to_netcdf(fname)
        fnames.append(fname)

    with cluster() as (s, [_a, _b]):
        with Client(s["address"]):
            with xr.open_mfdataset(
                fnames, parallel=parallel, concat_dim="time", combine="nested"
            ) as tf:
                assert_identical(tf["test"], da)


# TODO: move this to test_backends.py
@requires_cftime
@requires_netCDF4
@pytest.mark.parametrize("parallel", (True, False))
def test_open_mfdataset_multiple_files_parallel(parallel, tmp_path):
    if parallel:
        pytest.skip(
            "Flaky in CI. Would be a welcome contribution to make a similar test reliable."
        )
    lon = np.arange(100)
    time = xr.date_range("20010101", periods=100, calendar="360_day", use_cftime=True)
    data = np.random.random((time.size, lon.size))
    da = xr.DataArray(data, coords={"time": time, "lon": lon}, name="test")

    fnames = []
    for i in range(0, 100, 10):
        fname = tmp_path / f"test_{i}.nc"
        da.isel(time=slice(i, i + 10)).to_netcdf(fname)
        fnames.append(fname)

    for get in [dask.threaded.get, dask.multiprocessing.get, dask.local.get_sync, None]:
        with dask.config.set(scheduler=get):
            with xr.open_mfdataset(
                fnames, parallel=parallel, concat_dim="time", combine="nested"
            ) as tf:
                assert_identical(tf["test"], da)


@pytest.mark.parametrize("engine,nc_format", ENGINES_AND_FORMATS)
def test_dask_distributed_read_netcdf_integration_test(
    loop,  # noqa: F811
    tmp_netcdf_filename,
    engine,
    nc_format,
):
    if engine not in ENGINES:
        pytest.skip("engine not available")

    chunks = {"dim1": 4, "dim2": 3, "dim3": 6}

    with cluster() as (s, [_a, _b]):
        with Client(s["address"], loop=loop):
            original = create_test_data()
            original.to_netcdf(tmp_netcdf_filename, engine=engine, format=nc_format)

            with xr.open_dataset(
                tmp_netcdf_filename, chunks=chunks, engine=engine
            ) as restored:
                assert isinstance(restored.var1.data, da.Array)
                computed = restored.compute()
                assert_allclose(original, computed)


# fixture vendored from dask
# heads-up, this is using quite private zarr API
# https://github.com/dask/dask/blob/e04734b4d8959ba259801f2e2a490cb4ee8d891f/dask/tests/test_distributed.py#L338-L358
@pytest.fixture
def zarr(client):  # noqa: F811
    zarr_lib = pytest.importorskip("zarr")
    # Zarr-Python 3 lazily allocates a dedicated thread/IO loop
    # for to execute async tasks. To avoid having this thread
    # be picked up as a "leaked thread", we manually trigger it's
    # creation before using zarr
    try:
        _ = zarr_lib.core.sync._get_loop()
        _ = zarr_lib.core.sync._get_executor()
        yield zarr_lib
    except AttributeError:
        yield zarr_lib
    finally:
        # Zarr-Python 3 lazily allocates an IO thread, a thread pool executor, and
        # an IO loop. Here we clean up these resources to avoid leaking threads
        # In normal operations, this is done as by an atexit handler when Zarr
        # is shutting down.
        with contextlib.suppress(AttributeError):
            zarr_lib.core.sync.cleanup_resources()


@requires_zarr
@pytest.mark.parametrize("consolidated", [True, False])
@pytest.mark.parametrize("compute", [True, False])
def test_dask_distributed_zarr_integration_test(
    client,  # noqa: F811
    zarr,
    consolidated: bool,
    compute: bool,
) -> None:
    if consolidated:
        write_kwargs: dict[str, Any] = {"consolidated": True}
        read_kwargs: dict[str, Any] = {"backend_kwargs": {"consolidated": True}}
    else:
        write_kwargs = read_kwargs = {}
    chunks = {"dim1": 4, "dim2": 3, "dim3": 5}
    original = create_test_data().chunk(chunks)
    with create_tmp_file(allow_cleanup_failure=ON_WINDOWS, suffix=".zarrc") as filename:
        maybe_futures = original.to_zarr(  # type: ignore[call-overload]  #mypy bug?
            filename, compute=compute, **write_kwargs
        )
        if not compute:
            maybe_futures.compute()
        with xr.open_dataset(
            filename, chunks="auto", engine="zarr", **read_kwargs
        ) as restored:
            assert isinstance(restored.var1.data, da.Array)
            computed = restored.compute()
            assert_allclose(original, computed)


@gen_cluster(client=True)
async def test_async(c, s, a, b) -> None:
    x = create_test_data()
    assert not dask.is_dask_collection(x)
    y = x.chunk({"dim2": 4}) + 10
    assert dask.is_dask_collection(y)
    assert dask.is_dask_collection(y.var1)
    assert dask.is_dask_collection(y.var2)

    z = c.persist(y)
    assert str(z)

    assert dask.is_dask_collection(z)
    assert dask.is_dask_collection(z.var1)
    assert dask.is_dask_collection(z.var2)
    assert len(y.__dask_graph__()) > len(z.__dask_graph__())

    assert not futures_of(y)
    assert futures_of(z)

    future = c.compute(z)
    w = await future
    assert not dask.is_dask_collection(w)
    assert_allclose(x + 10, w)

    assert s.tasks


def test_hdf5_lock() -> None:
    assert isinstance(HDF5_LOCK, SerializableLock)


@gen_cluster(client=True)
async def test_serializable_locks(c, s, a, b) -> None:
    def f(x, lock=None):
        with lock:
            return x + 1

    # note, the creation of Lock needs to be done inside a cluster
    for lock in [
        HDF5_LOCK,
        Lock(),
        Lock("filename.nc"),
        CombinedLock([HDF5_LOCK]),
        CombinedLock([HDF5_LOCK, Lock("filename.nc")]),
    ]:
        futures = c.map(f, list(range(10)), lock=lock)
        await c.gather(futures)

        lock2 = pickle.loads(pickle.dumps(lock))
        assert type(lock) is type(lock2)
