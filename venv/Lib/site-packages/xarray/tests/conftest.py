from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import pytest

import xarray as xr
from xarray import DataArray, Dataset, DataTree
from xarray.tests import create_test_data, has_cftime, requires_dask


@pytest.fixture(autouse=True)
def handle_numpy_1_warnings():
    """Handle NumPy 1.x DeprecationWarnings for out-of-bound integer conversions.

    NumPy 1.x raises DeprecationWarning when converting out-of-bounds values
    (e.g., 255 to int8), while NumPy 2.x raises OverflowError. This fixture
    suppresses the warning in NumPy 1.x environments to allow tests to pass.
    """
    # Only apply for NumPy < 2.0
    if np.__version__.startswith("1."):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                "NumPy will stop allowing conversion of out-of-bound Python integers",
                DeprecationWarning,
            )
            yield
    else:
        yield


@pytest.fixture(params=["numpy", pytest.param("dask", marks=requires_dask)])
def backend(request):
    return request.param


@pytest.fixture(params=["numbagg", "bottleneck", None])
def compute_backend(request):
    if request.param is None:
        options = dict(use_bottleneck=False, use_numbagg=False)
    elif request.param == "bottleneck":
        options = dict(use_bottleneck=True, use_numbagg=False)
    elif request.param == "numbagg":
        options = dict(use_bottleneck=False, use_numbagg=True)
    else:
        raise ValueError

    with xr.set_options(**options):
        yield request.param


@pytest.fixture(params=[1])
def ds(request, backend):
    if request.param == 1:
        ds = Dataset(
            dict(
                z1=(["y", "x"], np.random.randn(2, 8)),
                z2=(["time", "y"], np.random.randn(10, 2)),
            ),
            dict(
                x=("x", np.linspace(0, 1.0, 8)),
                time=("time", np.linspace(0, 1.0, 10)),
                c=("y", ["a", "b"]),
                y=range(2),
            ),
        )
    elif request.param == 2:
        ds = Dataset(
            dict(
                z1=(["time", "y"], np.random.randn(10, 2)),
                z2=(["time"], np.random.randn(10)),
                z3=(["x", "time"], np.random.randn(8, 10)),
            ),
            dict(
                x=("x", np.linspace(0, 1.0, 8)),
                time=("time", np.linspace(0, 1.0, 10)),
                c=("y", ["a", "b"]),
                y=range(2),
            ),
        )
    elif request.param == 3:
        ds = create_test_data()
    else:
        raise ValueError

    if backend == "dask":
        return ds.chunk()

    return ds


@pytest.fixture(params=[1])
def da(request, backend):
    if request.param == 1:
        times = pd.date_range("2000-01-01", freq="1D", periods=21)
        da = DataArray(
            np.random.random((3, 21, 4)),
            dims=("a", "time", "x"),
            coords=dict(time=times),
        )

    if request.param == 2:
        da = DataArray([0, np.nan, 1, 2, np.nan, 3, 4, 5, np.nan, 6, 7], dims="time")

    if request.param == "repeating_ints":
        da = DataArray(
            np.tile(np.arange(12), 5).reshape(5, 4, 3),
            coords={"x": list("abc"), "y": list("defg")},
            dims=list("zyx"),
        )

    if backend == "dask":
        return da.chunk()
    elif backend == "numpy":
        return da
    else:
        raise ValueError


@pytest.fixture(
    params=[
        False,
        pytest.param(
            True, marks=pytest.mark.skipif(not has_cftime, reason="no cftime")
        ),
    ]
)
def use_cftime(request):
    return request.param


@pytest.fixture(params=[Dataset, DataArray])
def type(request):
    return request.param


@pytest.fixture(params=[1])
def d(request, backend, type) -> DataArray | Dataset:
    """
    For tests which can test either a DataArray or a Dataset.
    """
    result: DataArray | Dataset
    if request.param == 1:
        ds = Dataset(
            dict(
                a=(["x", "z"], np.arange(24).reshape(2, 12)),
                b=(["y", "z"], np.arange(100, 136).reshape(3, 12).astype(np.float64)),
            ),
            dict(
                x=("x", np.linspace(0, 1.0, 2)),
                y=range(3),
                z=("z", pd.date_range("2000-01-01", periods=12)),
                w=("x", ["a", "b"]),
            ),
        )
        if type == DataArray:
            result = ds["a"].assign_coords(w=ds.coords["w"])
        elif type == Dataset:
            result = ds
        else:
            raise ValueError
    else:
        raise ValueError

    if backend == "dask":
        return result.chunk()
    elif backend == "numpy":
        return result
    else:
        raise ValueError


@pytest.fixture
def byte_attrs_dataset():
    """For testing issue #9407"""
    null_byte = b"\x00"
    other_bytes = bytes(range(1, 256))
    ds = Dataset({"x": 1}, coords={"x_coord": [1]})
    ds["x"].attrs["null_byte"] = null_byte
    ds["x"].attrs["other_bytes"] = other_bytes

    expected = ds.copy()
    expected["x"].attrs["null_byte"] = ""
    expected["x"].attrs["other_bytes"] = other_bytes.decode(errors="replace")

    return {
        "input": ds,
        "expected": expected,
        "h5netcdf_error": r"Invalid value provided for attribute .*: .*\. Null characters .*",
    }


@pytest.fixture(scope="module")
def create_test_datatree():
    """
    Create a test datatree with this structure:

    <xarray.DataTree>
    Group: /
    │   Dimensions:  (y: 3, x: 2)
    │   Dimensions without coordinates: y, x
    │   Data variables:
    │       a        (y) int64 24B 6 7 8
    │       set0     (x) int64 16B 9 10
    ├── Group: /set1
    │   │   Dimensions:  ()
    │   │   Data variables:
    │   │       a        int64 8B 0
    │   │       b        int64 8B 1
    │   ├── Group: /set1/set1
    │   └── Group: /set1/set2
    ├── Group: /set2
    │   │   Dimensions:  (x: 2)
    │   │   Dimensions without coordinates: x
    │   │   Data variables:
    │   │       a        (x) int64 16B 2 3
    │   │       b        (x) float64 16B 0.1 0.2
    │   └── Group: /set2/set1
    └── Group: /set3

    The structure has deliberately repeated names of tags, variables, and
    dimensions in order to better check for bugs caused by name conflicts.
    """

    def _create_test_datatree(modify=lambda ds: ds):
        set1_data = modify(xr.Dataset({"a": 0, "b": 1}))
        set2_data = modify(xr.Dataset({"a": ("x", [2, 3]), "b": ("x", [0.1, 0.2])}))
        root_data = modify(xr.Dataset({"a": ("y", [6, 7, 8]), "set0": ("x", [9, 10])}))

        root = DataTree.from_dict(
            {
                "/": root_data,
                "/set1": set1_data,
                "/set1/set1": None,
                "/set1/set2": None,
                "/set2": set2_data,
                "/set2/set1": None,
                "/set3": None,
            }
        )

        return root

    return _create_test_datatree


@pytest.fixture(scope="module")
def simple_datatree(create_test_datatree):
    """
    Invoke create_test_datatree fixture (callback).

    Returns a DataTree.
    """
    return create_test_datatree()


@pytest.fixture(params=["s", "ms", "us", "ns"])
def time_unit(request):
    return request.param
