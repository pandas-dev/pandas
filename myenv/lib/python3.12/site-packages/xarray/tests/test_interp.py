from __future__ import annotations

from itertools import combinations, permutations
from typing import cast

import numpy as np
import pandas as pd
import pytest

import xarray as xr
from xarray.coding.cftimeindex import _parse_array_of_cftime_strings
from xarray.core.types import InterpOptions
from xarray.tests import (
    assert_allclose,
    assert_equal,
    assert_identical,
    has_dask,
    has_scipy,
    requires_cftime,
    requires_dask,
    requires_scipy,
)
from xarray.tests.test_dataset import create_test_data

try:
    import scipy
except ImportError:
    pass


def get_example_data(case: int) -> xr.DataArray:
    if case == 0:
        # 2D
        x = np.linspace(0, 1, 100)
        y = np.linspace(0, 0.1, 30)
        return xr.DataArray(
            np.sin(x[:, np.newaxis]) * np.cos(y),
            dims=["x", "y"],
            coords={"x": x, "y": y, "x2": ("x", x**2)},
        )
    elif case == 1:
        # 2D chunked single dim
        return get_example_data(0).chunk({"y": 3})
    elif case == 2:
        # 2D chunked both dims
        return get_example_data(0).chunk({"x": 25, "y": 3})
    elif case == 3:
        # 3D
        x = np.linspace(0, 1, 100)
        y = np.linspace(0, 0.1, 30)
        z = np.linspace(0.1, 0.2, 10)
        return xr.DataArray(
            np.sin(x[:, np.newaxis, np.newaxis]) * np.cos(y[:, np.newaxis]) * z,
            dims=["x", "y", "z"],
            coords={"x": x, "y": y, "x2": ("x", x**2), "z": z},
        )
    elif case == 4:
        # 3D chunked single dim
        return get_example_data(3).chunk({"z": 5})
    else:
        raise ValueError("case must be 1-4")


def test_keywargs():
    if not has_scipy:
        pytest.skip("scipy is not installed.")

    da = get_example_data(0)
    assert_equal(da.interp(x=[0.5, 0.8]), da.interp({"x": [0.5, 0.8]}))


@pytest.mark.parametrize("method", ["linear", "cubic"])
@pytest.mark.parametrize("dim", ["x", "y"])
@pytest.mark.parametrize(
    "case", [pytest.param(0, id="no_chunk"), pytest.param(1, id="chunk_y")]
)
def test_interpolate_1d(method: InterpOptions, dim: str, case: int) -> None:
    if not has_scipy:
        pytest.skip("scipy is not installed.")

    if not has_dask and case in [1]:
        pytest.skip("dask is not installed in the environment.")

    da = get_example_data(case)
    xdest = np.linspace(0.0, 0.9, 80)
    actual = da.interp(method=method, coords={dim: xdest})

    # scipy interpolation for the reference
    def func(obj, new_x):
        return scipy.interpolate.interp1d(
            da[dim],
            obj.data,
            axis=obj.get_axis_num(dim),
            bounds_error=False,
            fill_value=np.nan,
            kind=method,
        )(new_x)

    if dim == "x":
        coords = {"x": xdest, "y": da["y"], "x2": ("x", func(da["x2"], xdest))}
    else:  # y
        coords = {"x": da["x"], "y": xdest, "x2": da["x2"]}

    expected = xr.DataArray(func(da, xdest), dims=["x", "y"], coords=coords)
    assert_allclose(actual, expected)


@pytest.mark.parametrize("method", ["cubic", "zero"])
def test_interpolate_1d_methods(method: InterpOptions) -> None:
    if not has_scipy:
        pytest.skip("scipy is not installed.")

    da = get_example_data(0)
    dim = "x"
    xdest = np.linspace(0.0, 0.9, 80)

    actual = da.interp(method=method, coords={dim: xdest})

    # scipy interpolation for the reference
    def func(obj, new_x):
        return scipy.interpolate.interp1d(
            da[dim],
            obj.data,
            axis=obj.get_axis_num(dim),
            bounds_error=False,
            fill_value=np.nan,
            kind=method,
        )(new_x)

    coords = {"x": xdest, "y": da["y"], "x2": ("x", func(da["x2"], xdest))}
    expected = xr.DataArray(func(da, xdest), dims=["x", "y"], coords=coords)
    assert_allclose(actual, expected)


@pytest.mark.parametrize("use_dask", [False, True])
def test_interpolate_vectorize(use_dask: bool) -> None:
    if not has_scipy:
        pytest.skip("scipy is not installed.")

    if not has_dask and use_dask:
        pytest.skip("dask is not installed in the environment.")

    # scipy interpolation for the reference
    def func(obj, dim, new_x):
        shape = [s for i, s in enumerate(obj.shape) if i != obj.get_axis_num(dim)]
        for s in new_x.shape[::-1]:
            shape.insert(obj.get_axis_num(dim), s)

        return scipy.interpolate.interp1d(
            da[dim],
            obj.data,
            axis=obj.get_axis_num(dim),
            bounds_error=False,
            fill_value=np.nan,
        )(new_x).reshape(shape)

    da = get_example_data(0)
    if use_dask:
        da = da.chunk({"y": 5})

    # xdest is 1d but has different dimension
    xdest = xr.DataArray(
        np.linspace(0.1, 0.9, 30),
        dims="z",
        coords={"z": np.random.randn(30), "z2": ("z", np.random.randn(30))},
    )

    actual = da.interp(x=xdest, method="linear")

    expected = xr.DataArray(
        func(da, "x", xdest),
        dims=["z", "y"],
        coords={
            "z": xdest["z"],
            "z2": xdest["z2"],
            "y": da["y"],
            "x": ("z", xdest.values),
            "x2": ("z", func(da["x2"], "x", xdest)),
        },
    )
    assert_allclose(actual, expected.transpose("z", "y", transpose_coords=True))

    # xdest is 2d
    xdest = xr.DataArray(
        np.linspace(0.1, 0.9, 30).reshape(6, 5),
        dims=["z", "w"],
        coords={
            "z": np.random.randn(6),
            "w": np.random.randn(5),
            "z2": ("z", np.random.randn(6)),
        },
    )

    actual = da.interp(x=xdest, method="linear")

    expected = xr.DataArray(
        func(da, "x", xdest),
        dims=["z", "w", "y"],
        coords={
            "z": xdest["z"],
            "w": xdest["w"],
            "z2": xdest["z2"],
            "y": da["y"],
            "x": (("z", "w"), xdest.data),
            "x2": (("z", "w"), func(da["x2"], "x", xdest)),
        },
    )
    assert_allclose(actual, expected.transpose("z", "w", "y", transpose_coords=True))


@pytest.mark.parametrize(
    "case", [pytest.param(3, id="no_chunk"), pytest.param(4, id="chunked")]
)
def test_interpolate_nd(case: int) -> None:
    if not has_scipy:
        pytest.skip("scipy is not installed.")

    if not has_dask and case == 4:
        pytest.skip("dask is not installed in the environment.")

    da = get_example_data(case)

    # grid -> grid
    xdestnp = np.linspace(0.1, 1.0, 11)
    ydestnp = np.linspace(0.0, 0.2, 10)
    actual = da.interp(x=xdestnp, y=ydestnp, method="linear")

    # linear interpolation is separateable
    expected = da.interp(x=xdestnp, method="linear")
    expected = expected.interp(y=ydestnp, method="linear")
    assert_allclose(actual.transpose("x", "y", "z"), expected.transpose("x", "y", "z"))

    # grid -> 1d-sample
    xdest = xr.DataArray(np.linspace(0.1, 1.0, 11), dims="y")
    ydest = xr.DataArray(np.linspace(0.0, 0.2, 11), dims="y")
    actual = da.interp(x=xdest, y=ydest, method="linear")

    # linear interpolation is separateable
    expected_data = scipy.interpolate.RegularGridInterpolator(
        (da["x"], da["y"]),
        da.transpose("x", "y", "z").values,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )(np.stack([xdest, ydest], axis=-1))
    expected = xr.DataArray(
        expected_data,
        dims=["y", "z"],
        coords={
            "z": da["z"],
            "y": ydest,
            "x": ("y", xdest.values),
            "x2": da["x2"].interp(x=xdest),
        },
    )
    assert_allclose(actual.transpose("y", "z"), expected)

    # reversed order
    actual = da.interp(y=ydest, x=xdest, method="linear")
    assert_allclose(actual.transpose("y", "z"), expected)


@requires_scipy
def test_interpolate_nd_nd() -> None:
    """Interpolate nd array with an nd indexer sharing coordinates."""
    # Create original array
    a = [0, 2]
    x = [0, 1, 2]
    da = xr.DataArray(
        np.arange(6).reshape(2, 3), dims=("a", "x"), coords={"a": a, "x": x}
    )

    # Create indexer into `a` with dimensions (y, x)
    y = [10]
    c = {"x": x, "y": y}
    ia = xr.DataArray([[1, 2, 2]], dims=("y", "x"), coords=c)
    out = da.interp(a=ia)
    expected = xr.DataArray([[1.5, 4, 5]], dims=("y", "x"), coords=c)
    xr.testing.assert_allclose(out.drop_vars("a"), expected)

    # If the *shared* indexing coordinates do not match, interp should fail.
    with pytest.raises(ValueError):
        c = {"x": [1], "y": y}
        ia = xr.DataArray([[1]], dims=("y", "x"), coords=c)
        da.interp(a=ia)

    with pytest.raises(ValueError):
        c = {"x": [5, 6, 7], "y": y}
        ia = xr.DataArray([[1]], dims=("y", "x"), coords=c)
        da.interp(a=ia)


@requires_scipy
def test_interpolate_nd_with_nan() -> None:
    """Interpolate an array with an nd indexer and `NaN` values."""

    # Create indexer into `a` with dimensions (y, x)
    x = [0, 1, 2]
    y = [10, 20]
    c = {"x": x, "y": y}
    a = np.arange(6, dtype=float).reshape(2, 3)
    a[0, 1] = np.nan
    ia = xr.DataArray(a, dims=("y", "x"), coords=c)

    da = xr.DataArray([1, 2, 2], dims=("a"), coords={"a": [0, 2, 4]})
    out = da.interp(a=ia)
    expected = xr.DataArray(
        [[1.0, np.nan, 2.0], [2.0, 2.0, np.nan]], dims=("y", "x"), coords=c
    )
    xr.testing.assert_allclose(out.drop_vars("a"), expected)

    db = 2 * da
    ds = xr.Dataset({"da": da, "db": db})
    out2 = ds.interp(a=ia)
    expected_ds = xr.Dataset({"da": expected, "db": 2 * expected})
    xr.testing.assert_allclose(out2.drop_vars("a"), expected_ds)


@pytest.mark.parametrize("method", ["linear"])
@pytest.mark.parametrize(
    "case", [pytest.param(0, id="no_chunk"), pytest.param(1, id="chunk_y")]
)
def test_interpolate_scalar(method: InterpOptions, case: int) -> None:
    if not has_scipy:
        pytest.skip("scipy is not installed.")

    if not has_dask and case in [1]:
        pytest.skip("dask is not installed in the environment.")

    da = get_example_data(case)
    xdest = 0.4

    actual = da.interp(x=xdest, method=method)

    # scipy interpolation for the reference
    def func(obj, new_x):
        return scipy.interpolate.interp1d(
            da["x"],
            obj.data,
            axis=obj.get_axis_num("x"),
            bounds_error=False,
            fill_value=np.nan,
        )(new_x)

    coords = {"x": xdest, "y": da["y"], "x2": func(da["x2"], xdest)}
    expected = xr.DataArray(func(da, xdest), dims=["y"], coords=coords)
    assert_allclose(actual, expected)


@pytest.mark.parametrize("method", ["linear"])
@pytest.mark.parametrize(
    "case", [pytest.param(3, id="no_chunk"), pytest.param(4, id="chunked")]
)
def test_interpolate_nd_scalar(method: InterpOptions, case: int) -> None:
    if not has_scipy:
        pytest.skip("scipy is not installed.")

    if not has_dask and case in [4]:
        pytest.skip("dask is not installed in the environment.")

    da = get_example_data(case)
    xdest = 0.4
    ydest = 0.05

    actual = da.interp(x=xdest, y=ydest, method=method)
    # scipy interpolation for the reference
    expected_data = scipy.interpolate.RegularGridInterpolator(
        (da["x"], da["y"]),
        da.transpose("x", "y", "z").values,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )(np.stack([xdest, ydest], axis=-1))

    coords = {"x": xdest, "y": ydest, "x2": da["x2"].interp(x=xdest), "z": da["z"]}
    expected = xr.DataArray(expected_data[0], dims=["z"], coords=coords)
    assert_allclose(actual, expected)


@pytest.mark.parametrize("use_dask", [True, False])
def test_nans(use_dask: bool) -> None:
    if not has_scipy:
        pytest.skip("scipy is not installed.")

    da = xr.DataArray([0, 1, np.nan, 2], dims="x", coords={"x": range(4)})

    if not has_dask and use_dask:
        pytest.skip("dask is not installed in the environment.")
        da = da.chunk()

    actual = da.interp(x=[0.5, 1.5])
    # not all values are nan
    assert actual.count() > 0


@pytest.mark.parametrize("use_dask", [True, False])
def test_errors(use_dask: bool) -> None:
    if not has_scipy:
        pytest.skip("scipy is not installed.")

    # akima and spline are unavailable
    da = xr.DataArray([0, 1, np.nan, 2], dims="x", coords={"x": range(4)})
    if not has_dask and use_dask:
        pytest.skip("dask is not installed in the environment.")
        da = da.chunk()

    for method in ["akima", "spline"]:
        with pytest.raises(ValueError):
            da.interp(x=[0.5, 1.5], method=method)  # type: ignore

    # not sorted
    if use_dask:
        da = get_example_data(3)
    else:
        da = get_example_data(0)

    result = da.interp(x=[-1, 1, 3], kwargs={"fill_value": 0.0})
    assert not np.isnan(result.values).any()
    result = da.interp(x=[-1, 1, 3])
    assert np.isnan(result.values).any()

    # invalid method
    with pytest.raises(ValueError):
        da.interp(x=[2, 0], method="boo")  # type: ignore
    with pytest.raises(ValueError):
        da.interp(y=[2, 0], method="boo")  # type: ignore

    # object-type DataArray cannot be interpolated
    da = xr.DataArray(["a", "b", "c"], dims="x", coords={"x": [0, 1, 2]})
    with pytest.raises(TypeError):
        da.interp(x=0)


@requires_scipy
def test_dtype() -> None:
    data_vars = dict(
        a=("time", np.array([1, 1.25, 2])),
        b=("time", np.array([True, True, False], dtype=bool)),
        c=("time", np.array(["start", "start", "end"], dtype=str)),
    )
    time = np.array([0, 0.25, 1], dtype=float)
    expected = xr.Dataset(data_vars, coords=dict(time=time))
    actual = xr.Dataset(
        {k: (dim, arr[[0, -1]]) for k, (dim, arr) in data_vars.items()},
        coords=dict(time=time[[0, -1]]),
    )
    actual = actual.interp(time=time, method="linear")
    assert_identical(expected, actual)


@requires_scipy
def test_sorted() -> None:
    # unsorted non-uniform gridded data
    x = np.random.randn(100)
    y = np.random.randn(30)
    z = np.linspace(0.1, 0.2, 10) * 3.0
    da = xr.DataArray(
        np.cos(x[:, np.newaxis, np.newaxis]) * np.cos(y[:, np.newaxis]) * z,
        dims=["x", "y", "z"],
        coords={"x": x, "y": y, "x2": ("x", x**2), "z": z},
    )

    x_new = np.linspace(0, 1, 30)
    y_new = np.linspace(0, 1, 20)

    da_sorted = da.sortby("x")
    assert_allclose(da.interp(x=x_new), da_sorted.interp(x=x_new, assume_sorted=True))
    da_sorted = da.sortby(["x", "y"])
    assert_allclose(
        da.interp(x=x_new, y=y_new),
        da_sorted.interp(x=x_new, y=y_new, assume_sorted=True),
    )

    with pytest.raises(ValueError):
        da.interp(x=[0, 1, 2], assume_sorted=True)


@requires_scipy
def test_dimension_wo_coords() -> None:
    da = xr.DataArray(
        np.arange(12).reshape(3, 4), dims=["x", "y"], coords={"y": [0, 1, 2, 3]}
    )
    da_w_coord = da.copy()
    da_w_coord["x"] = np.arange(3)

    assert_equal(da.interp(x=[0.1, 0.2, 0.3]), da_w_coord.interp(x=[0.1, 0.2, 0.3]))
    assert_equal(
        da.interp(x=[0.1, 0.2, 0.3], y=[0.5]),
        da_w_coord.interp(x=[0.1, 0.2, 0.3], y=[0.5]),
    )


@requires_scipy
def test_dataset() -> None:
    ds = create_test_data()
    ds.attrs["foo"] = "var"
    ds["var1"].attrs["buz"] = "var2"
    new_dim2 = xr.DataArray([0.11, 0.21, 0.31], dims="z")
    interpolated = ds.interp(dim2=new_dim2)

    assert_allclose(interpolated["var1"], ds["var1"].interp(dim2=new_dim2))
    assert interpolated["var3"].equals(ds["var3"])

    # make sure modifying interpolated does not affect the original dataset
    interpolated["var1"][:, 1] = 1.0
    interpolated["var2"][:, 1] = 1.0
    interpolated["var3"][:, 1] = 1.0

    assert not interpolated["var1"].equals(ds["var1"])
    assert not interpolated["var2"].equals(ds["var2"])
    assert not interpolated["var3"].equals(ds["var3"])
    # attrs should be kept
    assert interpolated.attrs["foo"] == "var"
    assert interpolated["var1"].attrs["buz"] == "var2"


@pytest.mark.parametrize("case", [pytest.param(0, id="2D"), pytest.param(3, id="3D")])
def test_interpolate_dimorder(case: int) -> None:
    """Make sure the resultant dimension order is consistent with .sel()"""
    if not has_scipy:
        pytest.skip("scipy is not installed.")

    da = get_example_data(case)

    new_x = xr.DataArray([0, 1, 2], dims="x")
    assert da.interp(x=new_x).dims == da.sel(x=new_x, method="nearest").dims

    new_y = xr.DataArray([0, 1, 2], dims="y")
    actual = da.interp(x=new_x, y=new_y).dims
    expected = da.sel(x=new_x, y=new_y, method="nearest").dims
    assert actual == expected
    # reversed order
    actual = da.interp(y=new_y, x=new_x).dims
    expected = da.sel(y=new_y, x=new_x, method="nearest").dims
    assert actual == expected

    new_x = xr.DataArray([0, 1, 2], dims="a")
    assert da.interp(x=new_x).dims == da.sel(x=new_x, method="nearest").dims
    assert da.interp(y=new_x).dims == da.sel(y=new_x, method="nearest").dims
    new_y = xr.DataArray([0, 1, 2], dims="a")
    actual = da.interp(x=new_x, y=new_y).dims
    expected = da.sel(x=new_x, y=new_y, method="nearest").dims
    assert actual == expected

    new_x = xr.DataArray([[0], [1], [2]], dims=["a", "b"])
    assert da.interp(x=new_x).dims == da.sel(x=new_x, method="nearest").dims
    assert da.interp(y=new_x).dims == da.sel(y=new_x, method="nearest").dims

    if case == 3:
        new_x = xr.DataArray([[0], [1], [2]], dims=["a", "b"])
        new_z = xr.DataArray([[0], [1], [2]], dims=["a", "b"])
        actual = da.interp(x=new_x, z=new_z).dims
        expected = da.sel(x=new_x, z=new_z, method="nearest").dims
        assert actual == expected

        actual = da.interp(z=new_z, x=new_x).dims
        expected = da.sel(z=new_z, x=new_x, method="nearest").dims
        assert actual == expected

        actual = da.interp(x=0.5, z=new_z).dims
        expected = da.sel(x=0.5, z=new_z, method="nearest").dims
        assert actual == expected


@requires_scipy
def test_interp_like() -> None:
    ds = create_test_data()
    ds.attrs["foo"] = "var"
    ds["var1"].attrs["buz"] = "var2"

    other = xr.DataArray(np.random.randn(3), dims=["dim2"], coords={"dim2": [0, 1, 2]})
    interpolated = ds.interp_like(other)

    assert_allclose(interpolated["var1"], ds["var1"].interp(dim2=other["dim2"]))
    assert_allclose(interpolated["var1"], ds["var1"].interp_like(other))
    assert interpolated["var3"].equals(ds["var3"])

    # attrs should be kept
    assert interpolated.attrs["foo"] == "var"
    assert interpolated["var1"].attrs["buz"] == "var2"

    other = xr.DataArray(
        np.random.randn(3), dims=["dim3"], coords={"dim3": ["a", "b", "c"]}
    )

    actual = ds.interp_like(other)
    expected = ds.reindex_like(other)
    assert_allclose(actual, expected)


@requires_scipy
@pytest.mark.parametrize(
    "x_new, expected",
    [
        (pd.date_range("2000-01-02", periods=3), [1, 2, 3]),
        (
            np.array(
                [np.datetime64("2000-01-01T12:00"), np.datetime64("2000-01-02T12:00")]
            ),
            [0.5, 1.5],
        ),
        (["2000-01-01T12:00", "2000-01-02T12:00"], [0.5, 1.5]),
        (["2000-01-01T12:00", "2000-01-02T12:00", "NaT"], [0.5, 1.5, np.nan]),
        (["2000-01-01T12:00"], 0.5),
        pytest.param("2000-01-01T12:00", 0.5, marks=pytest.mark.xfail),
    ],
)
@pytest.mark.filterwarnings("ignore:Converting non-nanosecond")
def test_datetime(x_new, expected) -> None:
    da = xr.DataArray(
        np.arange(24),
        dims="time",
        coords={"time": pd.date_range("2000-01-01", periods=24)},
    )

    actual = da.interp(time=x_new)
    expected_da = xr.DataArray(
        np.atleast_1d(expected),
        dims=["time"],
        coords={"time": (np.atleast_1d(x_new).astype("datetime64[ns]"))},
    )

    assert_allclose(actual, expected_da)


@requires_scipy
def test_datetime_single_string() -> None:
    da = xr.DataArray(
        np.arange(24),
        dims="time",
        coords={"time": pd.date_range("2000-01-01", periods=24)},
    )
    actual = da.interp(time="2000-01-01T12:00")
    expected = xr.DataArray(0.5)

    assert_allclose(actual.drop_vars("time"), expected)


@requires_cftime
@requires_scipy
def test_cftime() -> None:
    times = xr.cftime_range("2000", periods=24, freq="D")
    da = xr.DataArray(np.arange(24), coords=[times], dims="time")

    times_new = xr.cftime_range("2000-01-01T12:00:00", periods=3, freq="D")
    actual = da.interp(time=times_new)
    expected = xr.DataArray([0.5, 1.5, 2.5], coords=[times_new], dims=["time"])

    assert_allclose(actual, expected)


@requires_cftime
@requires_scipy
def test_cftime_type_error() -> None:
    times = xr.cftime_range("2000", periods=24, freq="D")
    da = xr.DataArray(np.arange(24), coords=[times], dims="time")

    times_new = xr.cftime_range(
        "2000-01-01T12:00:00", periods=3, freq="D", calendar="noleap"
    )
    with pytest.raises(TypeError):
        da.interp(time=times_new)


@requires_cftime
@requires_scipy
def test_cftime_list_of_strings() -> None:
    from cftime import DatetimeProlepticGregorian

    times = xr.cftime_range(
        "2000", periods=24, freq="D", calendar="proleptic_gregorian"
    )
    da = xr.DataArray(np.arange(24), coords=[times], dims="time")

    times_new = ["2000-01-01T12:00", "2000-01-02T12:00", "2000-01-03T12:00"]
    actual = da.interp(time=times_new)

    times_new_array = _parse_array_of_cftime_strings(
        np.array(times_new), DatetimeProlepticGregorian
    )
    expected = xr.DataArray([0.5, 1.5, 2.5], coords=[times_new_array], dims=["time"])

    assert_allclose(actual, expected)


@requires_cftime
@requires_scipy
def test_cftime_single_string() -> None:
    from cftime import DatetimeProlepticGregorian

    times = xr.cftime_range(
        "2000", periods=24, freq="D", calendar="proleptic_gregorian"
    )
    da = xr.DataArray(np.arange(24), coords=[times], dims="time")

    times_new = "2000-01-01T12:00"
    actual = da.interp(time=times_new)

    times_new_array = _parse_array_of_cftime_strings(
        np.array(times_new), DatetimeProlepticGregorian
    )
    expected = xr.DataArray(0.5, coords={"time": times_new_array})

    assert_allclose(actual, expected)


@requires_scipy
def test_datetime_to_non_datetime_error() -> None:
    da = xr.DataArray(
        np.arange(24),
        dims="time",
        coords={"time": pd.date_range("2000-01-01", periods=24)},
    )
    with pytest.raises(TypeError):
        da.interp(time=0.5)


@requires_cftime
@requires_scipy
def test_cftime_to_non_cftime_error() -> None:
    times = xr.cftime_range("2000", periods=24, freq="D")
    da = xr.DataArray(np.arange(24), coords=[times], dims="time")

    with pytest.raises(TypeError):
        da.interp(time=0.5)


@requires_scipy
def test_datetime_interp_noerror() -> None:
    # GH:2667
    a = xr.DataArray(
        np.arange(21).reshape(3, 7),
        dims=["x", "time"],
        coords={
            "x": [1, 2, 3],
            "time": pd.date_range("01-01-2001", periods=7, freq="D"),
        },
    )
    xi = xr.DataArray(
        np.linspace(1, 3, 50),
        dims=["time"],
        coords={"time": pd.date_range("01-01-2001", periods=50, freq="h")},
    )
    a.interp(x=xi, time=xi.time)  # should not raise an error


@requires_cftime
@requires_scipy
def test_3641() -> None:
    times = xr.cftime_range("0001", periods=3, freq="500YE")
    da = xr.DataArray(range(3), dims=["time"], coords=[times])
    da.interp(time=["0002-05-01"])


@requires_scipy
@pytest.mark.parametrize("method", ["nearest", "linear"])
def test_decompose(method: InterpOptions) -> None:
    da = xr.DataArray(
        np.arange(6).reshape(3, 2),
        dims=["x", "y"],
        coords={"x": [0, 1, 2], "y": [-0.1, -0.3]},
    )
    x_new = xr.DataArray([0.5, 1.5, 2.5], dims=["x1"])
    y_new = xr.DataArray([-0.15, -0.25], dims=["y1"])
    x_broadcast, y_broadcast = xr.broadcast(x_new, y_new)
    assert x_broadcast.ndim == 2

    actual = da.interp(x=x_new, y=y_new, method=method).drop_vars(("x", "y"))
    expected = da.interp(x=x_broadcast, y=y_broadcast, method=method).drop_vars(
        ("x", "y")
    )
    assert_allclose(actual, expected)


@requires_scipy
@requires_dask
@pytest.mark.parametrize(
    "method", ["linear", "nearest", "zero", "slinear", "quadratic", "cubic"]
)
@pytest.mark.parametrize("chunked", [True, False])
@pytest.mark.parametrize(
    "data_ndim,interp_ndim,nscalar",
    [
        (data_ndim, interp_ndim, nscalar)
        for data_ndim in range(1, 4)
        for interp_ndim in range(1, data_ndim + 1)
        for nscalar in range(0, interp_ndim + 1)
    ],
)
def test_interpolate_chunk_1d(
    method: InterpOptions, data_ndim, interp_ndim, nscalar, chunked: bool
) -> None:
    """Interpolate nd array with multiple independent indexers

    It should do a series of 1d interpolation
    """

    # 3d non chunked data
    x = np.linspace(0, 1, 5)
    y = np.linspace(2, 4, 7)
    z = np.linspace(-0.5, 0.5, 11)
    da = xr.DataArray(
        data=np.sin(x[:, np.newaxis, np.newaxis])
        * np.cos(y[:, np.newaxis])
        * np.exp(z),
        coords=[("x", x), ("y", y), ("z", z)],
    )
    kwargs = {"fill_value": "extrapolate"}

    # choose the data dimensions
    for data_dims in permutations(da.dims, data_ndim):
        # select only data_ndim dim
        da = da.isel(  # take the middle line
            {dim: len(da.coords[dim]) // 2 for dim in da.dims if dim not in data_dims}
        )

        # chunk data
        da = da.chunk(chunks={dim: i + 1 for i, dim in enumerate(da.dims)})

        # choose the interpolation dimensions
        for interp_dims in permutations(da.dims, interp_ndim):
            # choose the scalar interpolation dimensions
            for scalar_dims in combinations(interp_dims, nscalar):
                dest = {}
                for dim in interp_dims:
                    if dim in scalar_dims:
                        # take the middle point
                        dest[dim] = 0.5 * (da.coords[dim][0] + da.coords[dim][-1])
                    else:
                        # pick some points, including outside the domain
                        before = 2 * da.coords[dim][0] - da.coords[dim][1]
                        after = 2 * da.coords[dim][-1] - da.coords[dim][-2]

                        dest[dim] = cast(
                            xr.DataArray,
                            np.linspace(
                                before.item(), after.item(), len(da.coords[dim]) * 13
                            ),
                        )
                        if chunked:
                            dest[dim] = xr.DataArray(data=dest[dim], dims=[dim])
                            dest[dim] = dest[dim].chunk(2)
                actual = da.interp(method=method, **dest, kwargs=kwargs)
                expected = da.compute().interp(method=method, **dest, kwargs=kwargs)

                assert_identical(actual, expected)

                # all the combinations are usually not necessary
                break
            break
        break


@requires_scipy
@requires_dask
@pytest.mark.parametrize("method", ["linear", "nearest"])
@pytest.mark.filterwarnings("ignore:Increasing number of chunks")
def test_interpolate_chunk_advanced(method: InterpOptions) -> None:
    """Interpolate nd array with an nd indexer sharing coordinates."""
    # Create original array
    x = np.linspace(-1, 1, 5)
    y = np.linspace(-1, 1, 7)
    z = np.linspace(-1, 1, 11)
    t = np.linspace(0, 1, 13)
    q = np.linspace(0, 1, 17)
    da = xr.DataArray(
        data=np.sin(x[:, np.newaxis, np.newaxis, np.newaxis, np.newaxis])
        * np.cos(y[:, np.newaxis, np.newaxis, np.newaxis])
        * np.exp(z[:, np.newaxis, np.newaxis])
        * t[:, np.newaxis]
        + q,
        dims=("x", "y", "z", "t", "q"),
        coords={"x": x, "y": y, "z": z, "t": t, "q": q, "label": "dummy_attr"},
    )

    # Create indexer into `da` with shared coordinate ("full-twist" Möbius strip)
    theta = np.linspace(0, 2 * np.pi, 5)
    w = np.linspace(-0.25, 0.25, 7)
    r = xr.DataArray(
        data=1 + w[:, np.newaxis] * np.cos(theta),
        coords=[("w", w), ("theta", theta)],
    )

    xda = r * np.cos(theta)
    yda = r * np.sin(theta)
    zda = xr.DataArray(
        data=w[:, np.newaxis] * np.sin(theta),
        coords=[("w", w), ("theta", theta)],
    )

    kwargs = {"fill_value": None}
    expected = da.interp(t=0.5, x=xda, y=yda, z=zda, kwargs=kwargs, method=method)

    da = da.chunk(2)
    xda = xda.chunk(1)
    zda = zda.chunk(3)
    actual = da.interp(t=0.5, x=xda, y=yda, z=zda, kwargs=kwargs, method=method)
    assert_identical(actual, expected)


@requires_scipy
def test_interp1d_bounds_error() -> None:
    """Ensure exception on bounds error is raised if requested"""
    da = xr.DataArray(
        np.sin(0.3 * np.arange(4)),
        [("time", np.arange(4))],
    )

    with pytest.raises(ValueError):
        da.interp(time=3.5, kwargs=dict(bounds_error=True))

    # default is to fill with nans, so this should pass
    da.interp(time=3.5)


@requires_scipy
@pytest.mark.parametrize(
    "x, expect_same_attrs",
    [
        (2.5, True),
        (np.array([2.5, 5]), True),
        (("x", np.array([0, 0.5, 1, 2]), dict(unit="s")), False),
    ],
)
def test_coord_attrs(x, expect_same_attrs: bool) -> None:
    base_attrs = dict(foo="bar")
    ds = xr.Dataset(
        data_vars=dict(a=2 * np.arange(5)),
        coords={"x": ("x", np.arange(5), base_attrs)},
    )

    has_same_attrs = ds.interp(x=x).x.attrs == base_attrs
    assert expect_same_attrs == has_same_attrs


@requires_scipy
def test_interp1d_complex_out_of_bounds() -> None:
    """Ensure complex nans are used by default"""
    da = xr.DataArray(
        np.exp(0.3j * np.arange(4)),
        [("time", np.arange(4))],
    )

    expected = da.interp(time=3.5, kwargs=dict(fill_value=np.nan + np.nan * 1j))
    actual = da.interp(time=3.5)
    assert_identical(actual, expected)
