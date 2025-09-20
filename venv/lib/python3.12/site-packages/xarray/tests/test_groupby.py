from __future__ import annotations

import datetime
import operator
import warnings
from itertools import pairwise
from typing import Literal, cast
from unittest import mock

import numpy as np
import pandas as pd
import pytest
from packaging.version import Version

import xarray as xr
from xarray import DataArray, Dataset, Variable, date_range
from xarray.core.groupby import _consolidate_slices
from xarray.core.types import InterpOptions, ResampleCompatible
from xarray.groupers import (
    BinGrouper,
    EncodedGroups,
    Grouper,
    SeasonGrouper,
    SeasonResampler,
    TimeResampler,
    UniqueGrouper,
    season_to_month_tuple,
)
from xarray.namedarray.pycompat import is_chunked_array
from xarray.structure.alignment import broadcast
from xarray.tests import (
    _ALL_CALENDARS,
    InaccessibleArray,
    assert_allclose,
    assert_equal,
    assert_identical,
    create_test_data,
    has_cftime,
    has_dask,
    has_dask_ge_2024_08_1,
    has_flox,
    has_pandas_ge_2_2,
    raise_if_dask_computes,
    requires_cftime,
    requires_dask,
    requires_dask_ge_2024_08_1,
    requires_flox,
    requires_flox_0_9_12,
    requires_pandas_ge_2_2,
    requires_scipy,
)


@pytest.fixture
def dataset() -> xr.Dataset:
    ds = xr.Dataset(
        {
            "foo": (("x", "y", "z"), np.random.randn(3, 4, 2)),
            "baz": ("x", ["e", "f", "g"]),
            "cat": ("y", pd.Categorical(["cat1", "cat2", "cat2", "cat1"])),
        },
        {"x": ("x", ["a", "b", "c"], {"name": "x"}), "y": [1, 2, 3, 4], "z": [1, 2]},
    )
    ds["boo"] = (("z", "y"), [["f", "g", "h", "j"]] * 2)

    return ds


@pytest.fixture
def array(dataset) -> xr.DataArray:
    return dataset["foo"]


def test_consolidate_slices() -> None:
    assert _consolidate_slices([slice(3), slice(3, 5)]) == [slice(5)]
    assert _consolidate_slices([slice(2, 3), slice(3, 6)]) == [slice(2, 6)]
    assert _consolidate_slices([slice(2, 3, 1), slice(3, 6, 1)]) == [slice(2, 6, 1)]

    slices = [slice(2, 3), slice(5, 6)]
    assert _consolidate_slices(slices) == slices

    # ignore type because we're checking for an error anyway
    with pytest.raises(ValueError):
        _consolidate_slices([slice(3), 4])  # type: ignore[list-item]


@pytest.mark.filterwarnings("ignore:return type")
def test_groupby_dims_property(dataset) -> None:
    with pytest.warns(FutureWarning, match="The return type of"):
        assert dataset.groupby("x").dims == dataset.isel(x=[1]).dims
    with pytest.warns(FutureWarning, match="The return type of"):
        assert dataset.groupby("y").dims == dataset.isel(y=[1]).dims

    assert tuple(dataset.groupby("x").dims) == tuple(dataset.isel(x=slice(1, 2)).dims)
    assert tuple(dataset.groupby("y").dims) == tuple(dataset.isel(y=slice(1, 2)).dims)

    dataset = dataset.drop_vars(["cat"])
    stacked = dataset.stack({"xy": ("x", "y")})
    assert tuple(stacked.groupby("xy").dims) == tuple(stacked.isel(xy=[0]).dims)


def test_groupby_sizes_property(dataset) -> None:
    assert dataset.groupby("x").sizes == dataset.isel(x=[1]).sizes
    assert dataset.groupby("y").sizes == dataset.isel(y=[1]).sizes
    dataset = dataset.drop_vars("cat")
    stacked = dataset.stack({"xy": ("x", "y")})
    assert stacked.groupby("xy").sizes == stacked.isel(xy=[0]).sizes


def test_multi_index_groupby_map(dataset) -> None:
    # regression test for GH873
    ds = dataset.isel(z=1, drop=True)[["foo"]]
    expected = 2 * ds
    actual = (
        ds.stack(space=["x", "y"])
        .groupby("space")
        .map(lambda x: 2 * x)
        .unstack("space")
    )
    assert_equal(expected, actual)


@pytest.mark.parametrize("grouper", [dict(group="x"), dict(x=UniqueGrouper())])
def test_reduce_numeric_only(dataset, grouper: dict) -> None:
    gb = dataset.groupby(**grouper)
    with xr.set_options(use_flox=False):
        expected = gb.sum()
    with xr.set_options(use_flox=True):
        actual = gb.sum()
    assert_identical(expected, actual)


def test_multi_index_groupby_sum() -> None:
    # regression test for GH873
    ds = xr.Dataset(
        {"foo": (("x", "y", "z"), np.ones((3, 4, 2)))},
        {"x": ["a", "b", "c"], "y": [1, 2, 3, 4]},
    )
    expected = ds.sum("z")
    actual = ds.stack(space=["x", "y"]).groupby("space").sum("z").unstack("space")
    assert_equal(expected, actual)

    with pytest.raises(NotImplementedError):
        actual = (
            ds.stack(space=["x", "y"])
            .groupby(space=UniqueGrouper(), z=UniqueGrouper())
            .sum("z")
            .unstack("space")
        )
        assert_equal(expected, ds)

    if not has_pandas_ge_2_2:
        # the next line triggers a mysterious multiindex error on pandas 2.0
        return

    actual = ds.stack(space=["x", "y"]).groupby("space").sum(...).unstack("space")
    assert_equal(expected, actual)


@requires_pandas_ge_2_2
def test_multi_index_propagation() -> None:
    # regression test for GH9648
    times = pd.date_range("2023-01-01", periods=4)
    locations = ["A", "B"]
    data = [[0.5, 0.7], [0.6, 0.5], [0.4, 0.6], [0.4, 0.9]]

    da = xr.DataArray(
        data, dims=["time", "location"], coords={"time": times, "location": locations}
    )
    da = da.stack(multiindex=["time", "location"])
    grouped = da.groupby("multiindex")

    with xr.set_options(use_flox=True):
        actual = grouped.sum()
    with xr.set_options(use_flox=False):
        expected = grouped.first()
    assert_identical(actual, expected)


def test_groupby_da_datetime() -> None:
    # test groupby with a DataArray of dtype datetime for GH1132
    # create test data
    times = pd.date_range("2000-01-01", periods=4)
    foo = xr.DataArray([1, 2, 3, 4], coords=dict(time=times), dims="time")
    # create test index
    reference_dates = [times[0], times[2]]
    labels = reference_dates[0:1] * 2 + reference_dates[1:2] * 2
    ind = xr.DataArray(
        labels, coords=dict(time=times), dims="time", name="reference_date"
    )
    g = foo.groupby(ind)
    actual = g.sum(dim="time")
    expected = xr.DataArray(
        [3, 7], coords=dict(reference_date=reference_dates), dims="reference_date"
    )
    assert_equal(expected, actual)


def test_groupby_duplicate_coordinate_labels() -> None:
    # fix for https://stackoverflow.com/questions/38065129
    array = xr.DataArray([1, 2, 3], [("x", [1, 1, 2])])
    expected = xr.DataArray([3, 3], [("x", [1, 2])])
    actual = array.groupby("x").sum()
    assert_equal(expected, actual)


def test_groupby_input_mutation() -> None:
    # regression test for GH2153
    array = xr.DataArray([1, 2, 3], [("x", [2, 2, 1])])
    array_copy = array.copy()
    expected = xr.DataArray([3, 3], [("x", [1, 2])])
    actual = array.groupby("x").sum()
    assert_identical(expected, actual)
    assert_identical(array, array_copy)  # should not modify inputs


@pytest.mark.parametrize("use_flox", [True, False])
def test_groupby_indexvariable(use_flox: bool) -> None:
    # regression test for GH7919
    array = xr.DataArray([1, 2, 3], [("x", [2, 2, 1])])
    iv = xr.IndexVariable(dims="x", data=pd.Index(array.x.values))
    with xr.set_options(use_flox=use_flox):
        actual = array.groupby(iv).sum()
    actual = array.groupby(iv).sum()
    expected = xr.DataArray([3, 3], [("x", [1, 2])])
    assert_identical(expected, actual)


@pytest.mark.parametrize(
    "obj",
    [
        xr.DataArray([1, 2, 3, 4, 5, 6], [("x", [1, 1, 1, 2, 2, 2])]),
        xr.Dataset({"foo": ("x", [1, 2, 3, 4, 5, 6])}, {"x": [1, 1, 1, 2, 2, 2]}),
    ],
)
def test_groupby_map_shrink_groups(obj) -> None:
    expected = obj.isel(x=[0, 1, 3, 4])
    actual = obj.groupby("x").map(lambda f: f.isel(x=[0, 1]))
    assert_identical(expected, actual)


@pytest.mark.parametrize(
    "obj",
    [
        xr.DataArray([1, 2, 3], [("x", [1, 2, 2])]),
        xr.Dataset({"foo": ("x", [1, 2, 3])}, {"x": [1, 2, 2]}),
    ],
)
def test_groupby_map_change_group_size(obj) -> None:
    def func(group):
        if group.sizes["x"] == 1:
            result = group.isel(x=[0, 0])
        else:
            result = group.isel(x=[0])
        return result

    expected = obj.isel(x=[0, 0, 1])
    actual = obj.groupby("x").map(func)
    assert_identical(expected, actual)


def test_da_groupby_map_func_args() -> None:
    def func(arg1, arg2, arg3=0):
        return arg1 + arg2 + arg3

    array = xr.DataArray([1, 1, 1], [("x", [1, 2, 3])])
    expected = xr.DataArray([3, 3, 3], [("x", [1, 2, 3])])
    actual = array.groupby("x").map(func, args=(1,), arg3=1)
    assert_identical(expected, actual)


def test_ds_groupby_map_func_args() -> None:
    def func(arg1, arg2, arg3=0):
        return arg1 + arg2 + arg3

    dataset = xr.Dataset({"foo": ("x", [1, 1, 1])}, {"x": [1, 2, 3]})
    expected = xr.Dataset({"foo": ("x", [3, 3, 3])}, {"x": [1, 2, 3]})
    actual = dataset.groupby("x").map(func, args=(1,), arg3=1)
    assert_identical(expected, actual)


def test_da_groupby_empty() -> None:
    empty_array = xr.DataArray([], dims="dim")

    with pytest.raises(ValueError):
        empty_array.groupby("dim")


@requires_dask
def test_dask_da_groupby_quantile() -> None:
    # Scalar quantile
    expected = xr.DataArray(
        data=[2, 5], coords={"x": [1, 2], "quantile": 0.5}, dims="x"
    )
    array = xr.DataArray(
        data=[1, 2, 3, 4, 5, 6], coords={"x": [1, 1, 1, 2, 2, 2]}, dims="x"
    )

    # will work blockwise with flox
    actual = array.chunk(x=3).groupby("x").quantile(0.5)
    assert_identical(expected, actual)

    # will work blockwise with flox
    actual = array.chunk(x=-1).groupby("x").quantile(0.5)
    assert_identical(expected, actual)


@requires_dask
def test_dask_da_groupby_median() -> None:
    expected = xr.DataArray(data=[2, 5], coords={"x": [1, 2]}, dims="x")
    array = xr.DataArray(
        data=[1, 2, 3, 4, 5, 6], coords={"x": [1, 1, 1, 2, 2, 2]}, dims="x"
    )
    with xr.set_options(use_flox=False):
        actual = array.chunk(x=1).groupby("x").median()
    assert_identical(expected, actual)

    with xr.set_options(use_flox=True):
        actual = array.chunk(x=1).groupby("x").median()
    assert_identical(expected, actual)

    # will work blockwise with flox
    actual = array.chunk(x=3).groupby("x").median()
    assert_identical(expected, actual)

    # will work blockwise with flox
    actual = array.chunk(x=-1).groupby("x").median()
    assert_identical(expected, actual)


@pytest.mark.parametrize("use_flox", [pytest.param(True, marks=requires_flox), False])
def test_da_groupby_quantile(use_flox: bool) -> None:
    array = xr.DataArray(
        data=[1, 2, 3, 4, 5, 6], coords={"x": [1, 1, 1, 2, 2, 2]}, dims="x"
    )

    # Scalar quantile
    expected = xr.DataArray(
        data=[2, 5], coords={"x": [1, 2], "quantile": 0.5}, dims="x"
    )

    with xr.set_options(use_flox=use_flox):
        actual = array.groupby("x").quantile(0.5)
        assert_identical(expected, actual)

    # Vector quantile
    expected = xr.DataArray(
        data=[[1, 3], [4, 6]],
        coords={"x": [1, 2], "quantile": [0, 1]},
        dims=("x", "quantile"),
    )
    with xr.set_options(use_flox=use_flox):
        actual = array.groupby("x").quantile([0, 1])
    assert_identical(expected, actual)

    array = xr.DataArray(
        data=[np.nan, 2, 3, 4, 5, 6], coords={"x": [1, 1, 1, 2, 2, 2]}, dims="x"
    )

    for skipna in (True, False, None):
        e = [np.nan, 5] if skipna is False else [2.5, 5]

        expected = xr.DataArray(data=e, coords={"x": [1, 2], "quantile": 0.5}, dims="x")
        with xr.set_options(use_flox=use_flox):
            actual = array.groupby("x").quantile(0.5, skipna=skipna)
        assert_identical(expected, actual)

    # Multiple dimensions
    array = xr.DataArray(
        data=[[1, 11, 26], [2, 12, 22], [3, 13, 23], [4, 16, 24], [5, 15, 25]],
        coords={"x": [1, 1, 1, 2, 2], "y": [0, 0, 1]},
        dims=("x", "y"),
    )

    actual_x = array.groupby("x").quantile(0, dim=...)
    expected_x = xr.DataArray(
        data=[1, 4], coords={"x": [1, 2], "quantile": 0}, dims="x"
    )
    assert_identical(expected_x, actual_x)

    actual_y = array.groupby("y").quantile(0, dim=...)
    expected_y = xr.DataArray(
        data=[1, 22], coords={"y": [0, 1], "quantile": 0}, dims="y"
    )
    assert_identical(expected_y, actual_y)

    actual_xx = array.groupby("x").quantile(0)
    expected_xx = xr.DataArray(
        data=[[1, 11, 22], [4, 15, 24]],
        coords={"x": [1, 2], "y": [0, 0, 1], "quantile": 0},
        dims=("x", "y"),
    )
    assert_identical(expected_xx, actual_xx)

    actual_yy = array.groupby("y").quantile(0)
    expected_yy = xr.DataArray(
        data=[[1, 26], [2, 22], [3, 23], [4, 24], [5, 25]],
        coords={"x": [1, 1, 1, 2, 2], "y": [0, 1], "quantile": 0},
        dims=("x", "y"),
    )
    assert_identical(expected_yy, actual_yy)

    times = pd.date_range("2000-01-01", periods=365)
    x = [0, 1]
    foo = xr.DataArray(
        np.reshape(np.arange(365 * 2), (365, 2)),
        coords={"time": times, "x": x},
        dims=("time", "x"),
    )
    g = foo.groupby(foo.time.dt.month)

    actual = g.quantile(0, dim=...)
    expected = xr.DataArray(
        data=[
            0.0,
            62.0,
            120.0,
            182.0,
            242.0,
            304.0,
            364.0,
            426.0,
            488.0,
            548.0,
            610.0,
            670.0,
        ],
        coords={"month": np.arange(1, 13), "quantile": 0},
        dims="month",
    )
    assert_identical(expected, actual)

    actual = g.quantile(0, dim="time")[:2]
    expected = xr.DataArray(
        data=[[0.0, 1], [62.0, 63]],
        coords={"month": [1, 2], "x": [0, 1], "quantile": 0},
        dims=("month", "x"),
    )
    assert_identical(expected, actual)

    # method keyword
    array = xr.DataArray(data=[1, 2, 3, 4], coords={"x": [1, 1, 2, 2]}, dims="x")

    expected = xr.DataArray(
        data=[1, 3], coords={"x": [1, 2], "quantile": 0.5}, dims="x"
    )
    actual = array.groupby("x").quantile(0.5, method="lower")
    assert_identical(expected, actual)


def test_ds_groupby_quantile() -> None:
    ds = xr.Dataset(
        data_vars={"a": ("x", [1, 2, 3, 4, 5, 6])}, coords={"x": [1, 1, 1, 2, 2, 2]}
    )

    # Scalar quantile
    expected = xr.Dataset(
        data_vars={"a": ("x", [2, 5])}, coords={"quantile": 0.5, "x": [1, 2]}
    )
    actual = ds.groupby("x").quantile(0.5)
    assert_identical(expected, actual)

    # Vector quantile
    expected = xr.Dataset(
        data_vars={"a": (("x", "quantile"), [[1, 3], [4, 6]])},
        coords={"x": [1, 2], "quantile": [0, 1]},
    )
    actual = ds.groupby("x").quantile([0, 1])
    assert_identical(expected, actual)

    ds = xr.Dataset(
        data_vars={"a": ("x", [np.nan, 2, 3, 4, 5, 6])},
        coords={"x": [1, 1, 1, 2, 2, 2]},
    )

    for skipna in (True, False, None):
        e = [np.nan, 5] if skipna is False else [2.5, 5]

        expected = xr.Dataset(
            data_vars={"a": ("x", e)}, coords={"quantile": 0.5, "x": [1, 2]}
        )
        actual = ds.groupby("x").quantile(0.5, skipna=skipna)
        assert_identical(expected, actual)

    # Multiple dimensions
    ds = xr.Dataset(
        data_vars={
            "a": (
                ("x", "y"),
                [[1, 11, 26], [2, 12, 22], [3, 13, 23], [4, 16, 24], [5, 15, 25]],
            )
        },
        coords={"x": [1, 1, 1, 2, 2], "y": [0, 0, 1]},
    )

    actual_x = ds.groupby("x").quantile(0, dim=...)
    expected_x = xr.Dataset({"a": ("x", [1, 4])}, coords={"x": [1, 2], "quantile": 0})
    assert_identical(expected_x, actual_x)

    actual_y = ds.groupby("y").quantile(0, dim=...)
    expected_y = xr.Dataset({"a": ("y", [1, 22])}, coords={"y": [0, 1], "quantile": 0})
    assert_identical(expected_y, actual_y)

    actual_xx = ds.groupby("x").quantile(0)
    expected_xx = xr.Dataset(
        {"a": (("x", "y"), [[1, 11, 22], [4, 15, 24]])},
        coords={"x": [1, 2], "y": [0, 0, 1], "quantile": 0},
    )
    assert_identical(expected_xx, actual_xx)

    actual_yy = ds.groupby("y").quantile(0)
    expected_yy = xr.Dataset(
        {"a": (("x", "y"), [[1, 26], [2, 22], [3, 23], [4, 24], [5, 25]])},
        coords={"x": [1, 1, 1, 2, 2], "y": [0, 1], "quantile": 0},
    ).transpose()
    assert_identical(expected_yy, actual_yy)

    times = pd.date_range("2000-01-01", periods=365)
    x = [0, 1]
    foo = xr.Dataset(
        {"a": (("time", "x"), np.reshape(np.arange(365 * 2), (365, 2)))},
        coords=dict(time=times, x=x),
    )
    g = foo.groupby(foo.time.dt.month)

    actual = g.quantile(0, dim=...)
    expected = xr.Dataset(
        {
            "a": (
                "month",
                [
                    0.0,
                    62.0,
                    120.0,
                    182.0,
                    242.0,
                    304.0,
                    364.0,
                    426.0,
                    488.0,
                    548.0,
                    610.0,
                    670.0,
                ],
            )
        },
        coords={"month": np.arange(1, 13), "quantile": 0},
    )
    assert_identical(expected, actual)

    actual = g.quantile(0, dim="time").isel(month=slice(None, 2))
    expected = xr.Dataset(
        data_vars={"a": (("month", "x"), [[0.0, 1], [62.0, 63]])},
        coords={"month": [1, 2], "x": [0, 1], "quantile": 0},
    )
    assert_identical(expected, actual)

    ds = xr.Dataset(data_vars={"a": ("x", [1, 2, 3, 4])}, coords={"x": [1, 1, 2, 2]})

    # method keyword
    expected = xr.Dataset(
        data_vars={"a": ("x", [1, 3])}, coords={"quantile": 0.5, "x": [1, 2]}
    )
    actual = ds.groupby("x").quantile(0.5, method="lower")
    assert_identical(expected, actual)


@pytest.mark.filterwarnings(
    "default:The `interpolation` argument to quantile was renamed to `method`:FutureWarning"
)
@pytest.mark.parametrize("as_dataset", [False, True])
def test_groupby_quantile_interpolation_deprecated(as_dataset: bool) -> None:
    array = xr.DataArray(data=[1, 2, 3, 4], coords={"x": [1, 1, 2, 2]}, dims="x")

    arr: xr.DataArray | xr.Dataset
    arr = array.to_dataset(name="name") if as_dataset else array

    with pytest.warns(
        FutureWarning,
        match="`interpolation` argument to quantile was renamed to `method`",
    ):
        actual = arr.quantile(0.5, interpolation="lower")

    expected = arr.quantile(0.5, method="lower")

    assert_identical(actual, expected)

    with warnings.catch_warnings(record=True):
        with pytest.raises(TypeError, match="interpolation and method keywords"):
            arr.quantile(0.5, method="lower", interpolation="lower")


def test_da_groupby_assign_coords() -> None:
    actual = xr.DataArray(
        [[3, 4, 5], [6, 7, 8]], dims=["y", "x"], coords={"y": range(2), "x": range(3)}
    )
    actual1 = actual.groupby("x").assign_coords({"y": [-1, -2]})
    actual2 = actual.groupby("x").assign_coords(y=[-1, -2])
    expected = xr.DataArray(
        [[3, 4, 5], [6, 7, 8]], dims=["y", "x"], coords={"y": [-1, -2], "x": range(3)}
    )
    assert_identical(expected, actual1)
    assert_identical(expected, actual2)


repr_da = xr.DataArray(
    np.random.randn(10, 20, 6, 24),
    dims=["x", "y", "z", "t"],
    coords={
        "z": ["a", "b", "c", "a", "b", "c"],
        "x": [1, 1, 1, 2, 2, 3, 4, 5, 3, 4],
        "t": xr.date_range("2001-01-01", freq="ME", periods=24, use_cftime=False),
        "month": ("t", list(range(1, 13)) * 2),
    },
)


@pytest.mark.parametrize("dim", ["x", "y", "z", "month"])
@pytest.mark.parametrize("obj", [repr_da, repr_da.to_dataset(name="a")])
def test_groupby_repr(obj, dim) -> None:
    actual = repr(obj.groupby(dim))
    N = len(np.unique(obj[dim]))
    expected = f"<{obj.__class__.__name__}GroupBy"
    expected += f", grouped over 1 grouper(s), {N} groups in total:"
    expected += f"\n    {dim!r}: UniqueGrouper({dim!r}), {N}/{N} groups with labels "
    if dim == "x":
        expected += "1, 2, 3, 4, 5>"
    elif dim == "y":
        expected += "0, 1, 2, 3, 4, 5, ..., 15, 16, 17, 18, 19>"
    elif dim == "z":
        expected += "'a', 'b', 'c'>"
    elif dim == "month":
        expected += "1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12>"
    assert actual == expected


@pytest.mark.parametrize("obj", [repr_da, repr_da.to_dataset(name="a")])
def test_groupby_repr_datetime(obj) -> None:
    actual = repr(obj.groupby("t.month"))
    expected = f"<{obj.__class__.__name__}GroupBy"
    expected += ", grouped over 1 grouper(s), 12 groups in total:\n"
    expected += "    'month': UniqueGrouper('month'), 12/12 groups with labels "
    expected += "1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12>"
    assert actual == expected


@pytest.mark.filterwarnings("ignore:No index created for dimension id:UserWarning")
@pytest.mark.filterwarnings("ignore:invalid value encountered in divide:RuntimeWarning")
@pytest.mark.parametrize("shuffle", [True, False])
@pytest.mark.parametrize(
    "chunk",
    [
        pytest.param(
            dict(lat=1), marks=pytest.mark.skipif(not has_dask, reason="no dask")
        ),
        pytest.param(
            dict(lat=2, lon=2), marks=pytest.mark.skipif(not has_dask, reason="no dask")
        ),
        False,
    ],
)
def test_groupby_drops_nans(shuffle: bool, chunk: Literal[False] | dict) -> None:
    if shuffle and chunk and not has_dask_ge_2024_08_1:
        pytest.skip()
    # GH2383
    # nan in 2D data variable (requires stacking)
    ds = xr.Dataset(
        {
            "variable": (("lat", "lon", "time"), np.arange(60.0).reshape((4, 3, 5))),
            "id": (("lat", "lon"), np.arange(12.0).reshape((4, 3))),
        },
        coords={"lat": np.arange(4), "lon": np.arange(3), "time": np.arange(5)},
    )

    ds["id"].values[0, 0] = np.nan
    ds["id"].values[3, 0] = np.nan
    ds["id"].values[-1, -1] = np.nan

    if chunk:
        ds["variable"] = ds["variable"].chunk(chunk)
    grouped = ds.groupby(ds.id)
    if shuffle:
        grouped = grouped.shuffle_to_chunks().groupby(ds.id)

    # non reduction operation
    expected1 = ds.copy()
    expected1.variable.data[0, 0, :] = np.nan
    expected1.variable.data[-1, -1, :] = np.nan
    expected1.variable.data[3, 0, :] = np.nan
    actual1 = grouped.map(lambda x: x).transpose(*ds.variable.dims)
    assert_identical(actual1, expected1)

    # reduction along grouped dimension
    actual2 = grouped.mean()
    stacked = ds.stack({"xy": ["lat", "lon"]})
    expected2 = (
        stacked.variable.where(stacked.id.notnull())
        .rename({"xy": "id"})
        .to_dataset()
        .reset_index("id", drop=True)
        .assign(id=stacked.id.values)
        .dropna("id")
        .transpose(*actual2.variable.dims)
    )
    assert_identical(actual2, expected2)

    # reduction operation along a different dimension
    actual3 = grouped.mean("time")
    expected3 = ds.mean("time").where(ds.id.notnull())
    assert_identical(actual3, expected3)

    # NaN in non-dimensional coordinate
    array = xr.DataArray([1, 2, 3], [("x", [1, 2, 3])])
    array["x1"] = ("x", [1, 1, np.nan])
    expected4 = xr.DataArray(3, [("x1", [1])])
    actual4 = array.groupby("x1").sum()
    assert_equal(expected4, actual4)

    # NaT in non-dimensional coordinate
    array["t"] = (
        "x",
        [
            np.datetime64("2001-01-01"),
            np.datetime64("2001-01-01"),
            np.datetime64("NaT"),
        ],
    )
    expected5 = xr.DataArray(3, [("t", [np.datetime64("2001-01-01")])])
    actual5 = array.groupby("t").sum()
    assert_equal(expected5, actual5)

    # test for repeated coordinate labels
    array = xr.DataArray([0, 1, 2, 4, 3, 4], [("x", [np.nan, 1, 1, np.nan, 2, np.nan])])
    expected6 = xr.DataArray([3, 3], [("x", [1, 2])])
    actual6 = array.groupby("x").sum()
    assert_equal(expected6, actual6)


def test_groupby_grouping_errors() -> None:
    dataset = xr.Dataset({"foo": ("x", [1, 1, 1])}, {"x": [1, 2, 3]})
    with pytest.raises(
        ValueError, match=r"None of the data falls within bins with edges"
    ):
        dataset.groupby_bins("x", bins=[0.1, 0.2, 0.3])

    with pytest.raises(
        ValueError, match=r"None of the data falls within bins with edges"
    ):
        dataset.to_dataarray().groupby_bins("x", bins=[0.1, 0.2, 0.3])

    with pytest.raises(ValueError, match=r"All bin edges are NaN."):
        dataset.groupby_bins("x", bins=[np.nan, np.nan, np.nan])

    with pytest.raises(ValueError, match=r"All bin edges are NaN."):
        dataset.to_dataarray().groupby_bins("x", bins=[np.nan, np.nan, np.nan])

    with pytest.raises(ValueError, match=r"Failed to group data."):
        dataset.groupby(dataset.foo * np.nan)

    with pytest.raises(ValueError, match=r"Failed to group data."):
        dataset.to_dataarray().groupby(dataset.foo * np.nan)

    with pytest.raises(TypeError, match=r"Cannot group by a Grouper object"):
        dataset.groupby(UniqueGrouper(labels=[1, 2, 3]))  # type: ignore[arg-type]

    with pytest.raises(TypeError, match=r"got multiple values for argument"):
        UniqueGrouper(dataset.x, labels=[1, 2, 3])  # type: ignore[misc]


def test_groupby_reduce_dimension_error(array) -> None:
    grouped = array.groupby("y")
    # assert_identical(array, grouped.mean())

    with pytest.raises(ValueError, match=r"cannot reduce over dimensions"):
        grouped.mean("huh")

    with pytest.raises(ValueError, match=r"cannot reduce over dimensions"):
        grouped.mean(("x", "y", "asd"))

    assert_identical(array.mean("x"), grouped.reduce(np.mean, "x"))
    assert_allclose(array.mean(["x", "z"]), grouped.reduce(np.mean, ["x", "z"]))

    grouped = array.groupby("y")
    assert_identical(array, grouped.mean())

    assert_identical(array.mean("x"), grouped.reduce(np.mean, "x"))
    assert_allclose(array.mean(["x", "z"]), grouped.reduce(np.mean, ["x", "z"]))


def test_groupby_multiple_string_args(array) -> None:
    with pytest.raises(TypeError):
        array.groupby("x", squeeze="y")


def test_groupby_bins_timeseries() -> None:
    ds = xr.Dataset()
    ds["time"] = xr.DataArray(
        pd.date_range("2010-08-01", "2010-08-15", freq="15min"), dims="time"
    )
    ds["val"] = xr.DataArray(np.ones(ds["time"].shape), dims="time")
    time_bins = pd.date_range(start="2010-08-01", end="2010-08-15", freq="24h")
    actual = ds.groupby_bins("time", time_bins).sum()
    expected = xr.DataArray(
        96 * np.ones((14,)),
        dims=["time_bins"],
        coords={"time_bins": pd.cut(time_bins, time_bins).categories},  # type: ignore[arg-type]
    ).to_dataset(name="val")
    assert_identical(actual, expected)


def test_groupby_none_group_name() -> None:
    # GH158
    # xarray should not fail if a DataArray's name attribute is None

    data = np.arange(10) + 10
    da = xr.DataArray(data)  # da.name = None
    key = xr.DataArray(np.floor_divide(data, 2))

    mean = da.groupby(key).mean()
    assert "group" in mean.dims


def test_groupby_getitem(dataset) -> None:
    assert_identical(dataset.sel(x=["a"]), dataset.groupby("x")["a"])
    assert_identical(dataset.sel(z=[1]), dataset.groupby("z")[1])
    assert_identical(dataset.foo.sel(x=["a"]), dataset.foo.groupby("x")["a"])
    assert_identical(dataset.foo.sel(z=[1]), dataset.foo.groupby("z")[1])
    assert_identical(dataset.cat.sel(y=[1]), dataset.cat.groupby("y")[1])

    with pytest.raises(
        NotImplementedError, match="Cannot broadcast 1d-only pandas extension array."
    ):
        dataset.groupby("boo")
    dataset = dataset.drop_vars(["cat"])
    actual = dataset.groupby("boo")["f"].unstack().transpose("x", "y", "z")
    expected = dataset.sel(y=[1], z=[1, 2]).transpose("x", "y", "z")
    assert_identical(expected, actual)


def test_groupby_dataset() -> None:
    data = Dataset(
        {"z": (["x", "y"], np.random.randn(3, 5))},
        {"x": ("x", list("abc")), "c": ("x", [0, 1, 0]), "y": range(5)},
    )
    groupby = data.groupby("x")
    assert len(groupby) == 3
    expected_groups = {"a": slice(0, 1), "b": slice(1, 2), "c": slice(2, 3)}
    assert groupby.groups == expected_groups
    expected_items = [
        ("a", data.isel(x=[0])),
        ("b", data.isel(x=[1])),
        ("c", data.isel(x=[2])),
    ]
    for actual1, expected1 in zip(groupby, expected_items, strict=True):
        assert actual1[0] == expected1[0]
        assert_equal(actual1[1], expected1[1])

    def identity(x):
        return x

    for k in ["x", "c", "y"]:
        actual2 = data.groupby(k).map(identity)
        assert_equal(data, actual2)


def test_groupby_dataset_returns_new_type() -> None:
    data = Dataset({"z": (["x", "y"], np.random.randn(3, 5))})

    actual1 = data.groupby("x").map(lambda ds: ds["z"])
    expected1 = data["z"]
    assert_identical(expected1, actual1)

    actual2 = data["z"].groupby("x").map(lambda x: x.to_dataset())
    expected2 = data
    assert_identical(expected2, actual2)


def test_groupby_dataset_iter() -> None:
    data = create_test_data()
    for n, (t, sub) in enumerate(list(data.groupby("dim1"))[:3]):
        assert data["dim1"][n] == t
        assert_equal(data["var1"][[n]], sub["var1"])
        assert_equal(data["var2"][[n]], sub["var2"])
        assert_equal(data["var3"][:, [n]], sub["var3"])


def test_groupby_dataset_errors() -> None:
    data = create_test_data()
    with pytest.raises(TypeError, match=r"`group` must be"):
        data.groupby(np.arange(10))  # type: ignore[arg-type,unused-ignore]
    with pytest.raises(ValueError, match=r"length does not match"):
        data.groupby(data["dim1"][:3])
    with pytest.raises(TypeError, match=r"`group` must be"):
        data.groupby(data.coords["dim1"].to_index())  # type: ignore[arg-type]


@pytest.mark.parametrize("use_flox", [True, False])
@pytest.mark.parametrize(
    "by_func",
    [
        pytest.param(lambda x: x, id="group-by-string"),
        pytest.param(lambda x: {x: UniqueGrouper()}, id="group-by-unique-grouper"),
    ],
)
@pytest.mark.parametrize("letters_as_coord", [True, False])
def test_groupby_dataset_reduce_ellipsis(
    by_func, use_flox: bool, letters_as_coord: bool
) -> None:
    data = Dataset(
        {
            "xy": (["x", "y"], np.random.randn(3, 4)),
            "xonly": ("x", np.random.randn(3)),
            "yonly": ("y", np.random.randn(4)),
            "letters": ("y", ["a", "a", "b", "b"]),
        }
    )

    if letters_as_coord:
        data = data.set_coords("letters")

    expected = data.mean("y")
    expected["yonly"] = expected["yonly"].variable.set_dims({"x": 3})
    gb = data.groupby(by_func("x"))
    with xr.set_options(use_flox=use_flox):
        actual = gb.mean(...)
    assert_allclose(expected, actual)

    with xr.set_options(use_flox=use_flox):
        actual = gb.mean("y")
    assert_allclose(expected, actual)

    letters = data["letters"]
    expected = Dataset(
        {
            "xy": data["xy"].groupby(letters).mean(...),
            "xonly": (data["xonly"].mean().variable.set_dims({"letters": 2})),
            "yonly": data["yonly"].groupby(letters).mean(),
        }
    )
    gb = data.groupby(by_func("letters"))
    with xr.set_options(use_flox=use_flox):
        actual = gb.mean(...)
    assert_allclose(expected, actual)


def test_groupby_dataset_math() -> None:
    def reorder_dims(x):
        return x.transpose("dim1", "dim2", "dim3", "time")

    ds = create_test_data()
    ds["dim1"] = ds["dim1"]
    grouped = ds.groupby("dim1")

    expected = reorder_dims(ds + ds.coords["dim1"])
    actual = grouped + ds.coords["dim1"]
    assert_identical(expected, reorder_dims(actual))

    actual = ds.coords["dim1"] + grouped
    assert_identical(expected, reorder_dims(actual))

    ds2 = 2 * ds
    expected = reorder_dims(ds + ds2)
    actual = grouped + ds2
    assert_identical(expected, reorder_dims(actual))

    actual = ds2 + grouped
    assert_identical(expected, reorder_dims(actual))


def test_groupby_math_more() -> None:
    ds = create_test_data()
    grouped = ds.groupby("numbers")
    zeros = DataArray([0, 0, 0, 0], [("numbers", range(4))])
    expected = (ds + Variable("dim3", np.zeros(10))).transpose(
        "dim3", "dim1", "dim2", "time"
    )
    actual = grouped + zeros
    assert_equal(expected, actual)

    actual = zeros + grouped
    assert_equal(expected, actual)

    with pytest.raises(ValueError, match=r"incompat.* grouped binary"):
        grouped + ds
    with pytest.raises(ValueError, match=r"incompat.* grouped binary"):
        ds + grouped
    with pytest.raises(TypeError, match=r"only support binary ops"):
        grouped + 1  # type: ignore[operator]
    with pytest.raises(TypeError, match=r"only support binary ops"):
        grouped + grouped  # type: ignore[operator]
    with pytest.raises(TypeError, match=r"in-place operations"):
        ds += grouped  # type: ignore[arg-type]

    ds = Dataset(
        {
            "x": ("time", np.arange(100)),
            "time": pd.date_range("2000-01-01", periods=100),
        }
    )
    with pytest.raises(ValueError, match=r"incompat.* grouped binary"):
        ds + ds.groupby("time.month")


def test_groupby_math_bitshift() -> None:
    # create new dataset of int's only
    ds = Dataset(
        {
            "x": ("index", np.ones(4, dtype=int)),
            "y": ("index", np.ones(4, dtype=int) * -1),
            "level": ("index", [1, 1, 2, 2]),
            "index": [0, 1, 2, 3],
        }
    )
    shift = DataArray([1, 2, 1], [("level", [1, 2, 8])])

    left_expected = Dataset(
        {
            "x": ("index", [2, 2, 4, 4]),
            "y": ("index", [-2, -2, -4, -4]),
            "level": ("index", [2, 2, 8, 8]),
            "index": [0, 1, 2, 3],
        }
    )

    left_manual = []
    for lev, group in ds.groupby("level"):
        shifter = shift.sel(level=lev)
        left_manual.append(group << shifter)
    left_actual = xr.concat(left_manual, dim="index").reset_coords(names="level")
    assert_equal(left_expected, left_actual)

    left_actual = (ds.groupby("level") << shift).reset_coords(names="level")
    assert_equal(left_expected, left_actual)

    right_expected = Dataset(
        {
            "x": ("index", [0, 0, 2, 2]),
            "y": ("index", [-1, -1, -2, -2]),
            "level": ("index", [0, 0, 4, 4]),
            "index": [0, 1, 2, 3],
        }
    )
    right_manual = []
    for lev, group in left_expected.groupby("level"):
        shifter = shift.sel(level=lev)
        right_manual.append(group >> shifter)
    right_actual = xr.concat(right_manual, dim="index").reset_coords(names="level")
    assert_equal(right_expected, right_actual)

    right_actual = (left_expected.groupby("level") >> shift).reset_coords(names="level")
    assert_equal(right_expected, right_actual)


@pytest.mark.parametrize(
    "x_bins", ((0, 2, 4, 6), pd.IntervalIndex.from_breaks((0, 2, 4, 6), closed="left"))
)
@pytest.mark.parametrize("use_flox", [True, False])
def test_groupby_bins_cut_kwargs(use_flox: bool, x_bins) -> None:
    da = xr.DataArray(np.arange(12).reshape(6, 2), dims=("x", "y"))

    with xr.set_options(use_flox=use_flox):
        actual = da.groupby_bins(
            "x", bins=x_bins, include_lowest=True, right=False
        ).mean()
    expected = xr.DataArray(
        np.array([[1.0, 2.0], [5.0, 6.0], [9.0, 10.0]]),
        dims=("x_bins", "y"),
        coords={
            "x_bins": (
                "x_bins",
                x_bins
                if isinstance(x_bins, pd.IntervalIndex)
                else pd.IntervalIndex.from_breaks(x_bins, closed="left"),
            )
        },
    )
    assert_identical(expected, actual)

    with xr.set_options(use_flox=use_flox):
        actual = da.groupby(
            x=BinGrouper(bins=x_bins, include_lowest=True, right=False),
        ).mean()
    assert_identical(expected, actual)

    with xr.set_options(use_flox=use_flox):
        labels = ["one", "two", "three"]
        actual = da.groupby(x=BinGrouper(bins=x_bins, labels=labels)).sum()
        assert actual.xindexes["x_bins"].index.equals(pd.Index(labels))  # type: ignore[attr-defined]


@pytest.mark.parametrize("indexed_coord", [True, False])
@pytest.mark.parametrize(
    ["groupby_method", "args"],
    (
        ("groupby_bins", ("x", np.arange(0, 8, 3))),
        ("groupby", ({"x": BinGrouper(bins=np.arange(0, 8, 3))},)),
    ),
)
def test_groupby_bins_math(groupby_method, args, indexed_coord) -> None:
    N = 7
    da = DataArray(np.random.random((N, N)), dims=("x", "y"))
    if indexed_coord:
        da["x"] = np.arange(N)
        da["y"] = np.arange(N)

    g = getattr(da, groupby_method)(*args)
    mean = g.mean()
    expected = da.isel(x=slice(1, None)) - mean.isel(x_bins=("x", [0, 0, 0, 1, 1, 1]))
    actual = g - mean
    assert_identical(expected, actual)


def test_groupby_math_nD_group() -> None:
    N = 40
    da = DataArray(
        np.random.random((N, N)),
        dims=("x", "y"),
        coords={
            "labels": (
                "x",
                np.repeat(["a", "b", "c", "d", "e", "f", "g", "h"], repeats=N // 8),
            ),
        },
    )
    da["labels2d"] = xr.broadcast(da.labels, da)[0]

    g = da.groupby("labels2d")
    mean = g.mean()
    expected = da - mean.sel(labels2d=da.labels2d)
    expected["labels"] = expected.labels.broadcast_like(expected.labels2d)
    actual = g - mean
    assert_identical(expected, actual)

    da["num"] = (
        "x",
        np.repeat([1, 2, 3, 4, 5, 6, 7, 8], repeats=N // 8),
    )
    da["num2d"] = xr.broadcast(da.num, da)[0]
    g = da.groupby_bins("num2d", bins=[0, 4, 6])
    mean = g.mean()
    idxr = np.digitize(da.num2d, bins=(0, 4, 6), right=True)[:30, :] - 1
    expanded_mean = mean.drop_vars("num2d_bins").isel(num2d_bins=(("x", "y"), idxr))
    expected = da.isel(x=slice(30)) - expanded_mean
    expected["labels"] = expected.labels.broadcast_like(expected.labels2d)
    expected["num"] = expected.num.broadcast_like(expected.num2d)
    # mean.num2d_bins.data is a pandas IntervalArray so needs to be put in `numpy` to allow indexing
    expected["num2d_bins"] = (("x", "y"), mean.num2d_bins.data.to_numpy()[idxr])
    actual = g - mean
    assert_identical(expected, actual)


def test_groupby_dataset_math_virtual() -> None:
    ds = Dataset({"x": ("t", [1, 2, 3])}, {"t": pd.date_range("20100101", periods=3)})
    grouped = ds.groupby("t.day")
    actual = grouped - grouped.mean(...)
    expected = Dataset({"x": ("t", [0, 0, 0])}, ds[["t", "t.day"]])
    assert_identical(actual, expected)


def test_groupby_math_dim_order() -> None:
    da = DataArray(
        np.ones((10, 10, 12)),
        dims=("x", "y", "time"),
        coords={"time": pd.date_range("2001-01-01", periods=12, freq="6h")},
    )
    grouped = da.groupby("time.day")
    result = grouped - grouped.mean()
    assert result.dims == da.dims


def test_groupby_dataset_nan() -> None:
    # nan should be excluded from groupby
    ds = Dataset({"foo": ("x", [1, 2, 3, 4])}, {"bar": ("x", [1, 1, 2, np.nan])})
    actual = ds.groupby("bar").mean(...)
    expected = Dataset({"foo": ("bar", [1.5, 3]), "bar": [1, 2]})
    assert_identical(actual, expected)


def test_groupby_dataset_order() -> None:
    # groupby should preserve variables order
    ds = Dataset()
    for vn in ["a", "b", "c"]:
        ds[vn] = DataArray(np.arange(10), dims=["t"])
    data_vars_ref = list(ds.data_vars.keys())
    ds = ds.groupby("t").mean(...)
    data_vars = list(ds.data_vars.keys())
    assert data_vars == data_vars_ref
    # coords are now at the end of the list, so the test below fails
    # all_vars = list(ds.variables.keys())
    # all_vars_ref = list(ds.variables.keys())
    # .assertEqual(all_vars, all_vars_ref)


def test_groupby_dataset_fillna() -> None:
    ds = Dataset({"a": ("x", [np.nan, 1, np.nan, 3])}, {"x": [0, 1, 2, 3]})
    expected = Dataset({"a": ("x", range(4))}, {"x": [0, 1, 2, 3]})
    for target in [ds, expected]:
        target.coords["b"] = ("x", [0, 0, 1, 1])
    actual = ds.groupby("b").fillna(DataArray([0, 2], dims="b"))
    assert_identical(expected, actual)

    actual = ds.groupby("b").fillna(Dataset({"a": ("b", [0, 2])}))
    assert_identical(expected, actual)

    # attrs with groupby
    ds.attrs["attr"] = "ds"
    ds.a.attrs["attr"] = "da"
    actual = ds.groupby("b").fillna(Dataset({"a": ("b", [0, 2])}))
    assert actual.attrs == ds.attrs
    assert actual.a.name == "a"
    assert actual.a.attrs == ds.a.attrs


def test_groupby_dataset_where() -> None:
    # groupby
    ds = Dataset({"a": ("x", range(5))}, {"c": ("x", [0, 0, 1, 1, 1])})
    cond = Dataset({"a": ("c", [True, False])})
    expected = ds.copy(deep=True)
    expected["a"].values = np.array([0, 1] + [np.nan] * 3)
    actual = ds.groupby("c").where(cond)
    assert_identical(expected, actual)

    # attrs with groupby
    ds.attrs["attr"] = "ds"
    ds.a.attrs["attr"] = "da"
    actual = ds.groupby("c").where(cond)
    assert actual.attrs == ds.attrs
    assert actual.a.name == "a"
    assert actual.a.attrs == ds.a.attrs


def test_groupby_dataset_assign() -> None:
    ds = Dataset({"a": ("x", range(3))}, {"b": ("x", ["A"] * 2 + ["B"])})
    actual = ds.groupby("b").assign(c=lambda ds: 2 * ds.a)
    expected = ds.merge({"c": ("x", [0, 2, 4])})
    assert_identical(actual, expected)

    actual = ds.groupby("b").assign(c=lambda ds: ds.a.sum())
    expected = ds.merge({"c": ("x", [1, 1, 2])})
    assert_identical(actual, expected)

    actual = ds.groupby("b").assign_coords(c=lambda ds: ds.a.sum())
    expected = expected.set_coords("c")
    assert_identical(actual, expected)


def test_groupby_dataset_map_dataarray_func() -> None:
    # regression GH6379
    ds = Dataset({"foo": ("x", [1, 2, 3, 4])}, coords={"x": [0, 0, 1, 1]})
    actual = ds.groupby("x").map(lambda grp: grp.foo.mean())
    expected = DataArray([1.5, 3.5], coords={"x": [0, 1]}, dims="x", name="foo")
    assert_identical(actual, expected)


def test_groupby_dataarray_map_dataset_func() -> None:
    # regression GH6379
    da = DataArray([1, 2, 3, 4], coords={"x": [0, 0, 1, 1]}, dims="x", name="foo")
    actual = da.groupby("x").map(lambda grp: grp.mean().to_dataset())
    expected = xr.Dataset({"foo": ("x", [1.5, 3.5])}, coords={"x": [0, 1]})
    assert_identical(actual, expected)


@requires_flox
@pytest.mark.parametrize("kwargs", [{"method": "map-reduce"}, {"engine": "numpy"}])
def test_groupby_flox_kwargs(kwargs) -> None:
    ds = Dataset({"a": ("x", range(5))}, {"c": ("x", [0, 0, 1, 1, 1])})
    with xr.set_options(use_flox=False):
        expected = ds.groupby("c").mean()
    with xr.set_options(use_flox=True):
        actual = ds.groupby("c").mean(**kwargs)
    assert_identical(expected, actual)


class TestDataArrayGroupBy:
    @pytest.fixture(autouse=True)
    def setup(self) -> None:
        self.attrs = {"attr1": "value1", "attr2": 2929}
        self.x = np.random.random((10, 20))
        self.v = Variable(["x", "y"], self.x)
        self.va = Variable(["x", "y"], self.x, self.attrs)
        self.ds = Dataset({"foo": self.v})
        self.dv = self.ds["foo"]

        self.mindex = pd.MultiIndex.from_product(
            [["a", "b"], [1, 2]], names=("level_1", "level_2")
        )
        self.mda = DataArray([0, 1, 2, 3], coords={"x": self.mindex}, dims="x")

        self.da = self.dv.copy()
        self.da.coords["abc"] = ("y", np.array(["a"] * 9 + ["c"] + ["b"] * 10))
        self.da.coords["y"] = 20 + 100 * self.da["y"]

    def test_stack_groupby_unsorted_coord(self) -> None:
        data = [[0, 1], [2, 3]]
        data_flat = [0, 1, 2, 3]
        dims = ["x", "y"]
        y_vals = [2, 3]

        arr = xr.DataArray(data, dims=dims, coords={"y": y_vals})
        actual1 = arr.stack(z=dims).groupby("z").first()
        midx1 = pd.MultiIndex.from_product([[0, 1], [2, 3]], names=dims)
        expected1 = xr.DataArray(data_flat, dims=["z"], coords={"z": midx1})
        assert_equal(actual1, expected1)

        # GH: 3287.  Note that y coord values are not in sorted order.
        arr = xr.DataArray(data, dims=dims, coords={"y": y_vals[::-1]})
        actual2 = arr.stack(z=dims).groupby("z").first()
        midx2 = pd.MultiIndex.from_product([[0, 1], [3, 2]], names=dims)
        expected2 = xr.DataArray(data_flat, dims=["z"], coords={"z": midx2})
        assert_equal(actual2, expected2)

    def test_groupby_iter(self) -> None:
        for (act_x, act_dv), (exp_x, exp_ds) in zip(
            self.dv.groupby("y"), self.ds.groupby("y"), strict=True
        ):
            assert exp_x == act_x
            assert_identical(exp_ds["foo"], act_dv)
            for (_, exp_dv), (_, act_dv) in zip(
                self.dv.groupby("x"), self.dv.groupby("x"), strict=True
            ):
                assert_identical(exp_dv, act_dv)

    def test_groupby_properties(self) -> None:
        grouped = self.da.groupby("abc")
        expected_groups = {"a": range(9), "c": [9], "b": range(10, 20)}
        assert expected_groups.keys() == grouped.groups.keys()
        for key, expected_group in expected_groups.items():
            actual_group = grouped.groups[key]

            # TODO: array_api doesn't allow slice:
            assert not isinstance(expected_group, slice)
            assert not isinstance(actual_group, slice)

            np.testing.assert_array_equal(expected_group, actual_group)
        assert 3 == len(grouped)

    @pytest.mark.parametrize(
        "by, use_da", [("x", False), ("y", False), ("y", True), ("abc", False)]
    )
    @pytest.mark.parametrize("shortcut", [True, False])
    def test_groupby_map_identity(self, by, use_da, shortcut) -> None:
        expected = self.da
        if use_da:
            by = expected.coords[by]

        def identity(x):
            return x

        grouped = expected.groupby(by)
        actual = grouped.map(identity, shortcut=shortcut)
        assert_identical(expected, actual)

    def test_groupby_sum(self) -> None:
        array = self.da
        grouped = array.groupby("abc")

        expected_sum_all = Dataset(
            {
                "foo": Variable(
                    ["abc"],
                    np.array(
                        [
                            self.x[:, :9].sum(),
                            self.x[:, 10:].sum(),
                            self.x[:, 9:10].sum(),
                        ]
                    ).T,
                ),
                "abc": Variable(["abc"], np.array(["a", "b", "c"])),
            }
        )["foo"]
        assert_allclose(expected_sum_all, grouped.reduce(np.sum, dim=...))
        assert_allclose(expected_sum_all, grouped.sum(...))

        expected = DataArray(
            [
                array["y"].values[idx].sum()
                for idx in [slice(9), slice(10, None), slice(9, 10)]
            ],
            [["a", "b", "c"]],
            ["abc"],
        )
        actual = array["y"].groupby("abc").map(np.sum)
        assert_allclose(expected, actual)
        actual = array["y"].groupby("abc").sum(...)
        assert_allclose(expected, actual)

        expected_sum_axis1 = Dataset(
            {
                "foo": (
                    ["x", "abc"],
                    np.array(
                        [
                            self.x[:, :9].sum(1),
                            self.x[:, 10:].sum(1),
                            self.x[:, 9:10].sum(1),
                        ]
                    ).T,
                ),
                "abc": Variable(["abc"], np.array(["a", "b", "c"])),
            }
        )["foo"]
        assert_allclose(expected_sum_axis1, grouped.reduce(np.sum, "y"))
        assert_allclose(expected_sum_axis1, grouped.sum("y"))

    @pytest.mark.parametrize("use_flox", [True, False])
    @pytest.mark.parametrize("shuffle", [True, False])
    @pytest.mark.parametrize(
        "chunk",
        [
            pytest.param(
                True, marks=pytest.mark.skipif(not has_dask, reason="no dask")
            ),
            False,
        ],
    )
    @pytest.mark.parametrize("method", ["sum", "mean", "median"])
    def test_groupby_reductions(
        self, use_flox: bool, method: str, shuffle: bool, chunk: bool
    ) -> None:
        if shuffle and chunk and not has_dask_ge_2024_08_1:
            pytest.skip()

        array = self.da
        if chunk:
            array.data = array.chunk({"y": 5}).data
        reduction = getattr(np, method)
        expected = Dataset(
            {
                "foo": Variable(
                    ["x", "abc"],
                    np.array(
                        [
                            reduction(self.x[:, :9], axis=-1),
                            reduction(self.x[:, 10:], axis=-1),
                            reduction(self.x[:, 9:10], axis=-1),
                        ]
                    ).T,
                ),
                "abc": Variable(["abc"], np.array(["a", "b", "c"])),
            }
        )["foo"]

        with raise_if_dask_computes():
            grouped = array.groupby("abc")
            if shuffle:
                grouped = grouped.shuffle_to_chunks().groupby("abc")

            with xr.set_options(use_flox=use_flox):
                actual = getattr(grouped, method)(dim="y")
        assert_allclose(expected, actual)

    def test_groupby_count(self) -> None:
        array = DataArray(
            [0, 0, np.nan, np.nan, 0, 0],
            coords={"cat": ("x", ["a", "b", "b", "c", "c", "c"])},
            dims="x",
        )
        actual = array.groupby("cat").count()
        expected = DataArray([1, 1, 2], coords=[("cat", ["a", "b", "c"])])
        assert_identical(actual, expected)

    @pytest.mark.parametrize("shortcut", [True, False])
    @pytest.mark.parametrize("keep_attrs", [None, True, False])
    def test_groupby_reduce_keep_attrs(
        self, shortcut: bool, keep_attrs: bool | None
    ) -> None:
        array = self.da
        array.attrs["foo"] = "bar"

        actual = array.groupby("abc").reduce(
            np.mean, keep_attrs=keep_attrs, shortcut=shortcut
        )
        with xr.set_options(use_flox=False):
            expected = array.groupby("abc").mean(keep_attrs=keep_attrs)
        assert_identical(expected, actual)

    @pytest.mark.parametrize("keep_attrs", [None, True, False])
    def test_groupby_keep_attrs(self, keep_attrs: bool | None) -> None:
        array = self.da
        array.attrs["foo"] = "bar"

        with xr.set_options(use_flox=False):
            expected = array.groupby("abc").mean(keep_attrs=keep_attrs)
        with xr.set_options(use_flox=True):
            actual = array.groupby("abc").mean(keep_attrs=keep_attrs)

        # values are tested elsewhere, here we just check data
        # TODO: add check_attrs kwarg to assert_allclose
        actual.data = expected.data
        assert_identical(expected, actual)

    def test_groupby_map_center(self) -> None:
        def center(x):
            return x - np.mean(x)

        array = self.da
        grouped = array.groupby("abc")

        expected_ds = array.to_dataset()
        exp_data = np.hstack(
            [center(self.x[:, :9]), center(self.x[:, 9:10]), center(self.x[:, 10:])]
        )
        expected_ds["foo"] = (["x", "y"], exp_data)
        expected_centered = expected_ds["foo"]
        assert_allclose(expected_centered, grouped.map(center))

    def test_groupby_map_ndarray(self) -> None:
        # regression test for #326
        array = self.da
        grouped = array.groupby("abc")
        actual = grouped.map(np.asarray)  # type: ignore[arg-type] # TODO: Not sure using np.asarray like this makes sense with array api
        assert_equal(array, actual)

    def test_groupby_map_changes_metadata(self) -> None:
        def change_metadata(x):
            x.coords["x"] = x.coords["x"] * 2
            x.attrs["fruit"] = "lemon"
            return x

        array = self.da
        grouped = array.groupby("abc")
        actual = grouped.map(change_metadata)
        expected = array.copy()
        expected = change_metadata(expected)
        assert_equal(expected, actual)

    def test_groupby_math_squeeze(self) -> None:
        array = self.da
        grouped = array.groupby("x")

        expected = array + array.coords["x"]
        actual = grouped + array.coords["x"]
        assert_identical(expected, actual)

        actual = array.coords["x"] + grouped
        assert_identical(expected, actual)

        ds = array.coords["x"].to_dataset(name="X")
        expected = array + ds
        actual = grouped + ds
        assert_identical(expected, actual)

        actual = ds + grouped
        assert_identical(expected, actual)

    def test_groupby_math(self) -> None:
        array = self.da
        grouped = array.groupby("abc")
        expected_agg = (grouped.mean(...) - np.arange(3)).rename(None)
        actual = grouped - DataArray(range(3), [("abc", ["a", "b", "c"])])
        actual_agg = actual.groupby("abc").mean(...)
        assert_allclose(expected_agg, actual_agg)

        with pytest.raises(TypeError, match=r"only support binary ops"):
            grouped + 1  # type: ignore[type-var]
        with pytest.raises(TypeError, match=r"only support binary ops"):
            grouped + grouped  # type: ignore[type-var]
        with pytest.raises(TypeError, match=r"in-place operations"):
            array += grouped  # type: ignore[arg-type]

    def test_groupby_math_not_aligned(self) -> None:
        array = DataArray(
            range(4), {"b": ("x", [0, 0, 1, 1]), "x": [0, 1, 2, 3]}, dims="x"
        )
        other = DataArray([10], coords={"b": [0]}, dims="b")
        actual = array.groupby("b") + other
        expected = DataArray([10, 11, np.nan, np.nan], array.coords)
        assert_identical(expected, actual)

        # regression test for #7797
        other = array.groupby("b").sum()
        actual = array.sel(x=[0, 1]).groupby("b") - other
        expected = DataArray([-1, 0], {"b": ("x", [0, 0]), "x": [0, 1]}, dims="x")
        assert_identical(expected, actual)

        other = DataArray([10], coords={"c": 123, "b": [0]}, dims="b")
        actual = array.groupby("b") + other
        expected = DataArray([10, 11, np.nan, np.nan], array.coords)
        expected.coords["c"] = (["x"], [123] * 2 + [np.nan] * 2)
        assert_identical(expected, actual)

        other_ds = Dataset({"a": ("b", [10])}, {"b": [0]})
        actual_ds = array.groupby("b") + other_ds
        expected_ds = Dataset({"a": ("x", [10, 11, np.nan, np.nan])}, array.coords)
        assert_identical(expected_ds, actual_ds)

    def test_groupby_restore_dim_order(self) -> None:
        array = DataArray(
            np.random.randn(5, 3),
            coords={"a": ("x", range(5)), "b": ("y", range(3))},
            dims=["x", "y"],
        )
        for by, expected_dims in [
            ("x", ("x", "y")),
            ("y", ("x", "y")),
            ("a", ("a", "y")),
            ("b", ("x", "b")),
        ]:
            result = array.groupby(by).map(lambda x: x.squeeze())
            assert result.dims == expected_dims

    def test_groupby_restore_coord_dims(self) -> None:
        array = DataArray(
            np.random.randn(5, 3),
            coords={
                "a": ("x", range(5)),
                "b": ("y", range(3)),
                "c": (("x", "y"), np.random.randn(5, 3)),
            },
            dims=["x", "y"],
        )

        for by, expected_dims in [
            ("x", ("x", "y")),
            ("y", ("x", "y")),
            ("a", ("a", "y")),
            ("b", ("x", "b")),
        ]:
            result = array.groupby(by, restore_coord_dims=True).map(
                lambda x: x.squeeze()
            )["c"]
            assert result.dims == expected_dims

    def test_groupby_first_and_last(self) -> None:
        array = DataArray([1, 2, 3, 4, 5], dims="x")
        by = DataArray(["a"] * 2 + ["b"] * 3, dims="x", name="ab")

        expected = DataArray([1, 3], [("ab", ["a", "b"])])
        actual = array.groupby(by).first()
        assert_identical(expected, actual)

        expected = DataArray([2, 5], [("ab", ["a", "b"])])
        actual = array.groupby(by).last()
        assert_identical(expected, actual)

        array = DataArray(np.random.randn(5, 3), dims=["x", "y"])
        expected = DataArray(array[[0, 2]], {"ab": ["a", "b"]}, ["ab", "y"])
        actual = array.groupby(by).first()
        assert_identical(expected, actual)

        actual = array.groupby("x").first()
        expected = array  # should be a no-op
        assert_identical(expected, actual)

        # TODO: groupby_bins too

    def make_groupby_multidim_example_array(self) -> DataArray:
        return DataArray(
            [[[0, 1], [2, 3]], [[5, 10], [15, 20]]],
            coords={
                "lon": (["ny", "nx"], [[30, 40], [40, 50]]),
                "lat": (["ny", "nx"], [[10, 10], [20, 20]]),
            },
            dims=["time", "ny", "nx"],
        )

    def test_groupby_multidim(self) -> None:
        array = self.make_groupby_multidim_example_array()
        for dim, expected_sum in [
            ("lon", DataArray([5, 28, 23], coords=[("lon", [30.0, 40.0, 50.0])])),
            ("lat", DataArray([16, 40], coords=[("lat", [10.0, 20.0])])),
        ]:
            actual_sum = array.groupby(dim).sum(...)
            assert_identical(expected_sum, actual_sum)

        if has_flox:
            # GH9803
            # reduce over one dim of a nD grouper
            array.coords["labels"] = (("ny", "nx"), np.array([["a", "b"], ["b", "a"]]))
            actual = array.groupby("labels").sum("nx")
            expected_np = np.array([[[0, 1], [3, 2]], [[5, 10], [20, 15]]])
            expected = xr.DataArray(
                expected_np,
                dims=("time", "ny", "labels"),
                coords={"labels": ["a", "b"]},
            )
            assert_identical(expected, actual)

    def test_groupby_multidim_map(self) -> None:
        array = self.make_groupby_multidim_example_array()
        actual = array.groupby("lon").map(lambda x: x - x.mean())
        expected = DataArray(
            [[[-2.5, -6.0], [-5.0, -8.5]], [[2.5, 3.0], [8.0, 8.5]]],
            coords=array.coords,
            dims=array.dims,
        )
        assert_identical(expected, actual)

    @pytest.mark.parametrize("use_flox", [True, False])
    @pytest.mark.parametrize("coords", [np.arange(4), np.arange(4)[::-1], [2, 0, 3, 1]])
    @pytest.mark.parametrize(
        "cut_kwargs",
        (
            {"labels": None, "include_lowest": True},
            {"labels": None, "include_lowest": False},
            {"labels": ["a", "b"]},
            {"labels": [1.2, 3.5]},
            {"labels": ["b", "a"]},
        ),
    )
    def test_groupby_bins(
        self,
        coords: np.typing.ArrayLike,
        use_flox: bool,
        cut_kwargs: dict,
    ) -> None:
        array = DataArray(
            np.arange(4), dims="dim_0", coords={"dim_0": coords}, name="a"
        )
        # the first value should not be part of any group ("right" binning)
        array[0] = 99
        # bins follow conventions for pandas.cut
        # https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.cut.html
        bins = [0, 1.5, 5]

        df = array.to_dataframe()
        df["dim_0_bins"] = pd.cut(array["dim_0"], bins, **cut_kwargs)  # type: ignore[call-overload]

        expected_df = df.groupby("dim_0_bins", observed=True).sum()
        expected = expected_df.to_xarray().assign_coords(
            dim_0_bins=cast(pd.CategoricalIndex, expected_df.index).categories
        )["a"]

        with xr.set_options(use_flox=use_flox):
            gb = array.groupby_bins("dim_0", bins=bins, **cut_kwargs)
            shuffled = gb.shuffle_to_chunks().groupby_bins(
                "dim_0", bins=bins, **cut_kwargs
            )
            actual = gb.sum()
            assert_identical(expected, actual)
            assert_identical(expected, shuffled.sum())

            actual = gb.map(lambda x: x.sum())
            assert_identical(expected, actual)
            assert_identical(expected, shuffled.map(lambda x: x.sum()))

            # make sure original array dims are unchanged
            assert len(array.dim_0) == 4

    def test_groupby_bins_ellipsis(self) -> None:
        da = xr.DataArray(np.ones((2, 3, 4)))
        bins = [-1, 0, 1, 2]
        with xr.set_options(use_flox=False):
            actual = da.groupby_bins("dim_0", bins).mean(...)
        with xr.set_options(use_flox=True):
            expected = da.groupby_bins("dim_0", bins).mean(...)
        assert_allclose(actual, expected)

    @pytest.mark.parametrize("use_flox", [True, False])
    def test_groupby_bins_gives_correct_subset(self, use_flox: bool) -> None:
        # GH7766
        rng = np.random.default_rng(42)
        coords = rng.normal(5, 5, 1000)
        bins = np.logspace(-4, 1, 10)
        labels = [
            "one",
            "two",
            "three",
            "four",
            "five",
            "six",
            "seven",
            "eight",
            "nine",
        ]
        # xArray
        # Make a mock dataarray
        darr = xr.DataArray(coords, coords=[coords], dims=["coords"])
        expected = xr.DataArray(
            [np.nan, np.nan, 1, 1, 1, 8, 31, 104, 542],
            dims="coords_bins",
            coords={"coords_bins": labels},
        )
        gb = darr.groupby_bins("coords", bins, labels=labels)
        with xr.set_options(use_flox=use_flox):
            actual = gb.count()
        assert_identical(actual, expected)

    def test_groupby_bins_empty(self) -> None:
        array = DataArray(np.arange(4), [("x", range(4))])
        # one of these bins will be empty
        bins = [0, 4, 5]
        bin_coords = pd.cut(array["x"], bins).categories  # type: ignore[call-overload]
        actual = array.groupby_bins("x", bins).sum()
        expected = DataArray([6, np.nan], dims="x_bins", coords={"x_bins": bin_coords})
        assert_identical(expected, actual)
        # make sure original array is unchanged
        # (was a problem in earlier versions)
        assert len(array.x) == 4

    def test_groupby_bins_multidim(self) -> None:
        array = self.make_groupby_multidim_example_array()
        bins = [0, 15, 20]
        bin_coords = pd.cut(array["lat"].values.flat, bins).categories  # type: ignore[call-overload]
        expected = DataArray([16, 40], dims="lat_bins", coords={"lat_bins": bin_coords})
        actual = array.groupby_bins("lat", bins).map(lambda x: x.sum())
        assert_identical(expected, actual)
        # modify the array coordinates to be non-monotonic after unstacking
        array["lat"].data = np.array([[10.0, 20.0], [20.0, 10.0]])
        expected = DataArray([28, 28], dims="lat_bins", coords={"lat_bins": bin_coords})
        actual = array.groupby_bins("lat", bins).map(lambda x: x.sum())
        assert_identical(expected, actual)

        bins = [-2, -1, 0, 1, 2]
        field = DataArray(np.ones((5, 3)), dims=("x", "y"))
        by = DataArray(
            np.array([[-1.5, -1.5, 0.5, 1.5, 1.5] * 3]).reshape(5, 3), dims=("x", "y")
        )
        actual = field.groupby_bins(by, bins=bins).count()

        bincoord = np.array(
            [
                pd.Interval(left, right, closed="right")
                for left, right in pairwise(bins)
            ],
            dtype=object,
        )
        expected = DataArray(
            np.array([6, np.nan, 3, 6]),
            dims="group_bins",
            coords={"group_bins": bincoord},
        )
        assert_identical(actual, expected)

    def test_groupby_bins_sort(self) -> None:
        data = xr.DataArray(
            np.arange(100), dims="x", coords={"x": np.linspace(-100, 100, num=100)}
        )
        binned_mean = data.groupby_bins("x", bins=11).mean()
        assert binned_mean.to_index().is_monotonic_increasing

        with xr.set_options(use_flox=True):
            actual = data.groupby_bins("x", bins=11).count()
        with xr.set_options(use_flox=False):
            expected = data.groupby_bins("x", bins=11).count()
        assert_identical(actual, expected)

    def test_groupby_assign_coords(self) -> None:
        array = DataArray([1, 2, 3, 4], {"c": ("x", [0, 0, 1, 1])}, dims="x")
        actual = array.groupby("c").assign_coords(d=lambda a: a.mean())
        expected = array.copy()
        expected.coords["d"] = ("x", [1.5, 1.5, 3.5, 3.5])
        assert_identical(actual, expected)

    def test_groupby_fillna(self) -> None:
        a = DataArray([np.nan, 1, np.nan, 3], coords={"x": range(4)}, dims="x")
        fill_value = DataArray([0, 1], dims="y")
        actual = a.fillna(fill_value)
        expected = DataArray(
            [[0, 1], [1, 1], [0, 1], [3, 3]], coords={"x": range(4)}, dims=("x", "y")
        )
        assert_identical(expected, actual)

        b = DataArray(range(4), coords={"x": range(4)}, dims="x")
        expected = b.copy()
        for target in [a, expected]:
            target.coords["b"] = ("x", [0, 0, 1, 1])
        actual = a.groupby("b").fillna(DataArray([0, 2], dims="b"))
        assert_identical(expected, actual)

    @pytest.mark.parametrize("use_flox", [True, False])
    def test_groupby_fastpath_for_monotonic(self, use_flox: bool) -> None:
        # Fixes https://github.com/pydata/xarray/issues/6220
        # Fixes https://github.com/pydata/xarray/issues/9279
        index = [1, 2, 3, 4, 7, 9, 10]
        array = DataArray(np.arange(len(index)), [("idx", index)])
        array_rev = array.copy().assign_coords({"idx": index[::-1]})
        fwd = array.groupby("idx", squeeze=False)
        rev = array_rev.groupby("idx", squeeze=False)

        for gb in [fwd, rev]:
            assert all(isinstance(elem, slice) for elem in gb.encoded.group_indices)

        with xr.set_options(use_flox=use_flox):
            assert_identical(fwd.sum(), array)
            assert_identical(rev.sum(), array_rev)


class TestDataArrayResample:
    @pytest.mark.parametrize("shuffle", [True, False])
    @pytest.mark.parametrize(
        "resample_freq",
        [
            "24h",
            "123456s",
            "1234567890us",
            pd.Timedelta(hours=2),
            pd.offsets.MonthBegin(),
            pd.offsets.Second(123456),
            datetime.timedelta(days=1, hours=6),
        ],
    )
    def test_resample(
        self, use_cftime: bool, shuffle: bool, resample_freq: ResampleCompatible
    ) -> None:
        if use_cftime and not has_cftime:
            pytest.skip()
        times = xr.date_range(
            "2000-01-01", freq="6h", periods=10, use_cftime=use_cftime
        )

        def resample_as_pandas(array, *args, **kwargs):
            array_ = array.copy(deep=True)
            if use_cftime:
                array_["time"] = times.to_datetimeindex(time_unit="ns")
            result = DataArray.from_series(
                array_.to_series().resample(*args, **kwargs).mean()
            )
            if use_cftime:
                result = result.convert_calendar(
                    calendar="standard", use_cftime=use_cftime
                )
            return result

        array = DataArray(np.arange(10), [("time", times)])

        rs = array.resample(time=resample_freq)
        shuffled = rs.shuffle_to_chunks().resample(time=resample_freq)
        actual = rs.mean()
        expected = resample_as_pandas(array, resample_freq)
        assert_identical(expected, actual)
        assert_identical(expected, shuffled.mean())

        assert_identical(expected, rs.reduce(np.mean))
        assert_identical(expected, shuffled.reduce(np.mean))

        rs = array.resample(time="24h", closed="right")
        actual = rs.mean()
        shuffled = rs.shuffle_to_chunks().resample(time="24h", closed="right")
        expected = resample_as_pandas(array, "24h", closed="right")
        assert_identical(expected, actual)
        assert_identical(expected, shuffled.mean())

        with pytest.raises(ValueError, match=r"Index must be monotonic"):
            array[[2, 0, 1]].resample(time=resample_freq)

        reverse = array.isel(time=slice(-1, None, -1))
        with pytest.raises(ValueError):
            reverse.resample(time=resample_freq).mean()

    def test_resample_doctest(self, use_cftime: bool) -> None:
        # run the doctest example here so we are not surprised
        da = xr.DataArray(
            np.array([1, 2, 3, 1, 2, np.nan]),
            dims="time",
            coords=dict(
                time=(
                    "time",
                    xr.date_range(
                        "2001-01-01", freq="ME", periods=6, use_cftime=use_cftime
                    ),
                ),
                labels=("time", np.array(["a", "b", "c", "c", "b", "a"])),
            ),
        )
        actual = da.resample(time="3ME").count()
        expected = DataArray(
            [1, 3, 1],
            dims="time",
            coords={
                "time": xr.date_range(
                    "2001-01-01", freq="3ME", periods=3, use_cftime=use_cftime
                )
            },
        )
        assert_identical(actual, expected)

    def test_da_resample_func_args(self) -> None:
        def func(arg1, arg2, arg3=0.0):
            return arg1.mean("time") + arg2 + arg3

        times = pd.date_range("2000", periods=3, freq="D")
        da = xr.DataArray([1.0, 1.0, 1.0], coords=[times], dims=["time"])
        expected = xr.DataArray([3.0, 3.0, 3.0], coords=[times], dims=["time"])
        actual = da.resample(time="D").map(func, args=(1.0,), arg3=1.0)
        assert_identical(actual, expected)

    def test_resample_first_last(self, use_cftime) -> None:
        times = xr.date_range(
            "2000-01-01", freq="6h", periods=10, use_cftime=use_cftime
        )
        array = DataArray(np.arange(10), [("time", times)])

        # resample to same frequency
        actual = array.resample(time="6h").first()
        assert_identical(array, actual)

        actual = array.resample(time="1D").first()
        expected = DataArray([0, 4, 8], [("time", times[::4])])
        assert_identical(expected, actual)

        # verify that labels don't use the first value
        actual = array.resample(time="24h").first()
        expected = array.isel(time=[0, 4, 8])
        assert_identical(expected, actual)

        # missing values
        array = array.astype(float)
        array[:2] = np.nan
        actual = array.resample(time="1D").first()
        expected = DataArray([2, 4, 8], [("time", times[::4])])
        assert_identical(expected, actual)

        actual = array.resample(time="1D").first(skipna=False)
        expected = DataArray([np.nan, 4, 8], [("time", times[::4])])
        assert_identical(expected, actual)

        # regression test for https://stackoverflow.com/questions/33158558/
        array = Dataset({"time": times})["time"]
        actual = array.resample(time="1D").last()
        expected = array.isel(time=[3, 7, 9]).assign_coords(time=times[::4])
        assert_identical(expected, actual)

        # missing periods, GH10169
        actual = array.isel(time=[0, 1, 2, 3, 8, 9]).resample(time="1D").last()
        expected = DataArray(
            np.array([times[3], np.datetime64("NaT"), times[9]]),
            dims="time",
            coords={"time": times[::4]},
            name="time",
        )
        assert_identical(expected, actual)

    def test_resample_bad_resample_dim(self) -> None:
        times = pd.date_range("2000-01-01", freq="6h", periods=10)
        array = DataArray(np.arange(10), [("__resample_dim__", times)])
        with pytest.raises(ValueError, match=r"Proxy resampling dimension"):
            array.resample(__resample_dim__="1D").first()

    @requires_scipy
    def test_resample_drop_nondim_coords(self) -> None:
        xs = np.arange(6)
        ys = np.arange(3)
        times = pd.date_range("2000-01-01", freq="6h", periods=5)
        data = np.tile(np.arange(5), (6, 3, 1))
        xx, yy = np.meshgrid(xs * 5, ys * 2.5)
        tt = np.arange(len(times), dtype=int)
        array = DataArray(data, {"time": times, "x": xs, "y": ys}, ("x", "y", "time"))
        xcoord = DataArray(xx.T, {"x": xs, "y": ys}, ("x", "y"))
        ycoord = DataArray(yy.T, {"x": xs, "y": ys}, ("x", "y"))
        tcoord = DataArray(tt, {"time": times}, ("time",))
        ds = Dataset({"data": array, "xc": xcoord, "yc": ycoord, "tc": tcoord})
        ds = ds.set_coords(["xc", "yc", "tc"])

        # Select the data now, with the auxiliary coordinates in place
        array = ds["data"]

        # Re-sample
        actual = array.resample(time="12h", restore_coord_dims=True).mean("time")
        assert "tc" not in actual.coords

        # Up-sample - filling
        actual = array.resample(time="1h", restore_coord_dims=True).ffill()
        assert "tc" not in actual.coords

        # Up-sample - interpolation
        actual = array.resample(time="1h", restore_coord_dims=True).interpolate(
            "linear"
        )
        assert "tc" not in actual.coords

    def test_resample_keep_attrs(self) -> None:
        times = pd.date_range("2000-01-01", freq="6h", periods=10)
        array = DataArray(np.ones(10), [("time", times)])
        array.attrs["meta"] = "data"

        result = array.resample(time="1D").mean(keep_attrs=True)
        expected = DataArray([1, 1, 1], [("time", times[::4])], attrs=array.attrs)
        assert_identical(result, expected)

    def test_resample_skipna(self) -> None:
        times = pd.date_range("2000-01-01", freq="6h", periods=10)
        array = DataArray(np.ones(10), [("time", times)])
        array[1] = np.nan

        result = array.resample(time="1D").mean(skipna=False)
        expected = DataArray([np.nan, 1, 1], [("time", times[::4])])
        assert_identical(result, expected)

    def test_upsample(self) -> None:
        times = pd.date_range("2000-01-01", freq="6h", periods=5)
        array = DataArray(np.arange(5), [("time", times)])

        # Forward-fill
        actual = array.resample(time="3h").ffill()
        expected = DataArray(array.to_series().resample("3h").ffill())
        assert_identical(expected, actual)

        # Backward-fill
        actual = array.resample(time="3h").bfill()
        expected = DataArray(array.to_series().resample("3h").bfill())
        assert_identical(expected, actual)

        # As frequency
        actual = array.resample(time="3h").asfreq()
        expected = DataArray(array.to_series().resample("3h").asfreq())
        assert_identical(expected, actual)

        # Pad
        actual = array.resample(time="3h").pad()
        expected = DataArray(array.to_series().resample("3h").ffill())
        assert_identical(expected, actual)

        # Nearest
        rs = array.resample(time="3h")
        actual = rs.nearest()
        new_times = rs.groupers[0].full_index
        expected = DataArray(array.reindex(time=new_times, method="nearest"))
        assert_identical(expected, actual)

    def test_upsample_nd(self) -> None:
        # Same as before, but now we try on multi-dimensional DataArrays.
        xs = np.arange(6)
        ys = np.arange(3)
        times = pd.date_range("2000-01-01", freq="6h", periods=5)
        data = np.tile(np.arange(5), (6, 3, 1))
        array = DataArray(data, {"time": times, "x": xs, "y": ys}, ("x", "y", "time"))

        # Forward-fill
        actual = array.resample(time="3h").ffill()
        expected_data = np.repeat(data, 2, axis=-1)
        expected_times = times.to_series().resample("3h").asfreq().index
        expected_data = expected_data[..., : len(expected_times)]
        expected = DataArray(
            expected_data,
            {"time": expected_times, "x": xs, "y": ys},
            ("x", "y", "time"),
        )
        assert_identical(expected, actual)

        # Backward-fill
        actual = array.resample(time="3h").ffill()
        expected_data = np.repeat(np.flipud(data.T).T, 2, axis=-1)
        expected_data = np.flipud(expected_data.T).T
        expected_times = times.to_series().resample("3h").asfreq().index
        expected_data = expected_data[..., : len(expected_times)]
        expected = DataArray(
            expected_data,
            {"time": expected_times, "x": xs, "y": ys},
            ("x", "y", "time"),
        )
        assert_identical(expected, actual)

        # As frequency
        actual = array.resample(time="3h").asfreq()
        expected_data = np.repeat(data, 2, axis=-1).astype(float)[..., :-1]
        expected_data[..., 1::2] = np.nan
        expected_times = times.to_series().resample("3h").asfreq().index
        expected = DataArray(
            expected_data,
            {"time": expected_times, "x": xs, "y": ys},
            ("x", "y", "time"),
        )
        assert_identical(expected, actual)

        # Pad
        actual = array.resample(time="3h").pad()
        expected_data = np.repeat(data, 2, axis=-1)
        expected_data[..., 1::2] = expected_data[..., ::2]
        expected_data = expected_data[..., :-1]
        expected_times = times.to_series().resample("3h").asfreq().index
        expected = DataArray(
            expected_data,
            {"time": expected_times, "x": xs, "y": ys},
            ("x", "y", "time"),
        )
        assert_identical(expected, actual)

    def test_upsample_tolerance(self) -> None:
        # Test tolerance keyword for upsample methods bfill, pad, nearest
        times = pd.date_range("2000-01-01", freq="1D", periods=2)
        times_upsampled = pd.date_range("2000-01-01", freq="6h", periods=5)
        array = DataArray(np.arange(2), [("time", times)])

        # Forward fill
        actual = array.resample(time="6h").ffill(tolerance="12h")
        expected = DataArray([0.0, 0.0, 0.0, np.nan, 1.0], [("time", times_upsampled)])
        assert_identical(expected, actual)

        # Backward fill
        actual = array.resample(time="6h").bfill(tolerance="12h")
        expected = DataArray([0.0, np.nan, 1.0, 1.0, 1.0], [("time", times_upsampled)])
        assert_identical(expected, actual)

        # Nearest
        actual = array.resample(time="6h").nearest(tolerance="6h")
        expected = DataArray([0, 0, np.nan, 1, 1], [("time", times_upsampled)])
        assert_identical(expected, actual)

    @requires_scipy
    def test_upsample_interpolate(self) -> None:
        from scipy.interpolate import interp1d

        xs = np.arange(6)
        ys = np.arange(3)
        times = pd.date_range("2000-01-01", freq="6h", periods=5)

        z = np.arange(5) ** 2
        data = np.tile(z, (6, 3, 1))
        array = DataArray(data, {"time": times, "x": xs, "y": ys}, ("x", "y", "time"))

        expected_times = times.to_series().resample("1h").asfreq().index
        # Split the times into equal sub-intervals to simulate the 6 hour
        # to 1 hour up-sampling
        new_times_idx = np.linspace(0, len(times) - 1, len(times) * 5)
        kinds: list[InterpOptions] = [
            "linear",
            "nearest",
            "zero",
            "slinear",
            "quadratic",
            "cubic",
            "polynomial",
        ]
        for kind in kinds:
            kwargs = {}
            if kind == "polynomial":
                kwargs["order"] = 1
            actual = array.resample(time="1h").interpolate(kind, **kwargs)
            # using interp1d, polynomial order is to set directly in kind using int
            f = interp1d(
                np.arange(len(times)),
                data,
                kind=kwargs["order"] if kind == "polynomial" else kind,  # type: ignore[arg-type,unused-ignore]
                axis=-1,
                bounds_error=True,
                assume_sorted=True,
            )
            expected_data = f(new_times_idx)
            expected = DataArray(
                expected_data,
                {"time": expected_times, "x": xs, "y": ys},
                ("x", "y", "time"),
            )
            # Use AllClose because there are some small differences in how
            # we upsample timeseries versus the integer indexing as I've
            # done here due to floating point arithmetic
            assert_allclose(expected, actual, rtol=1e-16)

    @requires_scipy
    def test_upsample_interpolate_bug_2197(self) -> None:
        dates = pd.date_range("2007-02-01", "2007-03-01", freq="D", unit="s")
        da = xr.DataArray(np.arange(len(dates)), [("time", dates)])
        result = da.resample(time="ME").interpolate("linear")
        expected_times = np.array(
            [np.datetime64("2007-02-28"), np.datetime64("2007-03-31")]
        )
        expected = xr.DataArray([27.0, np.nan], [("time", expected_times)])
        assert_equal(result, expected)

    @requires_scipy
    def test_upsample_interpolate_regression_1605(self) -> None:
        dates = pd.date_range("2016-01-01", "2016-03-31", freq="1D")
        expected = xr.DataArray(
            np.random.random((len(dates), 2, 3)),
            dims=("time", "x", "y"),
            coords={"time": dates},
        )
        actual = expected.resample(time="1D").interpolate("linear")
        assert_allclose(actual, expected, rtol=1e-16)

    @requires_dask
    @requires_scipy
    @pytest.mark.parametrize("chunked_time", [True, False])
    def test_upsample_interpolate_dask(self, chunked_time: bool) -> None:
        from scipy.interpolate import interp1d

        xs = np.arange(6)
        ys = np.arange(3)
        times = pd.date_range("2000-01-01", freq="6h", periods=5)

        z = np.arange(5) ** 2
        data = np.tile(z, (6, 3, 1))
        array = DataArray(data, {"time": times, "x": xs, "y": ys}, ("x", "y", "time"))
        chunks = {"x": 2, "y": 1}
        if chunked_time:
            chunks["time"] = 3

        expected_times = times.to_series().resample("1h").asfreq().index
        # Split the times into equal sub-intervals to simulate the 6 hour
        # to 1 hour up-sampling
        new_times_idx = np.linspace(0, len(times) - 1, len(times) * 5)
        kinds: list[InterpOptions] = [
            "linear",
            "nearest",
            "zero",
            "slinear",
            "quadratic",
            "cubic",
            "polynomial",
        ]
        for kind in kinds:
            kwargs = {}
            if kind == "polynomial":
                kwargs["order"] = 1
            actual = array.chunk(chunks).resample(time="1h").interpolate(kind, **kwargs)
            actual = actual.compute()
            # using interp1d, polynomial order is to set directly in kind using int
            f = interp1d(
                np.arange(len(times)),
                data,
                kind=kwargs["order"] if kind == "polynomial" else kind,  # type: ignore[arg-type,unused-ignore]
                axis=-1,
                bounds_error=True,
                assume_sorted=True,
            )
            expected_data = f(new_times_idx)
            expected = DataArray(
                expected_data,
                {"time": expected_times, "x": xs, "y": ys},
                ("x", "y", "time"),
            )
            # Use AllClose because there are some small differences in how
            # we upsample timeseries versus the integer indexing as I've
            # done here due to floating point arithmetic
            assert_allclose(expected, actual, rtol=1e-16)

    def test_resample_offset(self) -> None:
        times = pd.date_range("2000-01-01T02:03:01", freq="6h", periods=10)
        array = DataArray(np.arange(10), [("time", times)])

        offset = pd.Timedelta("11h")
        actual = array.resample(time="24h", offset=offset).mean()
        expected = DataArray(array.to_series().resample("24h", offset=offset).mean())
        assert_identical(expected, actual)

    def test_resample_origin(self) -> None:
        times = pd.date_range("2000-01-01T02:03:01", freq="6h", periods=10)
        array = DataArray(np.arange(10), [("time", times)])

        origin: Literal["start"] = "start"
        actual = array.resample(time="24h", origin=origin).mean()
        expected = DataArray(array.to_series().resample("24h", origin=origin).mean())
        assert_identical(expected, actual)


class TestDatasetResample:
    @pytest.mark.parametrize(
        "resample_freq",
        [
            "24h",
            "123456s",
            "1234567890us",
            pd.Timedelta(hours=2),
            pd.offsets.MonthBegin(),
            pd.offsets.Second(123456),
            datetime.timedelta(days=1, hours=6),
        ],
    )
    def test_resample(
        self, use_cftime: bool, resample_freq: ResampleCompatible
    ) -> None:
        if use_cftime and not has_cftime:
            pytest.skip()
        times = xr.date_range(
            "2000-01-01", freq="6h", periods=10, use_cftime=use_cftime
        )

        def resample_as_pandas(ds, *args, **kwargs):
            ds_ = ds.copy(deep=True)
            if use_cftime:
                ds_["time"] = times.to_datetimeindex(time_unit="ns")
            result = Dataset.from_dataframe(
                ds_.to_dataframe().resample(*args, **kwargs).mean()
            )
            if use_cftime:
                result = result.convert_calendar(
                    calendar="standard", use_cftime=use_cftime
                )
            return result

        ds = Dataset(
            {
                "foo": ("time", np.random.randint(1, 1000, 10)),
                "bar": ("time", np.random.randint(1, 1000, 10)),
                "time": times,
            }
        )

        actual = ds.resample(time=resample_freq).mean()
        expected = resample_as_pandas(ds, resample_freq)
        assert_identical(expected, actual)

        actual = ds.resample(time=resample_freq).reduce(np.mean)
        assert_identical(expected, actual)

        actual = ds.resample(time=resample_freq, closed="right").mean()
        expected = resample_as_pandas(ds, resample_freq, closed="right")
        assert_identical(expected, actual)

        with pytest.raises(ValueError, match=r"Index must be monotonic"):
            ds.isel(time=[2, 0, 1]).resample(time=resample_freq)

        reverse = ds.isel(time=slice(-1, None, -1))
        with pytest.raises(ValueError):
            reverse.resample(time=resample_freq).mean()

    def test_resample_and_first(self) -> None:
        times = pd.date_range("2000-01-01", freq="6h", periods=10)
        ds = Dataset(
            {
                "foo": (["time", "x", "y"], np.random.randn(10, 5, 3)),
                "bar": ("time", np.random.randn(10), {"meta": "data"}),
                "time": times,
            }
        )

        actual = ds.resample(time="1D").first(keep_attrs=True)
        expected = ds.isel(time=[0, 4, 8])
        assert_identical(expected, actual)

        # upsampling
        expected_time = pd.date_range("2000-01-01", freq="3h", periods=19)
        expected = ds.reindex(time=expected_time)
        rs = ds.resample(time="3h")
        for how in ["mean", "sum", "first", "last"]:
            method = getattr(rs, how)
            result = method()
            assert_equal(expected, result)
        for method in [np.mean]:
            result = rs.reduce(method)
            assert_equal(expected, result)

    def test_resample_min_count(self) -> None:
        times = pd.date_range("2000-01-01", freq="6h", periods=10)
        ds = Dataset(
            {
                "foo": (["time", "x", "y"], np.random.randn(10, 5, 3)),
                "bar": ("time", np.random.randn(10), {"meta": "data"}),
                "time": times,
            }
        )
        # inject nan
        ds["foo"] = xr.where(ds["foo"] > 2.0, np.nan, ds["foo"])

        actual = ds.resample(time="1D").sum(min_count=1)
        expected = xr.concat(
            [
                ds.isel(time=slice(i * 4, (i + 1) * 4)).sum("time", min_count=1)
                for i in range(3)
            ],
            dim=actual["time"],
            data_vars="all",
        )
        assert_allclose(expected, actual)

    def test_resample_by_mean_with_keep_attrs(self) -> None:
        times = pd.date_range("2000-01-01", freq="6h", periods=10)
        ds = Dataset(
            {
                "foo": (["time", "x", "y"], np.random.randn(10, 5, 3)),
                "bar": ("time", np.random.randn(10), {"meta": "data"}),
                "time": times,
            }
        )
        ds.attrs["dsmeta"] = "dsdata"

        resampled_ds = ds.resample(time="1D").mean(keep_attrs=True)
        actual = resampled_ds["bar"].attrs
        expected = ds["bar"].attrs
        assert expected == actual

        actual = resampled_ds.attrs
        expected = ds.attrs
        assert expected == actual

    def test_resample_by_mean_discarding_attrs(self) -> None:
        times = pd.date_range("2000-01-01", freq="6h", periods=10)
        ds = Dataset(
            {
                "foo": (["time", "x", "y"], np.random.randn(10, 5, 3)),
                "bar": ("time", np.random.randn(10), {"meta": "data"}),
                "time": times,
            }
        )
        ds.attrs["dsmeta"] = "dsdata"

        resampled_ds = ds.resample(time="1D").mean(keep_attrs=False)

        assert resampled_ds["bar"].attrs == {}
        assert resampled_ds.attrs == {}

    def test_resample_by_last_discarding_attrs(self) -> None:
        times = pd.date_range("2000-01-01", freq="6h", periods=10)
        ds = Dataset(
            {
                "foo": (["time", "x", "y"], np.random.randn(10, 5, 3)),
                "bar": ("time", np.random.randn(10), {"meta": "data"}),
                "time": times,
            }
        )
        ds.attrs["dsmeta"] = "dsdata"

        resampled_ds = ds.resample(time="1D").last(keep_attrs=False)

        assert resampled_ds["bar"].attrs == {}
        assert resampled_ds.attrs == {}

    @requires_scipy
    def test_resample_drop_nondim_coords(self) -> None:
        xs = np.arange(6)
        ys = np.arange(3)
        times = pd.date_range("2000-01-01", freq="6h", periods=5)
        data = np.tile(np.arange(5), (6, 3, 1))
        xx, yy = np.meshgrid(xs * 5, ys * 2.5)
        tt = np.arange(len(times), dtype=int)
        array = DataArray(data, {"time": times, "x": xs, "y": ys}, ("x", "y", "time"))
        xcoord = DataArray(xx.T, {"x": xs, "y": ys}, ("x", "y"))
        ycoord = DataArray(yy.T, {"x": xs, "y": ys}, ("x", "y"))
        tcoord = DataArray(tt, {"time": times}, ("time",))
        ds = Dataset({"data": array, "xc": xcoord, "yc": ycoord, "tc": tcoord})
        ds = ds.set_coords(["xc", "yc", "tc"])

        # Re-sample
        actual = ds.resample(time="12h").mean("time")
        assert "tc" not in actual.coords

        # Up-sample - filling
        actual = ds.resample(time="1h").ffill()
        assert "tc" not in actual.coords

        # Up-sample - interpolation
        actual = ds.resample(time="1h").interpolate("linear")
        assert "tc" not in actual.coords

    def test_resample_ds_da_are_the_same(self) -> None:
        time = pd.date_range("2000-01-01", freq="6h", periods=365 * 4)
        ds = xr.Dataset(
            {
                "foo": (("time", "x"), np.random.randn(365 * 4, 5)),
                "time": time,
                "x": np.arange(5),
            }
        )
        assert_allclose(
            ds.resample(time="ME").mean()["foo"], ds.foo.resample(time="ME").mean()
        )

    def test_ds_resample_apply_func_args(self) -> None:
        def func(arg1, arg2, arg3=0.0):
            return arg1.mean("time") + arg2 + arg3

        times = pd.date_range("2000", freq="D", periods=3)
        ds = xr.Dataset({"foo": ("time", [1.0, 1.0, 1.0]), "time": times})
        expected = xr.Dataset({"foo": ("time", [3.0, 3.0, 3.0]), "time": times})
        actual = ds.resample(time="D").map(func, args=(1.0,), arg3=1.0)
        assert_identical(expected, actual)


def test_groupby_cumsum() -> None:
    ds = xr.Dataset(
        {"foo": (("x",), [7, 3, 1, 1, 1, 1, 1])},
        coords={"x": [0, 1, 2, 3, 4, 5, 6], "group_id": ("x", [0, 0, 1, 1, 2, 2, 2])},
    )
    actual = ds.groupby("group_id").cumsum(dim="x")
    expected = xr.Dataset(
        {
            "foo": (("x",), [7, 10, 1, 2, 1, 2, 3]),
        },
        coords={
            "x": [0, 1, 2, 3, 4, 5, 6],
            "group_id": ds.group_id,
        },
    )
    # TODO: Remove drop_vars when GH6528 is fixed
    # when Dataset.cumsum propagates indexes, and the group variable?
    assert_identical(expected.drop_vars(["x", "group_id"]), actual)

    actual = ds.foo.groupby("group_id").cumsum(dim="x")
    expected.coords["group_id"] = ds.group_id
    expected.coords["x"] = np.arange(7)
    assert_identical(expected.foo, actual)


def test_groupby_cumprod() -> None:
    ds = xr.Dataset(
        {"foo": (("x",), [7, 3, 0, 1, 1, 2, 1])},
        coords={"x": [0, 1, 2, 3, 4, 5, 6], "group_id": ("x", [0, 0, 1, 1, 2, 2, 2])},
    )
    actual = ds.groupby("group_id").cumprod(dim="x")
    expected = xr.Dataset(
        {
            "foo": (("x",), [7, 21, 0, 0, 1, 2, 2]),
        },
        coords={
            "x": [0, 1, 2, 3, 4, 5, 6],
            "group_id": ds.group_id,
        },
    )
    # TODO: Remove drop_vars when GH6528 is fixed
    # when Dataset.cumsum propagates indexes, and the group variable?
    assert_identical(expected.drop_vars(["x", "group_id"]), actual)

    actual = ds.foo.groupby("group_id").cumprod(dim="x")
    expected.coords["group_id"] = ds.group_id
    expected.coords["x"] = np.arange(7)
    assert_identical(expected.foo, actual)


@pytest.mark.parametrize(
    "method, expected_array",
    [
        ("cumsum", [1.0, 2.0, 5.0, 6.0, 2.0, 2.0]),
        ("cumprod", [1.0, 2.0, 6.0, 6.0, 2.0, 2.0]),
    ],
)
def test_resample_cumsum(method: str, expected_array: list[float]) -> None:
    ds = xr.Dataset(
        {"foo": ("time", [1, 2, 3, 1, 2, np.nan])},
        coords={
            "time": xr.date_range("01-01-2001", freq="ME", periods=6, use_cftime=False),
        },
    )
    actual = getattr(ds.resample(time="3ME"), method)(dim="time")
    expected = xr.Dataset(
        {"foo": (("time",), expected_array)},
        coords={
            "time": xr.date_range("01-01-2001", freq="ME", periods=6, use_cftime=False),
        },
    )
    # TODO: Remove drop_vars when GH6528 is fixed
    # when Dataset.cumsum propagates indexes, and the group variable?
    assert_identical(expected.drop_vars(["time"]), actual)

    actual = getattr(ds.foo.resample(time="3ME"), method)(dim="time")
    expected.coords["time"] = ds.time
    assert_identical(expected.drop_vars(["time"]).foo, actual)


def test_groupby_binary_op_regression() -> None:
    # regression test for #7797
    # monthly timeseries that should return "zero anomalies" everywhere
    time = xr.date_range("2023-01-01", "2023-12-31", freq="MS")
    data = np.linspace(-1, 1, 12)
    x = xr.DataArray(data, coords={"time": time})
    clim = xr.DataArray(data, coords={"month": np.arange(1, 13, 1)})

    # seems to give the correct result if we use the full x, but not with a slice
    x_slice = x.sel(time=["2023-04-01"])

    # two typical ways of computing anomalies
    anom_gb = x_slice.groupby("time.month") - clim

    assert_identical(xr.zeros_like(anom_gb), anom_gb)


def test_groupby_multiindex_level() -> None:
    # GH6836
    midx = pd.MultiIndex.from_product([list("abc"), [0, 1]], names=("one", "two"))
    mda = xr.DataArray(np.random.rand(6, 3), [("x", midx), ("y", range(3))])
    groups = mda.groupby("one").groups
    assert groups == {"a": [0, 1], "b": [2, 3], "c": [4, 5]}


@requires_flox
@pytest.mark.parametrize("func", ["sum", "prod"])
@pytest.mark.parametrize("skipna", [True, False])
@pytest.mark.parametrize("min_count", [None, 1])
def test_min_count_vs_flox(func: str, min_count: int | None, skipna: bool) -> None:
    da = DataArray(
        data=np.array([np.nan, 1, 1, np.nan, 1, 1]),
        dims="x",
        coords={"labels": ("x", np.array([1, 2, 3, 1, 2, 3]))},
    )

    gb = da.groupby("labels")
    method = operator.methodcaller(func, min_count=min_count, skipna=skipna)
    with xr.set_options(use_flox=True):
        actual = method(gb)
    with xr.set_options(use_flox=False):
        expected = method(gb)
    assert_identical(actual, expected)


@pytest.mark.parametrize("use_flox", [True, False])
def test_min_count_error(use_flox: bool) -> None:
    if use_flox and not has_flox:
        pytest.skip()
    da = DataArray(
        data=np.array([np.nan, 1, 1, np.nan, 1, 1]),
        dims="x",
        coords={"labels": ("x", np.array([1, 2, 3, 1, 2, 3]))},
    )
    with xr.set_options(use_flox=use_flox):
        with pytest.raises(TypeError):
            da.groupby("labels").mean(min_count=1)


@requires_dask
def test_groupby_math_auto_chunk() -> None:
    da = xr.DataArray(
        [[1, 2, 3], [1, 2, 3], [1, 2, 3]],
        dims=("y", "x"),
        coords={"label": ("x", [2, 2, 1])},
    )
    sub = xr.DataArray(
        InaccessibleArray(np.array([1, 2])), dims="label", coords={"label": [1, 2]}
    )
    chunked = da.chunk(x=1, y=2)
    chunked.label.load()
    actual = chunked.groupby("label") - sub
    assert actual.chunksizes == {"x": (1, 1, 1), "y": (2, 1)}


@pytest.mark.parametrize("use_flox", [True, False])
def test_groupby_dim_no_dim_equal(use_flox: bool) -> None:
    # https://github.com/pydata/xarray/issues/8263
    da = DataArray(
        data=[1, 2, 3, 4], dims="lat", coords={"lat": np.linspace(0, 1.01, 4)}
    )
    with xr.set_options(use_flox=use_flox):
        actual1 = da.drop_vars("lat").groupby("lat").sum()
        actual2 = da.groupby("lat").sum()
    assert_identical(actual1, actual2.drop_vars("lat"))


@requires_flox
def test_default_flox_method() -> None:
    import flox.xarray

    da = xr.DataArray([1, 2, 3], dims="x", coords={"label": ("x", [2, 2, 1])})

    result = xr.DataArray([3, 3], dims="label", coords={"label": [1, 2]})
    with mock.patch("flox.xarray.xarray_reduce", return_value=result) as mocked_reduce:
        da.groupby("label").sum()

    kwargs = mocked_reduce.call_args.kwargs
    if Version(flox.__version__) < Version("0.9.0"):
        assert kwargs["method"] == "cohorts"
    else:
        assert "method" not in kwargs


@requires_cftime
@pytest.mark.filterwarnings("ignore")
def test_cftime_resample_gh_9108() -> None:
    import cftime

    ds = Dataset(
        {"pr": ("time", np.random.random((10,)))},
        coords={"time": xr.date_range("0001-01-01", periods=10, freq="D")},
    )
    actual = ds.resample(time="ME").mean()
    expected = ds.mean("time").expand_dims(
        time=[cftime.DatetimeGregorian(1, 1, 31, 0, 0, 0, 0, has_year_zero=False)]
    )
    assert actual.time.data[0].has_year_zero == ds.time.data[0].has_year_zero
    assert_equal(actual, expected)


def test_custom_grouper() -> None:
    class YearGrouper(Grouper):
        """
        An example re-implementation of ``.groupby("time.year")``.
        """

        def factorize(self, group) -> EncodedGroups:
            assert np.issubdtype(group.dtype, np.datetime64)
            year = group.dt.year.data
            codes_, uniques = pd.factorize(year)
            codes = group.copy(data=codes_).rename("year")
            return EncodedGroups(codes=codes, full_index=pd.Index(uniques))

        def reset(self):
            return type(self)()

    da = xr.DataArray(
        dims="time",
        data=np.arange(20),
        coords={"time": ("time", pd.date_range("2000-01-01", freq="3MS", periods=20))},
        name="foo",
    )
    ds = da.to_dataset()

    expected = ds.groupby("time.year").mean()
    actual = ds.groupby(time=YearGrouper()).mean()
    assert_identical(expected, actual)

    actual = ds.groupby({"time": YearGrouper()}).mean()
    assert_identical(expected, actual)

    expected = ds.foo.groupby("time.year").mean()
    actual = ds.foo.groupby(time=YearGrouper()).mean()
    assert_identical(expected, actual)

    actual = ds.foo.groupby({"time": YearGrouper()}).mean()
    assert_identical(expected, actual)

    for obj in [ds, ds.foo]:
        with pytest.raises(ValueError):
            obj.groupby("time.year", time=YearGrouper())
        with pytest.raises(ValueError):
            obj.groupby()


@pytest.mark.parametrize("use_flox", [True, False])
def test_weather_data_resample(use_flox):
    # from the docs
    times = pd.date_range("2000-01-01", "2001-12-31", name="time")
    annual_cycle = np.sin(2 * np.pi * (times.dayofyear.values / 365.25 - 0.28))

    base = 10 + 15 * annual_cycle.reshape(-1, 1)
    tmin_values = base + 3 * np.random.randn(annual_cycle.size, 3)
    tmax_values = base + 10 + 3 * np.random.randn(annual_cycle.size, 3)

    ds = xr.Dataset(
        {
            "tmin": (("time", "location"), tmin_values),
            "tmax": (("time", "location"), tmax_values),
        },
        {
            "time": ("time", times, {"time_key": "time_values"}),
            "location": ("location", ["IA", "IN", "IL"], {"loc_key": "loc_value"}),
        },
    )

    with xr.set_options(use_flox=use_flox):
        actual = ds.resample(time="1MS").mean()
    assert "location" in actual._indexes

    gb = ds.groupby(time=TimeResampler(freq="1MS"), location=UniqueGrouper())
    with xr.set_options(use_flox=use_flox):
        actual = gb.mean()
        expected = ds.resample(time="1MS").mean().sortby("location")
    assert_allclose(actual, expected)
    assert actual.time.attrs == ds.time.attrs
    assert actual.location.attrs == ds.location.attrs

    assert expected.time.attrs == ds.time.attrs
    assert expected.location.attrs == ds.location.attrs


@pytest.mark.parametrize("as_dataset", [True, False])
def test_multiple_groupers_string(as_dataset) -> None:
    obj = DataArray(
        np.array([1, 2, 3, 0, 2, np.nan]),
        dims="d",
        coords=dict(
            labels1=("d", np.array(["a", "b", "c", "c", "b", "a"])),
            labels2=("d", np.array(["x", "y", "z", "z", "y", "x"])),
        ),
        name="foo",
    )

    if as_dataset:
        obj = obj.to_dataset()  # type: ignore[assignment]

    expected = obj.groupby(labels1=UniqueGrouper(), labels2=UniqueGrouper()).mean()
    actual = obj.groupby(("labels1", "labels2")).mean()
    assert_identical(expected, actual)

    # Passes `"labels2"` to squeeze; will raise an error around kwargs rather than the
    # warning & type error in the future
    with pytest.warns(FutureWarning):
        with pytest.raises(TypeError):
            obj.groupby("labels1", "labels2")  # type: ignore[arg-type, misc]
    with pytest.raises(ValueError):
        obj.groupby("labels1", foo="bar")  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        obj.groupby("labels1", foo=UniqueGrouper())


@pytest.mark.parametrize("shuffle", [True, False])
@pytest.mark.parametrize("use_flox", [True, False])
def test_multiple_groupers(use_flox: bool, shuffle: bool) -> None:
    da = DataArray(
        np.array([1, 2, 3, 0, 2, np.nan]),
        dims="d",
        coords=dict(
            labels1=("d", np.array(["a", "b", "c", "c", "b", "a"])),
            labels2=("d", np.array(["x", "y", "z", "z", "y", "x"])),
        ),
        name="foo",
    )

    groupers: dict[str, Grouper]
    groupers = dict(labels1=UniqueGrouper(), labels2=UniqueGrouper())
    gb = da.groupby(groupers)
    if shuffle:
        gb = gb.shuffle_to_chunks().groupby(groupers)
    repr(gb)

    expected = DataArray(
        np.array([[1.0, np.nan, np.nan], [np.nan, 2.0, np.nan], [np.nan, np.nan, 1.5]]),
        dims=("labels1", "labels2"),
        coords={
            "labels1": np.array(["a", "b", "c"], dtype=object),
            "labels2": np.array(["x", "y", "z"], dtype=object),
        },
        name="foo",
    )
    with xr.set_options(use_flox=use_flox):
        actual = gb.mean()
    assert_identical(actual, expected)

    # -------
    coords = {"a": ("x", [0, 0, 1, 1]), "b": ("y", [0, 0, 1, 1])}
    square = DataArray(np.arange(16).reshape(4, 4), coords=coords, dims=["x", "y"])
    groupers = dict(a=UniqueGrouper(), b=UniqueGrouper())
    gb = square.groupby(groupers)
    if shuffle:
        gb = gb.shuffle_to_chunks().groupby(groupers)
    repr(gb)
    with xr.set_options(use_flox=use_flox):
        actual = gb.mean()
    expected = DataArray(
        np.array([[2.5, 4.5], [10.5, 12.5]]),
        dims=("a", "b"),
        coords={"a": [0, 1], "b": [0, 1]},
    )
    assert_identical(actual, expected)

    expected = square.astype(np.float64)
    expected["a"], expected["b"] = broadcast(square.a, square.b)
    with xr.set_options(use_flox=use_flox):
        assert_identical(
            square.groupby(x=UniqueGrouper(), y=UniqueGrouper()).mean(), expected
        )

    b = xr.DataArray(
        np.random.default_rng(0).random((2, 3, 4)),
        coords={"xy": (("x", "y"), [["a", "b", "c"], ["b", "c", "c"]], {"foo": "bar"})},
        dims=["x", "y", "z"],
    )
    groupers = dict(x=UniqueGrouper(), y=UniqueGrouper())
    gb = b.groupby(groupers)
    if shuffle:
        gb = gb.shuffle_to_chunks().groupby(groupers)
    repr(gb)
    with xr.set_options(use_flox=use_flox):
        assert_identical(gb.mean("z"), b.mean("z"))

    groupers = dict(x=UniqueGrouper(), xy=UniqueGrouper())
    gb = b.groupby(groupers)
    if shuffle:
        gb = gb.shuffle_to_chunks().groupby(groupers)
    repr(gb)
    with xr.set_options(use_flox=use_flox):
        actual = gb.mean()
    expected = b.drop_vars("xy").rename({"y": "xy"}).copy(deep=True)
    newval = b.isel(x=1, y=slice(1, None)).mean("y").data
    expected.loc[dict(x=1, xy=1)] = expected.sel(x=1, xy=0).data
    expected.loc[dict(x=1, xy=0)] = np.nan
    expected.loc[dict(x=1, xy=2)] = newval
    expected["xy"] = ("xy", ["a", "b", "c"], {"foo": "bar"})
    # TODO: is order of dims correct?
    assert_identical(actual, expected.transpose("z", "x", "xy"))

    if has_dask:
        b["xy"] = b["xy"].chunk()
        expected = xr.DataArray(
            [[[1, 1, 1], [np.nan, 1, 2]]] * 4,
            dims=("z", "x", "xy"),
            coords={"xy": ("xy", ["a", "b", "c"], {"foo": "bar"})},
        )
        with raise_if_dask_computes(max_computes=0):
            gb = b.groupby(x=UniqueGrouper(), xy=UniqueGrouper(labels=["a", "b", "c"]))
        assert is_chunked_array(gb.encoded.codes.data)
        assert not gb.encoded.group_indices
        if has_flox:
            with raise_if_dask_computes(max_computes=1):
                assert_identical(gb.count(), expected)
        else:
            with pytest.raises(ValueError, match="when lazily grouping"):
                gb.count()


@pytest.mark.parametrize("use_flox", [True, False])
@pytest.mark.parametrize("shuffle", [True, False])
def test_multiple_groupers_mixed(use_flox: bool, shuffle: bool) -> None:
    # This groupby has missing groups
    ds = xr.Dataset(
        {"foo": (("x", "y"), np.arange(12).reshape((4, 3)))},
        coords={"x": [10, 20, 30, 40], "letters": ("x", list("abba"))},
    )
    groupers: dict[str, Grouper] = dict(
        x=BinGrouper(bins=[5, 15, 25]), letters=UniqueGrouper()
    )
    gb = ds.groupby(groupers)
    if shuffle:
        gb = gb.shuffle_to_chunks().groupby(groupers)
    expected_data = np.array(
        [
            [[0.0, np.nan], [np.nan, 3.0]],
            [[1.0, np.nan], [np.nan, 4.0]],
            [[2.0, np.nan], [np.nan, 5.0]],
        ]
    )
    expected = xr.Dataset(
        {"foo": (("y", "x_bins", "letters"), expected_data)},
        coords={
            "x_bins": (
                "x_bins",
                np.array(
                    [
                        pd.Interval(5, 15, closed="right"),
                        pd.Interval(15, 25, closed="right"),
                    ],
                    dtype=object,
                ),
            ),
            "letters": ("letters", np.array(["a", "b"], dtype=object)),
        },
    )
    with xr.set_options(use_flox=use_flox):
        actual = gb.sum()
    assert_identical(actual, expected)

    # assert_identical(
    #     b.groupby(['x', 'y']).apply(lambda x: x - x.mean()),
    #     b - b.mean("z"),
    # )

    # gb = square.groupby(x=UniqueGrouper(), y=UniqueGrouper())
    # gb - gb.mean()

    # ------


@requires_flox_0_9_12
@pytest.mark.parametrize(
    "reduction", ["max", "min", "nanmax", "nanmin", "sum", "nansum", "prod", "nanprod"]
)
def test_groupby_preserve_dtype(reduction):
    # all groups are present, we should follow numpy exactly
    ds = xr.Dataset(
        {
            "test": (
                ["x", "y"],
                np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]], dtype="int16"),
            )
        },
        coords={"idx": ("x", [1, 2, 1])},
    )

    kwargs = {}
    if "nan" in reduction:
        kwargs["skipna"] = True
    # TODO: fix dtype with numbagg/bottleneck and use_flox=False
    with xr.set_options(use_numbagg=False, use_bottleneck=False):
        actual = getattr(ds.groupby("idx"), reduction.removeprefix("nan"))(
            **kwargs
        ).test.dtype
    expected = getattr(np, reduction)(ds.test.data, axis=0).dtype

    assert actual == expected


@requires_dask
@requires_flox_0_9_12
@pytest.mark.parametrize("reduction", ["any", "all", "count"])
def test_gappy_resample_reductions(reduction):
    # GH8090
    dates = (("1988-12-01", "1990-11-30"), ("2000-12-01", "2001-11-30"))
    times = [xr.date_range(*d, freq="D") for d in dates]

    da = xr.concat(
        [
            xr.DataArray(np.random.rand(len(t)), coords={"time": t}, dims="time")
            for t in times
        ],
        dim="time",
    ).chunk(time=100)

    rs = (da > 0.5).resample(time="YS-DEC")
    method = getattr(rs, reduction)
    with xr.set_options(use_flox=True):
        actual = method(dim="time")
    with xr.set_options(use_flox=False):
        expected = method(dim="time")
    assert_identical(expected, actual)


def test_groupby_transpose() -> None:
    # GH5361
    data = xr.DataArray(
        np.random.randn(4, 2),
        dims=["x", "z"],
        coords={"x": ["a", "b", "a", "c"], "y": ("x", [0, 1, 0, 2])},
    )
    first = data.T.groupby("x").sum()
    second = data.groupby("x").sum()

    assert_identical(first, second.transpose(*first.dims))


@requires_dask
@pytest.mark.parametrize(
    "grouper, expect_index",
    [
        [UniqueGrouper(labels=np.arange(1, 5)), pd.Index(np.arange(1, 5))],
        [UniqueGrouper(labels=np.arange(1, 5)[::-1]), pd.Index(np.arange(1, 5)[::-1])],
        [
            BinGrouper(bins=np.arange(1, 5)),
            pd.IntervalIndex.from_breaks(np.arange(1, 5)),
        ],
    ],
)
def test_lazy_grouping(grouper, expect_index):
    import dask.array

    data = DataArray(
        dims=("x", "y"),
        data=dask.array.arange(20, chunks=3).reshape((4, 5)),
        name="zoo",
    )
    with raise_if_dask_computes():
        encoded = grouper.factorize(data)
    assert encoded.codes.ndim == data.ndim
    pd.testing.assert_index_equal(encoded.full_index, expect_index)
    np.testing.assert_array_equal(encoded.unique_coord.values, np.array(expect_index))

    eager = (
        xr.Dataset({"foo": data}, coords={"zoo": data.compute()})
        .groupby(zoo=grouper)
        .count()
    )
    expected = Dataset(
        {"foo": (encoded.codes.name, np.ones(encoded.full_index.size))},
        coords={encoded.codes.name: expect_index},
    )
    assert_identical(eager, expected)

    if has_flox:
        lazy = (
            xr.Dataset({"foo": data}, coords={"zoo": data}).groupby(zoo=grouper).count()
        )
        assert_identical(eager, lazy)


@requires_dask
def test_lazy_grouping_errors() -> None:
    import dask.array

    data = DataArray(
        dims=("x",),
        data=dask.array.arange(20, chunks=3),
        name="foo",
        coords={"y": ("x", dask.array.arange(20, chunks=3))},
    )

    gb = data.groupby(y=UniqueGrouper(labels=np.arange(5, 10)))
    message = "not supported when lazily grouping by"
    with pytest.raises(ValueError, match=message):
        gb.map(lambda x: x)

    with pytest.raises(ValueError, match=message):
        gb.reduce(np.mean)

    with pytest.raises(ValueError, match=message):
        for _, _ in gb:
            pass


@requires_dask
def test_lazy_int_bins_error() -> None:
    import dask.array

    with pytest.raises(ValueError, match="Bin edges must be provided"):
        with raise_if_dask_computes():
            _ = BinGrouper(bins=4).factorize(DataArray(dask.array.arange(3)))


def test_time_grouping_seasons_specified() -> None:
    time = xr.date_range("2001-01-01", "2002-01-01", freq="D")
    ds = xr.Dataset({"foo": np.arange(time.size)}, coords={"time": ("time", time)})
    labels = ["DJF", "MAM", "JJA", "SON"]
    actual = ds.groupby({"time.season": UniqueGrouper(labels=labels)}).sum()
    expected = ds.groupby("time.season").sum()
    assert_identical(actual, expected.reindex(season=labels))


def test_multiple_grouper_unsorted_order() -> None:
    time = xr.date_range("2001-01-01", "2003-01-01", freq="MS")
    ds = xr.Dataset({"foo": np.arange(time.size)}, coords={"time": ("time", time)})
    labels = ["DJF", "MAM", "JJA", "SON"]
    actual = ds.groupby(
        {
            "time.season": UniqueGrouper(labels=labels),
            "time.year": UniqueGrouper(labels=[2002, 2001]),
        }
    ).sum()
    expected = (
        ds.groupby({"time.season": UniqueGrouper(), "time.year": UniqueGrouper()})
        .sum()
        .reindex(season=labels, year=[2002, 2001])
    )
    assert_identical(actual, expected.reindex(season=labels))

    b = xr.DataArray(
        np.random.default_rng(0).random((2, 3, 4)),
        coords={"x": [0, 1], "y": [0, 1, 2]},
        dims=["x", "y", "z"],
    )
    actual2 = b.groupby(
        x=UniqueGrouper(labels=[1, 0]), y=UniqueGrouper(labels=[2, 0, 1])
    ).sum()
    expected2 = b.reindex(x=[1, 0], y=[2, 0, 1]).transpose("z", ...)
    assert_identical(actual2, expected2)


def test_groupby_multiple_bin_grouper_missing_groups() -> None:
    from numpy import nan

    ds = xr.Dataset(
        {"foo": (("z"), np.arange(12))},
        coords={"x": ("z", np.arange(12)), "y": ("z", np.arange(12))},
    )

    actual = ds.groupby(
        x=BinGrouper(np.arange(0, 13, 4)), y=BinGrouper(bins=np.arange(0, 16, 2))
    ).count()
    expected = Dataset(
        {
            "foo": (
                ("x_bins", "y_bins"),
                np.array(
                    [
                        [2.0, 2.0, nan, nan, nan, nan, nan],
                        [nan, nan, 2.0, 2.0, nan, nan, nan],
                        [nan, nan, nan, nan, 2.0, 1.0, nan],
                    ]
                ),
            )
        },
        coords={
            "x_bins": ("x_bins", pd.IntervalIndex.from_breaks(np.arange(0, 13, 4))),
            "y_bins": ("y_bins", pd.IntervalIndex.from_breaks(np.arange(0, 16, 2))),
        },
    )
    assert_identical(actual, expected)


@requires_dask_ge_2024_08_1
def test_shuffle_simple() -> None:
    import dask

    da = xr.DataArray(
        dims="x",
        data=dask.array.from_array([1, 2, 3, 4, 5, 6], chunks=2),
        coords={"label": ("x", ["a", "b", "c", "a", "b", "c"])},
    )
    actual = da.groupby(label=UniqueGrouper()).shuffle_to_chunks()
    expected = da.isel(x=[0, 3, 1, 4, 2, 5])
    assert_identical(actual, expected)

    with pytest.raises(ValueError):
        da.chunk(x=2, eagerly_load_group=False).groupby("label").shuffle_to_chunks()


@requires_dask_ge_2024_08_1
@pytest.mark.parametrize(
    "chunks, expected_chunks",
    [
        ((1,), (1, 3, 3, 3)),
        ((10,), (10,)),
    ],
)
def test_shuffle_by(chunks, expected_chunks):
    import dask.array

    da = xr.DataArray(
        dims="x",
        data=dask.array.arange(10, chunks=chunks),
        coords={"x": [1, 2, 3, 1, 2, 3, 1, 2, 3, 0]},
        name="a",
    )
    ds = da.to_dataset()

    for obj in [ds, da]:
        actual = obj.groupby(x=UniqueGrouper()).shuffle_to_chunks()
        assert_identical(actual, obj.sortby("x"))
        assert actual.chunksizes["x"] == expected_chunks


@requires_dask
def test_groupby_dask_eager_load_warnings() -> None:
    ds = xr.Dataset(
        {"foo": (("z"), np.arange(12))},
        coords={"x": ("z", np.arange(12)), "y": ("z", np.arange(12))},
    ).chunk(z=6)

    with pytest.raises(ValueError, match="Please pass"):
        with pytest.warns(DeprecationWarning):
            ds.groupby("x", eagerly_compute_group=False)
    with pytest.raises(ValueError, match="Eagerly computing"):
        ds.groupby("x", eagerly_compute_group=True)  # type: ignore[arg-type]

    # This is technically fine but anyone iterating over the groupby object
    # will see an error, so let's warn and have them opt-in.
    ds.groupby(x=UniqueGrouper(labels=[1, 2, 3]))

    with pytest.warns(DeprecationWarning):
        ds.groupby(x=UniqueGrouper(labels=[1, 2, 3]), eagerly_compute_group=False)

    with pytest.raises(ValueError, match="Please pass"):
        with pytest.warns(DeprecationWarning):
            ds.groupby_bins("x", bins=3, eagerly_compute_group=False)
    with pytest.raises(ValueError, match="Eagerly computing"):
        ds.groupby_bins("x", bins=3, eagerly_compute_group=True)  # type: ignore[arg-type]
    ds.groupby_bins("x", bins=[1, 2, 3])
    with pytest.warns(DeprecationWarning):
        ds.groupby_bins("x", bins=[1, 2, 3], eagerly_compute_group=False)


class TestSeasonGrouperAndResampler:
    def test_season_to_month_tuple(self):
        assert season_to_month_tuple(["JF", "MAM", "JJAS", "OND"]) == (
            (1, 2),
            (3, 4, 5),
            (6, 7, 8, 9),
            (10, 11, 12),
        )
        assert season_to_month_tuple(["DJFM", "AM", "JJAS", "ON"]) == (
            (12, 1, 2, 3),
            (4, 5),
            (6, 7, 8, 9),
            (10, 11),
        )

    def test_season_grouper_raises_error_if_months_are_not_valid_or_not_continuous(
        self,
    ):
        calendar = "standard"
        time = date_range("2001-01-01", "2002-12-30", freq="D", calendar=calendar)
        da = DataArray(np.ones(time.size), dims="time", coords={"time": time})

        with pytest.raises(KeyError, match="IN"):
            da.groupby(time=SeasonGrouper(["INVALID_SEASON"]))

        with pytest.raises(KeyError, match="MD"):
            da.groupby(time=SeasonGrouper(["MDF"]))

    @pytest.mark.parametrize("calendar", _ALL_CALENDARS)
    def test_season_grouper_with_months_spanning_calendar_year_using_same_year(
        self, calendar
    ):
        time = date_range("2001-01-01", "2002-12-30", freq="MS", calendar=calendar)
        # fmt: off
        data = np.array(
            [
                1.0, 1.25, 1.5, 1.75, 2.0, 1.1, 1.35, 1.6, 1.85, 1.2, 1.45, 1.7,
                1.95, 1.05, 1.3, 1.55, 1.8, 1.15, 1.4, 1.65, 1.9, 1.25, 1.5, 1.75,
            ]

        )
        # fmt: on
        da = DataArray(data, dims="time", coords={"time": time})
        da["year"] = da.time.dt.year

        actual = da.groupby(
            year=UniqueGrouper(), time=SeasonGrouper(["NDJFM", "AMJ"])
        ).mean()

        # Expected if the same year "ND" is used for seasonal grouping
        expected = xr.DataArray(
            data=np.array([[1.38, 1.616667], [1.51, 1.5]]),
            dims=["year", "season"],
            coords={"year": [2001, 2002], "season": ["NDJFM", "AMJ"]},
        )

        assert_allclose(expected, actual)

    @pytest.mark.parametrize("calendar", _ALL_CALENDARS)
    def test_season_grouper_with_partial_years(self, calendar):
        time = date_range("2001-01-01", "2002-06-30", freq="MS", calendar=calendar)
        # fmt: off
        data = np.array(
            [
                1.0, 1.25, 1.5, 1.75, 2.0, 1.1, 1.35, 1.6, 1.85, 1.2, 1.45, 1.7,
                1.95, 1.05, 1.3, 1.55, 1.8, 1.15,
            ]
        )
        # fmt: on
        da = DataArray(data, dims="time", coords={"time": time})
        da["year"] = da.time.dt.year

        actual = da.groupby(
            year=UniqueGrouper(), time=SeasonGrouper(["NDJFM", "AMJ"])
        ).mean()

        # Expected if partial years are handled correctly
        expected = xr.DataArray(
            data=np.array([[1.38, 1.616667], [1.43333333, 1.5]]),
            dims=["year", "season"],
            coords={"year": [2001, 2002], "season": ["NDJFM", "AMJ"]},
        )

        assert_allclose(expected, actual)

    @pytest.mark.parametrize("calendar", ["standard"])
    def test_season_grouper_with_single_month_seasons(self, calendar):
        time = date_range("2001-01-01", "2002-12-30", freq="MS", calendar=calendar)
        # fmt: off
        data = np.array(
            [
                1.0, 1.25, 1.5, 1.75, 2.0, 1.1, 1.35, 1.6, 1.85, 1.2, 1.45, 1.7,
                1.95, 1.05, 1.3, 1.55, 1.8, 1.15, 1.4, 1.65, 1.9, 1.25, 1.5, 1.75,
            ]
        )
        # fmt: on
        da = DataArray(data, dims="time", coords={"time": time})
        da["year"] = da.time.dt.year

        # TODO: Consider supporting this if needed
        # It does not work without flox, because the group labels are not unique,
        # and so the stack/unstack approach does not work.
        with pytest.raises(ValueError):
            da.groupby(
                year=UniqueGrouper(),
                time=SeasonGrouper(
                    ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]
                ),
            ).mean()

        # Expected if single month seasons are handled correctly
        # expected = xr.DataArray(
        #     data=np.array(
        #         [
        #             [1.0, 1.25, 1.5, 1.75, 2.0, 1.1, 1.35, 1.6, 1.85, 1.2, 1.45, 1.7],
        #             [1.95, 1.05, 1.3, 1.55, 1.8, 1.15, 1.4, 1.65, 1.9, 1.25, 1.5, 1.75],
        #         ]
        #     ),
        #     dims=["year", "season"],
        #     coords={
        #         "year": [2001, 2002],
        #         "season": ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"],
        #     },
        # )
        # assert_allclose(expected, actual)

    @pytest.mark.parametrize("calendar", _ALL_CALENDARS)
    def test_season_grouper_with_months_spanning_calendar_year_using_previous_year(
        self, calendar
    ):
        time = date_range("2001-01-01", "2002-12-30", freq="MS", calendar=calendar)
        # fmt: off
        data = np.array(
            [
                1.0, 1.25, 1.5, 1.75, 2.0, 1.1, 1.35, 1.6, 1.85, 1.2, 1.45, 1.7,
                1.95, 1.05, 1.3, 1.55, 1.8, 1.15, 1.4, 1.65, 1.9, 1.25, 1.5, 1.75,
            ]
        )
        # fmt: on
        da = DataArray(data, dims="time", coords={"time": time})

        gb = da.resample(time=SeasonResampler(["NDJFM", "AMJ"], drop_incomplete=False))
        actual = gb.mean()

        # fmt: off
        new_time_da = xr.DataArray(
            dims="time",
            data=pd.DatetimeIndex(
                [
                    "2000-11-01", "2001-04-01", "2001-11-01", "2002-04-01", "2002-11-01"
                ]
            ),
        )
        # fmt: on
        if calendar != "standard":
            new_time_da = new_time_da.convert_calendar(
                calendar=calendar, align_on="date"
            )
        new_time = new_time_da.time.variable

        # Expected if the previous "ND" is used for seasonal grouping
        expected = xr.DataArray(
            data=np.array([1.25, 1.616667, 1.49, 1.5, 1.625]),
            dims="time",
            coords={"time": new_time},
        )
        assert_allclose(expected, actual)

    @pytest.mark.parametrize("calendar", _ALL_CALENDARS)
    def test_season_grouper_simple(self, calendar) -> None:
        time = date_range("2001-01-01", "2002-12-30", freq="D", calendar=calendar)
        da = DataArray(np.ones(time.size), dims="time", coords={"time": time})
        expected = da.groupby("time.season").mean()
        # note season order matches expected
        actual = da.groupby(
            time=SeasonGrouper(
                ["DJF", "JJA", "MAM", "SON"],  # drop_incomplete=False
            )
        ).mean()
        assert_identical(expected, actual)

    @pytest.mark.parametrize("seasons", [["JJA", "MAM", "SON", "DJF"]])
    def test_season_resampling_raises_unsorted_seasons(self, seasons):
        calendar = "standard"
        time = date_range("2001-01-01", "2002-12-30", freq="D", calendar=calendar)
        da = DataArray(np.ones(time.size), dims="time", coords={"time": time})
        with pytest.raises(ValueError, match="sort"):
            da.resample(time=SeasonResampler(seasons))

    @pytest.mark.parametrize(
        "use_cftime", [pytest.param(True, marks=requires_cftime), False]
    )
    @pytest.mark.parametrize("drop_incomplete", [True, False])
    @pytest.mark.parametrize(
        "seasons",
        [
            pytest.param(["DJF", "MAM", "JJA", "SON"], id="standard"),
            pytest.param(["NDJ", "FMA", "MJJ", "ASO"], id="nov-first"),
            pytest.param(["MAM", "JJA", "SON", "DJF"], id="standard-diff-order"),
            pytest.param(["JFM", "AMJ", "JAS", "OND"], id="december-same-year"),
            pytest.param(["DJF", "MAM", "JJA", "ON"], id="skip-september"),
            pytest.param(["JJAS"], id="jjas-only"),
        ],
    )
    def test_season_resampler(
        self, seasons: list[str], drop_incomplete: bool, use_cftime: bool
    ) -> None:
        calendar = "standard"
        time = date_range(
            "2001-01-01",
            "2002-12-30",
            freq="D",
            calendar=calendar,
            use_cftime=use_cftime,
        )
        da = DataArray(np.ones(time.size), dims="time", coords={"time": time})
        counts = da.resample(time="ME").count()

        seasons_as_ints = season_to_month_tuple(seasons)
        month = counts.time.dt.month.data
        year = counts.time.dt.year.data
        for season, as_ints in zip(seasons, seasons_as_ints, strict=True):
            if "DJ" in season:
                for imonth in as_ints[season.index("D") + 1 :]:
                    year[month == imonth] -= 1
        counts["time"] = (
            "time",
            [pd.Timestamp(f"{y}-{m}-01") for y, m in zip(year, month, strict=True)],
        )
        if has_cftime:
            counts = counts.convert_calendar(calendar, "time", align_on="date")

        expected_vals = []
        expected_time = []
        for year in [2001, 2002, 2003]:
            for season, as_ints in zip(seasons, seasons_as_ints, strict=True):
                out_year = year
                if "DJ" in season:
                    out_year = year - 1
                if out_year == 2003:
                    # this is a dummy year added to make sure we cover 2002-DJF
                    continue
                available = [
                    counts.sel(time=f"{out_year}-{month:02d}").data for month in as_ints
                ]
                if any(len(a) == 0 for a in available) and drop_incomplete:
                    continue
                output_label = pd.Timestamp(f"{out_year}-{as_ints[0]:02d}-01")
                expected_time.append(output_label)
                # use concatenate to handle empty array when dec value does not exist
                expected_vals.append(np.concatenate(available).sum())

        expected = (
            # we construct expected in the standard calendar
            xr.DataArray(expected_vals, dims="time", coords={"time": expected_time})
        )
        if has_cftime:
            # and then convert to the expected calendar,
            expected = expected.convert_calendar(
                calendar, align_on="date", use_cftime=use_cftime
            )
        # and finally sort since DJF will be out-of-order
        expected = expected.sortby("time")

        rs = SeasonResampler(seasons, drop_incomplete=drop_incomplete)
        # through resample
        actual = da.resample(time=rs).sum()
        assert_identical(actual, expected)

    @requires_cftime
    def test_season_resampler_errors(self):
        time = date_range("2001-01-01", "2002-12-30", freq="D", calendar="360_day")
        da = DataArray(np.ones(time.size), dims="time", coords={"time": time})

        # non-datetime array
        with pytest.raises(ValueError):
            DataArray(np.ones(5), dims="time").groupby(time=SeasonResampler(["DJF"]))

        # ndim > 1 array
        with pytest.raises(ValueError):
            DataArray(
                np.ones((5, 5)), dims=("t", "x"), coords={"x": np.arange(5)}
            ).groupby(x=SeasonResampler(["DJF"]))

        # overlapping seasons
        with pytest.raises(ValueError):
            da.groupby(time=SeasonResampler(["DJFM", "MAMJ", "JJAS", "SOND"])).sum()

    @requires_cftime
    def test_season_resampler_groupby_identical(self):
        time = date_range("2001-01-01", "2002-12-30", freq="D")
        da = DataArray(np.ones(time.size), dims="time", coords={"time": time})

        # through resample
        resampler = SeasonResampler(["DJF", "MAM", "JJA", "SON"])
        rs = da.resample(time=resampler).sum()

        # through groupby
        gb = da.groupby(time=resampler).sum()
        assert_identical(rs, gb)


# TODO: Possible property tests to add to this module
# 1. lambda x: x
# 2. grouped-reduce on unique coords is identical to array
# 3. group_over == groupby-reduce along other dimensions
# 4. result is equivalent for transposed input
