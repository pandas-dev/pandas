from __future__ import annotations

from copy import deepcopy
from typing import TYPE_CHECKING, Any, Callable

import numpy as np
import pandas as pd
import pytest

from xarray import DataArray, Dataset, Variable, concat
from xarray.core import dtypes, merge
from xarray.core.coordinates import Coordinates
from xarray.core.indexes import PandasIndex
from xarray.tests import (
    ConcatenatableArray,
    InaccessibleArray,
    UnexpectedDataAccess,
    assert_array_equal,
    assert_equal,
    assert_identical,
    requires_dask,
)
from xarray.tests.test_dataset import create_test_data

if TYPE_CHECKING:
    from xarray.core.types import CombineAttrsOptions, JoinOptions


# helper method to create multiple tests datasets to concat
def create_concat_datasets(
    num_datasets: int = 2, seed: int | None = None, include_day: bool = True
) -> list[Dataset]:
    rng = np.random.default_rng(seed)
    lat = rng.standard_normal(size=(1, 4))
    lon = rng.standard_normal(size=(1, 4))
    result = []
    variables = ["temperature", "pressure", "humidity", "precipitation", "cloud_cover"]
    for i in range(num_datasets):
        if include_day:
            data_tuple = (
                ["x", "y", "day"],
                rng.standard_normal(size=(1, 4, 2)),
            )
            data_vars = {v: data_tuple for v in variables}
            result.append(
                Dataset(
                    data_vars=data_vars,
                    coords={
                        "lat": (["x", "y"], lat),
                        "lon": (["x", "y"], lon),
                        "day": ["day" + str(i * 2 + 1), "day" + str(i * 2 + 2)],
                    },
                )
            )
        else:
            data_tuple = (
                ["x", "y"],
                rng.standard_normal(size=(1, 4)),
            )
            data_vars = {v: data_tuple for v in variables}
            result.append(
                Dataset(
                    data_vars=data_vars,
                    coords={"lat": (["x", "y"], lat), "lon": (["x", "y"], lon)},
                )
            )

    return result


# helper method to create multiple tests datasets to concat with specific types
def create_typed_datasets(
    num_datasets: int = 2, seed: int | None = None
) -> list[Dataset]:
    var_strings = ["a", "b", "c", "d", "e", "f", "g", "h"]
    result = []
    rng = np.random.default_rng(seed)
    lat = rng.standard_normal(size=(1, 4))
    lon = rng.standard_normal(size=(1, 4))
    for i in range(num_datasets):
        result.append(
            Dataset(
                data_vars={
                    "float": (["x", "y", "day"], rng.standard_normal(size=(1, 4, 2))),
                    "float2": (["x", "y", "day"], rng.standard_normal(size=(1, 4, 2))),
                    "string": (
                        ["x", "y", "day"],
                        rng.choice(var_strings, size=(1, 4, 2)),
                    ),
                    "int": (["x", "y", "day"], rng.integers(0, 10, size=(1, 4, 2))),
                    "datetime64": (
                        ["x", "y", "day"],
                        np.arange(
                            np.datetime64("2017-01-01"), np.datetime64("2017-01-09")
                        ).reshape(1, 4, 2),
                    ),
                    "timedelta64": (
                        ["x", "y", "day"],
                        np.reshape([pd.Timedelta(days=i) for i in range(8)], [1, 4, 2]),
                    ),
                },
                coords={
                    "lat": (["x", "y"], lat),
                    "lon": (["x", "y"], lon),
                    "day": ["day" + str(i * 2 + 1), "day" + str(i * 2 + 2)],
                },
            )
        )
    return result


def test_concat_compat() -> None:
    ds1 = Dataset(
        {
            "has_x_y": (("y", "x"), [[1, 2]]),
            "has_x": ("x", [1, 2]),
            "no_x_y": ("z", [1, 2]),
        },
        coords={"x": [0, 1], "y": [0], "z": [-1, -2]},
    )
    ds2 = Dataset(
        {
            "has_x_y": (("y", "x"), [[3, 4]]),
            "has_x": ("x", [1, 2]),
            "no_x_y": (("q", "z"), [[1, 2]]),
        },
        coords={"x": [0, 1], "y": [1], "z": [-1, -2], "q": [0]},
    )

    result = concat([ds1, ds2], dim="y", data_vars="minimal", compat="broadcast_equals")
    assert_equal(ds2.no_x_y, result.no_x_y.transpose())

    for var in ["has_x", "no_x_y"]:
        assert "y" not in result[var].dims and "y" not in result[var].coords
    with pytest.raises(ValueError, match=r"'q' not present in all datasets"):
        concat([ds1, ds2], dim="q")
    with pytest.raises(ValueError, match=r"'q' not present in all datasets"):
        concat([ds2, ds1], dim="q")


def test_concat_missing_var() -> None:
    datasets = create_concat_datasets(2, seed=123)
    expected = concat(datasets, dim="day")
    vars_to_drop = ["humidity", "precipitation", "cloud_cover"]

    expected = expected.drop_vars(vars_to_drop)
    expected["pressure"][..., 2:] = np.nan

    datasets[0] = datasets[0].drop_vars(vars_to_drop)
    datasets[1] = datasets[1].drop_vars(vars_to_drop + ["pressure"])
    actual = concat(datasets, dim="day")

    assert list(actual.data_vars.keys()) == ["temperature", "pressure"]
    assert_identical(actual, expected)


def test_concat_categorical() -> None:
    data1 = create_test_data(use_extension_array=True)
    data2 = create_test_data(use_extension_array=True)
    concatenated = concat([data1, data2], dim="dim1")
    assert (
        concatenated["var4"]
        == type(data2["var4"].variable.data.array)._concat_same_type(
            [
                data1["var4"].variable.data.array,
                data2["var4"].variable.data.array,
            ]
        )
    ).all()


def test_concat_missing_multiple_consecutive_var() -> None:
    datasets = create_concat_datasets(3, seed=123)
    expected = concat(datasets, dim="day")
    vars_to_drop = ["humidity", "pressure"]

    expected["pressure"][..., :4] = np.nan
    expected["humidity"][..., :4] = np.nan

    datasets[0] = datasets[0].drop_vars(vars_to_drop)
    datasets[1] = datasets[1].drop_vars(vars_to_drop)
    actual = concat(datasets, dim="day")

    assert list(actual.data_vars.keys()) == [
        "temperature",
        "precipitation",
        "cloud_cover",
        "pressure",
        "humidity",
    ]
    assert_identical(actual, expected)


def test_concat_all_empty() -> None:
    ds1 = Dataset()
    ds2 = Dataset()
    expected = Dataset()
    actual = concat([ds1, ds2], dim="new_dim")

    assert_identical(actual, expected)


def test_concat_second_empty() -> None:
    ds1 = Dataset(data_vars={"a": ("y", [0.1])}, coords={"x": 0.1})
    ds2 = Dataset(coords={"x": 0.1})

    expected = Dataset(data_vars={"a": ("y", [0.1, np.nan])}, coords={"x": 0.1})
    actual = concat([ds1, ds2], dim="y")
    assert_identical(actual, expected)

    expected = Dataset(
        data_vars={"a": ("y", [0.1, np.nan])}, coords={"x": ("y", [0.1, 0.1])}
    )
    actual = concat([ds1, ds2], dim="y", coords="all")
    assert_identical(actual, expected)

    # Check concatenating scalar data_var only present in ds1
    ds1["b"] = 0.1
    expected = Dataset(
        data_vars={"a": ("y", [0.1, np.nan]), "b": ("y", [0.1, np.nan])},
        coords={"x": ("y", [0.1, 0.1])},
    )
    actual = concat([ds1, ds2], dim="y", coords="all", data_vars="all")
    assert_identical(actual, expected)

    expected = Dataset(
        data_vars={"a": ("y", [0.1, np.nan]), "b": 0.1}, coords={"x": 0.1}
    )
    actual = concat([ds1, ds2], dim="y", coords="different", data_vars="different")
    assert_identical(actual, expected)


def test_concat_multiple_missing_variables() -> None:
    datasets = create_concat_datasets(2, seed=123)
    expected = concat(datasets, dim="day")
    vars_to_drop = ["pressure", "cloud_cover"]

    expected["pressure"][..., 2:] = np.nan
    expected["cloud_cover"][..., 2:] = np.nan

    datasets[1] = datasets[1].drop_vars(vars_to_drop)
    actual = concat(datasets, dim="day")

    # check the variables orders are the same
    assert list(actual.data_vars.keys()) == [
        "temperature",
        "pressure",
        "humidity",
        "precipitation",
        "cloud_cover",
    ]

    assert_identical(actual, expected)


@pytest.mark.parametrize("include_day", [True, False])
def test_concat_multiple_datasets_missing_vars(include_day: bool) -> None:
    vars_to_drop = [
        "temperature",
        "pressure",
        "humidity",
        "precipitation",
        "cloud_cover",
    ]

    datasets = create_concat_datasets(
        len(vars_to_drop), seed=123, include_day=include_day
    )
    expected = concat(datasets, dim="day")

    for i, name in enumerate(vars_to_drop):
        if include_day:
            expected[name][..., i * 2 : (i + 1) * 2] = np.nan
        else:
            expected[name][i : i + 1, ...] = np.nan

    # set up the test data
    datasets = [ds.drop_vars(varname) for ds, varname in zip(datasets, vars_to_drop)]

    actual = concat(datasets, dim="day")

    assert list(actual.data_vars.keys()) == [
        "pressure",
        "humidity",
        "precipitation",
        "cloud_cover",
        "temperature",
    ]
    assert_identical(actual, expected)


def test_concat_multiple_datasets_with_multiple_missing_variables() -> None:
    vars_to_drop_in_first = ["temperature", "pressure"]
    vars_to_drop_in_second = ["humidity", "precipitation", "cloud_cover"]
    datasets = create_concat_datasets(2, seed=123)
    expected = concat(datasets, dim="day")
    for name in vars_to_drop_in_first:
        expected[name][..., :2] = np.nan
    for name in vars_to_drop_in_second:
        expected[name][..., 2:] = np.nan

    # set up the test data
    datasets[0] = datasets[0].drop_vars(vars_to_drop_in_first)
    datasets[1] = datasets[1].drop_vars(vars_to_drop_in_second)

    actual = concat(datasets, dim="day")

    assert list(actual.data_vars.keys()) == [
        "humidity",
        "precipitation",
        "cloud_cover",
        "temperature",
        "pressure",
    ]
    assert_identical(actual, expected)


@pytest.mark.filterwarnings("ignore:Converting non-nanosecond")
def test_concat_type_of_missing_fill() -> None:
    datasets = create_typed_datasets(2, seed=123)
    expected1 = concat(datasets, dim="day", fill_value=dtypes.NA)
    expected2 = concat(datasets[::-1], dim="day", fill_value=dtypes.NA)
    vars = ["float", "float2", "string", "int", "datetime64", "timedelta64"]
    expected = [expected2, expected1]
    for i, exp in enumerate(expected):
        sl = slice(i * 2, (i + 1) * 2)
        exp["float2"][..., sl] = np.nan
        exp["datetime64"][..., sl] = np.nan
        exp["timedelta64"][..., sl] = np.nan
        var = exp["int"] * 1.0
        var[..., sl] = np.nan
        exp["int"] = var
        var = exp["string"].astype(object)
        var[..., sl] = np.nan
        exp["string"] = var

    # set up the test data
    datasets[1] = datasets[1].drop_vars(vars[1:])

    actual = concat(datasets, dim="day", fill_value=dtypes.NA)

    assert_identical(actual, expected[1])

    # reversed
    actual = concat(datasets[::-1], dim="day", fill_value=dtypes.NA)

    assert_identical(actual, expected[0])


def test_concat_order_when_filling_missing() -> None:
    vars_to_drop_in_first: list[str] = []
    # drop middle
    vars_to_drop_in_second = ["humidity"]
    datasets = create_concat_datasets(2, seed=123)
    expected1 = concat(datasets, dim="day")
    for name in vars_to_drop_in_second:
        expected1[name][..., 2:] = np.nan
    expected2 = concat(datasets[::-1], dim="day")
    for name in vars_to_drop_in_second:
        expected2[name][..., :2] = np.nan

    # set up the test data
    datasets[0] = datasets[0].drop_vars(vars_to_drop_in_first)
    datasets[1] = datasets[1].drop_vars(vars_to_drop_in_second)

    actual = concat(datasets, dim="day")

    assert list(actual.data_vars.keys()) == [
        "temperature",
        "pressure",
        "humidity",
        "precipitation",
        "cloud_cover",
    ]
    assert_identical(actual, expected1)

    actual = concat(datasets[::-1], dim="day")

    assert list(actual.data_vars.keys()) == [
        "temperature",
        "pressure",
        "precipitation",
        "cloud_cover",
        "humidity",
    ]
    assert_identical(actual, expected2)


@pytest.fixture
def concat_var_names() -> Callable:
    # create var names list with one missing value
    def get_varnames(var_cnt: int = 10, list_cnt: int = 10) -> list[list[str]]:
        orig = [f"d{i:02d}" for i in range(var_cnt)]
        var_names = []
        for i in range(0, list_cnt):
            l1 = orig.copy()
            var_names.append(l1)
        return var_names

    return get_varnames


@pytest.fixture
def create_concat_ds() -> Callable:
    def create_ds(
        var_names: list[list[str]],
        dim: bool = False,
        coord: bool = False,
        drop_idx: list[int] | None = None,
    ) -> list[Dataset]:
        out_ds = []
        ds = Dataset()
        ds = ds.assign_coords({"x": np.arange(2)})
        ds = ds.assign_coords({"y": np.arange(3)})
        ds = ds.assign_coords({"z": np.arange(4)})
        for i, dsl in enumerate(var_names):
            vlist = dsl.copy()
            if drop_idx is not None:
                vlist.pop(drop_idx[i])
            foo_data = np.arange(48, dtype=float).reshape(2, 2, 3, 4)
            dsi = ds.copy()
            if coord:
                dsi = ds.assign({"time": (["time"], [i * 2, i * 2 + 1])})
            for k in vlist:
                dsi = dsi.assign({k: (["time", "x", "y", "z"], foo_data.copy())})
            if not dim:
                dsi = dsi.isel(time=0)
            out_ds.append(dsi)
        return out_ds

    return create_ds


@pytest.mark.parametrize("dim", [True, False])
@pytest.mark.parametrize("coord", [True, False])
def test_concat_fill_missing_variables(
    concat_var_names, create_concat_ds, dim: bool, coord: bool
) -> None:
    var_names = concat_var_names()
    drop_idx = [0, 7, 6, 4, 4, 8, 0, 6, 2, 0]

    expected = concat(
        create_concat_ds(var_names, dim=dim, coord=coord), dim="time", data_vars="all"
    )
    for i, idx in enumerate(drop_idx):
        if dim:
            expected[var_names[0][idx]][i * 2 : i * 2 + 2] = np.nan
        else:
            expected[var_names[0][idx]][i] = np.nan

    concat_ds = create_concat_ds(var_names, dim=dim, coord=coord, drop_idx=drop_idx)
    actual = concat(concat_ds, dim="time", data_vars="all")

    assert list(actual.data_vars.keys()) == [
        "d01",
        "d02",
        "d03",
        "d04",
        "d05",
        "d06",
        "d07",
        "d08",
        "d09",
        "d00",
    ]
    assert_identical(actual, expected)


class TestConcatDataset:
    @pytest.fixture
    def data(self, request) -> Dataset:
        use_extension_array = request.param if hasattr(request, "param") else False
        return create_test_data(use_extension_array=use_extension_array).drop_dims(
            "dim3"
        )

    def rectify_dim_order(self, data, dataset) -> Dataset:
        # return a new dataset with all variable dimensions transposed into
        # the order in which they are found in `data`
        return Dataset(
            {k: v.transpose(*data[k].dims) for k, v in dataset.data_vars.items()},
            dataset.coords,
            attrs=dataset.attrs,
        )

    @pytest.mark.parametrize("coords", ["different", "minimal"])
    @pytest.mark.parametrize(
        "dim,data", [["dim1", True], ["dim2", False]], indirect=["data"]
    )
    def test_concat_simple(self, data, dim, coords) -> None:
        datasets = [g for _, g in data.groupby(dim, squeeze=False)]
        assert_identical(data, concat(datasets, dim, coords=coords))

    def test_concat_merge_variables_present_in_some_datasets(self, data) -> None:
        # coordinates present in some datasets but not others
        ds1 = Dataset(data_vars={"a": ("y", [0.1])}, coords={"x": 0.1})
        ds2 = Dataset(data_vars={"a": ("y", [0.2])}, coords={"z": 0.2})
        actual = concat([ds1, ds2], dim="y", coords="minimal")
        expected = Dataset({"a": ("y", [0.1, 0.2])}, coords={"x": 0.1, "z": 0.2})
        assert_identical(expected, actual)

        # data variables present in some datasets but not others
        split_data = [data.isel(dim1=slice(3)), data.isel(dim1=slice(3, None))]
        data0, data1 = deepcopy(split_data)
        data1["foo"] = ("bar", np.random.randn(10))
        actual = concat([data0, data1], "dim1", data_vars="minimal")
        expected = data.copy().assign(foo=data1.foo)
        assert_identical(expected, actual)

        # expand foo
        actual = concat([data0, data1], "dim1")
        foo = np.ones((8, 10), dtype=data1.foo.dtype) * np.nan
        foo[3:] = data1.foo.values[None, ...]
        expected = data.copy().assign(foo=(["dim1", "bar"], foo))
        assert_identical(expected, actual)

    @pytest.mark.parametrize("data", [False], indirect=["data"])
    def test_concat_2(self, data) -> None:
        dim = "dim2"
        datasets = [g.squeeze(dim) for _, g in data.groupby(dim, squeeze=False)]
        concat_over = [k for k, v in data.coords.items() if dim in v.dims and k != dim]
        actual = concat(datasets, data[dim], coords=concat_over)
        assert_identical(data, self.rectify_dim_order(data, actual))

    @pytest.mark.parametrize("coords", ["different", "minimal", "all"])
    @pytest.mark.parametrize("dim", ["dim1", "dim2"])
    def test_concat_coords_kwarg(self, data, dim, coords) -> None:
        data = data.copy(deep=True)
        # make sure the coords argument behaves as expected
        data.coords["extra"] = ("dim4", np.arange(3))
        datasets = [g.squeeze() for _, g in data.groupby(dim, squeeze=False)]

        actual = concat(datasets, data[dim], coords=coords)
        if coords == "all":
            expected = np.array([data["extra"].values for _ in range(data.sizes[dim])])
            assert_array_equal(actual["extra"].values, expected)

        else:
            assert_equal(data["extra"], actual["extra"])

    def test_concat(self, data) -> None:
        split_data = [
            data.isel(dim1=slice(3)),
            data.isel(dim1=3),
            data.isel(dim1=slice(4, None)),
        ]
        assert_identical(data, concat(split_data, "dim1"))

    def test_concat_dim_precedence(self, data) -> None:
        # verify that the dim argument takes precedence over
        # concatenating dataset variables of the same name
        dim = (2 * data["dim1"]).rename("dim1")
        datasets = [g for _, g in data.groupby("dim1", squeeze=False)]
        expected = data.copy()
        expected["dim1"] = dim
        assert_identical(expected, concat(datasets, dim))

    def test_concat_data_vars_typing(self) -> None:
        # Testing typing, can be removed if the next function works with annotations.
        data = Dataset({"foo": ("x", np.random.randn(10))})
        objs: list[Dataset] = [data.isel(x=slice(5)), data.isel(x=slice(5, None))]
        actual = concat(objs, dim="x", data_vars="minimal")
        assert_identical(data, actual)

    def test_concat_data_vars(self) -> None:
        data = Dataset({"foo": ("x", np.random.randn(10))})
        objs: list[Dataset] = [data.isel(x=slice(5)), data.isel(x=slice(5, None))]
        for data_vars in ["minimal", "different", "all", [], ["foo"]]:
            actual = concat(objs, dim="x", data_vars=data_vars)
            assert_identical(data, actual)

    def test_concat_coords(self):
        # TODO: annotating this func fails
        data = Dataset({"foo": ("x", np.random.randn(10))})
        expected = data.assign_coords(c=("x", [0] * 5 + [1] * 5))
        objs = [
            data.isel(x=slice(5)).assign_coords(c=0),
            data.isel(x=slice(5, None)).assign_coords(c=1),
        ]
        for coords in ["different", "all", ["c"]]:
            actual = concat(objs, dim="x", coords=coords)
            assert_identical(expected, actual)
        for coords in ["minimal", []]:
            with pytest.raises(merge.MergeError, match="conflicting values"):
                concat(objs, dim="x", coords=coords)

    def test_concat_constant_index(self):
        # TODO: annotating this func fails
        # GH425
        ds1 = Dataset({"foo": 1.5}, {"y": 1})
        ds2 = Dataset({"foo": 2.5}, {"y": 1})
        expected = Dataset({"foo": ("y", [1.5, 2.5]), "y": [1, 1]})
        for mode in ["different", "all", ["foo"]]:
            actual = concat([ds1, ds2], "y", data_vars=mode)
            assert_identical(expected, actual)
        with pytest.raises(merge.MergeError, match="conflicting values"):
            # previously dim="y", and raised error which makes no sense.
            # "foo" has dimension "y" so minimal should concatenate it?
            concat([ds1, ds2], "new_dim", data_vars="minimal")

    def test_concat_size0(self) -> None:
        data = create_test_data()
        split_data = [data.isel(dim1=slice(0, 0)), data]
        actual = concat(split_data, "dim1")
        assert_identical(data, actual)

        actual = concat(split_data[::-1], "dim1")
        assert_identical(data, actual)

    def test_concat_autoalign(self) -> None:
        ds1 = Dataset({"foo": DataArray([1, 2], coords=[("x", [1, 2])])})
        ds2 = Dataset({"foo": DataArray([1, 2], coords=[("x", [1, 3])])})
        actual = concat([ds1, ds2], "y")
        expected = Dataset(
            {
                "foo": DataArray(
                    [[1, 2, np.nan], [1, np.nan, 2]],
                    dims=["y", "x"],
                    coords={"x": [1, 2, 3]},
                )
            }
        )
        assert_identical(expected, actual)

    def test_concat_errors(self):
        # TODO: annotating this func fails
        data = create_test_data()
        split_data = [data.isel(dim1=slice(3)), data.isel(dim1=slice(3, None))]

        with pytest.raises(ValueError, match=r"must supply at least one"):
            concat([], "dim1")

        with pytest.raises(ValueError, match=r"Cannot specify both .*='different'"):
            concat(
                [data, data], dim="concat_dim", data_vars="different", compat="override"
            )

        with pytest.raises(ValueError, match=r"must supply at least one"):
            concat([], "dim1")

        with pytest.raises(ValueError, match=r"are not found in the coordinates"):
            concat([data, data], "new_dim", coords=["not_found"])

        with pytest.raises(ValueError, match=r"are not found in the data variables"):
            concat([data, data], "new_dim", data_vars=["not_found"])

        with pytest.raises(ValueError, match=r"global attributes not"):
            # call deepcopy separately to get unique attrs
            data0 = deepcopy(split_data[0])
            data1 = deepcopy(split_data[1])
            data1.attrs["foo"] = "bar"
            concat([data0, data1], "dim1", compat="identical")
        assert_identical(data, concat([data0, data1], "dim1", compat="equals"))

        with pytest.raises(ValueError, match=r"compat.* invalid"):
            concat(split_data, "dim1", compat="foobar")

        with pytest.raises(ValueError, match=r"unexpected value for"):
            concat([data, data], "new_dim", coords="foobar")

        with pytest.raises(
            ValueError, match=r"coordinate in some datasets but not others"
        ):
            concat([Dataset({"x": 0}), Dataset({"x": [1]})], dim="z")

        with pytest.raises(
            ValueError, match=r"coordinate in some datasets but not others"
        ):
            concat([Dataset({"x": 0}), Dataset({}, {"x": 1})], dim="z")

    def test_concat_join_kwarg(self) -> None:
        ds1 = Dataset({"a": (("x", "y"), [[0]])}, coords={"x": [0], "y": [0]})
        ds2 = Dataset({"a": (("x", "y"), [[0]])}, coords={"x": [1], "y": [0.0001]})

        expected: dict[JoinOptions, Any] = {}
        expected["outer"] = Dataset(
            {"a": (("x", "y"), [[0, np.nan], [np.nan, 0]])},
            {"x": [0, 1], "y": [0, 0.0001]},
        )
        expected["inner"] = Dataset(
            {"a": (("x", "y"), [[], []])}, {"x": [0, 1], "y": []}
        )
        expected["left"] = Dataset(
            {"a": (("x", "y"), np.array([0, np.nan], ndmin=2).T)},
            coords={"x": [0, 1], "y": [0]},
        )
        expected["right"] = Dataset(
            {"a": (("x", "y"), np.array([np.nan, 0], ndmin=2).T)},
            coords={"x": [0, 1], "y": [0.0001]},
        )
        expected["override"] = Dataset(
            {"a": (("x", "y"), np.array([0, 0], ndmin=2).T)},
            coords={"x": [0, 1], "y": [0]},
        )

        with pytest.raises(ValueError, match=r"cannot align.*exact.*dimensions.*'y'"):
            actual = concat([ds1, ds2], join="exact", dim="x")

        for join in expected:
            actual = concat([ds1, ds2], join=join, dim="x")
            assert_equal(actual, expected[join])

        # regression test for #3681
        actual = concat(
            [ds1.drop_vars("x"), ds2.drop_vars("x")], join="override", dim="y"
        )
        expected2 = Dataset(
            {"a": (("x", "y"), np.array([0, 0], ndmin=2))}, coords={"y": [0, 0.0001]}
        )
        assert_identical(actual, expected2)

    @pytest.mark.parametrize(
        "combine_attrs, var1_attrs, var2_attrs, expected_attrs, expect_exception",
        [
            (
                "no_conflicts",
                {"a": 1, "b": 2},
                {"a": 1, "c": 3},
                {"a": 1, "b": 2, "c": 3},
                False,
            ),
            ("no_conflicts", {"a": 1, "b": 2}, {}, {"a": 1, "b": 2}, False),
            ("no_conflicts", {}, {"a": 1, "c": 3}, {"a": 1, "c": 3}, False),
            (
                "no_conflicts",
                {"a": 1, "b": 2},
                {"a": 4, "c": 3},
                {"a": 1, "b": 2, "c": 3},
                True,
            ),
            ("drop", {"a": 1, "b": 2}, {"a": 1, "c": 3}, {}, False),
            ("identical", {"a": 1, "b": 2}, {"a": 1, "b": 2}, {"a": 1, "b": 2}, False),
            ("identical", {"a": 1, "b": 2}, {"a": 1, "c": 3}, {"a": 1, "b": 2}, True),
            (
                "override",
                {"a": 1, "b": 2},
                {"a": 4, "b": 5, "c": 3},
                {"a": 1, "b": 2},
                False,
            ),
            (
                "drop_conflicts",
                {"a": 41, "b": 42, "c": 43},
                {"b": 2, "c": 43, "d": 44},
                {"a": 41, "c": 43, "d": 44},
                False,
            ),
            (
                lambda attrs, context: {"a": -1, "b": 0, "c": 1} if any(attrs) else {},
                {"a": 41, "b": 42, "c": 43},
                {"b": 2, "c": 43, "d": 44},
                {"a": -1, "b": 0, "c": 1},
                False,
            ),
        ],
    )
    def test_concat_combine_attrs_kwarg(
        self, combine_attrs, var1_attrs, var2_attrs, expected_attrs, expect_exception
    ):
        ds1 = Dataset({"a": ("x", [0])}, coords={"x": [0]}, attrs=var1_attrs)
        ds2 = Dataset({"a": ("x", [0])}, coords={"x": [1]}, attrs=var2_attrs)

        if expect_exception:
            with pytest.raises(ValueError, match=f"combine_attrs='{combine_attrs}'"):
                concat([ds1, ds2], dim="x", combine_attrs=combine_attrs)
        else:
            actual = concat([ds1, ds2], dim="x", combine_attrs=combine_attrs)
            expected = Dataset(
                {"a": ("x", [0, 0])}, {"x": [0, 1]}, attrs=expected_attrs
            )

            assert_identical(actual, expected)

    @pytest.mark.parametrize(
        "combine_attrs, attrs1, attrs2, expected_attrs, expect_exception",
        [
            (
                "no_conflicts",
                {"a": 1, "b": 2},
                {"a": 1, "c": 3},
                {"a": 1, "b": 2, "c": 3},
                False,
            ),
            ("no_conflicts", {"a": 1, "b": 2}, {}, {"a": 1, "b": 2}, False),
            ("no_conflicts", {}, {"a": 1, "c": 3}, {"a": 1, "c": 3}, False),
            (
                "no_conflicts",
                {"a": 1, "b": 2},
                {"a": 4, "c": 3},
                {"a": 1, "b": 2, "c": 3},
                True,
            ),
            ("drop", {"a": 1, "b": 2}, {"a": 1, "c": 3}, {}, False),
            ("identical", {"a": 1, "b": 2}, {"a": 1, "b": 2}, {"a": 1, "b": 2}, False),
            ("identical", {"a": 1, "b": 2}, {"a": 1, "c": 3}, {"a": 1, "b": 2}, True),
            (
                "override",
                {"a": 1, "b": 2},
                {"a": 4, "b": 5, "c": 3},
                {"a": 1, "b": 2},
                False,
            ),
            (
                "drop_conflicts",
                {"a": 41, "b": 42, "c": 43},
                {"b": 2, "c": 43, "d": 44},
                {"a": 41, "c": 43, "d": 44},
                False,
            ),
            (
                lambda attrs, context: {"a": -1, "b": 0, "c": 1} if any(attrs) else {},
                {"a": 41, "b": 42, "c": 43},
                {"b": 2, "c": 43, "d": 44},
                {"a": -1, "b": 0, "c": 1},
                False,
            ),
        ],
    )
    def test_concat_combine_attrs_kwarg_variables(
        self, combine_attrs, attrs1, attrs2, expected_attrs, expect_exception
    ):
        """check that combine_attrs is used on data variables and coords"""
        ds1 = Dataset({"a": ("x", [0], attrs1)}, coords={"x": ("x", [0], attrs1)})
        ds2 = Dataset({"a": ("x", [0], attrs2)}, coords={"x": ("x", [1], attrs2)})

        if expect_exception:
            with pytest.raises(ValueError, match=f"combine_attrs='{combine_attrs}'"):
                concat([ds1, ds2], dim="x", combine_attrs=combine_attrs)
        else:
            actual = concat([ds1, ds2], dim="x", combine_attrs=combine_attrs)
            expected = Dataset(
                {"a": ("x", [0, 0], expected_attrs)},
                {"x": ("x", [0, 1], expected_attrs)},
            )

            assert_identical(actual, expected)

    def test_concat_promote_shape(self) -> None:
        # mixed dims within variables
        objs = [Dataset({}, {"x": 0}), Dataset({"x": [1]})]
        actual = concat(objs, "x")
        expected = Dataset({"x": [0, 1]})
        assert_identical(actual, expected)

        objs = [Dataset({"x": [0]}), Dataset({}, {"x": 1})]
        actual = concat(objs, "x")
        assert_identical(actual, expected)

        # mixed dims between variables
        objs = [Dataset({"x": [2], "y": 3}), Dataset({"x": [4], "y": 5})]
        actual = concat(objs, "x")
        expected = Dataset({"x": [2, 4], "y": ("x", [3, 5])})
        assert_identical(actual, expected)

        # mixed dims in coord variable
        objs = [Dataset({"x": [0]}, {"y": -1}), Dataset({"x": [1]}, {"y": ("x", [-2])})]
        actual = concat(objs, "x")
        expected = Dataset({"x": [0, 1]}, {"y": ("x", [-1, -2])})
        assert_identical(actual, expected)

        # scalars with mixed lengths along concat dim -- values should repeat
        objs = [Dataset({"x": [0]}, {"y": -1}), Dataset({"x": [1, 2]}, {"y": -2})]
        actual = concat(objs, "x")
        expected = Dataset({"x": [0, 1, 2]}, {"y": ("x", [-1, -2, -2])})
        assert_identical(actual, expected)

        # broadcast 1d x 1d -> 2d
        objs = [
            Dataset({"z": ("x", [-1])}, {"x": [0], "y": [0]}),
            Dataset({"z": ("y", [1])}, {"x": [1], "y": [0]}),
        ]
        actual = concat(objs, "x")
        expected = Dataset({"z": (("x", "y"), [[-1], [1]])}, {"x": [0, 1], "y": [0]})
        assert_identical(actual, expected)

        # regression GH6384
        objs = [
            Dataset({}, {"x": pd.Interval(-1, 0, closed="right")}),
            Dataset({"x": [pd.Interval(0, 1, closed="right")]}),
        ]
        actual = concat(objs, "x")
        expected = Dataset(
            {
                "x": [
                    pd.Interval(-1, 0, closed="right"),
                    pd.Interval(0, 1, closed="right"),
                ]
            }
        )
        assert_identical(actual, expected)

        # regression GH6416 (coord dtype) and GH6434
        time_data1 = np.array(["2022-01-01", "2022-02-01"], dtype="datetime64[ns]")
        time_data2 = np.array("2022-03-01", dtype="datetime64[ns]")
        time_expected = np.array(
            ["2022-01-01", "2022-02-01", "2022-03-01"], dtype="datetime64[ns]"
        )
        objs = [Dataset({}, {"time": time_data1}), Dataset({}, {"time": time_data2})]
        actual = concat(objs, "time")
        expected = Dataset({}, {"time": time_expected})
        assert_identical(actual, expected)
        assert isinstance(actual.indexes["time"], pd.DatetimeIndex)

    def test_concat_do_not_promote(self) -> None:
        # GH438
        objs = [
            Dataset({"y": ("t", [1])}, {"x": 1, "t": [0]}),
            Dataset({"y": ("t", [2])}, {"x": 1, "t": [0]}),
        ]
        expected = Dataset({"y": ("t", [1, 2])}, {"x": 1, "t": [0, 0]})
        actual = concat(objs, "t")
        assert_identical(expected, actual)

        objs = [
            Dataset({"y": ("t", [1])}, {"x": 1, "t": [0]}),
            Dataset({"y": ("t", [2])}, {"x": 2, "t": [0]}),
        ]
        with pytest.raises(ValueError):
            concat(objs, "t", coords="minimal")

    def test_concat_dim_is_variable(self) -> None:
        objs = [Dataset({"x": 0}), Dataset({"x": 1})]
        coord = Variable("y", [3, 4], attrs={"foo": "bar"})
        expected = Dataset({"x": ("y", [0, 1]), "y": coord})
        actual = concat(objs, coord)
        assert_identical(actual, expected)

    def test_concat_dim_is_dataarray(self) -> None:
        objs = [Dataset({"x": 0}), Dataset({"x": 1})]
        coord = DataArray([3, 4], dims="y", attrs={"foo": "bar"})
        expected = Dataset({"x": ("y", [0, 1]), "y": coord})
        actual = concat(objs, coord)
        assert_identical(actual, expected)

    def test_concat_multiindex(self) -> None:
        midx = pd.MultiIndex.from_product([[1, 2, 3], ["a", "b"]])
        midx_coords = Coordinates.from_pandas_multiindex(midx, "x")
        expected = Dataset(coords=midx_coords)
        actual = concat(
            [expected.isel(x=slice(2)), expected.isel(x=slice(2, None))], "x"
        )
        assert expected.equals(actual)
        assert isinstance(actual.x.to_index(), pd.MultiIndex)

    def test_concat_along_new_dim_multiindex(self) -> None:
        # see https://github.com/pydata/xarray/issues/6881
        level_names = ["x_level_0", "x_level_1"]
        midx = pd.MultiIndex.from_product([[1, 2, 3], ["a", "b"]], names=level_names)
        midx_coords = Coordinates.from_pandas_multiindex(midx, "x")
        ds = Dataset(coords=midx_coords)
        concatenated = concat([ds], "new")
        actual = list(concatenated.xindexes.get_all_coords("x"))
        expected = ["x"] + level_names
        assert actual == expected

    @pytest.mark.parametrize("fill_value", [dtypes.NA, 2, 2.0, {"a": 2, "b": 1}])
    def test_concat_fill_value(self, fill_value) -> None:
        datasets = [
            Dataset({"a": ("x", [2, 3]), "b": ("x", [-2, 1]), "x": [1, 2]}),
            Dataset({"a": ("x", [1, 2]), "b": ("x", [3, -1]), "x": [0, 1]}),
        ]
        if fill_value == dtypes.NA:
            # if we supply the default, we expect the missing value for a
            # float array
            fill_value_a = fill_value_b = np.nan
        elif isinstance(fill_value, dict):
            fill_value_a = fill_value["a"]
            fill_value_b = fill_value["b"]
        else:
            fill_value_a = fill_value_b = fill_value
        expected = Dataset(
            {
                "a": (("t", "x"), [[fill_value_a, 2, 3], [1, 2, fill_value_a]]),
                "b": (("t", "x"), [[fill_value_b, -2, 1], [3, -1, fill_value_b]]),
            },
            {"x": [0, 1, 2]},
        )
        actual = concat(datasets, dim="t", fill_value=fill_value)
        assert_identical(actual, expected)

    @pytest.mark.parametrize("dtype", [str, bytes])
    @pytest.mark.parametrize("dim", ["x1", "x2"])
    def test_concat_str_dtype(self, dtype, dim) -> None:
        data = np.arange(4).reshape([2, 2])

        da1 = Dataset(
            {
                "data": (["x1", "x2"], data),
                "x1": [0, 1],
                "x2": np.array(["a", "b"], dtype=dtype),
            }
        )
        da2 = Dataset(
            {
                "data": (["x1", "x2"], data),
                "x1": np.array([1, 2]),
                "x2": np.array(["c", "d"], dtype=dtype),
            }
        )
        actual = concat([da1, da2], dim=dim)

        assert np.issubdtype(actual.x2.dtype, dtype)

    def test_concat_avoids_index_auto_creation(self) -> None:
        # TODO once passing indexes={} directly to Dataset constructor is allowed then no need to create coords first
        coords = Coordinates(
            {"x": ConcatenatableArray(np.array([1, 2, 3]))}, indexes={}
        )
        datasets = [
            Dataset(
                {"a": (["x", "y"], ConcatenatableArray(np.zeros((3, 3))))},
                coords=coords,
            )
            for _ in range(2)
        ]
        # should not raise on concat
        combined = concat(datasets, dim="x")
        assert combined["a"].shape == (6, 3)
        assert combined["a"].dims == ("x", "y")

        # nor have auto-created any indexes
        assert combined.indexes == {}

        # should not raise on stack
        combined = concat(datasets, dim="z")
        assert combined["a"].shape == (2, 3, 3)
        assert combined["a"].dims == ("z", "x", "y")

        # nor have auto-created any indexes
        assert combined.indexes == {}

    def test_concat_avoids_index_auto_creation_new_1d_coord(self) -> None:
        # create 0D coordinates (without indexes)
        datasets = [
            Dataset(
                coords={"x": ConcatenatableArray(np.array(10))},
            )
            for _ in range(2)
        ]

        with pytest.raises(UnexpectedDataAccess):
            concat(datasets, dim="x", create_index_for_new_dim=True)

        # should not raise on concat iff create_index_for_new_dim=False
        combined = concat(datasets, dim="x", create_index_for_new_dim=False)
        assert combined["x"].shape == (2,)
        assert combined["x"].dims == ("x",)

        # nor have auto-created any indexes
        assert combined.indexes == {}

    def test_concat_promote_shape_without_creating_new_index(self) -> None:
        # different shapes but neither have indexes
        ds1 = Dataset(coords={"x": 0})
        ds2 = Dataset(data_vars={"x": [1]}).drop_indexes("x")
        actual = concat([ds1, ds2], dim="x", create_index_for_new_dim=False)
        expected = Dataset(data_vars={"x": [0, 1]}).drop_indexes("x")
        assert_identical(actual, expected, check_default_indexes=False)
        assert actual.indexes == {}


class TestConcatDataArray:
    def test_concat(self) -> None:
        ds = Dataset(
            {
                "foo": (["x", "y"], np.random.random((2, 3))),
                "bar": (["x", "y"], np.random.random((2, 3))),
            },
            {"x": [0, 1]},
        )
        foo = ds["foo"]
        bar = ds["bar"]

        # from dataset array:
        expected = DataArray(
            np.array([foo.values, bar.values]),
            dims=["w", "x", "y"],
            coords={"x": [0, 1]},
        )
        actual = concat([foo, bar], "w")
        assert_equal(expected, actual)
        # from iteration:
        grouped = [g.squeeze() for _, g in foo.groupby("x", squeeze=False)]
        stacked = concat(grouped, ds["x"])
        assert_identical(foo, stacked)
        # with an index as the 'dim' argument
        stacked = concat(grouped, pd.Index(ds["x"], name="x"))
        assert_identical(foo, stacked)

        actual2 = concat([foo[0], foo[1]], pd.Index([0, 1])).reset_coords(drop=True)
        expected = foo[:2].rename({"x": "concat_dim"})
        assert_identical(expected, actual2)

        actual3 = concat([foo[0], foo[1]], [0, 1]).reset_coords(drop=True)
        expected = foo[:2].rename({"x": "concat_dim"})
        assert_identical(expected, actual3)

        with pytest.raises(ValueError, match=r"not identical"):
            concat([foo, bar], dim="w", compat="identical")

        with pytest.raises(ValueError, match=r"not a valid argument"):
            concat([foo, bar], dim="w", data_vars="minimal")

    def test_concat_encoding(self) -> None:
        # Regression test for GH1297
        ds = Dataset(
            {
                "foo": (["x", "y"], np.random.random((2, 3))),
                "bar": (["x", "y"], np.random.random((2, 3))),
            },
            {"x": [0, 1]},
        )
        foo = ds["foo"]
        foo.encoding = {"complevel": 5}
        ds.encoding = {"unlimited_dims": "x"}
        assert concat([foo, foo], dim="x").encoding == foo.encoding
        assert concat([ds, ds], dim="x").encoding == ds.encoding

    @requires_dask
    def test_concat_lazy(self) -> None:
        import dask.array as da

        arrays = [
            DataArray(
                da.from_array(InaccessibleArray(np.zeros((3, 3))), 3), dims=["x", "y"]
            )
            for _ in range(2)
        ]
        # should not raise
        combined = concat(arrays, dim="z")
        assert combined.shape == (2, 3, 3)
        assert combined.dims == ("z", "x", "y")

    def test_concat_avoids_index_auto_creation(self) -> None:
        # TODO once passing indexes={} directly to DataArray constructor is allowed then no need to create coords first
        coords = Coordinates(
            {"x": ConcatenatableArray(np.array([1, 2, 3]))}, indexes={}
        )
        arrays = [
            DataArray(
                ConcatenatableArray(np.zeros((3, 3))),
                dims=["x", "y"],
                coords=coords,
            )
            for _ in range(2)
        ]
        # should not raise on concat
        combined = concat(arrays, dim="x")
        assert combined.shape == (6, 3)
        assert combined.dims == ("x", "y")

        # nor have auto-created any indexes
        assert combined.indexes == {}

        # should not raise on stack
        combined = concat(arrays, dim="z")
        assert combined.shape == (2, 3, 3)
        assert combined.dims == ("z", "x", "y")

        # nor have auto-created any indexes
        assert combined.indexes == {}

    @pytest.mark.parametrize("fill_value", [dtypes.NA, 2, 2.0])
    def test_concat_fill_value(self, fill_value) -> None:
        foo = DataArray([1, 2], coords=[("x", [1, 2])])
        bar = DataArray([1, 2], coords=[("x", [1, 3])])
        if fill_value == dtypes.NA:
            # if we supply the default, we expect the missing value for a
            # float array
            fill_value = np.nan
        expected = DataArray(
            [[1, 2, fill_value], [1, fill_value, 2]],
            dims=["y", "x"],
            coords={"x": [1, 2, 3]},
        )
        actual = concat((foo, bar), dim="y", fill_value=fill_value)
        assert_identical(actual, expected)

    def test_concat_join_kwarg(self) -> None:
        ds1 = Dataset(
            {"a": (("x", "y"), [[0]])}, coords={"x": [0], "y": [0]}
        ).to_dataarray()
        ds2 = Dataset(
            {"a": (("x", "y"), [[0]])}, coords={"x": [1], "y": [0.0001]}
        ).to_dataarray()

        expected: dict[JoinOptions, Any] = {}
        expected["outer"] = Dataset(
            {"a": (("x", "y"), [[0, np.nan], [np.nan, 0]])},
            {"x": [0, 1], "y": [0, 0.0001]},
        )
        expected["inner"] = Dataset(
            {"a": (("x", "y"), [[], []])}, {"x": [0, 1], "y": []}
        )
        expected["left"] = Dataset(
            {"a": (("x", "y"), np.array([0, np.nan], ndmin=2).T)},
            coords={"x": [0, 1], "y": [0]},
        )
        expected["right"] = Dataset(
            {"a": (("x", "y"), np.array([np.nan, 0], ndmin=2).T)},
            coords={"x": [0, 1], "y": [0.0001]},
        )
        expected["override"] = Dataset(
            {"a": (("x", "y"), np.array([0, 0], ndmin=2).T)},
            coords={"x": [0, 1], "y": [0]},
        )

        with pytest.raises(ValueError, match=r"cannot align.*exact.*dimensions.*'y'"):
            actual = concat([ds1, ds2], join="exact", dim="x")

        for join in expected:
            actual = concat([ds1, ds2], join=join, dim="x")
            assert_equal(actual, expected[join].to_dataarray())

    def test_concat_combine_attrs_kwarg(self) -> None:
        da1 = DataArray([0], coords=[("x", [0])], attrs={"b": 42})
        da2 = DataArray([0], coords=[("x", [1])], attrs={"b": 42, "c": 43})

        expected: dict[CombineAttrsOptions, Any] = {}
        expected["drop"] = DataArray([0, 0], coords=[("x", [0, 1])])
        expected["no_conflicts"] = DataArray(
            [0, 0], coords=[("x", [0, 1])], attrs={"b": 42, "c": 43}
        )
        expected["override"] = DataArray(
            [0, 0], coords=[("x", [0, 1])], attrs={"b": 42}
        )

        with pytest.raises(ValueError, match=r"combine_attrs='identical'"):
            actual = concat([da1, da2], dim="x", combine_attrs="identical")
        with pytest.raises(ValueError, match=r"combine_attrs='no_conflicts'"):
            da3 = da2.copy(deep=True)
            da3.attrs["b"] = 44
            actual = concat([da1, da3], dim="x", combine_attrs="no_conflicts")

        for combine_attrs in expected:
            actual = concat([da1, da2], dim="x", combine_attrs=combine_attrs)
            assert_identical(actual, expected[combine_attrs])

    @pytest.mark.parametrize("dtype", [str, bytes])
    @pytest.mark.parametrize("dim", ["x1", "x2"])
    def test_concat_str_dtype(self, dtype, dim) -> None:
        data = np.arange(4).reshape([2, 2])

        da1 = DataArray(
            data=data,
            dims=["x1", "x2"],
            coords={"x1": [0, 1], "x2": np.array(["a", "b"], dtype=dtype)},
        )
        da2 = DataArray(
            data=data,
            dims=["x1", "x2"],
            coords={"x1": np.array([1, 2]), "x2": np.array(["c", "d"], dtype=dtype)},
        )
        actual = concat([da1, da2], dim=dim)

        assert np.issubdtype(actual.x2.dtype, dtype)

    def test_concat_coord_name(self) -> None:
        da = DataArray([0], dims="a")
        da_concat = concat([da, da], dim=DataArray([0, 1], dims="b"))
        assert list(da_concat.coords) == ["b"]

        da_concat_std = concat([da, da], dim=DataArray([0, 1]))
        assert list(da_concat_std.coords) == ["dim_0"]


@pytest.mark.parametrize("attr1", ({"a": {"meta": [10, 20, 30]}}, {"a": [1, 2, 3]}, {}))
@pytest.mark.parametrize("attr2", ({"a": [1, 2, 3]}, {}))
def test_concat_attrs_first_variable(attr1, attr2) -> None:
    arrs = [
        DataArray([[1], [2]], dims=["x", "y"], attrs=attr1),
        DataArray([[3], [4]], dims=["x", "y"], attrs=attr2),
    ]

    concat_attrs = concat(arrs, "y").attrs
    assert concat_attrs == attr1


def test_concat_merge_single_non_dim_coord():
    # TODO: annotating this func fails
    da1 = DataArray([1, 2, 3], dims="x", coords={"x": [1, 2, 3], "y": 1})
    da2 = DataArray([4, 5, 6], dims="x", coords={"x": [4, 5, 6]})

    expected = DataArray(range(1, 7), dims="x", coords={"x": range(1, 7), "y": 1})

    for coords in ["different", "minimal"]:
        actual = concat([da1, da2], "x", coords=coords)
        assert_identical(actual, expected)

    with pytest.raises(ValueError, match=r"'y' not present in all datasets."):
        concat([da1, da2], dim="x", coords="all")

    da1 = DataArray([1, 2, 3], dims="x", coords={"x": [1, 2, 3], "y": 1})
    da2 = DataArray([4, 5, 6], dims="x", coords={"x": [4, 5, 6]})
    da3 = DataArray([7, 8, 9], dims="x", coords={"x": [7, 8, 9], "y": 1})
    for coords in ["different", "all"]:
        with pytest.raises(ValueError, match=r"'y' not present in all datasets"):
            concat([da1, da2, da3], dim="x", coords=coords)


def test_concat_preserve_coordinate_order() -> None:
    x = np.arange(0, 5)
    y = np.arange(0, 10)
    time = np.arange(0, 4)
    data = np.zeros((4, 10, 5), dtype=bool)

    ds1 = Dataset(
        {"data": (["time", "y", "x"], data[0:2])},
        coords={"time": time[0:2], "y": y, "x": x},
    )
    ds2 = Dataset(
        {"data": (["time", "y", "x"], data[2:4])},
        coords={"time": time[2:4], "y": y, "x": x},
    )

    expected = Dataset(
        {"data": (["time", "y", "x"], data)},
        coords={"time": time, "y": y, "x": x},
    )

    actual = concat([ds1, ds2], dim="time")

    # check dimension order
    for act, exp in zip(actual.dims, expected.dims):
        assert act == exp
        assert actual.sizes[act] == expected.sizes[exp]

    # check coordinate order
    for act, exp in zip(actual.coords, expected.coords):
        assert act == exp
        assert_identical(actual.coords[act], expected.coords[exp])


def test_concat_typing_check() -> None:
    ds = Dataset({"foo": 1}, {"bar": 2})
    da = Dataset({"foo": 3}, {"bar": 4}).to_dataarray(dim="foo")

    # concatenate a list of non-homogeneous types must raise TypeError
    with pytest.raises(
        TypeError,
        match="The elements in the input list need to be either all 'Dataset's or all 'DataArray's",
    ):
        concat([ds, da], dim="foo")  # type: ignore
    with pytest.raises(
        TypeError,
        match="The elements in the input list need to be either all 'Dataset's or all 'DataArray's",
    ):
        concat([da, ds], dim="foo")  # type: ignore


def test_concat_not_all_indexes() -> None:
    ds1 = Dataset(coords={"x": ("x", [1, 2])})
    # ds2.x has no default index
    ds2 = Dataset(coords={"x": ("y", [3, 4])})

    with pytest.raises(
        ValueError, match=r"'x' must have either an index or no index in all datasets.*"
    ):
        concat([ds1, ds2], dim="x")


def test_concat_index_not_same_dim() -> None:
    ds1 = Dataset(coords={"x": ("x", [1, 2])})
    ds2 = Dataset(coords={"x": ("y", [3, 4])})
    # TODO: use public API for setting a non-default index, when available
    ds2._indexes["x"] = PandasIndex([3, 4], "y")

    with pytest.raises(
        ValueError,
        match=r"Cannot concatenate along dimension 'x' indexes with dimensions.*",
    ):
        concat([ds1, ds2], dim="x")
