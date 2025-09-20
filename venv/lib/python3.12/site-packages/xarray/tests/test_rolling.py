from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

import xarray as xr
from xarray import DataArray, Dataset, set_options
from xarray.tests import (
    assert_allclose,
    assert_equal,
    assert_identical,
    has_dask,
    requires_dask,
    requires_dask_ge_2024_11_0,
    requires_numbagg,
)

pytestmark = [
    pytest.mark.filterwarnings("error:Mean of empty slice"),
    pytest.mark.filterwarnings("error:All-NaN (slice|axis) encountered"),
]


@pytest.mark.parametrize("func", ["mean", "sum"])
@pytest.mark.parametrize("min_periods", [1, 10])
def test_cumulative(d, func, min_periods) -> None:
    # One dim
    result = getattr(d.cumulative("z", min_periods=min_periods), func)()
    expected = getattr(d.rolling(z=d["z"].size, min_periods=min_periods), func)()
    assert_identical(result, expected)

    # Multiple dim
    result = getattr(d.cumulative(["z", "x"], min_periods=min_periods), func)()
    expected = getattr(
        d.rolling(z=d["z"].size, x=d["x"].size, min_periods=min_periods),
        func,
    )()
    assert_identical(result, expected)


def test_cumulative_vs_cum(d) -> None:
    result = d.cumulative("z").sum()
    expected = d.cumsum("z")
    # cumsum drops the coord of the dimension; cumulative doesn't
    expected = expected.assign_coords(z=result["z"])
    assert_identical(result, expected)


class TestDataArrayRolling:
    @pytest.mark.parametrize("da", (1, 2), indirect=True)
    @pytest.mark.parametrize("center", [True, False])
    @pytest.mark.parametrize("size", [1, 2, 3, 7])
    def test_rolling_iter(self, da: DataArray, center: bool, size: int) -> None:
        rolling_obj = da.rolling(time=size, center=center)
        rolling_obj_mean = rolling_obj.mean()

        assert len(rolling_obj.window_labels) == len(da["time"])
        assert_identical(rolling_obj.window_labels, da["time"])

        for i, (label, window_da) in enumerate(rolling_obj):
            assert label == da["time"].isel(time=i)

            actual = rolling_obj_mean.isel(time=i)
            expected = window_da.mean("time")

            np.testing.assert_allclose(actual.values, expected.values)

    @pytest.mark.parametrize("da", (1,), indirect=True)
    def test_rolling_repr(self, da) -> None:
        rolling_obj = da.rolling(time=7)
        assert repr(rolling_obj) == "DataArrayRolling [time->7]"
        rolling_obj = da.rolling(time=7, center=True)
        assert repr(rolling_obj) == "DataArrayRolling [time->7(center)]"
        rolling_obj = da.rolling(time=7, x=3, center=True)
        assert repr(rolling_obj) == "DataArrayRolling [time->7(center),x->3(center)]"

    @requires_dask
    def test_repeated_rolling_rechunks(self) -> None:
        # regression test for GH3277, GH2514
        dat = DataArray(np.random.rand(7653, 300), dims=("day", "item"))
        dat_chunk = dat.chunk({"item": 20})
        dat_chunk.rolling(day=10).mean().rolling(day=250).std()

    def test_rolling_doc(self, da) -> None:
        rolling_obj = da.rolling(time=7)

        # argument substitution worked
        assert "`mean`" in rolling_obj.mean.__doc__

    def test_rolling_properties(self, da) -> None:
        rolling_obj = da.rolling(time=4)

        assert rolling_obj.obj.get_axis_num("time") == 1

        # catching invalid args
        with pytest.raises(ValueError, match="window must be > 0"):
            da.rolling(time=-2)

        with pytest.raises(ValueError, match="min_periods must be greater than zero"):
            da.rolling(time=2, min_periods=0)

        with pytest.raises(
            KeyError,
            match=r"\('foo',\) not found in DataArray dimensions",
        ):
            da.rolling(foo=2)

    @requires_dask
    @pytest.mark.parametrize(
        "name", ("sum", "mean", "std", "min", "max", "median", "argmin", "argmax")
    )
    @pytest.mark.parametrize("center", (True, False, None))
    @pytest.mark.parametrize("min_periods", (1, None))
    @pytest.mark.parametrize("backend", ["numpy", "dask"], indirect=True)
    def test_rolling_wrapped_bottleneck(
        self, da, name, center, min_periods, compute_backend
    ) -> None:
        bn = pytest.importorskip("bottleneck", minversion="1.1")
        # Test all bottleneck functions
        rolling_obj = da.rolling(time=7, min_periods=min_periods)

        func_name = f"move_{name}"
        actual = getattr(rolling_obj, name)()
        window = 7
        expected = getattr(bn, func_name)(
            da.values, window=window, axis=1, min_count=min_periods
        )
        # index 0 is at the rightmost edge of the window
        # need to reverse index here
        # see GH #8541
        if func_name in ["move_argmin", "move_argmax"]:
            expected = window - 1 - expected

        # Using assert_allclose because we get tiny (1e-17) differences in numbagg.
        np.testing.assert_allclose(actual.values, expected)

        with pytest.warns(DeprecationWarning, match="Reductions are applied"):
            getattr(rolling_obj, name)(dim="time")

        # Test center
        rolling_obj = da.rolling(time=7, center=center)
        actual = getattr(rolling_obj, name)()["time"]
        # Using assert_allclose because we get tiny (1e-17) differences in numbagg.
        assert_allclose(actual, da["time"])

    @requires_dask
    @pytest.mark.parametrize("name", ("mean", "count"))
    @pytest.mark.parametrize("center", (True, False, None))
    @pytest.mark.parametrize("min_periods", (1, None))
    @pytest.mark.parametrize("window", (7, 8))
    @pytest.mark.parametrize("backend", ["dask"], indirect=True)
    def test_rolling_wrapped_dask(self, da, name, center, min_periods, window) -> None:
        # dask version
        rolling_obj = da.rolling(time=window, min_periods=min_periods, center=center)
        actual = getattr(rolling_obj, name)().load()
        if name != "count":
            with pytest.warns(DeprecationWarning, match="Reductions are applied"):
                getattr(rolling_obj, name)(dim="time")
        # numpy version
        rolling_obj = da.load().rolling(
            time=window, min_periods=min_periods, center=center
        )
        expected = getattr(rolling_obj, name)()

        # using all-close because rolling over ghost cells introduces some
        # precision errors
        assert_allclose(actual, expected)

        # with zero chunked array GH:2113
        rolling_obj = da.chunk().rolling(
            time=window, min_periods=min_periods, center=center
        )
        actual = getattr(rolling_obj, name)().load()
        assert_allclose(actual, expected)

    @pytest.mark.parametrize("center", (True, None))
    def test_rolling_wrapped_dask_nochunk(self, center) -> None:
        # GH:2113
        pytest.importorskip("dask.array")

        da_day_clim = xr.DataArray(
            np.arange(1, 367), coords=[np.arange(1, 367)], dims="dayofyear"
        )
        expected = da_day_clim.rolling(dayofyear=31, center=center).mean()
        actual = da_day_clim.chunk().rolling(dayofyear=31, center=center).mean()
        assert_allclose(actual, expected)

    @pytest.mark.parametrize("center", (True, False))
    @pytest.mark.parametrize("min_periods", (None, 1, 2, 3))
    @pytest.mark.parametrize("window", (1, 2, 3, 4))
    def test_rolling_pandas_compat(
        self, center, window, min_periods, compute_backend
    ) -> None:
        s = pd.Series(np.arange(10))
        da = DataArray.from_series(s)

        if min_periods is not None and window < min_periods:
            min_periods = window

        s_rolling = s.rolling(window, center=center, min_periods=min_periods).mean()
        da_rolling = da.rolling(
            index=window, center=center, min_periods=min_periods
        ).mean()
        da_rolling_np = da.rolling(
            index=window, center=center, min_periods=min_periods
        ).reduce(np.nanmean)

        np.testing.assert_allclose(np.asarray(s_rolling.values), da_rolling.values)
        np.testing.assert_allclose(s_rolling.index, da_rolling["index"])
        np.testing.assert_allclose(np.asarray(s_rolling.values), da_rolling_np.values)
        np.testing.assert_allclose(s_rolling.index, da_rolling_np["index"])

    @pytest.mark.parametrize("center", (True, False))
    @pytest.mark.parametrize("window", (1, 2, 3, 4))
    def test_rolling_construct(self, center: bool, window: int) -> None:
        s = pd.Series(np.arange(10))
        da = DataArray.from_series(s)

        s_rolling = s.rolling(window, center=center, min_periods=1).mean()
        da_rolling = da.rolling(index=window, center=center, min_periods=1)

        da_rolling_mean = da_rolling.construct("window").mean("window")
        np.testing.assert_allclose(np.asarray(s_rolling.values), da_rolling_mean.values)
        np.testing.assert_allclose(s_rolling.index, da_rolling_mean["index"])

        # with stride
        da_rolling_mean = da_rolling.construct("window", stride=2).mean("window")
        np.testing.assert_allclose(
            np.asarray(s_rolling.values[::2]), da_rolling_mean.values
        )
        np.testing.assert_allclose(s_rolling.index[::2], da_rolling_mean["index"])

        # with fill_value
        da_rolling_mean = da_rolling.construct("window", stride=2, fill_value=0.0).mean(
            "window"
        )
        assert da_rolling_mean.isnull().sum() == 0
        assert (da_rolling_mean == 0.0).sum() >= 0

    @pytest.mark.parametrize("da", (1, 2), indirect=True)
    @pytest.mark.parametrize("center", (True, False))
    @pytest.mark.parametrize("min_periods", (None, 1, 2, 3))
    @pytest.mark.parametrize("window", (1, 2, 3, 4))
    @pytest.mark.parametrize("name", ("sum", "mean", "std", "max"))
    def test_rolling_reduce(
        self, da, center, min_periods, window, name, compute_backend
    ) -> None:
        if min_periods is not None and window < min_periods:
            min_periods = window

        if da.isnull().sum() > 1 and window == 1:
            # this causes all nan slices
            window = 2

        rolling_obj = da.rolling(time=window, center=center, min_periods=min_periods)

        # add nan prefix to numpy methods to get similar # behavior as bottleneck
        actual = rolling_obj.reduce(getattr(np, f"nan{name}"))
        expected = getattr(rolling_obj, name)()
        assert_allclose(actual, expected)
        assert actual.sizes == expected.sizes

    @pytest.mark.parametrize("center", (True, False))
    @pytest.mark.parametrize("min_periods", (None, 1, 2, 3))
    @pytest.mark.parametrize("window", (1, 2, 3, 4))
    @pytest.mark.parametrize("name", ("sum", "max"))
    def test_rolling_reduce_nonnumeric(
        self, center, min_periods, window, name, compute_backend
    ) -> None:
        da = DataArray(
            [0, np.nan, 1, 2, np.nan, 3, 4, 5, np.nan, 6, 7], dims="time"
        ).isnull()

        if min_periods is not None and window < min_periods:
            min_periods = window

        rolling_obj = da.rolling(time=window, center=center, min_periods=min_periods)

        # add nan prefix to numpy methods to get similar behavior as bottleneck
        actual = rolling_obj.reduce(getattr(np, f"nan{name}"))
        expected = getattr(rolling_obj, name)()
        assert_allclose(actual, expected)
        assert actual.sizes == expected.sizes

    def test_rolling_count_correct(self, compute_backend) -> None:
        da = DataArray([0, np.nan, 1, 2, np.nan, 3, 4, 5, np.nan, 6, 7], dims="time")

        kwargs: list[dict[str, Any]] = [
            {"time": 11, "min_periods": 1},
            {"time": 11, "min_periods": None},
            {"time": 7, "min_periods": 2},
        ]
        expecteds = [
            DataArray([1, 1, 2, 3, 3, 4, 5, 6, 6, 7, 8], dims="time"),
            DataArray(
                [
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                    np.nan,
                ],
                dims="time",
            ),
            DataArray([np.nan, np.nan, 2, 3, 3, 4, 5, 5, 5, 5, 5], dims="time"),
        ]

        for kwarg, expected in zip(kwargs, expecteds, strict=True):
            result = da.rolling(**kwarg).count()
            assert_equal(result, expected)

            result = da.to_dataset(name="var1").rolling(**kwarg).count()["var1"]
            assert_equal(result, expected)

    @pytest.mark.parametrize("da", (1,), indirect=True)
    @pytest.mark.parametrize("center", (True, False))
    @pytest.mark.parametrize("min_periods", (None, 1))
    @pytest.mark.parametrize("name", ("sum", "mean", "max"))
    def test_ndrolling_reduce(
        self, da, center, min_periods, name, compute_backend
    ) -> None:
        rolling_obj = da.rolling(time=3, x=2, center=center, min_periods=min_periods)

        actual = getattr(rolling_obj, name)()
        expected = getattr(
            getattr(
                da.rolling(time=3, center=center, min_periods=min_periods), name
            )().rolling(x=2, center=center, min_periods=min_periods),
            name,
        )()

        assert_allclose(actual, expected)
        assert actual.sizes == expected.sizes

        if name == "mean":
            # test our reimplementation of nanmean using np.nanmean
            expected = getattr(rolling_obj.construct({"time": "tw", "x": "xw"}), name)(
                ["tw", "xw"]
            )
            count = rolling_obj.count()
            if min_periods is None:
                min_periods = 1
            assert_allclose(actual, expected.where(count >= min_periods))

    @pytest.mark.parametrize("center", (True, False, (True, False)))
    @pytest.mark.parametrize("fill_value", (np.nan, 0.0))
    def test_ndrolling_construct(self, center, fill_value) -> None:
        da = DataArray(
            np.arange(5 * 6 * 7).reshape(5, 6, 7).astype(float),
            dims=["x", "y", "z"],
            coords={"x": ["a", "b", "c", "d", "e"], "y": np.arange(6)},
        )
        actual = da.rolling(x=3, z=2, center=center).construct(
            x="x1", z="z1", fill_value=fill_value
        )
        if not isinstance(center, tuple):
            center = (center, center)
        expected = (
            da.rolling(x=3, center=center[0])
            .construct(x="x1", fill_value=fill_value)
            .rolling(z=2, center=center[1])
            .construct(z="z1", fill_value=fill_value)
        )
        assert_allclose(actual, expected)

    @pytest.mark.parametrize(
        "funcname, argument",
        [
            ("reduce", (np.mean,)),
            ("mean", ()),
            ("construct", ("window_dim",)),
            ("count", ()),
        ],
    )
    def test_rolling_keep_attrs(self, funcname, argument) -> None:
        attrs_da = {"da_attr": "test"}

        data = np.linspace(10, 15, 100)
        coords = np.linspace(1, 10, 100)

        da = DataArray(
            data, dims=("coord"), coords={"coord": coords}, attrs=attrs_da, name="name"
        )

        # attrs are now kept per default
        func = getattr(da.rolling(dim={"coord": 5}), funcname)
        result = func(*argument)
        assert result.attrs == attrs_da
        assert result.name == "name"

        # discard attrs
        func = getattr(da.rolling(dim={"coord": 5}), funcname)
        result = func(*argument, keep_attrs=False)
        assert result.attrs == {}
        assert result.name == "name"

        # test discard attrs using global option
        func = getattr(da.rolling(dim={"coord": 5}), funcname)
        with set_options(keep_attrs=False):
            result = func(*argument)
        assert result.attrs == {}
        assert result.name == "name"

        # keyword takes precedence over global option
        func = getattr(da.rolling(dim={"coord": 5}), funcname)
        with set_options(keep_attrs=False):
            result = func(*argument, keep_attrs=True)
        assert result.attrs == attrs_da
        assert result.name == "name"

        func = getattr(da.rolling(dim={"coord": 5}), funcname)
        with set_options(keep_attrs=True):
            result = func(*argument, keep_attrs=False)
        assert result.attrs == {}
        assert result.name == "name"

    @requires_dask
    @pytest.mark.parametrize("dtype", ["int", "float32", "float64"])
    def test_rolling_dask_dtype(self, dtype) -> None:
        data = DataArray(
            np.array([1, 2, 3], dtype=dtype), dims="x", coords={"x": [1, 2, 3]}
        )
        unchunked_result = data.rolling(x=3, min_periods=1).mean()
        chunked_result = data.chunk({"x": 1}).rolling(x=3, min_periods=1).mean()
        assert chunked_result.dtype == unchunked_result.dtype

    def test_rolling_mean_bool(self) -> None:
        bool_raster = DataArray(
            data=[0, 1, 1, 0, 1, 0],
            dims=("x"),
        ).astype(bool)

        expected = DataArray(
            data=[np.nan, 2 / 3, 2 / 3, 2 / 3, 1 / 3, np.nan],
            dims=("x"),
        )

        result = bool_raster.rolling(x=3, center=True).mean()
        assert_allclose(result, expected)


@requires_numbagg
class TestDataArrayRollingExp:
    @pytest.mark.parametrize("dim", ["time", "x"])
    @pytest.mark.parametrize(
        "window_type, window",
        [["span", 5], ["alpha", 0.5], ["com", 0.5], ["halflife", 5]],
    )
    @pytest.mark.parametrize("backend", ["numpy"], indirect=True)
    @pytest.mark.parametrize("func", ["mean", "sum", "var", "std"])
    def test_rolling_exp_runs(self, da, dim, window_type, window, func) -> None:
        da = da.where(da > 0.2)

        rolling_exp = da.rolling_exp(window_type=window_type, **{dim: window})
        result = getattr(rolling_exp, func)()
        assert isinstance(result, DataArray)

    @pytest.mark.parametrize("dim", ["time", "x"])
    @pytest.mark.parametrize(
        "window_type, window",
        [["span", 5], ["alpha", 0.5], ["com", 0.5], ["halflife", 5]],
    )
    @pytest.mark.parametrize("backend", ["numpy"], indirect=True)
    def test_rolling_exp_mean_pandas(self, da, dim, window_type, window) -> None:
        da = da.isel(a=0).where(lambda x: x > 0.2)

        result = da.rolling_exp(window_type=window_type, **{dim: window}).mean()
        assert isinstance(result, DataArray)

        pandas_array = da.to_pandas()
        assert pandas_array.index.name == "time"
        if dim == "x":
            pandas_array = pandas_array.T
        expected = xr.DataArray(
            pandas_array.ewm(**{window_type: window}).mean()
        ).transpose(*da.dims)

        assert_allclose(expected.variable, result.variable)

    @pytest.mark.parametrize("backend", ["numpy"], indirect=True)
    @pytest.mark.parametrize("func", ["mean", "sum"])
    def test_rolling_exp_keep_attrs(self, da, func) -> None:
        attrs = {"attrs": "da"}
        da.attrs = attrs

        # Equivalent of `da.rolling_exp(time=10).mean`
        rolling_exp_func = getattr(da.rolling_exp(time=10), func)

        # attrs are kept per default
        result = rolling_exp_func()
        assert result.attrs == attrs

        # discard attrs
        result = rolling_exp_func(keep_attrs=False)
        assert result.attrs == {}

        # test discard attrs using global option
        with set_options(keep_attrs=False):
            result = rolling_exp_func()
        assert result.attrs == {}

        # keyword takes precedence over global option
        with set_options(keep_attrs=False):
            result = rolling_exp_func(keep_attrs=True)
        assert result.attrs == attrs

        with set_options(keep_attrs=True):
            result = rolling_exp_func(keep_attrs=False)
        assert result.attrs == {}

        with pytest.warns(
            UserWarning,
            match="Passing ``keep_attrs`` to ``rolling_exp`` has no effect.",
        ):
            da.rolling_exp(time=10, keep_attrs=True)


class TestDatasetRolling:
    @pytest.mark.parametrize(
        "funcname, argument",
        [
            ("reduce", (np.mean,)),
            ("mean", ()),
            ("construct", ("window_dim",)),
            ("count", ()),
        ],
    )
    def test_rolling_keep_attrs(self, funcname, argument) -> None:
        global_attrs = {"units": "test", "long_name": "testing"}
        da_attrs = {"da_attr": "test"}
        da_not_rolled_attrs = {"da_not_rolled_attr": "test"}

        data = np.linspace(10, 15, 100)
        coords = np.linspace(1, 10, 100)

        ds = Dataset(
            data_vars={"da": ("coord", data), "da_not_rolled": ("no_coord", data)},
            coords={"coord": coords},
            attrs=global_attrs,
        )
        ds.da.attrs = da_attrs
        ds.da_not_rolled.attrs = da_not_rolled_attrs

        # attrs are now kept per default
        func = getattr(ds.rolling(dim={"coord": 5}), funcname)
        result = func(*argument)
        assert result.attrs == global_attrs
        assert result.da.attrs == da_attrs
        assert result.da_not_rolled.attrs == da_not_rolled_attrs
        assert result.da.name == "da"
        assert result.da_not_rolled.name == "da_not_rolled"

        # discard attrs
        func = getattr(ds.rolling(dim={"coord": 5}), funcname)
        result = func(*argument, keep_attrs=False)
        assert result.attrs == {}
        assert result.da.attrs == {}
        assert result.da_not_rolled.attrs == {}
        assert result.da.name == "da"
        assert result.da_not_rolled.name == "da_not_rolled"

        # test discard attrs using global option
        func = getattr(ds.rolling(dim={"coord": 5}), funcname)
        with set_options(keep_attrs=False):
            result = func(*argument)

        assert result.attrs == {}
        assert result.da.attrs == {}
        assert result.da_not_rolled.attrs == {}
        assert result.da.name == "da"
        assert result.da_not_rolled.name == "da_not_rolled"

        # keyword takes precedence over global option
        func = getattr(ds.rolling(dim={"coord": 5}), funcname)
        with set_options(keep_attrs=False):
            result = func(*argument, keep_attrs=True)

        assert result.attrs == global_attrs
        assert result.da.attrs == da_attrs
        assert result.da_not_rolled.attrs == da_not_rolled_attrs
        assert result.da.name == "da"
        assert result.da_not_rolled.name == "da_not_rolled"

        func = getattr(ds.rolling(dim={"coord": 5}), funcname)
        with set_options(keep_attrs=True):
            result = func(*argument, keep_attrs=False)

        assert result.attrs == {}
        assert result.da.attrs == {}
        assert result.da_not_rolled.attrs == {}
        assert result.da.name == "da"
        assert result.da_not_rolled.name == "da_not_rolled"

    def test_rolling_properties(self, ds) -> None:
        # catching invalid args
        with pytest.raises(ValueError, match="window must be > 0"):
            ds.rolling(time=-2)
        with pytest.raises(ValueError, match="min_periods must be greater than zero"):
            ds.rolling(time=2, min_periods=0)
        with pytest.raises(KeyError, match="time2"):
            ds.rolling(time2=2)
        with pytest.raises(
            KeyError,
            match=r"\('foo',\) not found in Dataset dimensions",
        ):
            ds.rolling(foo=2)

    @requires_dask_ge_2024_11_0
    def test_rolling_construct_automatic_rechunk(self):
        import dask

        # Construct dataset with chunk size of (400, 400, 1) or 1.22 MiB
        da = DataArray(
            dims=["latitude", "longitude", "time"],
            data=dask.array.random.random((400, 400, 400), chunks=(-1, -1, 1)),
            name="foo",
        )

        for obj in [da, da.to_dataset()]:
            # Dataset now has chunks of size (400, 400, 100 100) or 11.92 GiB
            rechunked = obj.rolling(time=100, center=True).construct(
                "window",
                sliding_window_view_kwargs=dict(
                    automatic_rechunk=True, writeable=False
                ),
            )
            not_rechunked = obj.rolling(time=100, center=True).construct(
                "window",
                sliding_window_view_kwargs=dict(
                    automatic_rechunk=False, writeable=True
                ),
            )
            assert rechunked.chunksizes != not_rechunked.chunksizes

            roller = obj.isel(time=slice(30)).rolling(time=10, center=True)
            one = roller.reduce(
                np.sum, sliding_window_view_kwargs=dict(automatic_rechunk=True)
            )
            two = roller.reduce(
                np.sum, sliding_window_view_kwargs=dict(automatic_rechunk=False)
            )
            assert_identical(one, two)

    @pytest.mark.parametrize(
        "name", ("sum", "mean", "std", "var", "min", "max", "median")
    )
    @pytest.mark.parametrize("center", (True, False, None))
    @pytest.mark.parametrize("min_periods", (1, None))
    @pytest.mark.parametrize("key", ("z1", "z2"))
    @pytest.mark.parametrize("backend", ["numpy"], indirect=True)
    def test_rolling_wrapped_bottleneck(
        self, ds, name, center, min_periods, key, compute_backend
    ) -> None:
        bn = pytest.importorskip("bottleneck", minversion="1.1")

        # Test all bottleneck functions
        rolling_obj = ds.rolling(time=7, min_periods=min_periods)

        func_name = f"move_{name}"
        actual = getattr(rolling_obj, name)()
        if key == "z1":  # z1 does not depend on 'Time' axis. Stored as it is.
            expected = ds[key]
        elif key == "z2":
            expected = getattr(bn, func_name)(
                ds[key].values, window=7, axis=0, min_count=min_periods
            )
        else:
            raise ValueError
        np.testing.assert_allclose(actual[key].values, expected)

        # Test center
        rolling_obj = ds.rolling(time=7, center=center)
        actual = getattr(rolling_obj, name)()["time"]
        assert_allclose(actual, ds["time"])

    @pytest.mark.parametrize("center", (True, False))
    @pytest.mark.parametrize("min_periods", (None, 1, 2, 3))
    @pytest.mark.parametrize("window", (1, 2, 3, 4))
    def test_rolling_pandas_compat(self, center, window, min_periods) -> None:
        df = pd.DataFrame(
            {
                "x": np.random.randn(20),
                "y": np.random.randn(20),
                "time": np.linspace(0, 1, 20),
            }
        )
        ds = Dataset.from_dataframe(df)

        if min_periods is not None and window < min_periods:
            min_periods = window

        df_rolling = df.rolling(window, center=center, min_periods=min_periods).mean()
        ds_rolling = ds.rolling(
            index=window, center=center, min_periods=min_periods
        ).mean()

        np.testing.assert_allclose(
            np.asarray(df_rolling["x"].values), ds_rolling["x"].values
        )
        np.testing.assert_allclose(df_rolling.index, ds_rolling["index"])

    @pytest.mark.parametrize("center", (True, False))
    @pytest.mark.parametrize("window", (1, 2, 3, 4))
    def test_rolling_construct(self, center: bool, window: int) -> None:
        df = pd.DataFrame(
            {
                "x": np.random.randn(20),
                "y": np.random.randn(20),
                "time": np.linspace(0, 1, 20),
            }
        )

        ds = Dataset.from_dataframe(df)
        df_rolling = df.rolling(window, center=center, min_periods=1).mean()
        ds_rolling = ds.rolling(index=window, center=center)

        ds_rolling_mean = ds_rolling.construct("window").mean("window")
        np.testing.assert_allclose(
            np.asarray(df_rolling["x"].values), ds_rolling_mean["x"].values
        )
        np.testing.assert_allclose(df_rolling.index, ds_rolling_mean["index"])

        # with fill_value
        ds_rolling_mean = ds_rolling.construct("window", stride=2, fill_value=0.0).mean(
            "window"
        )
        assert (ds_rolling_mean.isnull().sum() == 0).to_dataarray(dim="vars").all()
        assert (ds_rolling_mean["x"] == 0.0).sum() >= 0

    @pytest.mark.parametrize("center", (True, False))
    @pytest.mark.parametrize("window", (1, 2, 3, 4))
    def test_rolling_construct_stride(self, center: bool, window: int) -> None:
        df = pd.DataFrame(
            {
                "x": np.random.randn(20),
                "y": np.random.randn(20),
                "time": np.linspace(0, 1, 20),
            }
        )
        ds = Dataset.from_dataframe(df)
        df_rolling_mean = df.rolling(window, center=center, min_periods=1).mean()

        # With an index (dimension coordinate)
        ds_rolling = ds.rolling(index=window, center=center)
        ds_rolling_mean = ds_rolling.construct("w", stride=2).mean("w")
        np.testing.assert_allclose(
            np.asarray(df_rolling_mean["x"][::2].values), ds_rolling_mean["x"].values
        )
        np.testing.assert_allclose(df_rolling_mean.index[::2], ds_rolling_mean["index"])

        # Without index (https://github.com/pydata/xarray/issues/7021)
        ds2 = ds.drop_vars("index")
        ds2_rolling = ds2.rolling(index=window, center=center)
        ds2_rolling_mean = ds2_rolling.construct("w", stride=2).mean("w")
        np.testing.assert_allclose(
            np.asarray(df_rolling_mean["x"][::2].values), ds2_rolling_mean["x"].values
        )

        # Mixed coordinates, indexes and 2D coordinates
        ds3 = xr.Dataset(
            {"x": ("t", range(20)), "x2": ("y", range(5))},
            {
                "t": range(20),
                "y": ("y", range(5)),
                "t2": ("t", range(20)),
                "y2": ("y", range(5)),
                "yt": (["t", "y"], np.ones((20, 5))),
            },
        )
        ds3_rolling = ds3.rolling(t=window, center=center)
        ds3_rolling_mean = ds3_rolling.construct("w", stride=2).mean("w")
        for coord in ds3.coords:
            assert coord in ds3_rolling_mean.coords

    @pytest.mark.slow
    @pytest.mark.parametrize("ds", (1, 2), indirect=True)
    @pytest.mark.parametrize("center", (True, False))
    @pytest.mark.parametrize("min_periods", (None, 1, 2, 3))
    @pytest.mark.parametrize("window", (1, 2, 3, 4))
    @pytest.mark.parametrize(
        "name", ("sum", "mean", "std", "var", "min", "max", "median")
    )
    def test_rolling_reduce(self, ds, center, min_periods, window, name) -> None:
        if min_periods is not None and window < min_periods:
            min_periods = window

        if name == "std" and window == 1:
            pytest.skip("std with window == 1 is unstable in bottleneck")

        rolling_obj = ds.rolling(time=window, center=center, min_periods=min_periods)

        # add nan prefix to numpy methods to get similar behavior as bottleneck
        actual = rolling_obj.reduce(getattr(np, f"nan{name}"))
        expected = getattr(rolling_obj, name)()
        assert_allclose(actual, expected)
        assert ds.sizes == actual.sizes
        # make sure the order of data_var are not changed.
        assert list(ds.data_vars.keys()) == list(actual.data_vars.keys())

        # Make sure the dimension order is restored
        for key, src_var in ds.data_vars.items():
            assert src_var.dims == actual[key].dims

    @pytest.mark.parametrize("ds", (2,), indirect=True)
    @pytest.mark.parametrize("center", (True, False))
    @pytest.mark.parametrize("min_periods", (None, 1))
    @pytest.mark.parametrize("name", ("sum", "max"))
    @pytest.mark.parametrize("dask", (True, False))
    def test_ndrolling_reduce(self, ds, center, min_periods, name, dask) -> None:
        if dask and has_dask:
            ds = ds.chunk({"x": 4})

        rolling_obj = ds.rolling(time=4, x=3, center=center, min_periods=min_periods)

        actual = getattr(rolling_obj, name)()
        expected = getattr(
            getattr(
                ds.rolling(time=4, center=center, min_periods=min_periods), name
            )().rolling(x=3, center=center, min_periods=min_periods),
            name,
        )()
        assert_allclose(actual, expected)
        assert actual.sizes == expected.sizes

        # Do it in the opposite order
        expected = getattr(
            getattr(
                ds.rolling(x=3, center=center, min_periods=min_periods), name
            )().rolling(time=4, center=center, min_periods=min_periods),
            name,
        )()

        assert_allclose(actual, expected)
        assert actual.sizes == expected.sizes

    @pytest.mark.parametrize("center", (True, False, (True, False)))
    @pytest.mark.parametrize("fill_value", (np.nan, 0.0))
    @pytest.mark.parametrize("dask", (True, False))
    def test_ndrolling_construct(self, center, fill_value, dask) -> None:
        da = DataArray(
            np.arange(5 * 6 * 7).reshape(5, 6, 7).astype(float),
            dims=["x", "y", "z"],
            coords={"x": ["a", "b", "c", "d", "e"], "y": np.arange(6)},
        )
        ds = xr.Dataset({"da": da})
        if dask and has_dask:
            ds = ds.chunk({"x": 4})

        actual = ds.rolling(x=3, z=2, center=center).construct(
            x="x1", z="z1", fill_value=fill_value
        )
        if not isinstance(center, tuple):
            center = (center, center)
        expected = (
            ds.rolling(x=3, center=center[0])
            .construct(x="x1", fill_value=fill_value)
            .rolling(z=2, center=center[1])
            .construct(z="z1", fill_value=fill_value)
        )
        assert_allclose(actual, expected)

    @requires_dask
    @pytest.mark.filterwarnings("error")
    @pytest.mark.parametrize("ds", (2,), indirect=True)
    @pytest.mark.parametrize("name", ("mean", "max"))
    def test_raise_no_warning_dask_rolling_assert_close(self, ds, name) -> None:
        """
        This is a puzzle — I can't easily find the source of the warning. It
        requires `assert_allclose` to be run, for the `ds` param to be 2, and is
        different for `mean` and `max`. `sum` raises no warning.
        """

        ds = ds.chunk({"x": 4})

        rolling_obj = ds.rolling(time=4, x=3)

        actual = getattr(rolling_obj, name)()
        expected = getattr(getattr(ds.rolling(time=4), name)().rolling(x=3), name)()
        assert_allclose(actual, expected)


@requires_numbagg
class TestDatasetRollingExp:
    @pytest.mark.parametrize(
        "backend", ["numpy", pytest.param("dask", marks=requires_dask)], indirect=True
    )
    def test_rolling_exp(self, ds) -> None:
        result = ds.rolling_exp(time=10, window_type="span").mean()
        assert isinstance(result, Dataset)

    @pytest.mark.parametrize("backend", ["numpy"], indirect=True)
    def test_rolling_exp_keep_attrs(self, ds) -> None:
        attrs_global = {"attrs": "global"}
        attrs_z1 = {"attr": "z1"}

        ds.attrs = attrs_global
        ds.z1.attrs = attrs_z1

        # attrs are kept per default
        result = ds.rolling_exp(time=10).mean()
        assert result.attrs == attrs_global
        assert result.z1.attrs == attrs_z1

        # discard attrs
        result = ds.rolling_exp(time=10).mean(keep_attrs=False)
        assert result.attrs == {}
        # TODO: from #8114 — this arguably should be empty, but `apply_ufunc` doesn't do
        # that at the moment. We should change in `apply_func` rather than
        # special-case it here.
        #
        # assert result.z1.attrs == {}

        # test discard attrs using global option
        with set_options(keep_attrs=False):
            result = ds.rolling_exp(time=10).mean()
        assert result.attrs == {}
        # See above
        # assert result.z1.attrs == {}

        # keyword takes precedence over global option
        with set_options(keep_attrs=False):
            result = ds.rolling_exp(time=10).mean(keep_attrs=True)
        assert result.attrs == attrs_global
        assert result.z1.attrs == attrs_z1

        with set_options(keep_attrs=True):
            result = ds.rolling_exp(time=10).mean(keep_attrs=False)
        assert result.attrs == {}
        # See above
        # assert result.z1.attrs == {}

        with pytest.warns(
            UserWarning,
            match="Passing ``keep_attrs`` to ``rolling_exp`` has no effect.",
        ):
            ds.rolling_exp(time=10, keep_attrs=True)
