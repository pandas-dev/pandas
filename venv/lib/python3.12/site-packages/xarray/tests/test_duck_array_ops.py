from __future__ import annotations

import copy
import datetime as dt
import pickle
import warnings

import numpy as np
import pandas as pd
import pytest
from numpy import array, nan

from xarray import DataArray, Dataset, concat, date_range
from xarray.coding.times import _NS_PER_TIME_DELTA
from xarray.core import dtypes, duck_array_ops
from xarray.core.duck_array_ops import (
    array_notnull_equiv,
    concatenate,
    count,
    first,
    gradient,
    last,
    least_squares,
    mean,
    np_timedelta64_to_float,
    pd_timedelta_to_float,
    push,
    py_timedelta_to_float,
    stack,
    timedelta_to_numeric,
    where,
)
from xarray.core.extension_array import PandasExtensionArray
from xarray.core.types import NPDatetimeUnitOptions, PDDatetimeUnitOptions
from xarray.namedarray.pycompat import array_type
from xarray.testing import assert_allclose, assert_equal, assert_identical
from xarray.tests import (
    arm_xfail,
    assert_array_equal,
    has_dask,
    has_scipy,
    raise_if_dask_computes,
    requires_bottleneck,
    requires_cftime,
    requires_dask,
    requires_pyarrow,
)

dask_array_type = array_type("dask")


@pytest.fixture
def categorical1():
    return pd.Categorical(["cat1", "cat2", "cat2", "cat1", "cat2"])


@pytest.fixture
def categorical2():
    return pd.Categorical(["cat2", "cat1", "cat2", "cat3", "cat1"])


try:
    import pyarrow as pa

    @pytest.fixture
    def arrow1():
        return pd.arrays.ArrowExtensionArray(
            pa.array([{"x": 1, "y": True}, {"x": 2, "y": False}])
        )

    @pytest.fixture
    def arrow2():
        return pd.arrays.ArrowExtensionArray(
            pa.array([{"x": 3, "y": False}, {"x": 4, "y": True}])
        )

except ImportError:
    pass


@pytest.fixture
def int1():
    return pd.arrays.IntegerArray(
        np.array([1, 2, 3, 4, 5]), np.array([True, False, False, True, True])
    )


@pytest.fixture
def int2():
    return pd.arrays.IntegerArray(
        np.array([6, 7, 8, 9, 10]), np.array([True, True, False, True, False])
    )


class TestOps:
    @pytest.fixture(autouse=True)
    def setUp(self):
        self.x = array(
            [
                [
                    [nan, nan, 2.0, nan],
                    [nan, 5.0, 6.0, nan],
                    [8.0, 9.0, 10.0, nan],
                ],
                [
                    [nan, 13.0, 14.0, 15.0],
                    [nan, 17.0, 18.0, nan],
                    [nan, 21.0, nan, nan],
                ],
            ]
        )

    def test_first(self):
        expected_results = [
            array([[nan, 13, 2, 15], [nan, 5, 6, nan], [8, 9, 10, nan]]),
            array([[8, 5, 2, nan], [nan, 13, 14, 15]]),
            array([[2, 5, 8], [13, 17, 21]]),
        ]
        for axis, expected in zip(
            [0, 1, 2, -3, -2, -1], 2 * expected_results, strict=True
        ):
            actual = first(self.x, axis)
            assert_array_equal(expected, actual)

        expected = self.x[0]
        actual = first(self.x, axis=0, skipna=False)
        assert_array_equal(expected, actual)

        expected = self.x[..., 0]
        actual = first(self.x, axis=-1, skipna=False)
        assert_array_equal(expected, actual)

        with pytest.raises(IndexError, match=r"out of bounds"):
            first(self.x, 3)

    def test_last(self):
        expected_results = [
            array([[nan, 13, 14, 15], [nan, 17, 18, nan], [8, 21, 10, nan]]),
            array([[8, 9, 10, nan], [nan, 21, 18, 15]]),
            array([[2, 6, 10], [15, 18, 21]]),
        ]
        for axis, expected in zip(
            [0, 1, 2, -3, -2, -1], 2 * expected_results, strict=True
        ):
            actual = last(self.x, axis)
            assert_array_equal(expected, actual)

        expected = self.x[-1]
        actual = last(self.x, axis=0, skipna=False)
        assert_array_equal(expected, actual)

        expected = self.x[..., -1]
        actual = last(self.x, axis=-1, skipna=False)
        assert_array_equal(expected, actual)

        with pytest.raises(IndexError, match=r"out of bounds"):
            last(self.x, 3)

    def test_count(self):
        assert 12 == count(self.x)

        expected = array([[1, 2, 3], [3, 2, 1]])
        assert_array_equal(expected, count(self.x, axis=-1))

        assert 1 == count(np.datetime64("2000-01-01"))

    def test_where_type_promotion(self):
        result = where(np.array([True, False]), np.array([1, 2]), np.array(["a", "b"]))
        assert_array_equal(result, np.array([1, "b"], dtype=object))

        result = where([True, False], np.array([1, 2], np.float32), np.nan)
        assert result.dtype == np.float32
        assert_array_equal(result, np.array([1, np.nan], dtype=np.float32))

    def test_where_extension_duck_array(self, categorical1, categorical2):
        where_res = where(
            np.array([True, False, True, False, False]),
            PandasExtensionArray(categorical1),
            PandasExtensionArray(categorical2),
        )
        assert isinstance(where_res, PandasExtensionArray)
        assert (
            where_res == pd.Categorical(["cat1", "cat1", "cat2", "cat3", "cat1"])
        ).all()

    def test_concatenate_extension_duck_array(self, categorical1, categorical2):
        concate_res = concatenate(
            [PandasExtensionArray(categorical1), PandasExtensionArray(categorical2)]
        )
        assert isinstance(concate_res, PandasExtensionArray)
        assert (
            concate_res
            == type(categorical1)._concat_same_type((categorical1, categorical2))
        ).all()

    @requires_pyarrow
    def test_extension_array_pyarrow_concatenate(self, arrow1, arrow2):
        concatenated = concatenate(
            (PandasExtensionArray(arrow1), PandasExtensionArray(arrow2))
        )
        assert concatenated[2].array[0]["x"] == 3
        assert concatenated[3].array[0]["y"]

    @requires_pyarrow
    def test_extension_array_copy_arrow_type(self):
        arr = pd.array([pd.NA, 1, 2], dtype="int64[pyarrow]")
        # Relying on the `__getattr__` of `PandasExtensionArray` to do the deep copy
        # recursively only fails for `int64[pyarrow]` and similar types so this
        # test ensures that copying still works there.
        assert isinstance(
            copy.deepcopy(PandasExtensionArray(arr), memo=None).array, type(arr)
        )

    def test___getitem__extension_duck_array(self, categorical1):
        extension_duck_array = PandasExtensionArray(categorical1)
        assert (extension_duck_array[0:2] == categorical1[0:2]).all()
        assert isinstance(extension_duck_array[0:2], PandasExtensionArray)
        assert extension_duck_array[0] == categorical1[0]
        assert isinstance(extension_duck_array[0], PandasExtensionArray)
        mask = [True, False, True, False, True]
        assert (extension_duck_array[mask] == categorical1[mask]).all()

    def test__setitem__extension_duck_array(self, categorical1):
        extension_duck_array = PandasExtensionArray(categorical1)
        extension_duck_array[2] = "cat1"  # already existing category
        assert extension_duck_array[2] == "cat1"
        with pytest.raises(TypeError, match="Cannot setitem on a Categorical"):
            extension_duck_array[2] = "cat4"  # new category

    def test_stack_type_promotion(self):
        result = stack([1, "b"])
        assert_array_equal(result, np.array([1, "b"], dtype=object))

    def test_concatenate_type_promotion(self):
        result = concatenate([np.array([1]), np.array(["b"])])
        assert_array_equal(result, np.array([1, "b"], dtype=object))

    @pytest.mark.filterwarnings("error")
    def test_all_nan_arrays(self):
        assert np.isnan(mean([np.nan, np.nan]))


@requires_dask
class TestDaskOps(TestOps):
    @pytest.fixture(autouse=True)
    def setUp(self):
        import dask.array

        self.x = dask.array.from_array(
            [
                [
                    [nan, nan, 2.0, nan],
                    [nan, 5.0, 6.0, nan],
                    [8.0, 9.0, 10.0, nan],
                ],
                [
                    [nan, 13.0, 14.0, 15.0],
                    [nan, 17.0, 18.0, nan],
                    [nan, 21.0, nan, nan],
                ],
            ],
            chunks=(2, 1, 2),
        )


def test_cumsum_1d():
    inputs = np.array([0, 1, 2, 3])
    expected = np.array([0, 1, 3, 6])
    actual = duck_array_ops.cumsum(inputs)
    assert_array_equal(expected, actual)

    actual = duck_array_ops.cumsum(inputs, axis=0)
    assert_array_equal(expected, actual)

    actual = duck_array_ops.cumsum(inputs, axis=-1)
    assert_array_equal(expected, actual)

    actual = duck_array_ops.cumsum(inputs, axis=(0,))
    assert_array_equal(expected, actual)

    actual = duck_array_ops.cumsum(inputs, axis=())
    assert_array_equal(inputs, actual)


def test_cumsum_2d():
    inputs = np.array([[1, 2], [3, 4]])

    expected = np.array([[1, 3], [4, 10]])
    actual = duck_array_ops.cumsum(inputs)
    assert_array_equal(expected, actual)

    actual = duck_array_ops.cumsum(inputs, axis=(0, 1))
    assert_array_equal(expected, actual)

    actual = duck_array_ops.cumsum(inputs, axis=())
    assert_array_equal(inputs, actual)


def test_cumprod_2d():
    inputs = np.array([[1, 2], [3, 4]])

    expected = np.array([[1, 2], [3, 2 * 3 * 4]])
    actual = duck_array_ops.cumprod(inputs)
    assert_array_equal(expected, actual)

    actual = duck_array_ops.cumprod(inputs, axis=(0, 1))
    assert_array_equal(expected, actual)

    actual = duck_array_ops.cumprod(inputs, axis=())
    assert_array_equal(inputs, actual)


class TestArrayNotNullEquiv:
    @pytest.mark.parametrize(
        "arr1, arr2",
        [
            (np.array([1, 2, 3]), np.array([1, 2, 3])),
            (np.array([1, 2, np.nan]), np.array([1, np.nan, 3])),
            (np.array([np.nan, 2, np.nan]), np.array([1, np.nan, np.nan])),
        ],
    )
    def test_equal(self, arr1, arr2):
        assert array_notnull_equiv(arr1, arr2)

    def test_some_not_equal(self):
        a = np.array([1, 2, 4])
        b = np.array([1, np.nan, 3])
        assert not array_notnull_equiv(a, b)

    def test_wrong_shape(self):
        a = np.array([[1, np.nan, np.nan, 4]])
        b = np.array([[1, 2], [np.nan, 4]])
        assert not array_notnull_equiv(a, b)

    @pytest.mark.parametrize(
        "val1, val2, val3, null",
        [
            (
                np.datetime64("2000"),
                np.datetime64("2001"),
                np.datetime64("2002"),
                np.datetime64("NaT"),
            ),
            (1.0, 2.0, 3.0, np.nan),
            ("foo", "bar", "baz", None),
            ("foo", "bar", "baz", np.nan),
        ],
    )
    def test_types(self, val1, val2, val3, null):
        dtype = object if isinstance(val1, str) else None
        arr1 = np.array([val1, null, val3, null], dtype=dtype)
        arr2 = np.array([val1, val2, null, null], dtype=dtype)
        assert array_notnull_equiv(arr1, arr2)


def construct_dataarray(dim_num, dtype, contains_nan, dask):
    # dimnum <= 3
    rng = np.random.default_rng(0)
    shapes = [16, 8, 4][:dim_num]
    dims = ("x", "y", "z")[:dim_num]

    if np.issubdtype(dtype, np.floating):
        array = rng.random(shapes).astype(dtype)
    elif np.issubdtype(dtype, np.integer):
        array = rng.integers(0, 10, size=shapes).astype(dtype)
    elif np.issubdtype(dtype, np.bool_):
        array = rng.integers(0, 1, size=shapes).astype(dtype)
    elif dtype is str:
        array = rng.choice(["a", "b", "c", "d"], size=shapes)
    else:
        raise ValueError

    if contains_nan:
        inds = rng.choice(range(array.size), int(array.size * 0.2))
        dtype, fill_value = dtypes.maybe_promote(array.dtype)
        array = array.astype(dtype)
        array.flat[inds] = fill_value

    da = DataArray(array, dims=dims, coords={"x": np.arange(16)}, name="da")

    if dask and has_dask:
        chunks = dict.fromkeys(dims, 4)
        da = da.chunk(chunks)

    return da


def from_series_or_scalar(se):
    if isinstance(se, pd.Series):
        return DataArray.from_series(se)
    else:  # scalar case
        return DataArray(se)


def series_reduce(da, func, dim, **kwargs):
    """convert DataArray to pd.Series, apply pd.func, then convert back to
    a DataArray. Multiple dims cannot be specified."""

    # pd no longer accepts skipna=None https://github.com/pandas-dev/pandas/issues/44178
    if kwargs.get("skipna", True) is None:
        kwargs["skipna"] = True

    if dim is None or da.ndim == 1:
        se = da.to_series()
        return from_series_or_scalar(getattr(se, func)(**kwargs))
    else:
        dims = list(da.dims)
        dims.remove(dim)
        d = dims[0]
        da1 = [
            series_reduce(da.isel(**{d: i}), func, dim, **kwargs)
            for i in range(len(da[d]))
        ]

        if d in da.coords:
            return concat(da1, dim=da[d])
        return concat(da1, dim=d)


def assert_dask_array(da, dask):
    if dask and da.ndim > 0:
        assert isinstance(da.data, dask_array_type)


@arm_xfail
@pytest.mark.filterwarnings("ignore:All-NaN .* encountered:RuntimeWarning")
@pytest.mark.parametrize("dask", [False, True] if has_dask else [False])
def test_datetime_mean(dask: bool, time_unit: PDDatetimeUnitOptions) -> None:
    # Note: only testing numpy, as dask is broken upstream
    dtype = f"M8[{time_unit}]"
    da = DataArray(
        np.array(["2010-01-01", "NaT", "2010-01-03", "NaT", "NaT"], dtype=dtype),
        dims=["time"],
    )
    if dask:
        # Trigger use case where a chunk is full of NaT
        da = da.chunk({"time": 3})

    expect = DataArray(np.array("2010-01-02", dtype="M8[ns]"))
    expect_nat = DataArray(np.array("NaT", dtype="M8[ns]"))

    actual = da.mean()
    if dask:
        assert actual.chunks is not None
    assert_equal(actual, expect)

    actual = da.mean(skipna=False)
    if dask:
        assert actual.chunks is not None
    assert_equal(actual, expect_nat)

    # tests for 1d array full of NaT
    assert_equal(da[[1]].mean(), expect_nat)
    assert_equal(da[[1]].mean(skipna=False), expect_nat)

    # tests for a 0d array
    assert_equal(da[0].mean(), da[0])
    assert_equal(da[0].mean(skipna=False), da[0])
    assert_equal(da[1].mean(), expect_nat)
    assert_equal(da[1].mean(skipna=False), expect_nat)


@requires_cftime
@pytest.mark.parametrize("dask", [False, True])
def test_cftime_datetime_mean(dask):
    if dask and not has_dask:
        pytest.skip("requires dask")

    times = date_range("2000", periods=4, use_cftime=True)
    da = DataArray(times, dims=["time"])
    da_2d = DataArray(times.values.reshape(2, 2))

    if dask:
        da = da.chunk({"time": 2})
        da_2d = da_2d.chunk({"dim_0": 2})

    expected = da.isel(time=0)
    # one compute needed to check the array contains cftime datetimes
    with raise_if_dask_computes(max_computes=1):
        result = da.isel(time=0).mean()
    assert_dask_array(result, dask)
    assert_equal(result, expected)

    expected = DataArray(times.date_type(2000, 1, 2, 12))
    with raise_if_dask_computes(max_computes=1):
        result = da.mean()
    assert_dask_array(result, dask)
    assert_equal(result, expected)

    with raise_if_dask_computes(max_computes=1):
        result = da_2d.mean()
    assert_dask_array(result, dask)
    assert_equal(result, expected)


@pytest.mark.parametrize("dask", [False, True])
def test_mean_over_long_spanning_datetime64(dask) -> None:
    if dask and not has_dask:
        pytest.skip("requires dask")
    array = np.array(["1678-01-01", "NaT", "2260-01-01"], dtype="datetime64[ns]")
    da = DataArray(array, dims=["time"])
    if dask:
        da = da.chunk({"time": 2})
    expected = DataArray(np.array("1969-01-01", dtype="datetime64[ns]"))
    result = da.mean()
    assert_equal(result, expected)


@requires_cftime
@requires_dask
def test_mean_over_non_time_dim_of_dataset_with_dask_backed_cftime_data():
    # Regression test for part two of GH issue 5897: averaging over a non-time
    # dimension still fails if the time variable is dask-backed.
    ds = Dataset(
        {
            "var1": (
                ("time",),
                date_range("2021-10-31", periods=10, freq="D", use_cftime=True),
            ),
            "var2": (("x",), list(range(10))),
        }
    )
    expected = ds.mean("x")
    result = ds.chunk({}).mean("x")
    assert_equal(result, expected)


@requires_cftime
def test_cftime_datetime_mean_long_time_period():
    import cftime

    times = np.array(
        [
            [
                cftime.DatetimeNoLeap(400, 12, 31, 0, 0, 0, 0),
                cftime.DatetimeNoLeap(520, 12, 31, 0, 0, 0, 0),
            ],
            [
                cftime.DatetimeNoLeap(520, 12, 31, 0, 0, 0, 0),
                cftime.DatetimeNoLeap(640, 12, 31, 0, 0, 0, 0),
            ],
            [
                cftime.DatetimeNoLeap(640, 12, 31, 0, 0, 0, 0),
                cftime.DatetimeNoLeap(760, 12, 31, 0, 0, 0, 0),
            ],
        ]
    )

    da = DataArray(times, dims=["time", "d2"])
    result = da.mean("d2")
    expected = DataArray(
        [
            cftime.DatetimeNoLeap(460, 12, 31, 0, 0, 0, 0),
            cftime.DatetimeNoLeap(580, 12, 31, 0, 0, 0, 0),
            cftime.DatetimeNoLeap(700, 12, 31, 0, 0, 0, 0),
        ],
        dims=["time"],
    )
    assert_equal(result, expected)


def test_empty_axis_dtype():
    ds = Dataset()
    ds["pos"] = [1, 2, 3]
    ds["data"] = ("pos", "time"), [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
    ds["var"] = "pos", [2, 3, 4]
    assert_identical(ds.mean(dim="time")["var"], ds["var"])
    assert_identical(ds.max(dim="time")["var"], ds["var"])
    assert_identical(ds.min(dim="time")["var"], ds["var"])
    assert_identical(ds.sum(dim="time")["var"], ds["var"])


@pytest.mark.parametrize("dim_num", [1, 2])
@pytest.mark.parametrize("dtype", [float, int, np.float32, np.bool_])
@pytest.mark.parametrize("dask", [False, True])
@pytest.mark.parametrize("func", ["sum", "min", "max", "mean", "var"])
# TODO test cumsum, cumprod
@pytest.mark.parametrize("skipna", [False, True])
@pytest.mark.parametrize("aggdim", [None, "x"])
def test_reduce(dim_num, dtype, dask, func, skipna, aggdim):
    if aggdim == "y" and dim_num < 2:
        pytest.skip("dim not in this test")

    if dtype == np.bool_ and func == "mean":
        pytest.skip("numpy does not support this")

    if dask and not has_dask:
        pytest.skip("requires dask")

    if dask and skipna is False and dtype == np.bool_:
        pytest.skip("dask does not compute object-typed array")

    rtol = 1e-04 if dtype == np.float32 else 1e-05

    da = construct_dataarray(dim_num, dtype, contains_nan=True, dask=dask)
    axis = None if aggdim is None else da.get_axis_num(aggdim)

    # TODO: remove these after resolving
    # https://github.com/dask/dask/issues/3245
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "Mean of empty slice")
        warnings.filterwarnings("ignore", "All-NaN slice")
        warnings.filterwarnings("ignore", "invalid value encountered in")

        if da.dtype.kind == "O" and skipna:
            # Numpy < 1.13 does not handle object-type array.
            try:
                if skipna:
                    expected = getattr(np, f"nan{func}")(da.values, axis=axis)
                else:
                    expected = getattr(np, func)(da.values, axis=axis)

                actual = getattr(da, func)(skipna=skipna, dim=aggdim)
                assert_dask_array(actual, dask)
                np.testing.assert_allclose(
                    actual.values, np.array(expected), rtol=1.0e-4, equal_nan=True
                )
            except (TypeError, AttributeError, ZeroDivisionError):
                # TODO currently, numpy does not support some methods such as
                # nanmean for object dtype
                pass

        actual = getattr(da, func)(skipna=skipna, dim=aggdim)

        # for dask case, make sure the result is the same for numpy backend
        expected = getattr(da.compute(), func)(skipna=skipna, dim=aggdim)
        assert_allclose(actual, expected, rtol=rtol)

        # make sure the compatibility with pandas' results.
        if func in ["var", "std"]:
            expected = series_reduce(da, func, skipna=skipna, dim=aggdim, ddof=0)
            assert_allclose(actual, expected, rtol=rtol)
            # also check ddof!=0 case
            actual = getattr(da, func)(skipna=skipna, dim=aggdim, ddof=5)
            if dask:
                assert isinstance(da.data, dask_array_type)
            expected = series_reduce(da, func, skipna=skipna, dim=aggdim, ddof=5)
            assert_allclose(actual, expected, rtol=rtol)
        else:
            expected = series_reduce(da, func, skipna=skipna, dim=aggdim)
            assert_allclose(actual, expected, rtol=rtol)

        # make sure the dtype argument
        if func not in ["max", "min"]:
            actual = getattr(da, func)(skipna=skipna, dim=aggdim, dtype=float)
            assert_dask_array(actual, dask)
            assert actual.dtype == float

        # without nan
        da = construct_dataarray(dim_num, dtype, contains_nan=False, dask=dask)
        actual = getattr(da, func)(skipna=skipna)
        if dask:
            assert isinstance(da.data, dask_array_type)
        expected = getattr(np, f"nan{func}")(da.values)
        if actual.dtype == object:
            assert actual.values == np.array(expected)
        else:
            assert np.allclose(actual.values, np.array(expected), rtol=rtol)


@pytest.mark.parametrize("dim_num", [1, 2])
@pytest.mark.parametrize("dtype", [float, int, np.float32, np.bool_, str])
@pytest.mark.parametrize("contains_nan", [True, False])
@pytest.mark.parametrize("dask", [False, True])
@pytest.mark.parametrize("func", ["min", "max"])
@pytest.mark.parametrize("skipna", [False, True])
@pytest.mark.parametrize("aggdim", ["x", "y"])
def test_argmin_max(dim_num, dtype, contains_nan, dask, func, skipna, aggdim):
    # pandas-dev/pandas#16830, we do not check consistency with pandas but
    # just make sure da[da.argmin()] == da.min()

    if aggdim == "y" and dim_num < 2:
        pytest.skip("dim not in this test")

    if dask and not has_dask:
        pytest.skip("requires dask")

    if contains_nan:
        if not skipna:
            pytest.skip("numpy's argmin (not nanargmin) does not handle object-dtype")
        if skipna and np.dtype(dtype).kind in "iufc":
            pytest.skip("numpy's nanargmin raises ValueError for all nan axis")
    da = construct_dataarray(dim_num, dtype, contains_nan=contains_nan, dask=dask)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", "All-NaN slice")

        actual = da.isel(
            **{aggdim: getattr(da, "arg" + func)(dim=aggdim, skipna=skipna).compute()}
        )
        expected = getattr(da, func)(dim=aggdim, skipna=skipna)
        assert_allclose(
            actual.drop_vars(list(actual.coords)),
            expected.drop_vars(list(expected.coords)),
        )


def test_argmin_max_error():
    da = construct_dataarray(2, np.bool_, contains_nan=True, dask=False)
    da[0] = np.nan
    with pytest.raises(ValueError):
        da.argmin(dim="y")


@pytest.mark.parametrize(
    ["array", "expected"],
    [
        (
            np.array([np.datetime64("2000-01-01"), np.datetime64("NaT")]),
            np.array([False, True]),
        ),
        (
            np.array([np.timedelta64(1, "h"), np.timedelta64("NaT")]),
            np.array([False, True]),
        ),
        (
            np.array([0.0, np.nan]),
            np.array([False, True]),
        ),
        (
            np.array([1j, np.nan]),
            np.array([False, True]),
        ),
        (
            np.array(["foo", np.nan], dtype=object),
            np.array([False, True]),
        ),
        (
            np.array([1, 2], dtype=int),
            np.array([False, False]),
        ),
        (
            np.array([True, False], dtype=bool),
            np.array([False, False]),
        ),
    ],
)
def test_isnull(array, expected):
    actual = duck_array_ops.isnull(array)
    np.testing.assert_equal(expected, actual)


@requires_dask
def test_isnull_with_dask():
    da = construct_dataarray(2, np.float32, contains_nan=True, dask=True)
    assert isinstance(da.isnull().data, dask_array_type)
    assert_equal(da.isnull().load(), da.load().isnull())


@pytest.mark.skipif(not has_dask, reason="This is for dask.")
@pytest.mark.parametrize("axis", [0, -1, 1])
@pytest.mark.parametrize("edge_order", [1, 2])
def test_dask_gradient(axis, edge_order):
    import dask.array as da

    array = np.array(np.random.randn(100, 5, 40))
    x = np.exp(np.linspace(0, 1, array.shape[axis]))

    darray = da.from_array(array, chunks=[(6, 30, 30, 20, 14), 5, 8])
    expected = gradient(array, x, axis=axis, edge_order=edge_order)
    actual = gradient(darray, x, axis=axis, edge_order=edge_order)

    assert isinstance(actual, da.Array)
    assert_array_equal(actual, expected)


@pytest.mark.parametrize("dim_num", [1, 2])
@pytest.mark.parametrize("dtype", [float, int, np.float32, np.bool_])
@pytest.mark.parametrize("dask", [False, True])
@pytest.mark.parametrize("func", ["sum", "prod"])
@pytest.mark.parametrize("aggdim", [None, "x"])
@pytest.mark.parametrize("contains_nan", [True, False])
@pytest.mark.parametrize("skipna", [True, False, None])
def test_min_count(dim_num, dtype, dask, func, aggdim, contains_nan, skipna):
    if dask and not has_dask:
        pytest.skip("requires dask")

    da = construct_dataarray(dim_num, dtype, contains_nan=contains_nan, dask=dask)
    min_count = 3

    # If using Dask, the function call should be lazy.
    with raise_if_dask_computes():
        actual = getattr(da, func)(dim=aggdim, skipna=skipna, min_count=min_count)

    expected = series_reduce(da, func, skipna=skipna, dim=aggdim, min_count=min_count)
    assert_allclose(actual, expected)
    assert_dask_array(actual, dask)


@pytest.mark.parametrize("dtype", [float, int, np.float32, np.bool_])
@pytest.mark.parametrize("dask", [False, True])
@pytest.mark.parametrize("func", ["sum", "prod"])
def test_min_count_nd(dtype, dask, func):
    if dask and not has_dask:
        pytest.skip("requires dask")

    min_count = 3
    dim_num = 3
    da = construct_dataarray(dim_num, dtype, contains_nan=True, dask=dask)

    # If using Dask, the function call should be lazy.
    with raise_if_dask_computes():
        actual = getattr(da, func)(
            dim=["x", "y", "z"], skipna=True, min_count=min_count
        )

    # Supplying all dims is equivalent to supplying `...` or `None`
    expected = getattr(da, func)(dim=..., skipna=True, min_count=min_count)

    assert_allclose(actual, expected)
    assert_dask_array(actual, dask)


@pytest.mark.parametrize("dask", [False, True])
@pytest.mark.parametrize("func", ["sum", "prod"])
@pytest.mark.parametrize("dim", [None, "a", "b"])
def test_min_count_specific(dask, func, dim):
    if dask and not has_dask:
        pytest.skip("requires dask")

    # Simple array with four non-NaN values.
    da = DataArray(np.ones((6, 6), dtype=np.float64) * np.nan, dims=("a", "b"))
    da[0][0] = 2
    da[0][3] = 2
    da[3][0] = 2
    da[3][3] = 2
    if dask:
        da = da.chunk({"a": 3, "b": 3})

    # Expected result if we set min_count to the number of non-NaNs in a
    # row/column/the entire array.
    if dim:
        min_count = 2
        expected = DataArray(
            [4.0, np.nan, np.nan] * 2, dims=("a" if dim == "b" else "b",)
        )
    else:
        min_count = 4
        expected = DataArray(8.0 if func == "sum" else 16.0)

    # Check for that min_count.
    with raise_if_dask_computes():
        actual = getattr(da, func)(dim, skipna=True, min_count=min_count)
    assert_dask_array(actual, dask)
    assert_allclose(actual, expected)

    # With min_count being one higher, should get all NaN.
    min_count += 1
    expected *= np.nan
    with raise_if_dask_computes():
        actual = getattr(da, func)(dim, skipna=True, min_count=min_count)
    assert_dask_array(actual, dask)
    assert_allclose(actual, expected)


@pytest.mark.parametrize("func", ["sum", "prod"])
def test_min_count_dataset(func):
    da = construct_dataarray(2, dtype=float, contains_nan=True, dask=False)
    ds = Dataset({"var1": da}, coords={"scalar": 0})
    actual = getattr(ds, func)(dim="x", skipna=True, min_count=3)["var1"]
    expected = getattr(ds["var1"], func)(dim="x", skipna=True, min_count=3)
    assert_allclose(actual, expected)


@pytest.mark.parametrize("dtype", [float, int, np.float32, np.bool_])
@pytest.mark.parametrize("dask", [False, True])
@pytest.mark.parametrize("skipna", [False, True])
@pytest.mark.parametrize("func", ["sum", "prod"])
def test_multiple_dims(dtype, dask, skipna, func):
    if dask and not has_dask:
        pytest.skip("requires dask")
    da = construct_dataarray(3, dtype, contains_nan=True, dask=dask)

    actual = getattr(da, func)(("x", "y"), skipna=skipna)
    expected = getattr(getattr(da, func)("x", skipna=skipna), func)("y", skipna=skipna)
    assert_allclose(actual, expected)


@pytest.mark.parametrize("dask", [True, False])
def test_datetime_to_numeric_datetime64(dask, time_unit: PDDatetimeUnitOptions):
    if dask and not has_dask:
        pytest.skip("requires dask")

    times = pd.date_range("2000", periods=5, freq="7D").as_unit(time_unit).values
    if dask:
        import dask.array

        times = dask.array.from_array(times, chunks=-1)

    with raise_if_dask_computes():
        result = duck_array_ops.datetime_to_numeric(times, datetime_unit="h")
    expected = 24 * np.arange(0, 35, 7)
    np.testing.assert_array_equal(result, expected)

    offset = times[1]
    with raise_if_dask_computes():
        result = duck_array_ops.datetime_to_numeric(
            times, offset=offset, datetime_unit="h"
        )
    expected = 24 * np.arange(-7, 28, 7)
    np.testing.assert_array_equal(result, expected)

    dtype = np.float32
    with raise_if_dask_computes():
        result = duck_array_ops.datetime_to_numeric(
            times, datetime_unit="h", dtype=dtype
        )
    expected2 = 24 * np.arange(0, 35, 7).astype(dtype)
    np.testing.assert_array_equal(result, expected2)


@requires_cftime
@pytest.mark.parametrize("dask", [True, False])
def test_datetime_to_numeric_cftime(dask):
    if dask and not has_dask:
        pytest.skip("requires dask")

    times = date_range(
        "2000", periods=5, freq="7D", calendar="standard", use_cftime=True
    ).values
    if dask:
        import dask.array

        times = dask.array.from_array(times, chunks=-1)
    with raise_if_dask_computes():
        result = duck_array_ops.datetime_to_numeric(times, datetime_unit="h", dtype=int)
    expected = 24 * np.arange(0, 35, 7)
    np.testing.assert_array_equal(result, expected)

    offset = times[1]
    with raise_if_dask_computes():
        result = duck_array_ops.datetime_to_numeric(
            times, offset=offset, datetime_unit="h", dtype=int
        )
    expected = 24 * np.arange(-7, 28, 7)
    np.testing.assert_array_equal(result, expected)

    dtype = np.float32
    with raise_if_dask_computes():
        result = duck_array_ops.datetime_to_numeric(
            times, datetime_unit="h", dtype=dtype
        )
    expected = 24 * np.arange(0, 35, 7).astype(dtype)
    np.testing.assert_array_equal(result, expected)

    with raise_if_dask_computes():
        if dask:
            time = dask.array.asarray(times[1])
        else:
            time = np.asarray(times[1])
        result = duck_array_ops.datetime_to_numeric(
            time, offset=times[0], datetime_unit="h", dtype=int
        )
    expected = np.array(24 * 7).astype(int)
    np.testing.assert_array_equal(result, expected)


@requires_cftime
def test_datetime_to_numeric_potential_overflow(time_unit: PDDatetimeUnitOptions):
    import cftime

    if time_unit == "ns":
        pytest.skip("out-of-bounds datetime64 overflow")
    dtype = f"M8[{time_unit}]"
    times = pd.date_range("2000", periods=5, freq="7D").values.astype(dtype)
    cftimes = date_range(
        "2000", periods=5, freq="7D", calendar="proleptic_gregorian", use_cftime=True
    ).values

    offset = np.datetime64("0001-01-01", time_unit)
    cfoffset = cftime.DatetimeProlepticGregorian(1, 1, 1)

    result = duck_array_ops.datetime_to_numeric(
        times, offset=offset, datetime_unit="D", dtype=int
    )
    cfresult = duck_array_ops.datetime_to_numeric(
        cftimes, offset=cfoffset, datetime_unit="D", dtype=int
    )

    expected = 730119 + np.arange(0, 35, 7)

    np.testing.assert_array_equal(result, expected)
    np.testing.assert_array_equal(cfresult, expected)


def test_py_timedelta_to_float():
    assert py_timedelta_to_float(dt.timedelta(days=1), "ns") == 86400 * 1e9
    assert py_timedelta_to_float(dt.timedelta(days=1e6), "ps") == 86400 * 1e18
    assert py_timedelta_to_float(dt.timedelta(days=1e6), "ns") == 86400 * 1e15
    assert py_timedelta_to_float(dt.timedelta(days=1e6), "us") == 86400 * 1e12
    assert py_timedelta_to_float(dt.timedelta(days=1e6), "ms") == 86400 * 1e9
    assert py_timedelta_to_float(dt.timedelta(days=1e6), "s") == 86400 * 1e6
    assert py_timedelta_to_float(dt.timedelta(days=1e6), "D") == 1e6


@pytest.mark.parametrize("np_dt_unit", ["D", "h", "m", "s", "ms", "us", "ns"])
def test_np_timedelta64_to_float(
    np_dt_unit: NPDatetimeUnitOptions, time_unit: PDDatetimeUnitOptions
):
    # tests any combination of source np.timedelta64 (NPDatetimeUnitOptions) with
    # np_timedelta_to_float with dedicated target unit (PDDatetimeUnitOptions)
    td = np.timedelta64(1, np_dt_unit)
    expected = _NS_PER_TIME_DELTA[np_dt_unit] / _NS_PER_TIME_DELTA[time_unit]

    out = np_timedelta64_to_float(td, datetime_unit=time_unit)
    np.testing.assert_allclose(out, expected)
    assert isinstance(out, float)

    out = np_timedelta64_to_float(np.atleast_1d(td), datetime_unit=time_unit)
    np.testing.assert_allclose(out, expected)


@pytest.mark.parametrize("np_dt_unit", ["D", "h", "m", "s", "ms", "us", "ns"])
def test_pd_timedelta_to_float(
    np_dt_unit: NPDatetimeUnitOptions, time_unit: PDDatetimeUnitOptions
):
    # tests any combination of source pd.Timedelta (NPDatetimeUnitOptions) with
    # np_timedelta_to_float with dedicated target unit (PDDatetimeUnitOptions)
    td = pd.Timedelta(1, np_dt_unit)
    expected = _NS_PER_TIME_DELTA[np_dt_unit] / _NS_PER_TIME_DELTA[time_unit]

    out = pd_timedelta_to_float(td, datetime_unit=time_unit)
    np.testing.assert_allclose(out, expected)
    assert isinstance(out, float)


@pytest.mark.parametrize(
    "td", [dt.timedelta(days=1), np.timedelta64(1, "D"), pd.Timedelta(1, "D"), "1 day"]
)
def test_timedelta_to_numeric(td, time_unit: PDDatetimeUnitOptions):
    # Scalar input
    out = timedelta_to_numeric(td, time_unit)
    expected = _NS_PER_TIME_DELTA["D"] / _NS_PER_TIME_DELTA[time_unit]
    np.testing.assert_allclose(out, expected)
    assert isinstance(out, float)


@pytest.mark.parametrize("use_dask", [True, False])
@pytest.mark.parametrize("skipna", [True, False])
def test_least_squares(use_dask, skipna):
    if use_dask and (not has_dask or not has_scipy):
        pytest.skip("requires dask and scipy")
    lhs = np.array([[1, 2], [1, 2], [3, 2]])
    rhs = DataArray(np.array([3, 5, 7]), dims=("y",))

    if use_dask:
        rhs = rhs.chunk({"y": 1})

    coeffs, residuals = least_squares(lhs, rhs.data, skipna=skipna)

    np.testing.assert_allclose(coeffs, [1.5, 1.25])
    np.testing.assert_allclose(residuals, [2.0])


@requires_dask
@requires_bottleneck
@pytest.mark.parametrize("method", ["sequential", "blelloch"])
@pytest.mark.parametrize(
    "arr",
    [
        [np.nan, 1, 2, 3, np.nan, np.nan, np.nan, np.nan, 4, 5, np.nan, 6],
        [
            np.nan,
            np.nan,
            np.nan,
            2,
            np.nan,
            np.nan,
            np.nan,
            9,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        ],
    ],
)
def test_push_dask(method, arr):
    import bottleneck
    import dask.array as da

    arr = np.array(arr)
    chunks = list(range(1, 11)) + [(1, 2, 3, 2, 2, 1, 1)]

    for n in [None, 1, 2, 3, 4, 5, 11]:
        expected = bottleneck.push(arr, axis=0, n=n)
        for c in chunks:
            with raise_if_dask_computes():
                actual = push(da.from_array(arr, chunks=c), axis=0, n=n, method=method)
            np.testing.assert_equal(actual, expected)


def test_extension_array_equality(categorical1, int1):
    int_duck_array = PandasExtensionArray(int1)
    categorical_duck_array = PandasExtensionArray(categorical1)
    assert (int_duck_array != categorical_duck_array).all()
    assert (categorical_duck_array == categorical1).all()
    assert (int_duck_array[0:2] == int1[0:2]).all()


def test_extension_array_singleton_equality(categorical1):
    categorical_duck_array = PandasExtensionArray(categorical1)
    assert (categorical_duck_array != "cat3").all()


def test_extension_array_repr(int1):
    int_duck_array = PandasExtensionArray(int1)
    assert repr(int1) in repr(int_duck_array)


def test_extension_array_attr():
    array = pd.Categorical(["cat2", "cat1", "cat2", "cat3", "cat1"])
    wrapped = PandasExtensionArray(array)
    assert_array_equal(array.categories, wrapped.categories)
    assert array.nbytes == wrapped.nbytes

    roundtripped = pickle.loads(pickle.dumps(wrapped))
    assert isinstance(roundtripped, PandasExtensionArray)
    assert (roundtripped == wrapped).all()

    interval_array = pd.arrays.IntervalArray.from_breaks([0, 1, 2, 3], closed="right")
    wrapped = PandasExtensionArray(interval_array)
    assert_array_equal(wrapped.left, interval_array.left, strict=True)
    assert wrapped.closed == interval_array.closed
