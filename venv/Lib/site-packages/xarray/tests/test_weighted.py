from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import numpy as np
import pytest

import xarray as xr
from xarray import DataArray, Dataset
from xarray.tests import (
    assert_allclose,
    assert_equal,
    raise_if_dask_computes,
    requires_cftime,
    requires_dask,
)


@pytest.mark.parametrize("as_dataset", (True, False))
def test_weighted_non_DataArray_weights(as_dataset: bool) -> None:
    data: DataArray | Dataset = DataArray([1, 2])
    if as_dataset:
        data = data.to_dataset(name="data")

    with pytest.raises(ValueError, match=r"`weights` must be a DataArray"):
        data.weighted([1, 2])  # type: ignore[arg-type]


@pytest.mark.parametrize("as_dataset", (True, False))
@pytest.mark.parametrize("weights", ([np.nan, 2], [np.nan, np.nan]))
def test_weighted_weights_nan_raises(as_dataset: bool, weights: list[float]) -> None:
    data: DataArray | Dataset = DataArray([1, 2])
    if as_dataset:
        data = data.to_dataset(name="data")

    with pytest.raises(ValueError, match="`weights` cannot contain missing values."):
        data.weighted(DataArray(weights))


@requires_dask
@pytest.mark.parametrize("as_dataset", (True, False))
@pytest.mark.parametrize("weights", ([np.nan, 2], [np.nan, np.nan]))
def test_weighted_weights_nan_raises_dask(as_dataset, weights):
    data = DataArray([1, 2]).chunk({"dim_0": -1})
    if as_dataset:
        data = data.to_dataset(name="data")

    weights = DataArray(weights).chunk({"dim_0": -1})

    with raise_if_dask_computes():
        weighted = data.weighted(weights)

    with pytest.raises(ValueError, match="`weights` cannot contain missing values."):
        weighted.sum().load()


@requires_cftime
@requires_dask
@pytest.mark.parametrize("time_chunks", (1, 5))
@pytest.mark.parametrize("resample_spec", ("1YS", "5YS", "10YS"))
def test_weighted_lazy_resample(time_chunks, resample_spec):
    # https://github.com/pydata/xarray/issues/4625

    # simple customized weighted mean function
    def mean_func(ds):
        return ds.weighted(ds.weights).mean("time")

    # example dataset
    t = xr.cftime_range(start="2000", periods=20, freq="1YS")
    weights = xr.DataArray(np.random.rand(len(t)), dims=["time"], coords={"time": t})
    data = xr.DataArray(
        np.random.rand(len(t)), dims=["time"], coords={"time": t, "weights": weights}
    )
    ds = xr.Dataset({"data": data}).chunk({"time": time_chunks})

    with raise_if_dask_computes():
        ds.resample(time=resample_spec).map(mean_func)


@pytest.mark.parametrize(
    ("weights", "expected"),
    (([1, 2], 3), ([2, 0], 2), ([0, 0], np.nan), ([-1, 1], np.nan)),
)
def test_weighted_sum_of_weights_no_nan(weights, expected):
    da = DataArray([1, 2])
    weights = DataArray(weights)
    result = da.weighted(weights).sum_of_weights()

    expected = DataArray(expected)

    assert_equal(expected, result)


@pytest.mark.parametrize(
    ("weights", "expected"),
    (([1, 2], 2), ([2, 0], np.nan), ([0, 0], np.nan), ([-1, 1], 1)),
)
def test_weighted_sum_of_weights_nan(weights, expected):
    da = DataArray([np.nan, 2])
    weights = DataArray(weights)
    result = da.weighted(weights).sum_of_weights()

    expected = DataArray(expected)

    assert_equal(expected, result)


def test_weighted_sum_of_weights_bool():
    # https://github.com/pydata/xarray/issues/4074

    da = DataArray([1, 2])
    weights = DataArray([True, True])
    result = da.weighted(weights).sum_of_weights()

    expected = DataArray(2)

    assert_equal(expected, result)


@pytest.mark.parametrize("da", ([1.0, 2], [1, np.nan], [np.nan, np.nan]))
@pytest.mark.parametrize("factor", [0, 1, 3.14])
@pytest.mark.parametrize("skipna", (True, False))
def test_weighted_sum_equal_weights(da, factor, skipna):
    # if all weights are 'f'; weighted sum is f times the ordinary sum

    da = DataArray(da)
    weights = xr.full_like(da, factor)

    expected = da.sum(skipna=skipna) * factor
    result = da.weighted(weights).sum(skipna=skipna)

    assert_equal(expected, result)


@pytest.mark.parametrize(
    ("weights", "expected"), (([1, 2], 5), ([0, 2], 4), ([0, 0], 0))
)
def test_weighted_sum_no_nan(weights, expected):
    da = DataArray([1, 2])

    weights = DataArray(weights)
    result = da.weighted(weights).sum()
    expected = DataArray(expected)

    assert_equal(expected, result)


@pytest.mark.parametrize(
    ("weights", "expected"), (([1, 2], 4), ([0, 2], 4), ([1, 0], 0), ([0, 0], 0))
)
@pytest.mark.parametrize("skipna", (True, False))
def test_weighted_sum_nan(weights, expected, skipna):
    da = DataArray([np.nan, 2])

    weights = DataArray(weights)
    result = da.weighted(weights).sum(skipna=skipna)

    if skipna:
        expected = DataArray(expected)
    else:
        expected = DataArray(np.nan)

    assert_equal(expected, result)


@pytest.mark.filterwarnings("error")
@pytest.mark.parametrize("da", ([1.0, 2], [1, np.nan], [np.nan, np.nan]))
@pytest.mark.parametrize("skipna", (True, False))
@pytest.mark.parametrize("factor", [1, 2, 3.14])
def test_weighted_mean_equal_weights(da, skipna, factor):
    # if all weights are equal (!= 0), should yield the same result as mean

    da = DataArray(da)

    # all weights as 1.
    weights = xr.full_like(da, factor)

    expected = da.mean(skipna=skipna)
    result = da.weighted(weights).mean(skipna=skipna)

    assert_equal(expected, result)


@pytest.mark.parametrize(
    ("weights", "expected"), (([4, 6], 1.6), ([1, 0], 1.0), ([0, 0], np.nan))
)
def test_weighted_mean_no_nan(weights, expected):
    da = DataArray([1, 2])
    weights = DataArray(weights)
    expected = DataArray(expected)

    result = da.weighted(weights).mean()

    assert_equal(expected, result)


@pytest.mark.parametrize(
    ("weights", "expected"),
    (
        (
            [0.25, 0.05, 0.15, 0.25, 0.15, 0.1, 0.05],
            [1.554595, 2.463784, 3.000000, 3.518378],
        ),
        (
            [0.05, 0.05, 0.1, 0.15, 0.15, 0.25, 0.25],
            [2.840000, 3.632973, 4.076216, 4.523243],
        ),
    ),
)
def test_weighted_quantile_no_nan(weights, expected):
    # Expected values were calculated by running the reference implementation
    # proposed in https://aakinshin.net/posts/weighted-quantiles/

    da = DataArray([1, 1.9, 2.2, 3, 3.7, 4.1, 5])
    q = [0.2, 0.4, 0.6, 0.8]
    weights = DataArray(weights)

    expected = DataArray(expected, coords={"quantile": q})
    result = da.weighted(weights).quantile(q)

    assert_allclose(expected, result)


def test_weighted_quantile_zero_weights():
    da = DataArray([0, 1, 2, 3])
    weights = DataArray([1, 0, 1, 0])
    q = 0.75

    result = da.weighted(weights).quantile(q)
    expected = DataArray([0, 2]).quantile(0.75)

    assert_allclose(expected, result)


def test_weighted_quantile_simple():
    # Check that weighted quantiles return the same value as numpy quantiles
    da = DataArray([0, 1, 2, 3])
    w = DataArray([1, 0, 1, 0])

    w_eps = DataArray([1, 0.0001, 1, 0.0001])
    q = 0.75

    expected = DataArray(np.quantile([0, 2], q), coords={"quantile": q})  # 1.5

    assert_equal(expected, da.weighted(w).quantile(q))
    assert_allclose(expected, da.weighted(w_eps).quantile(q), rtol=0.001)


@pytest.mark.parametrize("skipna", (True, False))
def test_weighted_quantile_nan(skipna):
    # Check skipna behavior
    da = DataArray([0, 1, 2, 3, np.nan])
    w = DataArray([1, 0, 1, 0, 1])
    q = [0.5, 0.75]

    result = da.weighted(w).quantile(q, skipna=skipna)

    if skipna:
        expected = DataArray(np.quantile([0, 2], q), coords={"quantile": q})
    else:
        expected = DataArray(np.full(len(q), np.nan), coords={"quantile": q})

    assert_allclose(expected, result)


@pytest.mark.parametrize(
    "da",
    (
        pytest.param([1, 1.9, 2.2, 3, 3.7, 4.1, 5], id="nonan"),
        pytest.param([1, 1.9, 2.2, 3, 3.7, 4.1, np.nan], id="singlenan"),
        pytest.param(
            [np.nan, np.nan, np.nan],
            id="allnan",
            marks=pytest.mark.filterwarnings(
                "ignore:All-NaN slice encountered:RuntimeWarning"
            ),
        ),
    ),
)
@pytest.mark.parametrize("q", (0.5, (0.2, 0.8)))
@pytest.mark.parametrize("skipna", (True, False))
@pytest.mark.parametrize("factor", [1, 3.14])
def test_weighted_quantile_equal_weights(
    da: list[float], q: float | tuple[float, ...], skipna: bool, factor: float
) -> None:
    # if all weights are equal (!= 0), should yield the same result as quantile

    data = DataArray(da)
    weights = xr.full_like(data, factor)

    expected = data.quantile(q, skipna=skipna)
    result = data.weighted(weights).quantile(q, skipna=skipna)

    assert_allclose(expected, result)


@pytest.mark.skip(reason="`method` argument is not currently exposed")
@pytest.mark.parametrize(
    "da",
    (
        [1, 1.9, 2.2, 3, 3.7, 4.1, 5],
        [1, 1.9, 2.2, 3, 3.7, 4.1, np.nan],
        [np.nan, np.nan, np.nan],
    ),
)
@pytest.mark.parametrize("q", (0.5, (0.2, 0.8)))
@pytest.mark.parametrize("skipna", (True, False))
@pytest.mark.parametrize(
    "method",
    [
        "linear",
        "interpolated_inverted_cdf",
        "hazen",
        "weibull",
        "median_unbiased",
        "normal_unbiased2",
    ],
)
def test_weighted_quantile_equal_weights_all_methods(da, q, skipna, factor, method):
    # If all weights are equal (!= 0), should yield the same result as numpy quantile

    da = DataArray(da)
    weights = xr.full_like(da, 3.14)

    expected = da.quantile(q, skipna=skipna, method=method)
    result = da.weighted(weights).quantile(q, skipna=skipna, method=method)

    assert_allclose(expected, result)


def test_weighted_quantile_bool():
    # https://github.com/pydata/xarray/issues/4074
    da = DataArray([1, 1])
    weights = DataArray([True, True])
    q = 0.5

    expected = DataArray([1], coords={"quantile": [q]}).squeeze()
    result = da.weighted(weights).quantile(q)

    assert_equal(expected, result)


@pytest.mark.parametrize("q", (-1, 1.1, (0.5, 1.1), ((0.2, 0.4), (0.6, 0.8))))
def test_weighted_quantile_with_invalid_q(q):
    da = DataArray([1, 1.9, 2.2, 3, 3.7, 4.1, 5])
    q = np.asarray(q)
    weights = xr.ones_like(da)

    if q.ndim <= 1:
        with pytest.raises(ValueError, match="q values must be between 0 and 1"):
            da.weighted(weights).quantile(q)
    else:
        with pytest.raises(ValueError, match="q must be a scalar or 1d"):
            da.weighted(weights).quantile(q)


@pytest.mark.parametrize(
    ("weights", "expected"), (([4, 6], 2.0), ([1, 0], np.nan), ([0, 0], np.nan))
)
@pytest.mark.parametrize("skipna", (True, False))
def test_weighted_mean_nan(weights, expected, skipna):
    da = DataArray([np.nan, 2])
    weights = DataArray(weights)

    if skipna:
        expected = DataArray(expected)
    else:
        expected = DataArray(np.nan)

    result = da.weighted(weights).mean(skipna=skipna)

    assert_equal(expected, result)


def test_weighted_mean_bool():
    # https://github.com/pydata/xarray/issues/4074
    da = DataArray([1, 1])
    weights = DataArray([True, True])
    expected = DataArray(1)

    result = da.weighted(weights).mean()

    assert_equal(expected, result)


@pytest.mark.parametrize(
    ("weights", "expected"),
    (([1, 2], 2 / 3), ([2, 0], 0), ([0, 0], 0), ([-1, 1], 0)),
)
def test_weighted_sum_of_squares_no_nan(weights, expected):
    da = DataArray([1, 2])
    weights = DataArray(weights)
    result = da.weighted(weights).sum_of_squares()

    expected = DataArray(expected)

    assert_equal(expected, result)


@pytest.mark.parametrize(
    ("weights", "expected"),
    (([1, 2], 0), ([2, 0], 0), ([0, 0], 0), ([-1, 1], 0)),
)
def test_weighted_sum_of_squares_nan(weights, expected):
    da = DataArray([np.nan, 2])
    weights = DataArray(weights)
    result = da.weighted(weights).sum_of_squares()

    expected = DataArray(expected)

    assert_equal(expected, result)


@pytest.mark.filterwarnings("error")
@pytest.mark.parametrize("da", ([1.0, 2], [1, np.nan]))
@pytest.mark.parametrize("skipna", (True, False))
@pytest.mark.parametrize("factor", [1, 2, 3.14])
def test_weighted_var_equal_weights(da, skipna, factor):
    # if all weights are equal (!= 0), should yield the same result as var

    da = DataArray(da)

    # all weights as 1.
    weights = xr.full_like(da, factor)

    expected = da.var(skipna=skipna)
    result = da.weighted(weights).var(skipna=skipna)

    assert_equal(expected, result)


@pytest.mark.parametrize(
    ("weights", "expected"), (([4, 6], 0.24), ([1, 0], 0.0), ([0, 0], np.nan))
)
def test_weighted_var_no_nan(weights, expected):
    da = DataArray([1, 2])
    weights = DataArray(weights)
    expected = DataArray(expected)

    result = da.weighted(weights).var()

    assert_equal(expected, result)


@pytest.mark.parametrize(
    ("weights", "expected"), (([4, 6], 0), ([1, 0], np.nan), ([0, 0], np.nan))
)
def test_weighted_var_nan(weights, expected):
    da = DataArray([np.nan, 2])
    weights = DataArray(weights)
    expected = DataArray(expected)

    result = da.weighted(weights).var()

    assert_equal(expected, result)


def test_weighted_var_bool():
    # https://github.com/pydata/xarray/issues/4074
    da = DataArray([1, 1])
    weights = DataArray([True, True])
    expected = DataArray(0)

    result = da.weighted(weights).var()

    assert_equal(expected, result)


@pytest.mark.filterwarnings("error")
@pytest.mark.parametrize("da", ([1.0, 2], [1, np.nan]))
@pytest.mark.parametrize("skipna", (True, False))
@pytest.mark.parametrize("factor", [1, 2, 3.14])
def test_weighted_std_equal_weights(da, skipna, factor):
    # if all weights are equal (!= 0), should yield the same result as std

    da = DataArray(da)

    # all weights as 1.
    weights = xr.full_like(da, factor)

    expected = da.std(skipna=skipna)
    result = da.weighted(weights).std(skipna=skipna)

    assert_equal(expected, result)


@pytest.mark.parametrize(
    ("weights", "expected"), (([4, 6], np.sqrt(0.24)), ([1, 0], 0.0), ([0, 0], np.nan))
)
def test_weighted_std_no_nan(weights, expected):
    da = DataArray([1, 2])
    weights = DataArray(weights)
    expected = DataArray(expected)

    result = da.weighted(weights).std()

    assert_equal(expected, result)


@pytest.mark.parametrize(
    ("weights", "expected"), (([4, 6], 0), ([1, 0], np.nan), ([0, 0], np.nan))
)
def test_weighted_std_nan(weights, expected):
    da = DataArray([np.nan, 2])
    weights = DataArray(weights)
    expected = DataArray(expected)

    result = da.weighted(weights).std()

    assert_equal(expected, result)


def test_weighted_std_bool():
    # https://github.com/pydata/xarray/issues/4074
    da = DataArray([1, 1])
    weights = DataArray([True, True])
    expected = DataArray(0)

    result = da.weighted(weights).std()

    assert_equal(expected, result)


def expected_weighted(da, weights, dim, skipna, operation):
    """
    Generate expected result using ``*`` and ``sum``. This is checked against
    the result of da.weighted which uses ``dot``
    """

    weighted_sum = (da * weights).sum(dim=dim, skipna=skipna)

    if operation == "sum":
        return weighted_sum

    masked_weights = weights.where(da.notnull())
    sum_of_weights = masked_weights.sum(dim=dim, skipna=True)
    valid_weights = sum_of_weights != 0
    sum_of_weights = sum_of_weights.where(valid_weights)

    if operation == "sum_of_weights":
        return sum_of_weights

    weighted_mean = weighted_sum / sum_of_weights

    if operation == "mean":
        return weighted_mean

    demeaned = da - weighted_mean
    sum_of_squares = ((demeaned**2) * weights).sum(dim=dim, skipna=skipna)

    if operation == "sum_of_squares":
        return sum_of_squares

    var = sum_of_squares / sum_of_weights

    if operation == "var":
        return var

    if operation == "std":
        return np.sqrt(var)


def check_weighted_operations(data, weights, dim, skipna):
    # check sum of weights
    result = data.weighted(weights).sum_of_weights(dim)
    expected = expected_weighted(data, weights, dim, skipna, "sum_of_weights")
    assert_allclose(expected, result)

    # check weighted sum
    result = data.weighted(weights).sum(dim, skipna=skipna)
    expected = expected_weighted(data, weights, dim, skipna, "sum")
    assert_allclose(expected, result)

    # check weighted mean
    result = data.weighted(weights).mean(dim, skipna=skipna)
    expected = expected_weighted(data, weights, dim, skipna, "mean")
    assert_allclose(expected, result)

    # check weighted sum of squares
    result = data.weighted(weights).sum_of_squares(dim, skipna=skipna)
    expected = expected_weighted(data, weights, dim, skipna, "sum_of_squares")
    assert_allclose(expected, result)

    # check weighted var
    result = data.weighted(weights).var(dim, skipna=skipna)
    expected = expected_weighted(data, weights, dim, skipna, "var")
    assert_allclose(expected, result)

    # check weighted std
    result = data.weighted(weights).std(dim, skipna=skipna)
    expected = expected_weighted(data, weights, dim, skipna, "std")
    assert_allclose(expected, result)


@pytest.mark.parametrize("dim", ("a", "b", "c", ("a", "b"), ("a", "b", "c"), None))
@pytest.mark.parametrize("add_nans", (True, False))
@pytest.mark.parametrize("skipna", (None, True, False))
@pytest.mark.filterwarnings("ignore:invalid value encountered in sqrt")
def test_weighted_operations_3D(dim, add_nans, skipna):
    dims = ("a", "b", "c")
    coords = dict(a=[0, 1, 2, 3], b=[0, 1, 2, 3], c=[0, 1, 2, 3])

    weights = DataArray(np.random.randn(4, 4, 4), dims=dims, coords=coords)

    data = np.random.randn(4, 4, 4)

    # add approximately 25 % NaNs (https://stackoverflow.com/a/32182680/3010700)
    if add_nans:
        c = int(data.size * 0.25)
        data.ravel()[np.random.choice(data.size, c, replace=False)] = np.nan

    data = DataArray(data, dims=dims, coords=coords)

    check_weighted_operations(data, weights, dim, skipna)

    data = data.to_dataset(name="data")
    check_weighted_operations(data, weights, dim, skipna)


@pytest.mark.parametrize("dim", ("a", "b", "c", ("a", "b"), ("a", "b", "c"), None))
@pytest.mark.parametrize("q", (0.5, (0.1, 0.9), (0.2, 0.4, 0.6, 0.8)))
@pytest.mark.parametrize("add_nans", (True, False))
@pytest.mark.parametrize("skipna", (None, True, False))
def test_weighted_quantile_3D(dim, q, add_nans, skipna):
    dims = ("a", "b", "c")
    coords = dict(a=[0, 1, 2], b=[0, 1, 2, 3], c=[0, 1, 2, 3, 4])

    data = np.arange(60).reshape(3, 4, 5).astype(float)

    # add approximately 25 % NaNs (https://stackoverflow.com/a/32182680/3010700)
    if add_nans:
        c = int(data.size * 0.25)
        data.ravel()[np.random.choice(data.size, c, replace=False)] = np.nan

    da = DataArray(data, dims=dims, coords=coords)

    # Weights are all ones, because we will compare against DataArray.quantile (non-weighted)
    weights = xr.ones_like(da)

    result = da.weighted(weights).quantile(q, dim=dim, skipna=skipna)
    expected = da.quantile(q, dim=dim, skipna=skipna)

    assert_allclose(expected, result)

    ds = da.to_dataset(name="data")
    result2 = ds.weighted(weights).quantile(q, dim=dim, skipna=skipna)

    assert_allclose(expected, result2.data)


@pytest.mark.parametrize(
    "coords_weights, coords_data, expected_value_at_weighted_quantile",
    [
        ([0, 1, 2, 3], [1, 2, 3, 4], 2.5),  # no weights for coord a == 4
        ([0, 1, 2, 3], [2, 3, 4, 5], 1.8),  # no weights for coord a == 4 or 5
        ([2, 3, 4, 5], [0, 1, 2, 3], 3.8),  # no weights for coord a == 0 or 1
    ],
)
def test_weighted_operations_nonequal_coords(
    coords_weights: Iterable[Any],
    coords_data: Iterable[Any],
    expected_value_at_weighted_quantile: float,
) -> None:
    """Check that weighted operations work with unequal coords.


    Parameters
    ----------
    coords_weights : Iterable[Any]
        The coords for the weights.
    coords_data : Iterable[Any]
        The coords for the data.
    expected_value_at_weighted_quantile : float
        The expected value for the quantile of the weighted data.
    """
    da_weights = DataArray(
        [0.5, 1.0, 1.0, 2.0], dims=("a",), coords=dict(a=coords_weights)
    )
    da_data = DataArray([1, 2, 3, 4], dims=("a",), coords=dict(a=coords_data))
    check_weighted_operations(da_data, da_weights, dim="a", skipna=None)

    quantile = 0.5
    da_actual = da_data.weighted(da_weights).quantile(quantile, dim="a")
    da_expected = DataArray(
        [expected_value_at_weighted_quantile], coords={"quantile": [quantile]}
    ).squeeze()
    assert_allclose(da_actual, da_expected)

    ds_data = da_data.to_dataset(name="data")
    check_weighted_operations(ds_data, da_weights, dim="a", skipna=None)

    ds_actual = ds_data.weighted(da_weights).quantile(quantile, dim="a")
    assert_allclose(ds_actual, da_expected.to_dataset(name="data"))


@pytest.mark.parametrize("shape_data", ((4,), (4, 4), (4, 4, 4)))
@pytest.mark.parametrize("shape_weights", ((4,), (4, 4), (4, 4, 4)))
@pytest.mark.parametrize("add_nans", (True, False))
@pytest.mark.parametrize("skipna", (None, True, False))
@pytest.mark.filterwarnings("ignore:invalid value encountered in sqrt")
def test_weighted_operations_different_shapes(
    shape_data, shape_weights, add_nans, skipna
):
    weights = DataArray(np.random.randn(*shape_weights))

    data = np.random.randn(*shape_data)

    # add approximately 25 % NaNs
    if add_nans:
        c = int(data.size * 0.25)
        data.ravel()[np.random.choice(data.size, c, replace=False)] = np.nan

    data = DataArray(data)

    check_weighted_operations(data, weights, "dim_0", skipna)
    check_weighted_operations(data, weights, None, skipna)

    data = data.to_dataset(name="data")
    check_weighted_operations(data, weights, "dim_0", skipna)
    check_weighted_operations(data, weights, None, skipna)


@pytest.mark.parametrize(
    "operation",
    ("sum_of_weights", "sum", "mean", "sum_of_squares", "var", "std", "quantile"),
)
@pytest.mark.parametrize("as_dataset", (True, False))
@pytest.mark.parametrize("keep_attrs", (True, False, None))
def test_weighted_operations_keep_attr(operation, as_dataset, keep_attrs):
    weights = DataArray(np.random.randn(2, 2), attrs=dict(attr="weights"))
    data = DataArray(np.random.randn(2, 2))

    if as_dataset:
        data = data.to_dataset(name="data")

    data.attrs = dict(attr="weights")

    kwargs = {"keep_attrs": keep_attrs}
    if operation == "quantile":
        kwargs["q"] = 0.5

    result = getattr(data.weighted(weights), operation)(**kwargs)

    if operation == "sum_of_weights":
        assert result.attrs == (weights.attrs if keep_attrs else {})
        assert result.attrs == (weights.attrs if keep_attrs else {})
    else:
        assert result.attrs == (weights.attrs if keep_attrs else {})
        assert result.attrs == (data.attrs if keep_attrs else {})


@pytest.mark.parametrize(
    "operation",
    ("sum_of_weights", "sum", "mean", "sum_of_squares", "var", "std", "quantile"),
)
def test_weighted_operations_keep_attr_da_in_ds(operation):
    # GH #3595

    weights = DataArray(np.random.randn(2, 2))
    data = DataArray(np.random.randn(2, 2), attrs=dict(attr="data"))
    data = data.to_dataset(name="a")

    kwargs = {"keep_attrs": True}
    if operation == "quantile":
        kwargs["q"] = 0.5

    result = getattr(data.weighted(weights), operation)(**kwargs)

    assert data.a.attrs == result.a.attrs


@pytest.mark.parametrize("operation", ("sum_of_weights", "sum", "mean", "quantile"))
@pytest.mark.parametrize("as_dataset", (True, False))
def test_weighted_bad_dim(operation, as_dataset):
    data = DataArray(np.random.randn(2, 2))
    weights = xr.ones_like(data)
    if as_dataset:
        data = data.to_dataset(name="data")

    kwargs = {"dim": "bad_dim"}
    if operation == "quantile":
        kwargs["q"] = 0.5

    with pytest.raises(
        ValueError,
        match=(
            f"Dimensions \\('bad_dim',\\) not found in {data.__class__.__name__}Weighted "
            # the order of (dim_0, dim_1) varies
            "dimensions \\(('dim_0', 'dim_1'|'dim_1', 'dim_0')\\)"
        ),
    ):
        getattr(data.weighted(weights), operation)(**kwargs)
