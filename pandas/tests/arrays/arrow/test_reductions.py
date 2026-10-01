from __future__ import annotations

import pytest

from pandas.compat import (
    pa_version_under20p0,
)

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")


@pytest.mark.parametrize("skipna", [True, False])
def test_boolean_reduce_series_all_null(all_boolean_reductions, skipna):
    # GH51624
    ser = pd.Series([None], dtype="float64[pyarrow]")
    result = getattr(ser, all_boolean_reductions)(skipna=skipna)
    if skipna:
        expected = all_boolean_reductions == "all"
    else:
        expected = pd.NA
    assert result is expected


@pytest.mark.parametrize(
    "method",
    ["sum", "prod", "mean", "median", "std", "var", "sem", "skew", "min", "max"],
)
@pytest.mark.parametrize("axis", [0, None, -1])
def test_reduction_axis_valid(method, axis):
    # axis equivalent to 0 (the only axis of a 1-D array) or None is accepted
    if method == "skew" and pa_version_under20p0:
        pytest.skip("pyarrow.compute.skew added in pyarrow 20.0.0")
    arr = pd.array([1, 2, 3], dtype="int64[pyarrow]")
    result = getattr(arr, method)(axis=axis)
    expected = getattr(arr, method)()
    assert result == expected or (pd.isna(result) and pd.isna(expected))


@pytest.mark.parametrize(
    "method",
    ["sum", "prod", "mean", "median", "std", "var", "sem", "skew", "min", "max"],
)
@pytest.mark.parametrize("axis", [1, 2])
def test_reduction_axis_out_of_bounds(method, axis):
    # axis beyond the array's dimensions is rejected rather than silently ignored
    arr = pd.array([1, 2, 3], dtype="int64[pyarrow]")
    msg = "`axis` must be fewer than the number of dimensions"
    with pytest.raises(ValueError, match=msg):
        getattr(arr, method)(axis=axis)


@pytest.mark.parametrize(
    "interpolation", ["linear", "lower", "higher", "nearest", "midpoint"]
)
@pytest.mark.parametrize("quantile", [0.5, [0.5, 0.5]])
def test_quantile(data, interpolation, quantile, request):
    pa_dtype = data.dtype.pyarrow_dtype

    data = data.take([0, 0, 0])
    ser = pd.Series(data)

    if (
        pa.types.is_string(pa_dtype)
        or pa.types.is_binary(pa_dtype)
        or pa.types.is_boolean(pa_dtype)
    ):
        # For string, bytes, and bool, we don't *expect* to have quantile work
        # Note this matches the non-pyarrow behavior
        msg = r"Function 'quantile' has no kernel matching input types \(.*\)"
        with pytest.raises(pa.ArrowNotImplementedError, match=msg):
            ser.quantile(q=quantile, interpolation=interpolation)
        return

    if (
        pa.types.is_integer(pa_dtype)
        or pa.types.is_floating(pa_dtype)
        or pa.types.is_decimal(pa_dtype)
    ):
        pass
    elif pa.types.is_temporal(data._pa_array.type):
        pass
    else:
        request.applymarker(
            pytest.mark.xfail(
                raises=pa.ArrowNotImplementedError,
                reason=f"quantile not supported by pyarrow for {pa_dtype}",
            )
        )
    data = data.take([0, 0, 0])
    ser = pd.Series(data)
    result = ser.quantile(q=quantile, interpolation=interpolation)

    if pa.types.is_timestamp(pa_dtype) and interpolation not in ["lower", "higher"]:
        # rounding error will make the check below fail
        #  (e.g. '2020-01-01 01:01:01.000001' vs '2020-01-01 01:01:01.000001024'),
        #  so we'll check for now that we match the numpy analogue
        if pa_dtype.tz:
            pd_dtype = f"M8[{pa_dtype.unit}, {pa_dtype.tz}]"
        else:
            pd_dtype = f"M8[{pa_dtype.unit}]"
        ser_np = ser.astype(pd_dtype)

        expected = ser_np.quantile(q=quantile, interpolation=interpolation)
        if quantile == 0.5:
            if pa_dtype.unit == "us":
                expected = expected.to_pydatetime(warn=False)
            assert result == expected
        else:
            if pa_dtype.unit == "us":
                expected = expected.dt.floor("us")
            tm.assert_series_equal(result, expected.astype(data.dtype))
        return

    if quantile == 0.5:
        assert result == data[0]
    else:
        # Just check the values
        expected = pd.Series(data.take([0, 0]), index=[0.5, 0.5])
        if (
            pa.types.is_integer(pa_dtype)
            or pa.types.is_floating(pa_dtype)
            or pa.types.is_decimal(pa_dtype)
        ):
            expected = expected.astype("float64[pyarrow]")
            result = result.astype("float64[pyarrow]")
        tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "take_idx, exp_idx",
    [[[0, 0, 2, 2, 4, 4], [4, 0]], [[0, 0, 0, 2, 4, 4], [0]]],
    ids=["multi_mode", "single_mode"],
)
def test_mode_dropna_true(data_for_grouping, take_idx, exp_idx):
    data = data_for_grouping.take(take_idx)
    ser = pd.Series(data)
    result = ser.mode(dropna=True)
    expected = pd.Series(data_for_grouping.take(exp_idx))
    tm.assert_series_equal(result, expected)


def test_mode_dropna_false_mode_na(data):
    # GH 50982
    more_nans = pd.Series([None, None, data[0]], dtype=data.dtype)
    result = more_nans.mode(dropna=False)
    expected = pd.Series([None], dtype=data.dtype)
    tm.assert_series_equal(result, expected)

    expected = pd.Series([data[0], None], dtype=data.dtype)
    result = expected.mode(dropna=False)
    tm.assert_series_equal(result, expected)
