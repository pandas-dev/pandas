import re

import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm
from pandas.core import missing


class TestDataFrameInterpolate:
    def test_interpolate_complex(self):
        # GH#53635
        ser = pd.Series([complex("1+1j"), float("nan"), complex("2+2j")])
        assert ser.dtype.kind == "c"

        res = ser.interpolate()
        expected = pd.Series([ser[0], ser[0] * 1.5, ser[2]])
        tm.assert_series_equal(res, expected)

        df = ser.to_frame()
        res = df.interpolate()
        expected = expected.to_frame()
        tm.assert_frame_equal(res, expected)

    def test_interpolate_datetimelike_values(self, frame_or_series):
        # GH#11312, GH#51005
        orig = pd.Series(pd.date_range("2012-01-01", periods=5))
        ser = orig.copy()
        ser[2] = pd.NaT

        res = frame_or_series(ser).interpolate()
        expected = frame_or_series(orig)
        tm.assert_equal(res, expected)

        # datetime64tz cast
        ser_tz = ser.dt.tz_localize("US/Pacific")
        res_tz = frame_or_series(ser_tz).interpolate()
        expected_tz = frame_or_series(orig.dt.tz_localize("US/Pacific"))
        tm.assert_equal(res_tz, expected_tz)

        # timedelta64 cast
        ser_td = ser - ser[0]
        res_td = frame_or_series(ser_td).interpolate()
        expected_td = frame_or_series(orig - orig[0])
        tm.assert_equal(res_td, expected_td)

    def test_interpolate_inplace(self, frame_or_series, request):
        # GH#44749
        obj = frame_or_series([1, np.nan, 2])
        orig = obj.values

        result = obj.interpolate(inplace=True)
        assert result is obj
        expected = frame_or_series([1, 1.5, 2])
        tm.assert_equal(obj, expected)

        # check we operated *actually* inplace
        assert np.shares_memory(orig, obj.values)
        assert orig.squeeze()[1] == 1.5

    @pytest.mark.parametrize("dtype", ["str", "object"])
    def test_interp_with_non_numeric(self, dtype, using_infer_string):
        if not using_infer_string:
            dtype = "object"
        df = pd.DataFrame(
            {
                "A": [1, 2, np.nan, 4],
                "B": [1, 4, 9, np.nan],
                "C": [1, 2, 3, 5],
                "D": list("abcd"),
            }
        )
        df["D"] = df["D"].astype(dtype)
        df_orig = df.copy()
        expected = pd.DataFrame(
            {
                "A": [1.0, 2.0, 3.0, 4.0],
                "B": [1.0, 4.0, 9.0, 9.0],
                "C": [1, 2, 3, 5],
                "D": list("abcd"),
            }
        )
        expected["D"] = expected["D"].astype(dtype)

        # the EA path (str dtype) and the object-block path word the error
        # differently; pin each to its dtype rather than hedging the case
        if dtype == "object":
            msg = "cannot interpolate with object dtype"
        else:
            msg = f"Cannot interpolate with {dtype} dtype"
        with pytest.raises(TypeError, match=msg):
            df.interpolate()
        tm.assert_frame_equal(df, df_orig)

        with pytest.raises(TypeError, match=msg):
            df.interpolate(inplace=True)
        # columns A and B already get interpolated before we hit the error
        tm.assert_frame_equal(df, expected)

    def test_interp_basic(self):
        df = pd.DataFrame(
            {
                "A": [1, 2, np.nan, 4],
                "B": [1, 4, 9, np.nan],
                "C": [1, 2, 3, 5],
            }
        )
        df_orig = df.copy()
        expected = pd.DataFrame(
            {
                "A": [1.0, 2.0, 3.0, 4.0],
                "B": [1.0, 4.0, 9.0, 9.0],
                "C": [1, 2, 3, 5],
            }
        )
        result = df.interpolate()
        tm.assert_frame_equal(result, expected)

        # check we didn't operate inplace GH#45791
        tm.assert_frame_equal(df, df_orig)
        bvalues = df["B"]._values
        cvalues = df["C"]._values
        assert not tm.shares_memory(bvalues, result["B"]._values)
        assert tm.shares_memory(cvalues, result["C"]._values)

        res = df.interpolate(inplace=True)
        assert res is df
        tm.assert_frame_equal(df, expected)

        # check we DID operate inplace
        assert tm.shares_memory(df["B"]._values, bvalues)
        assert tm.shares_memory(df["C"]._values, cvalues)

    def test_interp_basic_with_non_range_index(self):
        df = pd.DataFrame(
            {
                "A": [1, 2, np.nan, 4],
                "B": [1, 4, 9, np.nan],
                "C": [1, 2, 3, 5],
            }
        )
        result = df.set_index("C").interpolate()
        expected = df.set_index("C")
        expected.loc[3, "A"] = 3
        expected.loc[5, "B"] = 9
        tm.assert_frame_equal(result, expected)

    def test_interp_empty(self):
        # https://github.com/pandas-dev/pandas/issues/35598
        df = pd.DataFrame()
        result = df.interpolate()
        assert result is not df
        expected = df
        tm.assert_frame_equal(result, expected)

    def test_interp_bad_method(self):
        df = pd.DataFrame(
            {
                "A": [1, 2, np.nan, 4],
                "B": [1, 4, 9, np.nan],
                "C": [1, 2, 3, 5],
            }
        )
        msg = "Can not interpolate with method=not_a_method"
        with pytest.raises(ValueError, match=msg):
            df.interpolate(method="not_a_method")

    def test_interp_combo(self):
        df = pd.DataFrame(
            {
                "A": [1.0, 2.0, np.nan, 4.0],
                "B": [1, 4, 9, np.nan],
                "C": [1, 2, 3, 5],
                "D": list("abcd"),
            }
        )

        result = df["A"].interpolate()
        expected = pd.Series([1.0, 2.0, 3.0, 4.0], name="A")
        tm.assert_series_equal(result, expected)

    def test_interp_nan_idx(self):
        df = pd.DataFrame({"A": [1, 2, np.nan, 4], "B": [np.nan, 2, 3, 4]})
        df = df.set_index("A")
        msg = (
            "Interpolation with NaNs in the index has not been implemented. "
            "Try filling those NaNs before interpolating."
        )
        with pytest.raises(NotImplementedError, match=msg):
            df.interpolate(method="values")

    def test_interp_various(self):
        pytest.importorskip("scipy")
        df = pd.DataFrame(
            {"A": [1, 2, np.nan, 4, 5, np.nan, 7], "C": [1, 2, 3, 5, 8, 13, 21]}
        )
        df = df.set_index("C")
        expected = df.copy()
        result = df.interpolate(method="polynomial", order=1)

        expected.loc[3, "A"] = 2.66666667
        expected.loc[13, "A"] = 5.76923076
        tm.assert_frame_equal(result, expected)

        result = df.interpolate(method="cubic")
        # GH #15662.
        expected.loc[3, "A"] = 2.81547781
        expected.loc[13, "A"] = 5.52964175
        tm.assert_frame_equal(result, expected)

        result = df.interpolate(method="nearest")
        expected.loc[3, "A"] = 2
        expected.loc[13, "A"] = 5
        tm.assert_frame_equal(result, expected, check_dtype=False)

        result = df.interpolate(method="quadratic")
        expected.loc[3, "A"] = 2.82150771
        expected.loc[13, "A"] = 6.12648668
        tm.assert_frame_equal(result, expected)

        result = df.interpolate(method="slinear")
        expected.loc[3, "A"] = 2.66666667
        expected.loc[13, "A"] = 5.76923077
        tm.assert_frame_equal(result, expected)

        result = df.interpolate(method="zero")
        expected.loc[3, "A"] = 2.0
        expected.loc[13, "A"] = 5
        tm.assert_frame_equal(result, expected, check_dtype=False)

    def test_interp_nearest_single_valid_value(self):
        # GH#46230
        pytest.importorskip("scipy")
        df = pd.DataFrame({"a": [np.nan, 1]})
        result = df.interpolate(method="nearest")
        tm.assert_frame_equal(result, df)

    def test_interp_alt_scipy(self):
        pytest.importorskip("scipy")
        df = pd.DataFrame(
            {"A": [1, 2, np.nan, 4, 5, np.nan, 7], "C": [1, 2, 3, 5, 8, 13, 21]}
        )
        result = df.interpolate(method="barycentric")
        expected = df.copy()
        expected.loc[2, "A"] = 3
        expected.loc[5, "A"] = 6
        tm.assert_frame_equal(result, expected)

        result = df.interpolate(method="krogh")
        expectedk = df.copy()
        expectedk["A"] = expected["A"]
        tm.assert_frame_equal(result, expectedk)

        result = df.interpolate(method="pchip")
        expected.loc[2, "A"] = 3
        expected.loc[5, "A"] = 6.0

        tm.assert_frame_equal(result, expected)

    def test_interp_rowwise(self):
        df = pd.DataFrame(
            {
                0: [1, 2, np.nan, 4],
                1: [2, 3, 4, np.nan],
                2: [np.nan, 4, 5, 6],
                3: [4, np.nan, 6, 7],
                4: [1, 2, 3, 4],
            }
        )
        result = df.interpolate(axis=1)
        expected = df.copy()
        expected.loc[3, 1] = 5
        expected.loc[0, 2] = 3
        expected.loc[1, 3] = 3
        expected[4] = expected[4].astype(np.float64)
        tm.assert_frame_equal(result, expected)

        result = df.interpolate(axis=1, method="values")
        tm.assert_frame_equal(result, expected)

        result = df.interpolate(axis=0)
        expected = df.interpolate()
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "axis_name, axis_number",
        [
            pytest.param("rows", 0, id="rows_0"),
            pytest.param("index", 0, id="index_0"),
            pytest.param("columns", 1, id="columns_1"),
        ],
    )
    def test_interp_axis_names(self, axis_name, axis_number):
        # GH 29132: test axis names
        data = {0: [0, np.nan, 6], 1: [1, np.nan, 7], 2: [2, 5, 8]}

        df = pd.DataFrame(data, dtype=np.float64)
        result = df.interpolate(axis=axis_name, method="linear")
        expected = df.interpolate(axis=axis_number, method="linear")
        tm.assert_frame_equal(result, expected)

    def test_rowwise_alt(self):
        df = pd.DataFrame(
            {
                0: [0, 0.5, 1.0, np.nan, 4, 8, np.nan, np.nan, 64],
                1: [1, 2, 3, 4, 3, 2, 1, 0, -1],
            }
        )
        df.interpolate(axis=0)
        # TODO: assert something?

    @pytest.mark.parametrize(
        "check_scipy", [False, pytest.param(True, marks=td.skip_if_no("scipy"))]
    )
    def test_interp_leading_nans(self, check_scipy):
        df = pd.DataFrame(
            {"A": [np.nan, np.nan, 0.5, 0.25, 0], "B": [np.nan, -3, -3.5, np.nan, -4]}
        )
        result = df.interpolate()
        expected = df.copy()
        expected.loc[3, "B"] = -3.75
        tm.assert_frame_equal(result, expected)

        if check_scipy:
            result = df.interpolate(method="polynomial", order=1)
            tm.assert_frame_equal(result, expected)

    def test_interp_raise_on_only_mixed(self, axis):
        df = pd.DataFrame(
            {
                "A": [1, 2, np.nan, 4],
                "B": ["a", "b", "c", "d"],
                "C": [np.nan, 2, 5, 7],
                "D": [np.nan, np.nan, 9, 9],
                "E": [1, 2, 3, 4],
            }
        )
        msg = "DataFrame cannot interpolate with object dtype"
        with pytest.raises(TypeError, match=msg):
            df.astype("object").interpolate(axis=axis)

    def test_interp_raise_on_all_object_dtype(self):
        # GH 22985
        df = pd.DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}, dtype="object")
        msg = "DataFrame cannot interpolate with object dtype"
        with pytest.raises(TypeError, match=msg):
            df.interpolate()

    def test_interp_inplace(self):
        df = pd.DataFrame({"a": [1.0, 2.0, np.nan, 4.0]})
        df_orig = df.copy()
        expected = df.copy().interpolate()

        with tm.raises_chained_assignment_error():
            result = df["a"].interpolate(inplace=True)
        tm.assert_series_equal(result, expected["a"])
        tm.assert_frame_equal(df, df_orig)

    def test_interp_inplace_row(self):
        # GH 10395
        df = pd.DataFrame(
            {"a": [1.0, 2.0, 3.0, 4.0], "b": [np.nan, 2.0, 3.0, 4.0], "c": [3, 2, 2, 2]}
        )
        expected = df.interpolate(method="linear", axis=1, inplace=False)
        result = df.interpolate(method="linear", axis=1, inplace=True)
        assert result is df
        tm.assert_frame_equal(result, expected)

    def test_interp_ignore_all_good(self):
        # GH#6290
        df = pd.DataFrame(
            {
                "A": [1, 2, np.nan, 4],
                "B": [1, 2, 3, 4],
                "C": [1.0, 2.0, np.nan, 4.0],
                "D": [1.0, 2.0, 3.0, 4.0],
            }
        )
        expected = pd.DataFrame(
            {
                "A": np.array([1, 2, 3, 4], dtype="float64"),
                "B": np.array([1, 2, 3, 4], dtype="int64"),
                "C": np.array([1.0, 2.0, 3, 4.0], dtype="float64"),
                "D": np.array([1.0, 2.0, 3.0, 4.0], dtype="float64"),
            }
        )
        result = df.interpolate()
        tm.assert_frame_equal(result, expected)

        # all good
        result = df[["B", "D"]].interpolate()
        tm.assert_frame_equal(result, df[["B", "D"]])

    def test_interp_time_inplace_axis(self):
        # GH 9687
        periods = 5
        idx = pd.date_range(start="2014-01-01", periods=periods)
        data = np.random.default_rng(2).random((periods, periods))
        data[data < 0.5] = np.nan
        df = pd.DataFrame(index=idx, columns=idx, data=data)

        expected = df.interpolate(axis=0, method="time")
        result = df.interpolate(axis=0, method="time", inplace=True)
        assert result is df
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("axis_name, axis_number", [("index", 0), ("columns", 1)])
    def test_interp_string_axis(self, axis_name, axis_number):
        # https://github.com/pandas-dev/pandas/issues/25190
        x = np.linspace(0, 100, 3)
        y = np.sin(x)
        df = pd.DataFrame(
            data=np.tile(y, (10, 1)), index=np.arange(10), columns=x
        ).reindex(columns=x * 1.005)
        result = df.interpolate(method="linear", axis=axis_name)
        expected = df.interpolate(method="linear", axis=axis_number)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("multiblock", [True, False])
    @pytest.mark.parametrize("method", ["ffill", "bfill", "pad"])
    def test_interp_fillna_methods(self, axis, multiblock, method):
        # GH 12918
        df = pd.DataFrame(
            {
                "A": [1.0, 2.0, 3.0, 4.0, np.nan, 5.0],
                "B": [2.0, 4.0, 6.0, np.nan, 8.0, 10.0],
                "C": [3.0, 6.0, 9.0, np.nan, np.nan, 30.0],
            }
        )
        if multiblock:
            df["D"] = np.nan
            df["E"] = 1.0

        msg = f"Can not interpolate with method={method}"
        with pytest.raises(ValueError, match=msg):
            df.interpolate(method=method, axis=axis)

    def test_interpolate_empty_df(self):
        # GH#53199
        df = pd.DataFrame()
        expected = df.copy()
        result = df.interpolate(inplace=True)
        assert result is df
        tm.assert_frame_equal(result, expected)

    def test_interpolate_ea(self, any_int_ea_dtype):
        # GH#55347
        df = pd.DataFrame({"a": [1, None, None, None, 3]}, dtype=any_int_ea_dtype)
        orig = df.copy()
        result = df.interpolate(limit=2)
        expected = pd.DataFrame({"a": [1, 1.5, 2.0, None, 3]}, dtype="Float64")
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(df, orig)

    @pytest.mark.parametrize(
        "dtype",
        [
            "Float64",
            "Float32",
            pytest.param("float32[pyarrow]", marks=td.skip_if_no("pyarrow")),
            pytest.param("float64[pyarrow]", marks=td.skip_if_no("pyarrow")),
        ],
    )
    def test_interpolate_ea_float(self, dtype):
        # GH#55347
        df = pd.DataFrame({"a": [1, None, None, None, 3]}, dtype=dtype)
        orig = df.copy()
        result = df.interpolate(limit=2)
        expected = pd.DataFrame({"a": [1, 1.5, 2.0, None, 3]}, dtype=dtype)
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(df, orig)

    @pytest.mark.parametrize(
        "dtype",
        ["int64", "uint64", "int32", "int16", "int8", "uint32", "uint16", "uint8"],
    )
    def test_interpolate_arrow(self, dtype):
        # GH#55347
        pytest.importorskip("pyarrow")
        df = pd.DataFrame({"a": [1, None, None, None, 3]}, dtype=dtype + "[pyarrow]")
        result = df.interpolate(limit=2)
        expected = pd.DataFrame({"a": [1, 1.5, 2.0, None, 3]}, dtype="float64[pyarrow]")
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype",
        [
            "Int64",
            "Float64",
        ],
    )
    def test_interpolate_time_nullable_int_float(self, dtype):
        # GH#40252
        idx = pd.date_range("1970-01-02", periods=3, freq="D")

        df = pd.DataFrame({"a": [1, None, 2]}, index=idx, dtype=dtype)
        result = df.interpolate(method="time")
        expected = pd.DataFrame({"a": [1.0, 1.5, 2.0]}, index=idx, dtype="Float64")
        tm.assert_frame_equal(result, expected)

    def test_interpolate_fill_value_extrapolate_time(self):
        # GH#31949
        index = pd.to_datetime(["2020-01-01", "2020-01-02", "2020-01-04", "2020-01-08"])
        df = pd.DataFrame({"a": [np.nan, 1.0, 3.0, np.nan]}, index=index)
        result = df.interpolate(
            method="time", fill_value="extrapolate", limit_direction="both"
        )
        expected = pd.DataFrame({"a": [0.0, 1.0, 3.0, 7.0]}, index=index)
        tm.assert_frame_equal(result, expected)


def _spy_on_scipy_calls(monkeypatch):
    calls = []
    real = missing._interpolate_scipy_wrapper

    def spy(*args, **kwargs):
        calls.append(1)
        return real(*args, **kwargs)

    monkeypatch.setattr(missing, "_interpolate_scipy_wrapper", spy)
    return calls


def _interpolate_each_column(df, **kwargs):
    # Series always take the per-slice path
    return pd.DataFrame(
        {name: df[name].interpolate(**kwargs) for name in df}, index=df.index
    )


@pytest.mark.parametrize(
    "kwargs",
    [
        {"method": "nearest"},
        {"method": "zero"},
        {"method": "slinear"},
        {"method": "quadratic"},
        {"method": "cubic"},
        {"method": "polynomial", "order": 3},
        {"method": "nearest", "limit": 1},
        {"method": "nearest", "limit_direction": "both"},
        {"method": "nearest", "limit_area": "inside"},
        {"method": "nearest", "limit_area": "outside"},
        {"method": "slinear", "fill_value": "extrapolate"},
        {"method": "nearest", "fill_value": 7.5},
    ],
    ids=str,
)
@pytest.mark.parametrize(
    "index_values", [None, [0.0, 1.5, 2.0, 6.0, 6.5, 9.0, 20.0], [5, 1, 3, 0, 6, 2, 4]]
)
@pytest.mark.parametrize("dtype", ["float64", "float32", "complex128"])
def test_interpolate_frame_same_nans_matches_series(kwargs, index_values, dtype):
    # GH#48236 columns with their NaNs in the same places are interpolated together
    pytest.importorskip("scipy")
    first = [np.nan, 1.0, np.nan, 3.5, 5.0, np.nan, -2.0]
    second = [np.nan, 7.0, np.nan, -3.0, 0.25, np.nan, 8.0]
    third = [4.0, 1.0, np.nan, np.nan, 5.0, -1.0, 2.0]
    data = {
        "a": first,
        "b": second,
        "c": third,
        "d": first,
        "e": third,
        "has_inf": [4.0, np.inf, np.nan, np.nan, 5.0, -1.0, 2.0],
        "f": [9.0, 8.0, 6.0, 5.5, np.nan, 2.0, 1.0],
        "all_nan": [np.nan] * 7,
        "no_nan": [1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0],
    }
    df = pd.DataFrame(data, index=index_values, dtype=dtype)
    expected = _interpolate_each_column(df, **kwargs)
    tm.assert_frame_equal(df.interpolate(**kwargs), expected, check_exact=True)
    tm.assert_frame_equal(
        df.T.interpolate(axis=1, **kwargs), expected.T, check_exact=True
    )

    result = df.copy()
    result.interpolate(inplace=True, **kwargs)
    tm.assert_frame_equal(result, expected, check_exact=True)


def test_interpolate_frame_nans_differ_after_eighth_row():
    # GH#48236 these must not end up in the same group
    pytest.importorskip("scipy")
    first = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, np.nan, 10.0]
    second = [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, np.nan]
    df = pd.DataFrame({"a": first, "b": second})
    result = df.interpolate(method="nearest")
    expected = _interpolate_each_column(df, method="nearest")
    tm.assert_frame_equal(result, expected, check_exact=True)


def test_interpolate_frame_same_nans_one_valid_value():
    # GH#48236
    pytest.importorskip("scipy")
    df = pd.DataFrame({"a": [np.nan, 1.0, np.nan], "b": [np.nan, 2.0, np.nan]})
    result = df.interpolate(method="nearest")
    expected = _interpolate_each_column(df, method="nearest")
    tm.assert_frame_equal(result, expected, check_exact=True)


def test_interpolate_frame_same_nans_one_scipy_call_per_pattern(monkeypatch):
    # GH#48236 the equivalence tests can't tell how often SciPy was called
    pytest.importorskip("scipy")
    calls = _spy_on_scipy_calls(monkeypatch)
    first = [1.0, np.nan, 3.0, 4.0]
    second = [1.0, 2.0, np.nan, 4.0]
    df = pd.DataFrame({"a": first, "b": first, "c": second, "d": first})
    df.interpolate(method="nearest")
    assert len(calls) == 2


@pytest.mark.parametrize("chunk_size, n_calls", [(8, 3), (1, 5)])
def test_interpolate_frame_same_nans_chunks(chunk_size, n_calls, monkeypatch):
    # GH#48236 a group split over several SciPy calls must come out the same
    pytest.importorskip("scipy")
    values = [1.0, np.nan, 3.0, 4.0]
    df = pd.DataFrame(dict.fromkeys("abcde", values)) * [1.0, 2.0, 3.0, 4.0, 5.0]
    expected = _interpolate_each_column(df, method="nearest")
    calls = _spy_on_scipy_calls(monkeypatch)
    monkeypatch.setattr(missing, "_INTERP1D_CHUNK_SIZE", chunk_size)
    tm.assert_frame_equal(df.interpolate(method="nearest"), expected, check_exact=True)
    assert len(calls) == n_calls


def test_interpolate_frame_same_nans_list_fill_value_raises():
    # GH#48236 must not start broadcasting against the columns
    pytest.importorskip("scipy")
    values = [1.0, np.nan, 3.0]
    df = pd.DataFrame({"a": values, "b": values})
    with pytest.raises(ValueError, match="broadcast"):
        df.interpolate(method="nearest", fill_value=[1.0, 2.0])


def test_interpolate_frame_longdouble_matches_series():
    # GH#48236 SciPy can reject a 1-d longdouble y where it takes a 2-d one
    pytest.importorskip("scipy")
    values = np.array([1.0, np.nan, 3.0, 4.0, 5.0, 7.0], dtype=np.longdouble)
    df = pd.DataFrame({"a": values, "b": values})
    try:
        expected = df["a"].interpolate(method="cubic")
    except ValueError as err:
        with pytest.raises(ValueError, match=re.escape(str(err))):
            df.interpolate(method="cubic")
    else:
        result = df.interpolate(method="cubic")
        tm.assert_series_equal(result["a"], expected)


def test_interpolate_frame_same_nans_first_error_in_column_order():
    # GH#48236
    pytest.importorskip("scipy")
    df = pd.DataFrame(
        {"a": [1.0, 2.0, np.nan, np.nan], "b": [np.nan, np.nan, np.nan, 2.0]}
    )
    with pytest.raises(ValueError, match="derivatives") as expected:
        df["a"].interpolate(method="cubic")
    with pytest.raises(ValueError, match="derivatives") as result:
        df.interpolate(method="cubic")
    assert str(result.value) == str(expected.value)


def test_interpolate_frame_no_partial_columns_without_scipy(monkeypatch):
    # GH#48236 a frame whose columns are all valid or all NaN never needs SciPy
    real_import = missing.import_optional_dependency

    def no_scipy(name, *args, **kwargs):
        if name == "scipy":
            raise ImportError("Missing optional dependency 'scipy'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(missing, "import_optional_dependency", no_scipy)
    df = pd.DataFrame({"a": [1.0, 2.0, 3.0], "b": [np.nan, np.nan, np.nan]})
    tm.assert_frame_equal(df.interpolate(method="nearest"), df)


@pytest.mark.parametrize(
    "values",
    [[[1.0, 2.0], [3.0, 4.0]], [[1.0, np.nan, np.nan], [3.0, np.nan, np.nan]]],
    ids=["no_nan", "all_nan_columns"],
)
def test_interpolate_frame_inplace_readonly_nothing_to_fill(values):
    # GH#48236 a column with nothing to fill is not written to
    pytest.importorskip("scipy")
    arr = np.array(values)
    arr.flags.writeable = False
    df = pd.DataFrame(arr, copy=False)
    expected = df.copy()
    df.interpolate(method="nearest", inplace=True)
    tm.assert_frame_equal(df, expected)
