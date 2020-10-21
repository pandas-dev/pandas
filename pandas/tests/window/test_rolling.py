from datetime import datetime, timedelta

import numpy as np
import pytest

from pandas.errors import UnsupportedFunctionCall
import pandas.util._test_decorators as td

import pandas as pd
from pandas import DataFrame, Series, date_range
import pandas._testing as tm
from pandas.core.window import Rolling


def test_doc_string():

    df = DataFrame({"B": [0, 1, 2, np.nan, 4]})
    df
    df.rolling(2).sum()
    df.rolling(2, min_periods=1).sum()


def test_constructor(which):
    # GH 12669

    c = which.rolling

    # valid
    c(0)
    c(window=2)
    c(window=2, min_periods=1)
    c(window=2, min_periods=1, center=True)
    c(window=2, min_periods=1, center=False)

    # GH 13383

    msg = "window must be non-negative"

    with pytest.raises(ValueError, match=msg):
        c(-1)


@pytest.mark.parametrize("w", [2.0, "foo", np.array([2])])
def test_invalid_constructor(which, w):
    # not valid

    c = which.rolling

    msg = (
        "window must be an integer|"
        "passed window foo is not compatible with a datetimelike index"
    )
    with pytest.raises(ValueError, match=msg):
        c(window=w)

    msg = "min_periods must be an integer"
    with pytest.raises(ValueError, match=msg):
        c(window=2, min_periods=w)

    msg = "center must be a boolean"
    with pytest.raises(ValueError, match=msg):
        c(window=2, min_periods=1, center=w)


@td.skip_if_no_scipy
def test_constructor_with_win_type(which):
    # GH 13383
    c = which.rolling

    msg = "window must be > 0"

    with pytest.raises(ValueError, match=msg):
        c(-1, win_type="boxcar")


@pytest.mark.parametrize("window", [timedelta(days=3), pd.Timedelta(days=3)])
def test_constructor_with_timedelta_window(window):
    # GH 15440
    n = 10
    df = DataFrame(
        {"value": np.arange(n)}, index=pd.date_range("2015-12-24", periods=n, freq="D")
    )
    expected_data = np.append([0.0, 1.0], np.arange(3.0, 27.0, 3))

    result = df.rolling(window=window).sum()
    expected = DataFrame(
        {"value": expected_data},
        index=pd.date_range("2015-12-24", periods=n, freq="D"),
    )
    tm.assert_frame_equal(result, expected)
    expected = df.rolling("3D").sum()
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("window", [timedelta(days=3), pd.Timedelta(days=3), "3D"])
def test_constructor_timedelta_window_and_minperiods(window, raw):
    # GH 15305
    n = 10
    df = DataFrame(
        {"value": np.arange(n)}, index=pd.date_range("2017-08-08", periods=n, freq="D")
    )
    expected = DataFrame(
        {"value": np.append([np.NaN, 1.0], np.arange(3.0, 27.0, 3))},
        index=pd.date_range("2017-08-08", periods=n, freq="D"),
    )
    result_roll_sum = df.rolling(window=window, min_periods=2).sum()
    result_roll_generic = df.rolling(window=window, min_periods=2).apply(sum, raw=raw)
    tm.assert_frame_equal(result_roll_sum, expected)
    tm.assert_frame_equal(result_roll_generic, expected)


@pytest.mark.parametrize("method", ["std", "mean", "sum", "max", "min", "var"])
def test_numpy_compat(method):
    # see gh-12811
    r = Rolling(Series([2, 4, 6]), window=2)

    msg = "numpy operations are not valid with window objects"

    with pytest.raises(UnsupportedFunctionCall, match=msg):
        getattr(r, method)(1, 2, 3)
    with pytest.raises(UnsupportedFunctionCall, match=msg):
        getattr(r, method)(dtype=np.float64)


@pytest.mark.parametrize("closed", ["left", "right", "both", "neither"])
def test_closed_fixed(closed, arithmetic_win_operators):
    # GH 34315
    func_name = arithmetic_win_operators
    df_fixed = DataFrame({"A": [0, 1, 2, 3, 4]})
    df_time = DataFrame({"A": [0, 1, 2, 3, 4]}, index=date_range("2020", periods=5))

    result = getattr(df_fixed.rolling(2, closed=closed, min_periods=1), func_name)()
    expected = getattr(df_time.rolling("2D", closed=closed), func_name)().reset_index(
        drop=True
    )

    tm.assert_frame_equal(result, expected)


def test_closed_fixed_binary_col():
    # GH 34315
    data = [0, 1, 1, 0, 0, 1, 0, 1]
    df = DataFrame(
        {"binary_col": data},
        index=pd.date_range(start="2020-01-01", freq="min", periods=len(data)),
    )

    rolling = df.rolling(window=len(df), closed="left", min_periods=1)
    result = rolling.mean()
    expected = DataFrame(
        [np.nan, 0, 0.5, 2 / 3, 0.5, 0.4, 0.5, 0.428571],
        columns=["binary_col"],
        index=pd.date_range(start="2020-01-01", freq="min", periods=len(data)),
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("closed", ["neither", "left"])
def test_closed_empty(closed, arithmetic_win_operators):
    # GH 26005
    func_name = arithmetic_win_operators
    ser = Series(data=np.arange(5), index=pd.date_range("2000", periods=5, freq="2D"))
    roll = ser.rolling("1D", closed=closed)

    result = getattr(roll, func_name)()
    expected = Series([np.nan] * 5, index=ser.index)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("func", ["min", "max"])
def test_closed_one_entry(func):
    # GH24718
    ser = Series(data=[2], index=pd.date_range("2000", periods=1))
    result = getattr(ser.rolling("10D", closed="left"), func)()
    tm.assert_series_equal(result, Series([np.nan], index=ser.index))


@pytest.mark.parametrize("func", ["min", "max"])
def test_closed_one_entry_groupby(func):
    # GH24718
    ser = DataFrame(
        data={"A": [1, 1, 2], "B": [3, 2, 1]}, index=pd.date_range("2000", periods=3)
    )
    result = getattr(
        ser.groupby("A", sort=False)["B"].rolling("10D", closed="left"), func
    )()
    exp_idx = pd.MultiIndex.from_arrays(
        arrays=[[1, 1, 2], ser.index], names=("A", None)
    )
    expected = Series(data=[np.nan, 3, np.nan], index=exp_idx, name="B")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("input_dtype", ["int", "float"])
@pytest.mark.parametrize(
    "func,closed,expected",
    [
        ("min", "right", [0.0, 0, 0, 1, 2, 3, 4, 5, 6, 7]),
        ("min", "both", [0.0, 0, 0, 0, 1, 2, 3, 4, 5, 6]),
        ("min", "neither", [np.nan, 0, 0, 1, 2, 3, 4, 5, 6, 7]),
        ("min", "left", [np.nan, 0, 0, 0, 1, 2, 3, 4, 5, 6]),
        ("max", "right", [0.0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
        ("max", "both", [0.0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
        ("max", "neither", [np.nan, 0, 1, 2, 3, 4, 5, 6, 7, 8]),
        ("max", "left", [np.nan, 0, 1, 2, 3, 4, 5, 6, 7, 8]),
    ],
)
def test_closed_min_max_datetime(input_dtype, func, closed, expected):
    # see gh-21704
    ser = Series(
        data=np.arange(10).astype(input_dtype), index=pd.date_range("2000", periods=10)
    )

    result = getattr(ser.rolling("3D", closed=closed), func)()
    expected = Series(expected, index=ser.index)
    tm.assert_series_equal(result, expected)


def test_closed_uneven():
    # see gh-21704
    ser = Series(data=np.arange(10), index=pd.date_range("2000", periods=10))

    # uneven
    ser = ser.drop(index=ser.index[[1, 5]])
    result = ser.rolling("3D", closed="left").min()
    expected = Series([np.nan, 0, 0, 2, 3, 4, 6, 6], index=ser.index)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "func,closed,expected",
    [
        ("min", "right", [np.nan, 0, 0, 1, 2, 3, 4, 5, np.nan, np.nan]),
        ("min", "both", [np.nan, 0, 0, 0, 1, 2, 3, 4, 5, np.nan]),
        ("min", "neither", [np.nan, np.nan, 0, 1, 2, 3, 4, 5, np.nan, np.nan]),
        ("min", "left", [np.nan, np.nan, 0, 0, 1, 2, 3, 4, 5, np.nan]),
        ("max", "right", [np.nan, 1, 2, 3, 4, 5, 6, 6, np.nan, np.nan]),
        ("max", "both", [np.nan, 1, 2, 3, 4, 5, 6, 6, 6, np.nan]),
        ("max", "neither", [np.nan, np.nan, 1, 2, 3, 4, 5, 6, np.nan, np.nan]),
        ("max", "left", [np.nan, np.nan, 1, 2, 3, 4, 5, 6, 6, np.nan]),
    ],
)
def test_closed_min_max_minp(func, closed, expected):
    # see gh-21704
    ser = Series(data=np.arange(10), index=pd.date_range("2000", periods=10))
    ser[ser.index[-3:]] = np.nan
    result = getattr(ser.rolling("3D", min_periods=2, closed=closed), func)()
    expected = Series(expected, index=ser.index)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "closed,expected",
    [
        ("right", [0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8]),
        ("both", [0, 0.5, 1, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]),
        ("neither", [np.nan, 0, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]),
        ("left", [np.nan, 0, 0.5, 1, 2, 3, 4, 5, 6, 7]),
    ],
)
def test_closed_median_quantile(closed, expected):
    # GH 26005
    ser = Series(data=np.arange(10), index=pd.date_range("2000", periods=10))
    roll = ser.rolling("3D", closed=closed)
    expected = Series(expected, index=ser.index)

    result = roll.median()
    tm.assert_series_equal(result, expected)

    result = roll.quantile(0.5)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("roller", ["1s", 1])
def tests_empty_df_rolling(roller):
    # GH 15819 Verifies that datetime and integer rolling windows can be
    # applied to empty DataFrames
    expected = DataFrame()
    result = DataFrame().rolling(roller).sum()
    tm.assert_frame_equal(result, expected)

    # Verifies that datetime and integer rolling windows can be applied to
    # empty DataFrames with datetime index
    expected = DataFrame(index=pd.DatetimeIndex([]))
    result = DataFrame(index=pd.DatetimeIndex([])).rolling(roller).sum()
    tm.assert_frame_equal(result, expected)


def test_empty_window_median_quantile():
    # GH 26005
    expected = Series([np.nan, np.nan, np.nan])
    roll = Series(np.arange(3)).rolling(0)

    result = roll.median()
    tm.assert_series_equal(result, expected)

    result = roll.quantile(0.1)
    tm.assert_series_equal(result, expected)


def test_missing_minp_zero():
    # https://github.com/pandas-dev/pandas/pull/18921
    # minp=0
    x = Series([np.nan])
    result = x.rolling(1, min_periods=0).sum()
    expected = Series([0.0])
    tm.assert_series_equal(result, expected)

    # minp=1
    result = x.rolling(1, min_periods=1).sum()
    expected = Series([np.nan])
    tm.assert_series_equal(result, expected)


def test_missing_minp_zero_variable():
    # https://github.com/pandas-dev/pandas/pull/18921
    x = Series(
        [np.nan] * 4,
        index=pd.DatetimeIndex(
            ["2017-01-01", "2017-01-04", "2017-01-06", "2017-01-07"]
        ),
    )
    result = x.rolling(pd.Timedelta("2d"), min_periods=0).sum()
    expected = Series(0.0, index=x.index)
    tm.assert_series_equal(result, expected)


def test_multi_index_names():

    # GH 16789, 16825
    cols = pd.MultiIndex.from_product([["A", "B"], ["C", "D", "E"]], names=["1", "2"])
    df = DataFrame(np.ones((10, 6)), columns=cols)
    result = df.rolling(3).cov()

    tm.assert_index_equal(result.columns, df.columns)
    assert result.index.names == [None, "1", "2"]


def test_rolling_axis_sum(axis_frame):
    # see gh-23372.
    df = DataFrame(np.ones((10, 20)))
    axis = df._get_axis_number(axis_frame)

    if axis == 0:
        expected = DataFrame({i: [np.nan] * 2 + [3.0] * 8 for i in range(20)})
    else:
        # axis == 1
        expected = DataFrame([[np.nan] * 2 + [3.0] * 18] * 10)

    result = df.rolling(3, axis=axis_frame).sum()
    tm.assert_frame_equal(result, expected)


def test_rolling_axis_count(axis_frame):
    # see gh-26055
    df = DataFrame({"x": range(3), "y": range(3)})

    axis = df._get_axis_number(axis_frame)

    if axis in [0, "index"]:
        expected = DataFrame({"x": [1.0, 2.0, 2.0], "y": [1.0, 2.0, 2.0]})
    else:
        expected = DataFrame({"x": [1.0, 1.0, 1.0], "y": [2.0, 2.0, 2.0]})

    result = df.rolling(2, axis=axis_frame, min_periods=0).count()
    tm.assert_frame_equal(result, expected)


def test_readonly_array():
    # GH-27766
    arr = np.array([1, 3, np.nan, 3, 5])
    arr.setflags(write=False)
    result = Series(arr).rolling(2).mean()
    expected = Series([np.nan, 2, np.nan, np.nan, 4])
    tm.assert_series_equal(result, expected)


def test_rolling_datetime(axis_frame, tz_naive_fixture):
    # GH-28192
    tz = tz_naive_fixture
    df = DataFrame(
        {i: [1] * 2 for i in pd.date_range("2019-8-01", "2019-08-03", freq="D", tz=tz)}
    )
    if axis_frame in [0, "index"]:
        result = df.T.rolling("2D", axis=axis_frame).sum().T
    else:
        result = df.rolling("2D", axis=axis_frame).sum()
    expected = DataFrame(
        {
            **{
                i: [1.0] * 2
                for i in pd.date_range("2019-8-01", periods=1, freq="D", tz=tz)
            },
            **{
                i: [2.0] * 2
                for i in pd.date_range("2019-8-02", "2019-8-03", freq="D", tz=tz)
            },
        }
    )
    tm.assert_frame_equal(result, expected)


def test_rolling_window_as_string():
    # see gh-22590
    date_today = datetime.now()
    days = pd.date_range(date_today, date_today + timedelta(365), freq="D")

    npr = np.random.RandomState(seed=421)

    data = npr.randint(1, high=100, size=len(days))
    df = DataFrame({"DateCol": days, "metric": data})

    df.set_index("DateCol", inplace=True)
    result = df.rolling(window="21D", min_periods=2, closed="left")["metric"].agg("max")

    expData = (
        [np.nan] * 2
        + [88.0] * 16
        + [97.0] * 9
        + [98.0]
        + [99.0] * 21
        + [95.0] * 16
        + [93.0] * 5
        + [89.0] * 5
        + [96.0] * 21
        + [94.0] * 14
        + [90.0] * 13
        + [88.0] * 2
        + [90.0] * 9
        + [96.0] * 21
        + [95.0] * 6
        + [91.0]
        + [87.0] * 6
        + [92.0] * 21
        + [83.0] * 2
        + [86.0] * 10
        + [87.0] * 5
        + [98.0] * 21
        + [97.0] * 14
        + [93.0] * 7
        + [87.0] * 4
        + [86.0] * 4
        + [95.0] * 21
        + [85.0] * 14
        + [83.0] * 2
        + [76.0] * 5
        + [81.0] * 2
        + [98.0] * 21
        + [95.0] * 14
        + [91.0] * 7
        + [86.0]
        + [93.0] * 3
        + [95.0] * 20
    )

    expected = Series(
        expData, index=days.rename("DateCol")._with_freq(None), name="metric"
    )
    tm.assert_series_equal(result, expected)


def test_min_periods1():
    # GH#6795
    df = DataFrame([0, 1, 2, 1, 0], columns=["a"])
    result = df["a"].rolling(3, center=True, min_periods=1).max()
    expected = Series([1.0, 2.0, 2.0, 2.0, 1.0], name="a")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("constructor", [Series, DataFrame])
def test_rolling_count_with_min_periods(constructor):
    # GH 26996
    result = constructor(range(5)).rolling(3, min_periods=3).count()
    expected = constructor([np.nan, np.nan, 3.0, 3.0, 3.0])
    tm.assert_equal(result, expected)


@pytest.mark.parametrize("constructor", [Series, DataFrame])
def test_rolling_count_default_min_periods_with_null_values(constructor):
    # GH 26996
    values = [1, 2, 3, np.nan, 4, 5, 6]
    expected_counts = [1.0, 2.0, 3.0, 2.0, 2.0, 2.0, 3.0]

    # GH 31302
    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        result = constructor(values).rolling(3).count()
    expected = constructor(expected_counts)
    tm.assert_equal(result, expected)


@pytest.mark.parametrize(
    "df,expected,window,min_periods",
    [
        (
            DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}),
            [
                ({"A": [1], "B": [4]}, [0]),
                ({"A": [1, 2], "B": [4, 5]}, [0, 1]),
                ({"A": [1, 2, 3], "B": [4, 5, 6]}, [0, 1, 2]),
            ],
            3,
            None,
        ),
        (
            DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}),
            [
                ({"A": [1], "B": [4]}, [0]),
                ({"A": [1, 2], "B": [4, 5]}, [0, 1]),
                ({"A": [2, 3], "B": [5, 6]}, [1, 2]),
            ],
            2,
            1,
        ),
        (
            DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}),
            [
                ({"A": [1], "B": [4]}, [0]),
                ({"A": [1, 2], "B": [4, 5]}, [0, 1]),
                ({"A": [2, 3], "B": [5, 6]}, [1, 2]),
            ],
            2,
            2,
        ),
        (
            DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}),
            [
                ({"A": [1], "B": [4]}, [0]),
                ({"A": [2], "B": [5]}, [1]),
                ({"A": [3], "B": [6]}, [2]),
            ],
            1,
            1,
        ),
        (
            DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}),
            [
                ({"A": [1], "B": [4]}, [0]),
                ({"A": [2], "B": [5]}, [1]),
                ({"A": [3], "B": [6]}, [2]),
            ],
            1,
            0,
        ),
        (DataFrame({"A": [1], "B": [4]}), [], 2, None),
        (DataFrame({"A": [1], "B": [4]}), [], 2, 1),
        (DataFrame(), [({}, [])], 2, None),
        (
            DataFrame({"A": [1, np.nan, 3], "B": [np.nan, 5, 6]}),
            [
                ({"A": [1.0], "B": [np.nan]}, [0]),
                ({"A": [1, np.nan], "B": [np.nan, 5]}, [0, 1]),
                ({"A": [1, np.nan, 3], "B": [np.nan, 5, 6]}, [0, 1, 2]),
            ],
            3,
            2,
        ),
    ],
)
def test_iter_rolling_dataframe(df, expected, window, min_periods):
    # GH 11704
    expected = [DataFrame(values, index=index) for (values, index) in expected]

    for (expected, actual) in zip(
        expected, df.rolling(window, min_periods=min_periods)
    ):
        tm.assert_frame_equal(actual, expected)


@pytest.mark.parametrize(
    "expected,window",
    [
        (
            [
                ({"A": [1], "B": [4]}, [0]),
                ({"A": [1, 2], "B": [4, 5]}, [0, 1]),
                ({"A": [2, 3], "B": [5, 6]}, [1, 2]),
            ],
            "2D",
        ),
        (
            [
                ({"A": [1], "B": [4]}, [0]),
                ({"A": [1, 2], "B": [4, 5]}, [0, 1]),
                ({"A": [1, 2, 3], "B": [4, 5, 6]}, [0, 1, 2]),
            ],
            "3D",
        ),
        (
            [
                ({"A": [1], "B": [4]}, [0]),
                ({"A": [2], "B": [5]}, [1]),
                ({"A": [3], "B": [6]}, [2]),
            ],
            "1D",
        ),
    ],
)
def test_iter_rolling_on_dataframe(expected, window):
    # GH 11704
    df = DataFrame(
        {
            "A": [1, 2, 3, 4, 5],
            "B": [4, 5, 6, 7, 8],
            "C": date_range(start="2016-01-01", periods=5, freq="D"),
        }
    )

    expected = [DataFrame(values, index=index) for (values, index) in expected]
    for (expected, actual) in zip(expected, df.rolling(window, on="C")):
        tm.assert_frame_equal(actual, expected)


@pytest.mark.parametrize(
    "ser,expected,window, min_periods",
    [
        (
            Series([1, 2, 3]),
            [([1], [0]), ([1, 2], [0, 1]), ([1, 2, 3], [0, 1, 2])],
            3,
            None,
        ),
        (
            Series([1, 2, 3]),
            [([1], [0]), ([1, 2], [0, 1]), ([1, 2, 3], [0, 1, 2])],
            3,
            1,
        ),
        (Series([1, 2, 3]), [([1], [0]), ([1, 2], [0, 1]), ([2, 3], [1, 2])], 2, 1),
        (Series([1, 2, 3]), [([1], [0]), ([1, 2], [0, 1]), ([2, 3], [1, 2])], 2, 2),
        (Series([1, 2, 3]), [([1], [0]), ([2], [1]), ([3], [2])], 1, 0),
        (Series([1, 2, 3]), [([1], [0]), ([2], [1]), ([3], [2])], 1, 1),
        (Series([1, 2]), [([1], [0]), ([1, 2], [0, 1])], 2, 0),
        (Series([], dtype="int64"), [], 2, 1),
    ],
)
def test_iter_rolling_series(ser, expected, window, min_periods):
    # GH 11704
    expected = [Series(values, index=index) for (values, index) in expected]

    for (expected, actual) in zip(
        expected, ser.rolling(window, min_periods=min_periods)
    ):
        tm.assert_series_equal(actual, expected)


@pytest.mark.parametrize(
    "expected,expected_index,window",
    [
        (
            [[0], [1], [2], [3], [4]],
            [
                date_range("2020-01-01", periods=1, freq="D"),
                date_range("2020-01-02", periods=1, freq="D"),
                date_range("2020-01-03", periods=1, freq="D"),
                date_range("2020-01-04", periods=1, freq="D"),
                date_range("2020-01-05", periods=1, freq="D"),
            ],
            "1D",
        ),
        (
            [[0], [0, 1], [1, 2], [2, 3], [3, 4]],
            [
                date_range("2020-01-01", periods=1, freq="D"),
                date_range("2020-01-01", periods=2, freq="D"),
                date_range("2020-01-02", periods=2, freq="D"),
                date_range("2020-01-03", periods=2, freq="D"),
                date_range("2020-01-04", periods=2, freq="D"),
            ],
            "2D",
        ),
        (
            [[0], [0, 1], [0, 1, 2], [1, 2, 3], [2, 3, 4]],
            [
                date_range("2020-01-01", periods=1, freq="D"),
                date_range("2020-01-01", periods=2, freq="D"),
                date_range("2020-01-01", periods=3, freq="D"),
                date_range("2020-01-02", periods=3, freq="D"),
                date_range("2020-01-03", periods=3, freq="D"),
            ],
            "3D",
        ),
    ],
)
def test_iter_rolling_datetime(expected, expected_index, window):
    # GH 11704
    ser = Series(range(5), index=date_range(start="2020-01-01", periods=5, freq="D"))

    expected = [
        Series(values, index=idx) for (values, idx) in zip(expected, expected_index)
    ]

    for (expected, actual) in zip(expected, ser.rolling(window)):
        tm.assert_series_equal(actual, expected)


@pytest.mark.parametrize(
    "grouping,_index",
    [
        (
            {"level": 0},
            pd.MultiIndex.from_tuples(
                [(0, 0), (0, 0), (1, 1), (1, 1), (1, 1)], names=[None, None]
            ),
        ),
        (
            {"by": "X"},
            pd.MultiIndex.from_tuples(
                [(0, 0), (1, 0), (2, 1), (3, 1), (4, 1)], names=["X", None]
            ),
        ),
    ],
)
def test_rolling_positional_argument(grouping, _index, raw):
    # GH 34605

    def scaled_sum(*args):
        if len(args) < 2:
            raise ValueError("The function needs two arguments")
        array, scale = args
        return array.sum() / scale

    df = DataFrame(data={"X": range(5)}, index=[0, 0, 1, 1, 1])

    expected = DataFrame(data={"X": [0.0, 0.5, 1.0, 1.5, 2.0]}, index=_index)
    result = df.groupby(**grouping).rolling(1).apply(scaled_sum, raw=raw, args=(2,))
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("add", [0.0, 2.0])
def test_rolling_numerical_accuracy_kahan_mean(add):
    # GH: 36031 implementing kahan summation
    df = DataFrame(
        {"A": [3002399751580331.0 + add, -0.0, -0.0]},
        index=[
            pd.Timestamp("19700101 09:00:00"),
            pd.Timestamp("19700101 09:00:03"),
            pd.Timestamp("19700101 09:00:06"),
        ],
    )
    result = (
        df.resample("1s").ffill().rolling("3s", closed="left", min_periods=3).mean()
    )
    dates = pd.date_range("19700101 09:00:00", periods=7, freq="S")
    expected = DataFrame(
        {
            "A": [
                np.nan,
                np.nan,
                np.nan,
                3002399751580330.5,
                2001599834386887.25,
                1000799917193443.625,
                0.0,
            ]
        },
        index=dates,
    )
    tm.assert_frame_equal(result, expected)


def test_rolling_numerical_accuracy_kahan_sum():
    # GH: 13254
    df = DataFrame([2.186, -1.647, 0.0, 0.0, 0.0, 0.0], columns=["x"])
    result = df["x"].rolling(3).sum()
    expected = Series([np.nan, np.nan, 0.539, -1.647, 0.0, 0.0], name="x")
    tm.assert_series_equal(result, expected)


def test_rolling_numerical_accuracy_jump():
    # GH: 32761
    index = pd.date_range(start="2020-01-01", end="2020-01-02", freq="60s").append(
        pd.DatetimeIndex(["2020-01-03"])
    )
    data = np.random.rand(len(index))

    df = DataFrame({"data": data}, index=index)
    result = df.rolling("60s").mean()
    tm.assert_frame_equal(result, df[["data"]])


def test_rolling_numerical_accuracy_small_values():
    # GH: 10319
    s = Series(
        data=[0.00012456, 0.0003, -0.0, -0.0],
        index=date_range("1999-02-03", "1999-02-06"),
    )
    result = s.rolling(1).mean()
    tm.assert_series_equal(result, s)


def test_rolling_numerical_too_large_numbers():
    # GH: 11645
    dates = pd.date_range("2015-01-01", periods=10, freq="D")
    ds = Series(data=range(10), index=dates, dtype=np.float64)
    ds[2] = -9e33
    result = ds.rolling(5).mean()
    expected = Series(
        [np.nan, np.nan, np.nan, np.nan, -1.8e33, -1.8e33, -1.8e33, 5.0, 6.0, 7.0],
        index=dates,
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    ("func", "value"),
    [("sum", 2.0), ("max", 1.0), ("min", 1.0), ("mean", 1.0), ("median", 1.0)],
)
def test_rolling_mixed_dtypes_axis_1(func, value):
    # GH: 20649
    df = DataFrame(1, index=[1, 2], columns=["a", "b", "c"])
    df["c"] = 1.0
    result = getattr(df.rolling(window=2, min_periods=1, axis=1), func)()
    expected = DataFrame(
        {"a": [1.0, 1.0], "b": [value, value], "c": [value, value]}, index=[1, 2]
    )
    tm.assert_frame_equal(result, expected)


def test_rolling_axis_one_with_nan():
    # GH: 35596
    df = DataFrame(
        [
            [0, 1, 2, 4, np.nan, np.nan, np.nan],
            [0, 1, 2, np.nan, np.nan, np.nan, np.nan],
            [0, 2, 2, np.nan, 2, np.nan, 1],
        ]
    )
    result = df.rolling(window=7, min_periods=1, axis="columns").sum()
    expected = DataFrame(
        [
            [0.0, 1.0, 3.0, 7.0, 7.0, 7.0, 7.0],
            [0.0, 1.0, 3.0, 3.0, 3.0, 3.0, 3.0],
            [0.0, 2.0, 4.0, 4.0, 6.0, 6.0, 7.0],
        ]
    )
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "value",
    ["test", pd.to_datetime("2019-12-31"), pd.to_timedelta("1 days 06:05:01.00003")],
)
def test_rolling_axis_1_non_numeric_dtypes(value):
    # GH: 20649
    df = DataFrame({"a": [1, 2]})
    df["b"] = value
    result = df.rolling(window=2, min_periods=1, axis=1).sum()
    expected = DataFrame({"a": [1.0, 2.0]})
    tm.assert_frame_equal(result, expected)


def test_rolling_on_df_transposed():
    # GH: 32724
    df = DataFrame({"A": [1, None], "B": [4, 5], "C": [7, 8]})
    expected = DataFrame({"A": [1.0, np.nan], "B": [5.0, 5.0], "C": [11.0, 13.0]})
    result = df.rolling(min_periods=1, window=2, axis=1).sum()
    tm.assert_frame_equal(result, expected)

    result = df.T.rolling(min_periods=1, window=2).sum().T
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    ("index", "window"),
    [
        (
            pd.period_range(start="2020-01-01 08:00", end="2020-01-01 08:08", freq="T"),
            "2T",
        ),
        (
            pd.period_range(
                start="2020-01-01 08:00", end="2020-01-01 12:00", freq="30T"
            ),
            "1h",
        ),
    ],
)
@pytest.mark.parametrize(
    ("func", "values"),
    [
        ("min", [np.nan, 0, 0, 1, 2, 3, 4, 5, 6]),
        ("max", [np.nan, 0, 1, 2, 3, 4, 5, 6, 7]),
        ("sum", [np.nan, 0, 1, 3, 5, 7, 9, 11, 13]),
    ],
)
def test_rolling_period_index(index, window, func, values):
    # GH: 34225
    ds = Series([0, 1, 2, 3, 4, 5, 6, 7, 8], index=index)
    result = getattr(ds.rolling(window, closed="left"), func)()
    expected = Series(values, index=index)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("constructor", ["DataFrame", "Series"])
def test_rolling_sem(constructor):
    # GH: 26476
    obj = getattr(pd, constructor)([0, 1, 2])
    result = obj.rolling(2, min_periods=1).sem()
    if isinstance(result, DataFrame):
        result = Series(result[0].values)
    expected = Series([np.nan] + [0.707107] * 2)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    ("func", "third_value", "values"),
    [
        ("var", 1, [5e33, 0, 0.5, 0.5, 2, 0]),
        ("std", 1, [7.071068e16, 0, 0.7071068, 0.7071068, 1.414214, 0]),
        ("var", 2, [5e33, 0.5, 0, 0.5, 2, 0]),
        ("std", 2, [7.071068e16, 0.7071068, 0, 0.7071068, 1.414214, 0]),
    ],
)
def test_rolling_var_numerical_issues(func, third_value, values):
    # GH: 37051
    ds = Series([99999999999999999, 1, third_value, 2, 3, 1, 1])
    result = getattr(ds.rolling(2), func)()
    expected = Series([np.nan] + values)
    tm.assert_series_equal(result, expected)
