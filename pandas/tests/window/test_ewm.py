import numpy as np
import pytest

from pandas.errors import UnsupportedFunctionCall

from pandas import (
    DataFrame,
    DatetimeIndex,
    Series,
    date_range,
)
import pandas._testing as tm
from pandas.core.window import ExponentialMovingWindow


def test_doc_string():

    df = DataFrame({"B": [0, 1, 2, np.nan, 4]})
    df
    df.ewm(com=0.5).mean()


def test_constructor(frame_or_series):

    c = frame_or_series(range(5)).ewm

    # valid
    c(com=0.5)
    c(span=1.5)
    c(alpha=0.5)
    c(halflife=0.75)
    c(com=0.5, span=None)
    c(alpha=0.5, com=None)
    c(halflife=0.75, alpha=None)

    # not valid: mutually exclusive
    msg = "comass, span, halflife, and alpha are mutually exclusive"
    with pytest.raises(ValueError, match=msg):
        c(com=0.5, alpha=0.5)
    with pytest.raises(ValueError, match=msg):
        c(span=1.5, halflife=0.75)
    with pytest.raises(ValueError, match=msg):
        c(alpha=0.5, span=1.5)

    # not valid: com < 0
    msg = "comass must satisfy: comass >= 0"
    with pytest.raises(ValueError, match=msg):
        c(com=-0.5)

    # not valid: span < 1
    msg = "span must satisfy: span >= 1"
    with pytest.raises(ValueError, match=msg):
        c(span=0.5)

    # not valid: halflife <= 0
    msg = "halflife must satisfy: halflife > 0"
    with pytest.raises(ValueError, match=msg):
        c(halflife=0)

    # not valid: alpha <= 0 or alpha > 1
    msg = "alpha must satisfy: 0 < alpha <= 1"
    for alpha in (-0.5, 1.5):
        with pytest.raises(ValueError, match=msg):
            c(alpha=alpha)


@pytest.mark.parametrize("method", ["std", "mean", "var"])
def test_numpy_compat(method):
    # see gh-12811
    e = ExponentialMovingWindow(Series([2, 4, 6]), alpha=0.5)

    msg = "numpy operations are not valid with window objects"

    with pytest.raises(UnsupportedFunctionCall, match=msg):
        getattr(e, method)(1, 2, 3)
    with pytest.raises(UnsupportedFunctionCall, match=msg):
        getattr(e, method)(dtype=np.float64)


def test_ewma_times_not_datetime_type():
    msg = r"times must be datetime64\[ns\] dtype."
    with pytest.raises(ValueError, match=msg):
        Series(range(5)).ewm(times=np.arange(5))


def test_ewma_times_not_same_length():
    msg = "times must be the same length as the object."
    with pytest.raises(ValueError, match=msg):
        Series(range(5)).ewm(times=np.arange(4).astype("datetime64[ns]"))


def test_ewma_halflife_not_correct_type():
    msg = "halflife must be a string or datetime.timedelta object"
    with pytest.raises(ValueError, match=msg):
        Series(range(5)).ewm(halflife=1, times=np.arange(5).astype("datetime64[ns]"))


def test_ewma_halflife_without_times(halflife_with_times):
    msg = "halflife can only be a timedelta convertible argument if times is not None."
    with pytest.raises(ValueError, match=msg):
        Series(range(5)).ewm(halflife=halflife_with_times)


@pytest.mark.parametrize(
    "times",
    [
        np.arange(10).astype("datetime64[D]").astype("datetime64[ns]"),
        date_range("2000", freq="D", periods=10),
        date_range("2000", freq="D", periods=10).tz_localize("UTC"),
        "time_col",
    ],
)
@pytest.mark.parametrize("min_periods", [0, 2])
def test_ewma_with_times_equal_spacing(halflife_with_times, times, min_periods):
    halflife = halflife_with_times
    data = np.arange(10.0)
    data[::2] = np.nan
    df = DataFrame({"A": data, "time_col": date_range("2000", freq="D", periods=10)})
    with tm.assert_produces_warning(FutureWarning, match="nuisance columns"):
        # GH#42738
        result = df.ewm(halflife=halflife, min_periods=min_periods, times=times).mean()
        expected = df.ewm(halflife=1.0, min_periods=min_periods).mean()
    tm.assert_frame_equal(result, expected)


def test_ewma_with_times_variable_spacing(tz_aware_fixture):
    tz = tz_aware_fixture
    halflife = "23 days"
    times = DatetimeIndex(
        ["2020-01-01", "2020-01-10T00:04:05", "2020-02-23T05:00:23"]
    ).tz_localize(tz)
    data = np.arange(3)
    df = DataFrame(data)
    result = df.ewm(halflife=halflife, times=times).mean()
    expected = DataFrame([0.0, 0.5674161888241773, 1.545239952073459])
    tm.assert_frame_equal(result, expected)


def test_ewm_with_nat_raises(halflife_with_times):
    # GH#38535
    ser = Series(range(1))
    times = DatetimeIndex(["NaT"])
    with pytest.raises(ValueError, match="Cannot convert NaT values to integer"):
        ser.ewm(com=0.1, halflife=halflife_with_times, times=times)


def test_ewm_with_times_getitem(halflife_with_times):
    # GH 40164
    halflife = halflife_with_times
    data = np.arange(10.0)
    data[::2] = np.nan
    times = date_range("2000", freq="D", periods=10)
    df = DataFrame({"A": data, "B": data})
    result = df.ewm(halflife=halflife, times=times)["A"].mean()
    expected = df.ewm(halflife=1.0)["A"].mean()
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("arg", ["com", "halflife", "span", "alpha"])
def test_ewm_getitem_attributes_retained(arg, adjust, ignore_na):
    # GH 40164
    kwargs = {arg: 1, "adjust": adjust, "ignore_na": ignore_na}
    ewm = DataFrame({"A": range(1), "B": range(1)}).ewm(**kwargs)
    expected = {attr: getattr(ewm, attr) for attr in ewm._attributes}
    ewm_slice = ewm["A"]
    result = {attr: getattr(ewm, attr) for attr in ewm_slice._attributes}
    assert result == expected


def test_ewm_vol_deprecated():
    ser = Series(range(1))
    with tm.assert_produces_warning(FutureWarning):
        result = ser.ewm(com=0.1).vol()
    expected = ser.ewm(com=0.1).std()
    tm.assert_series_equal(result, expected)


def test_ewma_times_adjust_false_raises():
    # GH 40098
    with pytest.raises(
        NotImplementedError, match="times is not supported with adjust=False."
    ):
        Series(range(1)).ewm(
            0.1, adjust=False, times=date_range("2000", freq="D", periods=1)
        )


@pytest.mark.parametrize(
    "func, expected",
    [
        [
            "mean",
            DataFrame(
                {
                    0: range(5),
                    1: range(4, 9),
                    2: [7.428571, 9, 10.571429, 12.142857, 13.714286],
                },
                dtype=float,
            ),
        ],
        [
            "std",
            DataFrame(
                {
                    0: [np.nan] * 5,
                    1: [4.242641] * 5,
                    2: [4.6291, 5.196152, 5.781745, 6.380775, 6.989788],
                }
            ),
        ],
        [
            "var",
            DataFrame(
                {
                    0: [np.nan] * 5,
                    1: [18.0] * 5,
                    2: [21.428571, 27, 33.428571, 40.714286, 48.857143],
                }
            ),
        ],
    ],
)
def test_float_dtype_ewma(func, expected, float_numpy_dtype):
    # GH#42452

    df = DataFrame(
        {0: range(5), 1: range(6, 11), 2: range(10, 20, 2)}, dtype=float_numpy_dtype
    )
    e = df.ewm(alpha=0.5, axis=1)
    result = getattr(e, func)()

    tm.assert_frame_equal(result, expected)
