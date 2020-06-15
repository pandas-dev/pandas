import numpy as np
import pytest

from pandas.errors import UnsupportedFunctionCall

from pandas import DataFrame, Series, date_range
from pandas.core.window import EWM
import pandas._testing as tm


def test_doc_string():

    df = DataFrame({"B": [0, 1, 2, np.nan, 4]})
    df
    df.ewm(com=0.5).mean()


def test_constructor(which):

    c = which.ewm

    # valid
    c(com=0.5)
    c(span=1.5)
    c(alpha=0.5)
    c(halflife=0.75)
    c(com=0.5, span=None)
    c(alpha=0.5, com=None)
    c(halflife=0.75, alpha=None)

    # not valid: mutually exclusive
    with pytest.raises(ValueError):
        c(com=0.5, alpha=0.5)
    with pytest.raises(ValueError):
        c(span=1.5, halflife=0.75)
    with pytest.raises(ValueError):
        c(alpha=0.5, span=1.5)

    # not valid: com < 0
    with pytest.raises(ValueError):
        c(com=-0.5)

    # not valid: span < 1
    with pytest.raises(ValueError):
        c(span=0.5)

    # not valid: halflife <= 0
    with pytest.raises(ValueError):
        c(halflife=0)

    # not valid: alpha <= 0 or alpha > 1
    for alpha in (-0.5, 1.5):
        with pytest.raises(ValueError):
            c(alpha=alpha)


@pytest.mark.parametrize("method", ["std", "mean", "var"])
def test_numpy_compat(method):
    # see gh-12811
    e = EWM(Series([2, 4, 6]), alpha=0.5)

    msg = "numpy operations are not valid with window objects"

    with pytest.raises(UnsupportedFunctionCall, match=msg):
        getattr(e, method)(1, 2, 3)
    with pytest.raises(UnsupportedFunctionCall, match=msg):
        getattr(e, method)(dtype=np.float64)


def test_ewma_times_not_datetime_type():
    msg = "times must be datetime64\[ns\] dtype."
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


@pytest.mark.parametrize(
    "times",
    [
        np.arange(10).astype("datetime64[ns]"),
        date_range("2000", freq="ns", periods=10),
        date_range("2000", freq="ns", periods=10).tz_localize("UTC"),
        "time_col",
    ],
)
def test_ewma_with_times(halflife_with_times, times):
    halflife = halflife_with_times
    df = DataFrame(
        {"A": range(10), "time_col": date_range("2000", freq="ns", periods=10)}
    )
    result = df.ewm(halflife=halflife, times=times).mean()
    expected = None
    tm.assert_frame_equal(result, expected)
