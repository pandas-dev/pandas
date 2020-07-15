import numpy as np
import pytest

from pandas.errors import UnsupportedFunctionCall

import pandas as pd
from pandas import DataFrame, Series
import pandas._testing as tm
from pandas.core.window import Expanding


def test_doc_string():

    df = DataFrame({"B": [0, 1, 2, np.nan, 4]})
    df
    df.expanding(2).sum()


@pytest.mark.filterwarnings(
    "ignore:The `center` argument on `expanding` will be removed in the future"
)
def test_constructor(which):
    # GH 12669

    c = which.expanding

    # valid
    c(min_periods=1)
    c(min_periods=1, center=True)
    c(min_periods=1, center=False)

    # not valid
    for w in [2.0, "foo", np.array([2])]:
        msg = "min_periods must be an integer"
        with pytest.raises(ValueError, match=msg):
            c(min_periods=w)

        msg = "center must be a boolean"
        with pytest.raises(ValueError, match=msg):
            c(min_periods=1, center=w)


@pytest.mark.parametrize("method", ["std", "mean", "sum", "max", "min", "var"])
def test_numpy_compat(method):
    # see gh-12811
    e = Expanding(Series([2, 4, 6]), window=2)

    msg = "numpy operations are not valid with window objects"

    with pytest.raises(UnsupportedFunctionCall, match=msg):
        getattr(e, method)(1, 2, 3)
    with pytest.raises(UnsupportedFunctionCall, match=msg):
        getattr(e, method)(dtype=np.float64)


@pytest.mark.parametrize(
    "expander",
    [
        1,
        pytest.param(
            "ls",
            marks=pytest.mark.xfail(
                reason="GH#16425 expanding with offset not supported"
            ),
        ),
    ],
)
def test_empty_df_expanding(expander):
    # GH 15819 Verifies that datetime and integer expanding windows can be
    # applied to empty DataFrames

    expected = DataFrame()
    result = DataFrame().expanding(expander).sum()
    tm.assert_frame_equal(result, expected)

    # Verifies that datetime and integer expanding windows can be applied
    # to empty DataFrames with datetime index
    expected = DataFrame(index=pd.DatetimeIndex([]))
    result = DataFrame(index=pd.DatetimeIndex([])).expanding(expander).sum()
    tm.assert_frame_equal(result, expected)


def test_missing_minp_zero():
    # https://github.com/pandas-dev/pandas/pull/18921
    # minp=0
    x = pd.Series([np.nan])
    result = x.expanding(min_periods=0).sum()
    expected = pd.Series([0.0])
    tm.assert_series_equal(result, expected)

    # minp=1
    result = x.expanding(min_periods=1).sum()
    expected = pd.Series([np.nan])
    tm.assert_series_equal(result, expected)


def test_expanding_axis(axis_frame):
    # see gh-23372.
    df = DataFrame(np.ones((10, 20)))
    axis = df._get_axis_number(axis_frame)

    if axis == 0:
        expected = DataFrame(
            {i: [np.nan] * 2 + [float(j) for j in range(3, 11)] for i in range(20)}
        )
    else:
        # axis == 1
        expected = DataFrame([[np.nan] * 2 + [float(i) for i in range(3, 21)]] * 10)

    result = df.expanding(3, axis=axis_frame).sum()
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("constructor", [Series, DataFrame])
def test_expanding_count_with_min_periods(constructor):
    # GH 26996
    result = constructor(range(5)).expanding(min_periods=3).count()
    expected = constructor([np.nan, np.nan, 3.0, 4.0, 5.0])
    tm.assert_equal(result, expected)


@pytest.mark.parametrize("constructor", [Series, DataFrame])
def test_expanding_count_default_min_periods_with_null_values(constructor):
    # GH 26996
    values = [1, 2, 3, np.nan, 4, 5, 6]
    expected_counts = [1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 6.0]

    result = constructor(values).expanding().count()
    expected = constructor(expected_counts)
    tm.assert_equal(result, expected)


@pytest.mark.parametrize(
    "df,expected,min_periods",
    [
        (
            DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}),
            [
                ({"A": [1], "B": [4]}, [0]),
                ({"A": [1, 2], "B": [4, 5]}, [0, 1]),
                ({"A": [1, 2, 3], "B": [4, 5, 6]}, [0, 1, 2]),
            ],
            3,
        ),
        (
            DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}),
            [
                ({"A": [1], "B": [4]}, [0]),
                ({"A": [1, 2], "B": [4, 5]}, [0, 1]),
                ({"A": [1, 2, 3], "B": [4, 5, 6]}, [0, 1, 2]),
            ],
            2,
        ),
        (
            DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}),
            [
                ({"A": [1], "B": [4]}, [0]),
                ({"A": [1, 2], "B": [4, 5]}, [0, 1]),
                ({"A": [1, 2, 3], "B": [4, 5, 6]}, [0, 1, 2]),
            ],
            1,
        ),
        (DataFrame({"A": [1], "B": [4]}), [], 2),
        (DataFrame(), [({}, [])], 1),
        (
            DataFrame({"A": [1, np.nan, 3], "B": [np.nan, 5, 6]}),
            [
                ({"A": [1.0], "B": [np.nan]}, [0]),
                ({"A": [1, np.nan], "B": [np.nan, 5]}, [0, 1]),
                ({"A": [1, np.nan, 3], "B": [np.nan, 5, 6]}, [0, 1, 2]),
            ],
            3,
        ),
        (
            DataFrame({"A": [1, np.nan, 3], "B": [np.nan, 5, 6]}),
            [
                ({"A": [1.0], "B": [np.nan]}, [0]),
                ({"A": [1, np.nan], "B": [np.nan, 5]}, [0, 1]),
                ({"A": [1, np.nan, 3], "B": [np.nan, 5, 6]}, [0, 1, 2]),
            ],
            2,
        ),
        (
            DataFrame({"A": [1, np.nan, 3], "B": [np.nan, 5, 6]}),
            [
                ({"A": [1.0], "B": [np.nan]}, [0]),
                ({"A": [1, np.nan], "B": [np.nan, 5]}, [0, 1]),
                ({"A": [1, np.nan, 3], "B": [np.nan, 5, 6]}, [0, 1, 2]),
            ],
            1,
        ),
    ],
)
def test_iter_expanding_dataframe(df, expected, min_periods):
    # GH 11704
    expected = [DataFrame(values, index=index) for (values, index) in expected]

    for (expected, actual) in zip(expected, df.expanding(min_periods)):
        tm.assert_frame_equal(actual, expected)


@pytest.mark.parametrize(
    "ser,expected,min_periods",
    [
        (Series([1, 2, 3]), [([1], [0]), ([1, 2], [0, 1]), ([1, 2, 3], [0, 1, 2])], 3),
        (Series([1, 2, 3]), [([1], [0]), ([1, 2], [0, 1]), ([1, 2, 3], [0, 1, 2])], 2),
        (Series([1, 2, 3]), [([1], [0]), ([1, 2], [0, 1]), ([1, 2, 3], [0, 1, 2])], 1),
        (Series([1, 2]), [([1], [0]), ([1, 2], [0, 1])], 2),
        (Series([np.nan, 2]), [([np.nan], [0]), ([np.nan, 2], [0, 1])], 2),
        (Series([], dtype="int64"), [], 2),
    ],
)
def test_iter_expanding_series(ser, expected, min_periods):
    # GH 11704
    expected = [Series(values, index=index) for (values, index) in expected]

    for (expected, actual) in zip(expected, ser.expanding(min_periods)):
        tm.assert_series_equal(actual, expected)


def test_center_deprecate_warning():
    # GH 20647
    df = pd.DataFrame()
    with tm.assert_produces_warning(FutureWarning):
        df.expanding(center=True)

    with tm.assert_produces_warning(FutureWarning):
        df.expanding(center=False)

    with tm.assert_produces_warning(None):
        df.expanding()
