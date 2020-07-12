import pytest

import pandas.core.dtypes.concat as _concat

import pandas as pd
from pandas import DatetimeIndex, Period, PeriodIndex, Series, TimedeltaIndex
import pandas._testing as tm


@pytest.mark.parametrize(
    "to_concat, expected",
    [
        # int/float/str
        ([["a"], [1, 2]], ["i", "object"]),
        ([[3, 4], [1, 2]], ["i"]),
        ([[3, 4], [1, 2.1]], ["i", "f"]),
        # datetimelike
        ([DatetimeIndex(["2011-01-01"]), DatetimeIndex(["2011-01-02"])], ["datetime"]),
        ([TimedeltaIndex(["1 days"]), TimedeltaIndex(["2 days"])], ["timedelta"]),
        # datetimelike object
        (
            [
                DatetimeIndex(["2011-01-01"]),
                DatetimeIndex(["2011-01-02"], tz="US/Eastern"),
            ],
            ["datetime", "datetime64[ns, US/Eastern]"],
        ),
        (
            [
                DatetimeIndex(["2011-01-01"], tz="Asia/Tokyo"),
                DatetimeIndex(["2011-01-02"], tz="US/Eastern"),
            ],
            ["datetime64[ns, Asia/Tokyo]", "datetime64[ns, US/Eastern]"],
        ),
        ([TimedeltaIndex(["1 days"]), TimedeltaIndex(["2 hours"])], ["timedelta"]),
        (
            [
                DatetimeIndex(["2011-01-01"], tz="Asia/Tokyo"),
                TimedeltaIndex(["1 days"]),
            ],
            ["datetime64[ns, Asia/Tokyo]", "timedelta"],
        ),
    ],
)
def test_get_dtype_kinds(index_or_series, to_concat, expected):
    to_concat_klass = [index_or_series(c) for c in to_concat]
    result = _concat.get_dtype_kinds(to_concat_klass)
    assert result == set(expected)


@pytest.mark.parametrize(
    "to_concat, expected",
    [
        (
            [PeriodIndex(["2011-01"], freq="M"), PeriodIndex(["2011-01"], freq="M")],
            ["period[M]"],
        ),
        (
            [
                Series([Period("2011-01", freq="M")]),
                Series([Period("2011-02", freq="M")]),
            ],
            ["period[M]"],
        ),
        (
            [PeriodIndex(["2011-01"], freq="M"), PeriodIndex(["2011-01"], freq="D")],
            ["period[M]", "period[D]"],
        ),
        (
            [
                Series([Period("2011-01", freq="M")]),
                Series([Period("2011-02", freq="D")]),
            ],
            ["period[M]", "period[D]"],
        ),
    ],
)
def test_get_dtype_kinds_period(to_concat, expected):
    result = _concat.get_dtype_kinds(to_concat)
    assert result == set(expected)


def test_concat_mismatched_categoricals_with_empty():
    # concat_compat behavior on series._values should match pd.concat on series
    ser1 = Series(["a", "b", "c"], dtype="category")
    ser2 = Series([], dtype="category")

    result = _concat.concat_compat([ser1._values, ser2._values])
    expected = pd.concat([ser1, ser2])._values
    tm.assert_categorical_equal(result, expected)


@pytest.mark.parametrize("copy", [True, False])
def test_concat_single_dataframe_tz_aware(copy):
    # https://github.com/pandas-dev/pandas/issues/25257
    df = pd.DataFrame(
        {"timestamp": [pd.Timestamp("2020-04-08 09:00:00.709949+0000", tz="UTC")]}
    )
    expected = df.copy()
    result = pd.concat([df], copy=copy)
    tm.assert_frame_equal(result, expected)
