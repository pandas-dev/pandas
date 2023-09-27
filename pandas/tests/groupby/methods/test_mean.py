import pandas as pd
from pandas import (
    DataFrame,
    Index,
    Series,
)
import pandas._testing as tm


def test_groupby_mean_no_overflow():
    # Regression test for (#22487)
    df = DataFrame(
        {
            "user": ["A", "A", "A", "A", "A"],
            "connections": [4970, 4749, 4719, 4704, 18446744073699999744],
        }
    )
    assert df.groupby("user")["connections"].mean()["A"] == 3689348814740003840


def test_mean_on_timedelta():
    # GH 17382
    df = DataFrame({"time": pd.to_timedelta(range(10)), "cat": ["A", "B"] * 5})
    result = df.groupby("cat")["time"].mean()
    expected = Series(
        pd.to_timedelta([4, 5]), name="time", index=Index(["A", "B"], name="cat")
    )
    tm.assert_series_equal(result, expected)
