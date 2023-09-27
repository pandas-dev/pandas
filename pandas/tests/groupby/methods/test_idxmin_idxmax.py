import numpy as np
import pytest

import pandas as pd
from pandas import (
    DataFrame,
    Index,
    date_range,
)
import pandas._testing as tm


@pytest.mark.parametrize(
    "func, values",
    [
        ("idxmin", {"c_int": [0, 2], "c_float": [1, 3], "c_date": [1, 2]}),
        ("idxmax", {"c_int": [1, 3], "c_float": [0, 2], "c_date": [0, 3]}),
    ],
)
@pytest.mark.parametrize("numeric_only", [True, False])
def test_idxmin_idxmax_returns_int_types(func, values, numeric_only):
    # GH 25444
    df = DataFrame(
        {
            "name": ["A", "A", "B", "B"],
            "c_int": [1, 2, 3, 4],
            "c_float": [4.02, 3.03, 2.04, 1.05],
            "c_date": ["2019", "2018", "2016", "2017"],
        }
    )
    df["c_date"] = pd.to_datetime(df["c_date"])
    df["c_date_tz"] = df["c_date"].dt.tz_localize("US/Pacific")
    df["c_timedelta"] = df["c_date"] - df["c_date"].iloc[0]
    df["c_period"] = df["c_date"].dt.to_period("W")
    df["c_Integer"] = df["c_int"].astype("Int64")
    df["c_Floating"] = df["c_float"].astype("Float64")

    result = getattr(df.groupby("name"), func)(numeric_only=numeric_only)

    expected = DataFrame(values, index=Index(["A", "B"], name="name"))
    if numeric_only:
        expected = expected.drop(columns=["c_date"])
    else:
        expected["c_date_tz"] = expected["c_date"]
        expected["c_timedelta"] = expected["c_date"]
        expected["c_period"] = expected["c_date"]
    expected["c_Integer"] = expected["c_int"]
    expected["c_Floating"] = expected["c_float"]

    tm.assert_frame_equal(result, expected)


def test_idxmin_idxmax_axis1():
    df = DataFrame(
        np.random.default_rng(2).standard_normal((10, 4)), columns=["A", "B", "C", "D"]
    )
    df["A"] = [1, 2, 3, 1, 2, 3, 1, 2, 3, 4]

    gb = df.groupby("A")

    warn_msg = "DataFrameGroupBy.idxmax with axis=1 is deprecated"
    with tm.assert_produces_warning(FutureWarning, match=warn_msg):
        res = gb.idxmax(axis=1)

    alt = df.iloc[:, 1:].idxmax(axis=1)
    indexer = res.index.get_level_values(1)

    tm.assert_series_equal(alt[indexer], res.droplevel("A"))

    df["E"] = date_range("2016-01-01", periods=10)
    gb2 = df.groupby("A")

    msg = "'>' not supported between instances of 'Timestamp' and 'float'"
    with pytest.raises(TypeError, match=msg):
        with tm.assert_produces_warning(FutureWarning, match=warn_msg):
            gb2.idxmax(axis=1)
