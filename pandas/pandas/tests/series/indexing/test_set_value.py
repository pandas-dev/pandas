from datetime import datetime

import numpy as np

from pandas import (
    DatetimeIndex,
    Index,
    Series,
)
import pandas._testing as tm


def test_series_set_value(using_infer_string):
    # GH#1561, GH#51363 as of 3.0 we do not do inference in Index.insert

    dates = [datetime(2001, 1, 1), datetime(2001, 1, 2)]
    if using_infer_string:
        index = DatetimeIndex(dates)
    else:
        index = Index(dates, dtype=object)

    s = Series(dtype=object)
    s._set_value(dates[0], 1.0)
    s._set_value(dates[1], np.nan)

    expected = Series([1.0, np.nan], index=index)

    tm.assert_series_equal(s, expected)


def test_set_value_dt64(datetime_series):
    idx = datetime_series.index[10]
    res = datetime_series._set_value(idx, 0)
    assert res is None
    assert datetime_series[idx] == 0


def test_set_value_str_index(string_series):
    # equiv
    ser = string_series.copy()
    res = ser._set_value("foobar", 0)
    assert res is None
    assert ser.index[-1] == "foobar"
    assert ser["foobar"] == 0

    ser2 = string_series.copy()
    ser2.loc["foobar"] = 0
    assert ser2.index[-1] == "foobar"
    assert ser2["foobar"] == 0
