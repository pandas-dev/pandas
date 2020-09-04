import numpy as np

import pandas as pd


def test_xs_datetimelike_wrapping():
    # GH#31630 a case where we shouldn't wrap datetime64 in Timestamp
    arr = pd.date_range("2016-01-01", periods=3)._data._data

    ser = pd.Series(arr, dtype=object)
    for i in range(len(ser)):
        ser.iloc[i] = arr[i]
    assert ser.dtype == object
    assert isinstance(ser[0], np.datetime64)

    result = ser.xs(0)
    assert isinstance(result, np.datetime64)
