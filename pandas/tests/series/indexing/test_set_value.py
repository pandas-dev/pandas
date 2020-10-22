from datetime import datetime

import numpy as np

from pandas import DatetimeIndex, Series
import pandas._testing as tm


def test_series_set_value():
    # GH#1561

    dates = [datetime(2001, 1, 1), datetime(2001, 1, 2)]
    index = DatetimeIndex(dates)

    s = Series(dtype=object)
    s._set_value(dates[0], 1.0)
    s._set_value(dates[1], np.nan)

    expected = Series([1.0, np.nan], index=index)

    tm.assert_series_equal(s, expected)
