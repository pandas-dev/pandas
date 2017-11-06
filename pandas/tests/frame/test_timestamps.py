""" test DataFrame-boxed versions of the scalar Timestamp """
from datetime import datetime

import numpy as np

import pandas.util.testing as tm
from pandas.tseries import offsets

from pandas.compat import lrange
from pandas import (Timestamp, date_range,
                    Series, DataFrame, DatetimeIndex)


class TestDataFrameTimestamps(object):
    def test_series_map_box_timestamps(self):
        # GH#2689, GH#2627
        s = Series(date_range('1/1/2000', periods=10))

        def f(x):
            return (x.hour, x.day, x.month)

        # it works!
        s.map(f)
        s.apply(f)
        DataFrame(s).applymap(f)

    def test_frame_setitem_timestamp(self):
        # GH#2155
        columns = DatetimeIndex(start='1/1/2012', end='2/1/2012',
                                freq=offsets.BDay())
        index = lrange(10)
        data = DataFrame(columns=columns, index=index)
        t = datetime(2012, 11, 1)
        ts = Timestamp(t)
        data[ts] = np.nan  # works

    def test_to_html_timestamp(self):
        rng = date_range('2000-01-01', periods=10)
        df = DataFrame(np.random.randn(10, 4), index=rng)

        result = df.to_html()
        assert '2000-01-01' in result

    def test_compare_invalid(self):
        # GH 8058
        df = DataFrame(np.random.randn(5, 2))
        a = df[0]
        b = Series(np.random.randn(5))
        b.name = Timestamp('2000-01-01')
        tm.assert_series_equal(a / b, 1 / (b / a))
