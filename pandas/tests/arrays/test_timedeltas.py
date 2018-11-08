# -*- coding: utf-8 -*-

import pandas as pd
from pandas.core.arrays import TimedeltaArrayMixin as TimedeltaArray


class TestTimedeltaArray(object):
    def test_repr(self):
        tdi = pd.timedelta_range('1D', periods=10, freq='D')
        arr = TimedeltaArray(tdi)

        # non- truncated
        expected = (
            "<TimedeltaArrayMixin>\n"
            "['1 days 00:00:00', '2 days 00:00:00']\n"
            "Length: 2, dtype: timedelta64[ns], freq: None")
        result = repr(arr[:2])
        assert result == expected

        # truncated
        expected = (
            "<TimedeltaArrayMixin>\n"
            "["
            "'1 days 00:00:00', "
            "'2 days 00:00:00', "
            "'3 days 00:00:00', "
            "'...', "
            "'8 days 00:00:00', "
            "'9 days 00:00:00', "
            "'10 days 00:00:00"
            "']\n"
            "Length: 10, dtype: timedelta64[ns], freq: None")
        result = repr(arr)
        assert result == expected
