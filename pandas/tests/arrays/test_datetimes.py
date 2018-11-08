# -*- coding: utf-8 -*-

import pandas as pd
from pandas.core.arrays import DatetimeArrayMixin as DatetimeArray


class TestDatetimeArray(object):
    def test_repr(self):
        dti = pd.date_range('1994-07-01', periods=9, freq='W', tz='US/Central')
        arr = DatetimeArray(dti)

        # non-truncated
        expected = (
            "<DatetimeArrayMixin>\n"
            "['1994-07-03 00:00:00-05:00', '1994-07-10 00:00:00-05:00']\n"
            "Length: 2, dtype: datetime64[ns, US/Central], freq: W-SUN")
        result = repr(arr[:2])
        assert result == expected

        # truncated
        expected = (
            "<DatetimeArrayMixin>\n"
            "["
            "'1994-07-03 00:00:00-05:00', "
            "'1994-07-10 00:00:00-05:00', "
            "'1994-07-17 00:00:00-05:00', "
            "'...', "
            "'1994-08-14 00:00:00-05:00', "
            "'1994-08-21 00:00:00-05:00', "
            "'1994-08-28 00:00:00-05:00'"
            "]\n"
            "Length: 9, dtype: datetime64[ns, US/Central], freq: W-SUN")
        result = repr(arr)
        assert result == expected
