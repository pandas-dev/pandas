# -*- coding: utf-8 -*-
import pytest

import pandas as pd
from pandas.core.arrays import TimedeltaArrayMixin as TimedeltaArray


class TestTimedeltaArray(object):

    @pytest.mark.xfail(reason='DatetimeArray', strict=True)
    def test_repr(self):
        tdi = pd.timedelta_range('1D', periods=2, freq='D')
        arr = TimedeltaArray(tdi)

        # non- truncated
        expected = (
            "<TimedeltaArray>\n"
            "['1 days 00:00:00', '2 days 00:00:00']\n"
            "Length: 2, dtype: timedelta64[ns], freq: None")
        result = repr(arr)
        assert result == expected

    @pytest.mark.xfail(reason='DatetimeArray', strict=True)
    def test_repr_truncated(self):
        # truncated
        tdi = pd.timedelta_range('1D', periods=1000, freq='D')
        arr = TimedeltaArray(tdi)

        expected = (
            "<TimedeltaArray>\n"
            "['1 days 00:00:00', '2 days 00:00:00', "
            " '3 days 00:00:00', "
            " '...', "
            " '8 days 00:00:00', "
            " '9 days 00:00:00', "
            " '10 days 00:00:00"
            " ']\n"
            "Length: 10, dtype: timedelta64[ns], freq: None")
        result = repr(arr)
        assert result == expected
