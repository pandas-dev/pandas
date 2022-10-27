"""
Tests for yearandweek in isocalendar function in arrays\datetime.py
"""

from pandas.testing import assert_series_equal
import pandas as pd
from pandas import DataFrame

class TestIsoCalender:
    def test_isocalendar(self):
        # GH #48947

        idx = pd.date_range(start='2019-12-29', freq='D',periods=1)
        # isocalendar().yearandweek return value - type: <class 'pandas.core.series.Series'>
        resType = type(idx.isocalendar().yearandweek)

        # get the value of isocalendar().week using list for following value comparison
        res1 = idx.isocalendar().yearandweek.values.tolist()

        # expected value
        # res = DataFrame(data=[[201952]], index = idx, columns=['yearandweek'],dtype="UInt32")
        expextedValue = [201952]
        falseFormat = [2019052]

        # assert_frame_equal(res1,res)
        # assert_series_equal(res1,res)
        assert res1 == expextedValue   # test for return values
        assert resType is pd.Series  # test for return type
        assert not res1 == falseFormat # test for return format