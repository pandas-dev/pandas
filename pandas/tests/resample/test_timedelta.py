import numpy as np

import pandas as pd
from pandas import DataFrame
from pandas.core.indexes.timedeltas import timedelta_range
import pandas.util.testing as tm
from pandas.util.testing import assert_frame_equal


class TestTimedeltaIndex(object):
    def test_asfreq_bug(self):
        import datetime as dt
        df = DataFrame(data=[1, 3],
                       index=[dt.timedelta(), dt.timedelta(minutes=3)])
        result = df.resample('1T').asfreq()
        expected = DataFrame(data=[1, np.nan, np.nan, 3],
                             index=timedelta_range('0 day',
                                                   periods=4,
                                                   freq='1T'))
        assert_frame_equal(result, expected)

    def test_resample_with_nat(self):
        # GH 13223
        index = pd.to_timedelta(['0s', pd.NaT, '2s'])
        result = DataFrame({'value': [2, 3, 5]}, index).resample('1s').mean()
        expected = DataFrame({'value': [2.5, np.nan, 5.0]},
                             index=timedelta_range('0 day',
                                                   periods=3,
                                                   freq='1S'))
        assert_frame_equal(result, expected)

    def test_resample_as_freq_with_subperiod(self):
        # GH 13022
        index = timedelta_range('00:00:00', '00:10:00', freq='5T')
        df = DataFrame(data={'value': [1, 5, 10]}, index=index)
        result = df.resample('2T').asfreq()
        expected_data = {'value': [1, np.nan, np.nan, np.nan, np.nan, 10]}
        expected = DataFrame(data=expected_data,
                             index=timedelta_range('00:00:00',
                                                   '00:10:00', freq='2T'))
        tm.assert_frame_equal(result, expected)
