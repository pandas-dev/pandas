# pylint: disable=E1101

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series
from pandas.core.indexes.timedeltas import timedelta_range
from pandas.tests.resample.test_base import Base
import pandas.util.testing as tm
from pandas.util.testing import assert_frame_equal


class TestTimedeltaIndex(Base):
    _index_factory = lambda x: timedelta_range

    @pytest.fixture
    def _index_start(self):
        return '1 day'

    @pytest.fixture
    def _index_end(self):
        return '10 day'

    @pytest.fixture
    def _series_name(self):
        return 'tdi'

    def create_series(self):
        i = timedelta_range('1 day',
                            '10 day', freq='D')

        return Series(np.arange(len(i)), index=i, name='tdi')

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
