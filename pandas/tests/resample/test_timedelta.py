# pylint: disable=E1101

import numpy as np
import pytest

from pandas import DataFrame, Series
from pandas.core.indexes.timedeltas import timedelta_range
from pandas.util.testing import assert_frame_equal

from pandas.tests.resample.base import Base


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
