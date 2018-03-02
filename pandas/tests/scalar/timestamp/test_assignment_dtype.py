import pandas as pd
from pandas.core.dtypes.dtypes import DatetimeTZDtype


# referencing #19843: scalar assignment of a tz-aware is object dtype

class TestTimestampProperties(object):

    def test_scalar_assignment(self):
        df = pd.DataFrame(index=(0, 1, 2))

        df['now'] = pd.Timestamp('20130101', tz='UTC')

        assert isinstance(df.dtypes[0], DatetimeTZDtype)

    def test_datetime_index_assignment(self):
        df = pd.DataFrame(index=(0, 1, 2))

        di = pd.DatetimeIndex(
            [pd.Timestamp('20130101', tz='UTC')]).repeat(len(df))
        df['now'] = di

        assert isinstance(df.dtypes[0], DatetimeTZDtype)
