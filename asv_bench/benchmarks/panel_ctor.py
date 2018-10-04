import warnings
from datetime import datetime, timedelta

from pandas import DataFrame, Panel, DatetimeIndex, date_range

from .pandas_vb_common import setup  # noqa


class DifferentIndexes(object):
    goal_time = 0.2

    def setup(self):
        self.data_frames = {}
        start = datetime(1990, 1, 1)
        end = datetime(2012, 1, 1)
        for x in range(100):
            end += timedelta(days=1)
            idx = date_range(start, end)
            df = DataFrame({'a': 0, 'b': 1, 'c': 2}, index=idx)
            self.data_frames[x] = df

    def time_from_dict(self):
        with warnings.catch_warnings(record=True):
            Panel.from_dict(self.data_frames)


class SameIndexes(object):

    goal_time = 0.2

    def setup(self):
        idx = DatetimeIndex(start=datetime(1990, 1, 1),
                            end=datetime(2012, 1, 1),
                            freq='D')
        df = DataFrame({'a': 0, 'b': 1, 'c': 2}, index=idx)
        self.data_frames = dict(enumerate([df] * 100))

    def time_from_dict(self):
        with warnings.catch_warnings(record=True):
            Panel.from_dict(self.data_frames)


class TwoIndexes(object):

    goal_time = 0.2

    def setup(self):
        start = datetime(1990, 1, 1)
        end = datetime(2012, 1, 1)
        df1 = DataFrame({'a': 0, 'b': 1, 'c': 2},
                        index=DatetimeIndex(start=start, end=end, freq='D'))
        end += timedelta(days=1)
        df2 = DataFrame({'a': 0, 'b': 1, 'c': 2},
                        index=DatetimeIndex(start=start, end=end, freq='D'))
        dfs = [df1] * 50 + [df2] * 50
        self.data_frames = dict(enumerate(dfs))

    def time_from_dict(self):
        with warnings.catch_warnings(record=True):
            Panel.from_dict(self.data_frames)
