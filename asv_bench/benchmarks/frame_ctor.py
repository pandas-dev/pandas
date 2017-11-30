import numpy as np
import pandas.util.testing as tm
from pandas import DataFrame, Series, MultiIndex, Timestamp, date_range
try:
    from pandas.tseries import offsets
except:
    from pandas.core.datetools import * # noqa

from .pandas_vb_common import setup # noqa


class FromDicts(object):

    goal_time = 0.2

    def setup(self):
        N, K = 5000, 50
        index = tm.makeStringIndex(N)
        columns = tm.makeStringIndex(K)
        frame = DataFrame(np.random.randn(N, K), index=index, columns=columns)
        self.data = frame.to_dict()
        self.some_dict = list(self.data.values())[0]
        self.dict_list = frame.to_dict(orient='records')
        self.data2 = {i: {j: float(j) for j in range(100)}
                      for i in range(2000)}

    def time_frame_ctor_list_of_dict(self):
        DataFrame(self.dict_list)

    def time_frame_ctor_nested_dict(self):
        DataFrame(self.data)

    def time_series_ctor_from_dict(self):
        Series(self.some_dict)

    def time_frame_ctor_nested_dict_int64(self):
        # nested dict, integer indexes, regression described in #621
        DataFrame(self.data2)


class FromSeries(object):

    goal_time = 0.2

    def setup(self):
        mi = MultiIndex.from_product([range(100), range(100)])
        self.s = Series(np.random.randn(10000), index=mi)

    def time_frame_from_mi_series(self):
        DataFrame(self.s)

# ----------------------------------------------------------------------
# From dict with DatetimeIndex with all offsets

# dynamically generate benchmarks for every offset
#
# get_period_count & get_index_for_offset are there because blindly taking each
# offset times 1000 can easily go out of Timestamp bounds and raise errors.


def get_period_count(start_date, off):
    ten_offsets_in_days = ((start_date + (off * 10)) - start_date).days
    if (ten_offsets_in_days == 0):
        return 1000
    else:
        periods = 9 * (Timestamp.max - start_date).days // ten_offsets_in_days
        return min(periods, 1000)


def get_index_for_offset(off):
    start_date = Timestamp('1/1/1900')
    return date_range(start_date,
                      periods=get_period_count(start_date, off),
                      freq=off)


all_offsets = offsets.__all__
# extra cases
for off in ['FY5253', 'FY5253Quarter']:
    all_offsets.pop(all_offsets.index(off))
    all_offsets.extend([off + '_1', off + '_2'])


class FromDictwithTimestampOffsets(object):

    params = [all_offsets, [1, 2]]
    param_names = ['offset', 'n_steps']

    offset_kwargs = {'WeekOfMonth': {'weekday': 1, 'week': 1},
                     'LastWeekOfMonth': {'weekday': 1, 'week': 1},
                     'FY5253': {'startingMonth': 1, 'weekday': 1},
                     'FY5253Quarter': {'qtr_with_extra_week': 1,
                                       'startingMonth': 1,
                                       'weekday': 1}}

    offset_extra_cases = {'FY5253': {'variation': ['nearest', 'last']},
                          'FY5253Quarter': {'variation': ['nearest', 'last']}}

    def setup(self, offset, n_steps):
        np.random.seed(1234)
        extra = False
        if offset.endswith("_", None, -1):
            extra = int(offset[-1])
            offset = offset[:-2]

        kwargs = {}
        if offset in self.offset_kwargs:
            kwargs = self.offset_kwargs[offset]

        if extra:
            extras = self.offset_extra_cases[offset]
            for extra_arg in extras:
                kwargs[extra_arg] = extras[extra_arg][extra - 1]

        offset = getattr(offsets, offset)
        self.idx = get_index_for_offset(offset(n_steps, **kwargs))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = self.df.to_dict()

    def time_frame_ctor(self, offset, n_steps):
        DataFrame(self.d)


class FromRecords(object):

    goal_time = 0.2
    params = [None, 1000]
    param_names = ['nrows']

    def setup(self, nrows):
        N = 100000
        self.gen = ((x, (x * 20), (x * 100)) for x in range(N))

    def time_frame_from_records_generator(self, nrows):
        # issue-6700
        self.df = DataFrame.from_records(self.gen, nrows=nrows)
