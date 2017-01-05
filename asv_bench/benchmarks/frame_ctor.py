from .pandas_vb_common import *
try:
    from pandas.tseries.offsets import *
except:
    from pandas.core.datetools import *


#----------------------------------------------------------------------
# Creation from nested dict

class FromDicts(object):
    goal_time = 0.2

    def setup(self):
        (N, K) = (5000, 50)
        self.index = tm.makeStringIndex(N)
        self.columns = tm.makeStringIndex(K)
        self.frame = DataFrame(np.random.randn(N, K), index=self.index, columns=self.columns)
        try:
            self.data = self.frame.to_dict()
        except:
            self.data = self.frame.toDict()
        self.some_dict = self.data.values()[0]
        self.dict_list = [dict(zip(self.columns, row)) for row in self.frame.values]

        self.data2 = dict(
            ((i, dict(((j, float(j)) for j in range(100)))) for i in
             xrange(2000)))

    def time_frame_ctor_list_of_dict(self):
        DataFrame(self.dict_list)

    def time_frame_ctor_nested_dict(self):
        DataFrame(self.data)

    def time_series_ctor_from_dict(self):
        Series(self.some_dict)

    def time_frame_ctor_nested_dict_int64(self):
        # nested dict, integer indexes, regression described in #621
        DataFrame(self.data)


# from a mi-series

class frame_from_series(object):
    goal_time = 0.2

    def setup(self):
        self.mi = MultiIndex.from_tuples([(x, y) for x in range(100) for y in range(100)])
        self.s = Series(randn(10000), index=self.mi)

    def time_frame_from_mi_series(self):
        DataFrame(self.s)


#----------------------------------------------------------------------
# get_numeric_data

class frame_get_numeric_data(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 25))
        self.df['foo'] = 'bar'
        self.df['bar'] = 'baz'
        self.df = self.df.consolidate()

    def time_frame_get_numeric_data(self):
        self.df._get_numeric_data()


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
        return min((9 * ((Timestamp.max - start_date).days // ten_offsets_in_days)), 1000)


def get_index_for_offset(off):
    start_date = Timestamp('1/1/1900')
    return date_range(start_date, periods=min(1000, get_period_count(
        start_date, off)), freq=off)


all_offsets = offsets.__all__
# extra cases
for off in ['FY5253', 'FY5253Quarter']:
    all_offsets.pop(all_offsets.index(off))
    all_offsets.extend([off + '_1', off + '_2'])


class FrameConstructorDTIndexFromOffsets(object):

    params = [all_offsets, [1, 2]]
    param_names = ['offset', 'n_steps']

    offset_kwargs = {'WeekOfMonth': {'weekday': 1, 'week': 1},
                     'LastWeekOfMonth': {'weekday': 1, 'week': 1},
                     'FY5253': {'startingMonth': 1, 'weekday': 1},
                     'FY5253Quarter': {'qtr_with_extra_week': 1, 'startingMonth': 1, 'weekday': 1}}

    offset_extra_cases = {'FY5253': {'variation': ['nearest', 'last']},
                          'FY5253Quarter': {'variation': ['nearest', 'last']}}

    def setup(self, offset, n_steps):

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
                kwargs[extra_arg] = extras[extra_arg][extra -1]

        offset = getattr(offsets, offset)
        self.idx = get_index_for_offset(offset(n_steps, **kwargs))
        self.df = DataFrame(np.random.randn(len(self.idx), 10), index=self.idx)
        self.d = dict([(col, self.df[col]) for col in self.df.columns])

    def time_frame_ctor(self, offset, n_steps):
        DataFrame(self.d)
