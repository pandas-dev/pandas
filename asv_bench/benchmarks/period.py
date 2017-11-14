import pandas as pd
from pandas import Series, Period, PeriodIndex, date_range


class PeriodProperties(object):
    params = ['M', 'min']
    param_names = ['freq']

    def setup(self, freq):
        self.per = Period('2012-06-01', freq=freq)

    def time_year(self):
        self.per.year

    def time_month(self):
        self.per.month

    def time_day(self):
        self.per.day

    def time_hour(self):
        self.per.hour

    def time_minute(self):
        self.per.minute

    def time_second(self):
        self.per.second

    def time_is_leap_year(self):
        self.per.is_leap_year

    def time_quarter(self):
        self.per.quarter

    def time_qyear(self):
        self.per.qyear

    def time_week(self):
        self.per.week

    def time_daysinmonth(self):
        self.per.daysinmonth

    def time_dayofweek(self):
        self.per.dayofweek

    def time_dayofyear(self):
        self.per.dayofyear

    def time_start_time(self):
        self.per.start_time

    def time_end_time(self):
        self.per.end_time


class PeriodUnaryMethods(object):
    params = ['M', 'min']
    param_names = ['freq']

    def setup(self, freq):
        self.per = Period('2012-06-01', freq=freq)

    def time_to_timestamp(self):
        self.per.to_timestamp()

    def time_now(self):
        self.per.now()

    def time_asfreq(self):
        self.per.asfreq('A')


class PeriodIndexConstructor(object):
    goal_time = 0.2

    params = ['D']
    param_names = ['freq']

    def setup(self, freq):
        self.freq = freq
        self.rng = date_range('1985', periods=1000)
        self.rng2 = date_range('1985', periods=1000).to_pydatetime()

    def time_from_date_range(self):
        PeriodIndex(self.rng, freq=self.freq)

    def time_from_pydatetime(self):
        PeriodIndex(self.rng2, freq=self.freq)


class DataFramePeriodColumn(object):
    goal_time = 0.2

    def setup_cache(self):
        rng = pd.period_range(start='1/1/1990', freq='S', periods=20000)
        df = pd.DataFrame(index=range(len(rng)))
        return rng, df

    def time_setitem_period_column(self, tup):
        rng, df = tup
        df['col'] = rng


class Algorithms(object):
    goal_time = 0.2

    params = [PeriodIndex, Series]
    param_names = ['box_cls']

    def setup(self, box_cls):
        data = [Period('2011-01', freq='M'), Period('2011-02', freq='M'),
                Period('2011-03', freq='M'), Period('2011-04', freq='M')]

        if box_cls is PeriodIndex:
            self.vector = PeriodIndex(data, freq='M')
        else:
            self.vector = Series(data * 1000)

    def time_drop_duplicates(self):
        self.vector.drop_duplicates()

    def time_value_counts(self):
        self.vector.value_counts()


class PeriodStandardIndexing(object):
    goal_time = 0.2

    def setup(self):
        self.index = PeriodIndex(start='1985', periods=1000, freq='D')
        self.series = Series(range(1000), index=self.index)
        self.period = self.index[500]

    def time_get_loc(self):
        self.index.get_loc(self.period)

    def time_shape(self):
        self.index.shape

    def time_shallow_copy(self):
        self.index._shallow_copy()

    def time_series_loc(self):
        self.series.loc[self.period]

    def time_align(self):
        pd.DataFrame({'a': self.series, 'b': self.series[:500]})

    def time_intersection(self):
        self.index[:750].intersection(self.index[250:])
