from pandas import (DataFrame, Series, Period, PeriodIndex, date_range,
                    period_range)


class PeriodProperties(object):

    params = (['M', 'min'],
              ['year', 'month', 'day', 'hour', 'minute', 'second',
               'is_leap_year', 'quarter', 'qyear', 'week', 'daysinmonth',
               'dayofweek', 'dayofyear', 'start_time', 'end_time'])
    param_names = ['freq', 'attr']

    def setup(self, freq, attr):
        self.per = Period('2012-06-01', freq=freq)

    def time_property(self, freq, attr):
        getattr(self.per, attr)


class PeriodUnaryMethods(object):

    params = ['M', 'min']
    param_names = ['freq']

    def setup(self, freq):
        self.per = Period('2012-06-01', freq=freq)

    def time_to_timestamp(self, freq):
        self.per.to_timestamp()

    def time_now(self, freq):
        self.per.now(freq)

    def time_asfreq(self, freq):
        self.per.asfreq('A')


class PeriodIndexConstructor(object):

    goal_time = 0.2

    params = ['D']
    param_names = ['freq']

    def setup(self, freq):
        self.rng = date_range('1985', periods=1000)
        self.rng2 = date_range('1985', periods=1000).to_pydatetime()

    def time_from_date_range(self, freq):
        PeriodIndex(self.rng, freq=freq)

    def time_from_pydatetime(self, freq):
        PeriodIndex(self.rng2, freq=freq)


class DataFramePeriodColumn(object):

    goal_time = 0.2

    def setup(self):
        self.rng = period_range(start='1/1/1990', freq='S', periods=20000)
        self.df = DataFrame(index=range(len(self.rng)))

    def time_setitem_period_column(self):
        self.df['col'] = self.rng

    def time_set_index(self):
        # GH#21582 limited by comparisons of Period objects
        self.df['col2'] = self.rng
        self.df.set_index('col2', append=True)


class Algorithms(object):

    goal_time = 0.2

    params = ['index', 'series']
    param_names = ['typ']

    def setup(self, typ):
        data = [Period('2011-01', freq='M'), Period('2011-02', freq='M'),
                Period('2011-03', freq='M'), Period('2011-04', freq='M')]

        if typ == 'index':
            self.vector = PeriodIndex(data * 1000, freq='M')
        elif typ == 'series':
            self.vector = Series(data * 1000)

    def time_drop_duplicates(self, typ):
        self.vector.drop_duplicates()

    def time_value_counts(self, typ):
        self.vector.value_counts()


class Indexing(object):

    goal_time = 0.2
    pi_with_nas = PeriodIndex(['1985Q1', 'NaT', '1985Q2'] * 1000, freq='Q')
    pi_diff = PeriodIndex(['1985Q1', 'NaT', '1985Q3'] * 1000, freq='Q')

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
        DataFrame({'a': self.series, 'b': self.series[:500]})

    def time_intersection(self):
        self.index[:750].intersection(self.index[250:])

    def time_unique(self):
        self.index.unique()

    def time_dropna(self):
        self.pi_with_nas.dropna()

    def time_difference(self):
        self.pi_with_nas.difference(self.pi_diff)

    def time_symmetric_difference(self):
        self.pi_with_nas.symmetric_difference(self.pi_diff)
