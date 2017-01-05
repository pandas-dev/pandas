import pandas as pd
from pandas import Series, Period, PeriodIndex, date_range


class Constructor(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1985', periods=1000)
        self.rng2 = date_range('1985', periods=1000).to_pydatetime()

    def time_from_date_range(self):
        PeriodIndex(self.rng, freq='D')

    def time_from_pydatetime(self):
        PeriodIndex(self.rng2, freq='D')


class DataFrame(object):
    goal_time = 0.2

    def setup(self):
        self.rng = pd.period_range(start='1/1/1990', freq='S', periods=20000)
        self.df = pd.DataFrame(index=range(len(self.rng)))

    def time_setitem_period_column(self):
        self.df['col'] = self.rng


class Algorithms(object):
    goal_time = 0.2

    def setup(self):
        data = [Period('2011-01', freq='M'), Period('2011-02', freq='M'),
                Period('2011-03', freq='M'), Period('2011-04', freq='M')]
        self.s = Series(data * 1000)
        self.i = PeriodIndex(data, freq='M')

    def time_drop_duplicates_pseries(self):
        self.s.drop_duplicates()

    def time_drop_duplicates_pindex(self):
        self.i.drop_duplicates()

    def time_value_counts_pseries(self):
        self.s.value_counts()

    def time_value_counts_pindex(self):
        self.i.value_counts()


class period_standard_indexing(object):
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
