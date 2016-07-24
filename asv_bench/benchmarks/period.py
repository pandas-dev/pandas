from pandas import Series, Period, PeriodIndex, date_range


class create_period_index_from_date_range(object):
    goal_time = 0.2

    def time_period_index(self):
        # Simulate irregular PeriodIndex
        PeriodIndex(date_range('1985', periods=1000).to_pydatetime(), freq='D')


class period_setitem(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.rng = date_range(start='1/1/2000', periods=self.N, freq='T')
        if hasattr(Series, 'convert'):
            Series.resample = Series.convert
        self.ts = Series(np.random.randn(self.N), index=self.rng)
        self.rng = period_range(start='1/1/1990', freq='S', periods=20000)
        self.df = DataFrame(index=range(len(self.rng)))

    def time_period_setitem(self):
        self.df['col'] = self.rng


class period_algorithm(object):
    goal_time = 0.2

    def setup(self):
        data = [Period('2011-01', freq='M'), Period('2011-02', freq='M'),
                Period('2011-03', freq='M'), Period('2011-04', freq='M')]
        self.s = Series(data * 1000)
        self.i = PeriodIndex(data, freq='M')

    def time_period_series_drop_duplicates(self):
        self.s.drop_duplicates()

    def time_period_index_drop_duplicates(self):
        self.i.drop_duplicates()

    def time_period_series_value_counts(self):
        self.s.value_counts()

    def time_period_index_value_counts(self):
        self.i.value_counts()


