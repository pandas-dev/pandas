from .pandas_vb_common import *
import pandas.computation.expressions as expr


class Ops(object):
    goal_time = 0.2

    params = [[True, False], ['default', 1]]
    param_names = ['use_numexpr', 'threads']

    def setup(self, use_numexpr, threads):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))

        if threads != 'default':
            expr.set_numexpr_threads(threads)
        if not use_numexpr:
            expr.set_use_numexpr(False)


    def time_frame_add(self, use_numexpr, threads):
        (self.df + self.df2)

    def time_frame_mult(self, use_numexpr, threads):
        (self.df * self.df2)

    def time_frame_multi_and(self, use_numexpr, threads):
        self.df[((self.df > 0) & (self.df2 > 0))]

    def time_frame_comparison(self, use_numexpr, threads):
        (self.df > self.df2)

    def teardown(self, use_numexpr, threads):
        expr.set_use_numexpr(True)
        expr.set_numexpr_threads()


class Ops2(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 1000))
        self.df2 = DataFrame(np.random.randn(1000, 1000))

        self.df_int = DataFrame(
            np.random.random_integers(np.iinfo(np.int16).min,
                                      np.iinfo(np.int16).max,
                                      size=(1000, 1000)))
        self.df2_int = DataFrame(
            np.random.random_integers(np.iinfo(np.int16).min,
                                      np.iinfo(np.int16).max,
                                      size=(1000, 1000)))

    ## Division

    def time_frame_float_div(self):
        (self.df // self.df2)

    def time_frame_float_div_by_zero(self):
        (self.df / 0)

    def time_frame_float_floor_by_zero(self):
        (self.df // 0)

    def time_frame_int_div_by_zero(self):
        (self.df_int / 0)

    ## Modulo

    def time_frame_int_mod(self):
        (self.df / self.df2)

    def time_frame_float_mod(self):
        (self.df / self.df2)


class Timeseries(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.halfway = ((self.N // 2) - 1)
        self.s = Series(date_range('20010101', periods=self.N, freq='T'))
        self.ts = self.s[self.halfway]

        self.s2 = Series(date_range('20010101', periods=self.N, freq='s'))

    def time_series_timestamp_compare(self):
        (self.s <= self.ts)

    def time_timestamp_series_compare(self):
        (self.ts >= self.s)

    def time_timestamp_ops_diff1(self):
        self.s2.diff()

    def time_timestamp_ops_diff2(self):
        (self.s - self.s.shift())



class TimeseriesTZ(Timeseries):

    def setup(self):
        self.N = 1000000
        self.halfway = ((self.N // 2) - 1)
        self.s = Series(date_range('20010101', periods=self.N, freq='T', tz='US/Eastern'))
        self.ts = self.s[self.halfway]

        self.s2 = Series(date_range('20010101', periods=self.N, freq='s',  tz='US/Eastern'))