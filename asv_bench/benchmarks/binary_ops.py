import numpy as np
from pandas import DataFrame, Series, date_range
from pandas.core.algorithms import checked_add_with_arr
try:
    import pandas.core.computation.expressions as expr
except ImportError:
    import pandas.computation.expressions as expr

from .pandas_vb_common import setup # noqa


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
        self.df + self.df2

    def time_frame_mult(self, use_numexpr, threads):
        self.df * self.df2

    def time_frame_multi_and(self, use_numexpr, threads):
        self.df[(self.df > 0) & (self.df2 > 0)]

    def time_frame_comparison(self, use_numexpr, threads):
        self.df > self.df2

    def teardown(self, use_numexpr, threads):
        expr.set_use_numexpr(True)
        expr.set_numexpr_threads()


class Ops2(object):

    goal_time = 0.2

    def setup(self):
        N = 10**3
        self.df = DataFrame(np.random.randn(N, N))
        self.df2 = DataFrame(np.random.randn(N, N))

        self.df_int = DataFrame(np.random.randint(np.iinfo(np.int16).min,
                                                  np.iinfo(np.int16).max,
                                                  size=(N, N)))
        self.df2_int = DataFrame(np.random.randint(np.iinfo(np.int16).min,
                                                   np.iinfo(np.int16).max,
                                                   size=(N, N)))

    # Division

    def time_frame_float_div(self):
        self.df // self.df2

    def time_frame_float_div_by_zero(self):
        self.df / 0

    def time_frame_float_floor_by_zero(self):
        self.df // 0

    def time_frame_int_div_by_zero(self):
        self.df_int / 0

    # Modulo

    def time_frame_int_mod(self):
        self.df_int % self.df2_int

    def time_frame_float_mod(self):
        self.df % self.df2


class Timeseries(object):

    goal_time = 0.2

    params = [None, 'US/Eastern']
    param_names = ['tz']

    def setup(self, tz):
        N = 10**6
        halfway = (N // 2) - 1
        self.s = Series(date_range('20010101', periods=N, freq='T', tz=tz))
        self.ts = self.s[halfway]

        self.s2 = Series(date_range('20010101', periods=N, freq='s', tz=tz))

    def time_series_timestamp_compare(self, tz):
        self.s <= self.ts

    def time_timestamp_series_compare(self, tz):
        self.ts >= self.s

    def time_timestamp_ops_diff(self, tz):
        self.s2.diff()

    def time_timestamp_ops_diff_with_shift(self, tz):
        self.s - self.s.shift()


class AddOverflowScalar(object):

    goal_time = 0.2

    params = [1, -1, 0]
    param_names = ['scalar']

    def setup(self, scalar):
        N = 10**6
        self.arr = np.arange(N)

    def time_add_overflow_scalar(self, scalar):
        checked_add_with_arr(self.arr, scalar)


class AddOverflowArray(object):

    goal_time = 0.2

    def setup(self):
        N = 10**6
        self.arr = np.arange(N)
        self.arr_rev = np.arange(-N, 0)
        self.arr_mixed = np.array([1, -1]).repeat(N / 2)
        self.arr_nan_1 = np.random.choice([True, False], size=N)
        self.arr_nan_2 = np.random.choice([True, False], size=N)

    def time_add_overflow_arr_rev(self):
        checked_add_with_arr(self.arr, self.arr_rev)

    def time_add_overflow_arr_mask_nan(self):
        checked_add_with_arr(self.arr, self.arr_mixed, arr_mask=self.arr_nan_1)

    def time_add_overflow_b_mask_nan(self):
        checked_add_with_arr(self.arr, self.arr_mixed,
                             b_mask=self.arr_nan_1)

    def time_add_overflow_both_arg_nan(self):
        checked_add_with_arr(self.arr, self.arr_mixed, arr_mask=self.arr_nan_1,
                             b_mask=self.arr_nan_2)
