from pandas_vb_common import *
import pandas.computation.expressions as expr
import pandas as pd


class eval_frame_add_all_threads(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))

    def time_eval_frame_add_all_threads(self):
        pd.eval('df + df2 + df3 + df4')


class eval_frame_add_one_thread(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))
        expr.set_numexpr_threads(1)

    def time_eval_frame_add_one_thread(self):
        pd.eval('df + df2 + df3 + df4')


class eval_frame_add_python(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))

    def time_eval_frame_add_python(self):
        pd.eval('df + df2 + df3 + df4', engine='python')


class eval_frame_add_python_one_thread(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))
        expr.set_numexpr_threads(1)

    def time_eval_frame_add_python_one_thread(self):
        pd.eval('df + df2 + df3 + df4', engine='python')


class eval_frame_and_all_threads(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))

    def time_eval_frame_and_all_threads(self):
        pd.eval('(df > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)')


class eval_frame_and_python_one_thread(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))
        expr.set_numexpr_threads(1)

    def time_eval_frame_and_python_one_thread(self):
        pd.eval('(df > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)', engine='python')


class eval_frame_and_python(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))

    def time_eval_frame_and_python(self):
        pd.eval('(df > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)', engine='python')


class eval_frame_chained_cmp_all_threads(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))

    def time_eval_frame_chained_cmp_all_threads(self):
        pd.eval('df < df2 < df3 < df4')


class eval_frame_chained_cmp_python_one_thread(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))
        expr.set_numexpr_threads(1)

    def time_eval_frame_chained_cmp_python_one_thread(self):
        pd.eval('df < df2 < df3 < df4', engine='python')


class eval_frame_chained_cmp_python(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))

    def time_eval_frame_chained_cmp_python(self):
        pd.eval('df < df2 < df3 < df4', engine='python')


class eval_frame_mult_all_threads(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))

    def time_eval_frame_mult_all_threads(self):
        pd.eval('df * df2 * df3 * df4')


class eval_frame_mult_one_thread(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))
        expr.set_numexpr_threads(1)

    def time_eval_frame_mult_one_thread(self):
        pd.eval('df * df2 * df3 * df4')


class eval_frame_mult_python(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))

    def time_eval_frame_mult_python(self):
        pd.eval('df * df2 * df3 * df4', engine='python')


class eval_frame_mult_python_one_thread(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))
        expr.set_numexpr_threads(1)

    def time_eval_frame_mult_python_one_thread(self):
        pd.eval('df * df2 * df3 * df4', engine='python')


class query_datetime_index(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.halfway = ((self.N // 2) - 1)
        self.index = date_range('20010101', periods=self.N, freq='T')
        self.s = Series(self.index)
        self.ts = self.s.iloc[self.halfway]
        self.df = DataFrame({'a': np.random.randn(self.N), }, index=self.index)

    def time_query_datetime_index(self):
        self.df.query('index < @ts')


class query_datetime_series(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.halfway = ((self.N // 2) - 1)
        self.index = date_range('20010101', periods=self.N, freq='T')
        self.s = Series(self.index)
        self.ts = self.s.iloc[self.halfway]
        self.df = DataFrame({'dates': self.s.values, })

    def time_query_datetime_series(self):
        self.df.query('dates < @ts')


class query_with_boolean_selection(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.halfway = ((self.N // 2) - 1)
        self.index = date_range('20010101', periods=self.N, freq='T')
        self.s = Series(self.index)
        self.ts = self.s.iloc[self.halfway]
        self.N = 1000000
        self.df = DataFrame({'a': np.random.randn(self.N), })
        self.min_val = self.df['a'].min()
        self.max_val = self.df['a'].max()

    def time_query_with_boolean_selection(self):
        self.df.query('(a >= @min_val) & (a <= @max_val)')