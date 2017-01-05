from .pandas_vb_common import *
import pandas as pd
import pandas.computation.expressions as expr


class Eval(object):
    goal_time = 0.2

    params = [['numexpr', 'python'], [1, 'all']]
    param_names = ['engine', 'threads']

    def setup(self, engine, threads):
        self.df = DataFrame(np.random.randn(20000, 100))
        self.df2 = DataFrame(np.random.randn(20000, 100))
        self.df3 = DataFrame(np.random.randn(20000, 100))
        self.df4 = DataFrame(np.random.randn(20000, 100))

        if threads == 1:
            expr.set_numexpr_threads(1)

    def time_add(self, engine, threads):
        df, df2, df3, df4 = self.df, self.df2, self.df3, self.df4
        pd.eval('df + df2 + df3 + df4', engine=engine)

    def time_and(self, engine, threads):
        df, df2, df3, df4 = self.df, self.df2, self.df3, self.df4
        pd.eval('(df > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)', engine=engine)

    def time_chained_cmp(self, engine, threads):
        df, df2, df3, df4 = self.df, self.df2, self.df3, self.df4
        pd.eval('df < df2 < df3 < df4', engine=engine)

    def time_mult(self, engine, threads):
        df, df2, df3, df4 = self.df, self.df2, self.df3, self.df4
        pd.eval('df * df2 * df3 * df4', engine=engine)

    def teardown(self, engine, threads):
        expr.set_numexpr_threads()


class Query(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.halfway = ((self.N // 2) - 1)
        self.index = date_range('20010101', periods=self.N, freq='T')
        self.s = Series(self.index)
        self.ts = self.s.iloc[self.halfway]
        self.df = DataFrame({'a': np.random.randn(self.N), }, index=self.index)
        self.df2 = DataFrame({'dates': self.s.values,})

        self.df3 = DataFrame({'a': np.random.randn(self.N),})
        self.min_val = self.df3['a'].min()
        self.max_val = self.df3['a'].max()

    def time_query_datetime_index(self):
        ts = self.ts
        self.df.query('index < @ts')

    def time_query_datetime_series(self):
        ts = self.ts
        self.df2.query('dates < @ts')

    def time_query_with_boolean_selection(self):
        min_val, max_val = self.min_val, self.max_val
        self.df.query('(a >= @min_val) & (a <= @max_val)')
