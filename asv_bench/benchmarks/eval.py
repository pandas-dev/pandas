import numpy as np
import pandas as pd
try:
    import pandas.core.computation.expressions as expr
except ImportError:
    import pandas.computation.expressions as expr


class Eval(object):

    goal_time = 0.2

    params = [['numexpr', 'python'], [1, 'all']]
    param_names = ['engine', 'threads']

    def setup(self, engine, threads):
        self.df = pd.DataFrame(np.random.randn(20000, 100))
        self.df2 = pd.DataFrame(np.random.randn(20000, 100))
        self.df3 = pd.DataFrame(np.random.randn(20000, 100))
        self.df4 = pd.DataFrame(np.random.randn(20000, 100))

        if threads == 1:
            expr.set_numexpr_threads(1)

    def time_add(self, engine, threads):
        pd.eval('self.df + self.df2 + self.df3 + self.df4', engine=engine)

    def time_and(self, engine, threads):
        pd.eval('(self.df > 0) & (self.df2 > 0) & '
                '(self.df3 > 0) & (self.df4 > 0)', engine=engine)

    def time_chained_cmp(self, engine, threads):
        pd.eval('self.df < self.df2 < self.df3 < self.df4', engine=engine)

    def time_mult(self, engine, threads):
        pd.eval('self.df * self.df2 * self.df3 * self.df4', engine=engine)

    def teardown(self, engine, threads):
        expr.set_numexpr_threads()


class Query(object):

    goal_time = 0.2

    def setup(self):
        N = 10**6
        halfway = (N // 2) - 1
        index = pd.date_range('20010101', periods=N, freq='T')
        s = pd.Series(index)
        self.ts = s.iloc[halfway]
        self.df = pd.DataFrame({'a': np.random.randn(N), 'dates': s},
                               index=index)
        data = np.random.randn(N)
        self.min_val = data.min()
        self.max_val = data.max()

    def time_query_datetime_index(self):
        self.df.query('index < @self.ts')

    def time_query_datetime_column(self):
        self.df.query('dates < @self.ts')

    def time_query_with_boolean_selection(self):
        self.df.query('(a >= @self.min_val) & (a <= @self.max_val)')


from .pandas_vb_common import setup  # noqa: F401
