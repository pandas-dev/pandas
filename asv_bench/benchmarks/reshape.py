import string
from itertools import product

import numpy as np
from pandas import DataFrame, MultiIndex, date_range, melt, wide_to_long
import pandas as pd


class Melt(object):

    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(10000, 3), columns=['A', 'B', 'C'])
        self.df['id1'] = np.random.randint(0, 10, 10000)
        self.df['id2'] = np.random.randint(100, 1000, 10000)

    def time_melt_dataframe(self):
        melt(self.df, id_vars=['id1', 'id2'])


class Pivot(object):

    goal_time = 0.2

    def setup(self):
        N = 10000
        index = date_range('1/1/2000', periods=N, freq='h')
        data = {'value': np.random.randn(N * 50),
                'variable': np.arange(50).repeat(N),
                'date': np.tile(index.values, 50)}
        self.df = DataFrame(data)

    def time_reshape_pivot_time_series(self):
        self.df.pivot('date', 'variable', 'value')


class SimpleReshape(object):

    goal_time = 0.2

    def setup(self):
        arrays = [np.arange(100).repeat(100),
                  np.roll(np.tile(np.arange(100), 100), 25)]
        index = MultiIndex.from_arrays(arrays)
        self.df = DataFrame(np.random.randn(10000, 4), index=index)
        self.udf = self.df.unstack(1)

    def time_stack(self):
        self.udf.stack()

    def time_unstack(self):
        self.df.unstack(1)


class Unstack(object):

    goal_time = 0.2

    def setup(self):
        m = 100
        n = 1000

        levels = np.arange(m)
        index = MultiIndex.from_product([levels] * 2)
        columns = np.arange(n)
        values = np.arange(m * m * n).reshape(m * m, n)
        self.df = DataFrame(values, index, columns)
        self.df2 = self.df.iloc[:-1]

    def time_full_product(self):
        self.df.unstack()

    def time_without_last_row(self):
        self.df2.unstack()


class SparseIndex(object):

    goal_time = 0.2

    def setup(self):
        NUM_ROWS = 1000
        self.df = DataFrame({'A': np.random.randint(50, size=NUM_ROWS),
                             'B': np.random.randint(50, size=NUM_ROWS),
                             'C': np.random.randint(-10, 10, size=NUM_ROWS),
                             'D': np.random.randint(-10, 10, size=NUM_ROWS),
                             'E': np.random.randint(10, size=NUM_ROWS),
                             'F': np.random.randn(NUM_ROWS)})
        self.df = self.df.set_index(['A', 'B', 'C', 'D', 'E'])

    def time_unstack(self):
        self.df.unstack()


class WideToLong(object):

    goal_time = 0.2

    def setup(self):
        nyrs = 20
        nidvars = 20
        N = 5000
        self.letters = list('ABCD')
        yrvars = [l + str(num)
                  for l, num in product(self.letters, range(1, nyrs + 1))]
        columns = [str(i) for i in range(nidvars)] + yrvars
        self.df = DataFrame(np.random.randn(N, nidvars + len(yrvars)),
                            columns=columns)
        self.df['id'] = self.df.index

    def time_wide_to_long_big(self):
        wide_to_long(self.df, self.letters, i='id', j='year')


class PivotTable(object):

    goal_time = 0.2

    def setup(self):
        N = 100000
        fac1 = np.array(['A', 'B', 'C'], dtype='O')
        fac2 = np.array(['one', 'two'], dtype='O')
        ind1 = np.random.randint(0, 3, size=N)
        ind2 = np.random.randint(0, 2, size=N)
        self.df = DataFrame({'key1': fac1.take(ind1),
                             'key2': fac2.take(ind2),
                             'key3': fac2.take(ind2),
                             'value1': np.random.randn(N),
                             'value2': np.random.randn(N),
                             'value3': np.random.randn(N)})

    def time_pivot_table(self):
        self.df.pivot_table(index='key1', columns=['key2', 'key3'])


class GetDummies(object):
    goal_time = 0.2

    def setup(self):
        categories = list(string.ascii_letters[:12])
        s = pd.Series(np.random.choice(categories, size=1000000),
                      dtype=pd.api.types.CategoricalDtype(categories))
        self.s = s

    def time_get_dummies_1d(self):
        pd.get_dummies(self.s, sparse=False)

    def time_get_dummies_1d_sparse(self):
        pd.get_dummies(self.s, sparse=True)


from .pandas_vb_common import setup  # noqa: F401
