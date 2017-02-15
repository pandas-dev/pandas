from .pandas_vb_common import *
try:
    import pandas.computation.expressions as expr
except:
    expr = None


class Int64Indexing(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_getitem_scalar(self):
        self.s[800000]

    def time_getitem_slice(self):
        self.s[:800000]

    def time_getitem_list_like(self):
        self.s[[800000]]

    def time_getitem_array(self):
        self.s[np.arange(10000)]

    def time_iloc_array(self):
        self.s.iloc[np.arange(10000)]

    def time_iloc_list_like(self):
        self.s.iloc[[800000]]

    def time_iloc_scalar(self):
        self.s.iloc[800000]

    def time_iloc_slice(self):
        self.s.iloc[:800000]

    def time_ix_array(self):
        self.s.ix[np.arange(10000)]

    def time_ix_list_like(self):
        self.s.ix[[800000]]

    def time_ix_scalar(self):
        self.s.ix[800000]

    def time_ix_slice(self):
        self.s.ix[:800000]

    def time_loc_array(self):
        self.s.loc[np.arange(10000)]

    def time_loc_list_like(self):
        self.s.loc[[800000]]

    def time_loc_scalar(self):
        self.s.loc[800000]

    def time_loc_slice(self):
        self.s.loc[:800000]


class StringIndexing(object):
    goal_time = 0.2

    def setup(self):
        self.index = tm.makeStringIndex(1000000)
        self.s = Series(np.random.rand(1000000), index=self.index)
        self.lbl = self.s.index[800000]

    def time_getitem_label_slice(self):
        self.s[:self.lbl]

    def time_getitem_pos_slice(self):
        self.s[:800000]

    def time_get_value(self):
        self.s.get_value(self.lbl)


class DatetimeIndexing(object):
    goal_time = 0.2

    def setup(self):
        tm.N = 1000
        self.ts = tm.makeTimeSeries()
        self.dt = self.ts.index[500]

    def time_getitem_scalar(self):
        self.ts[self.dt]


class DataFrameIndexing(object):
    goal_time = 0.2

    def setup(self):
        self.index = tm.makeStringIndex(1000)
        self.columns = tm.makeStringIndex(30)
        self.df = DataFrame(np.random.randn(1000, 30), index=self.index,
                            columns=self.columns)
        self.idx = self.index[100]
        self.col = self.columns[10]

        self.df2 = DataFrame(np.random.randn(10000, 4),
                             columns=['A', 'B', 'C', 'D'])
        self.indexer = (self.df2['B'] > 0)
        self.obj_indexer = self.indexer.astype('O')

        # duptes
        self.idx_dupe = (np.array(range(30)) * 99)
        self.df3 = DataFrame({'A': ([0.1] * 1000), 'B': ([1] * 1000),})
        self.df3 = concat([self.df3, (2 * self.df3), (3 * self.df3)])

        self.df_big = DataFrame(dict(A=(['foo'] * 1000000)))

    def time_get_value(self):
        self.df.get_value(self.idx, self.col)

    def time_get_value_ix(self):
        self.df.ix[(self.idx, self.col)]

    def time_getitem_scalar(self):
        self.df[self.col][self.idx]

    def time_boolean_rows(self):
        self.df2[self.indexer]

    def time_boolean_rows_object(self):
        self.df2[self.obj_indexer]

    def time_iloc_dups(self):
        self.df3.iloc[self.idx_dupe]

    def time_loc_dups(self):
        self.df3.loc[self.idx_dupe]

    def time_iloc_big(self):
        self.df_big.iloc[:100, 0]


class IndexingMethods(object):
    # GH 13166
    goal_time = 0.2

    def setup(self):
        a = np.arange(100000)
        self.ind = pd.Float64Index(a * 4.8000000418824129e-08)

        self.s = Series(np.random.rand(100000))
        self.ts = Series(np.random.rand(100000),
                         index=date_range('2011-01-01', freq='S', periods=100000))
        self.indexer = ([True, False, True, True, False] * 20000)

    def time_get_loc_float(self):
        self.ind.get_loc(0)

    def time_take_dtindex(self):
        self.ts.take(self.indexer)

    def time_take_intindex(self):
        self.s.take(self.indexer)


class MultiIndexing(object):
    goal_time = 0.2

    def setup(self):
        self.mi = MultiIndex.from_tuples([(x, y) for x in range(1000) for y in range(1000)])
        self.s = Series(np.random.randn(1000000), index=self.mi)
        self.df = DataFrame(self.s)

        # slicers
        np.random.seed(1234)
        self.idx = pd.IndexSlice
        self.n = 100000
        self.mdt = pandas.DataFrame()
        self.mdt['A'] = np.random.choice(range(10000, 45000, 1000), self.n)
        self.mdt['B'] = np.random.choice(range(10, 400), self.n)
        self.mdt['C'] = np.random.choice(range(1, 150), self.n)
        self.mdt['D'] = np.random.choice(range(10000, 45000), self.n)
        self.mdt['x'] = np.random.choice(range(400), self.n)
        self.mdt['y'] = np.random.choice(range(25), self.n)
        self.test_A = 25000
        self.test_B = 25
        self.test_C = 40
        self.test_D = 35000
        self.eps_A = 5000
        self.eps_B = 5
        self.eps_C = 5
        self.eps_D = 5000
        self.mdt2 = self.mdt.set_index(['A', 'B', 'C', 'D']).sortlevel()
        self.miint = MultiIndex.from_product(
            [np.arange(1000),
             np.arange(1000)], names=['one', 'two'])

        import string
        self.mistring = MultiIndex.from_product(
            [np.arange(1000),
             np.arange(20), list(string.ascii_letters)],
            names=['one', 'two', 'three'])

    def time_series_xs_mi_ix(self):
        self.s.ix[999]

    def time_frame_xs_mi_ix(self):
        self.df.ix[999]

    def time_multiindex_slicers(self):
        self.mdt2.loc[self.idx[
            (self.test_A - self.eps_A):(self.test_A + self.eps_A),
            (self.test_B - self.eps_B):(self.test_B + self.eps_B),
            (self.test_C - self.eps_C):(self.test_C + self.eps_C),
            (self.test_D - self.eps_D):(self.test_D + self.eps_D)], :]

    def time_multiindex_get_indexer(self):
        self.miint.get_indexer(
            np.array([(0, 10), (0, 11), (0, 12),
                      (0, 13), (0, 14), (0, 15),
                      (0, 16), (0, 17), (0, 18),
                      (0, 19)], dtype=object))

    def time_multiindex_string_get_loc(self):
        self.mistring.get_loc((999, 19, 'Z'))

    def time_is_monotonic(self):
        self.miint.is_monotonic


class PanelIndexing(object):
    goal_time = 0.2

    def setup(self):
        self.p = Panel(np.random.randn(100, 100, 100))
        self.inds = range(0, 100, 10)

    def time_subset(self):
        self.p.ix[(self.inds, self.inds, self.inds)]
