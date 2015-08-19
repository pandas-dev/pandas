from pandas_vb_common import *
import pandas.computation.expressions as expr


class dataframe_getitem_scalar(object):
    goal_time = 0.2

    def setup(self):
        self.index = tm.makeStringIndex(1000)
        self.columns = tm.makeStringIndex(30)
        self.df = DataFrame(np.random.rand(1000, 30), index=self.index, columns=self.columns)
        self.idx = self.index[100]
        self.col = self.columns[10]

    def time_dataframe_getitem_scalar(self):
        self.df[self.col][self.idx]


class datamatrix_getitem_scalar(object):
    goal_time = 0.2

    def setup(self):
        try:
            self.klass = DataMatrix
        except:
            self.klass = DataFrame
        self.index = tm.makeStringIndex(1000)
        self.columns = tm.makeStringIndex(30)
        self.df = self.klass(np.random.rand(1000, 30), index=self.index, columns=self.columns)
        self.idx = self.index[100]
        self.col = self.columns[10]

    def time_datamatrix_getitem_scalar(self):
        self.df[self.col][self.idx]


class series_get_value(object):
    goal_time = 0.2

    def setup(self):
        self.index = tm.makeStringIndex(1000)
        self.s = Series(np.random.rand(1000), index=self.index)
        self.idx = self.index[100]

    def time_series_get_value(self):
        self.s.get_value(self.idx)


class time_series_getitem_scalar(object):
    goal_time = 0.2

    def setup(self):
        tm.N = 1000
        self.ts = tm.makeTimeSeries()
        self.dt = self.ts.index[500]

    def time_time_series_getitem_scalar(self):
        self.ts[self.dt]


class frame_iloc_big(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(dict(A=(['foo'] * 1000000)))

    def time_frame_iloc_big(self):
        self.df.iloc[:100, 0]


class frame_iloc_dups(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'A': ([0.1] * 3000), 'B': ([1] * 3000), })
        self.idx = (np.array(range(30)) * 99)
        self.df2 = DataFrame({'A': ([0.1] * 1000), 'B': ([1] * 1000), })
        self.df2 = concat([self.df2, (2 * self.df2), (3 * self.df2)])

    def time_frame_iloc_dups(self):
        self.df2.iloc[self.idx]


class frame_loc_dups(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'A': ([0.1] * 3000), 'B': ([1] * 3000), })
        self.idx = (np.array(range(30)) * 99)
        self.df2 = DataFrame({'A': ([0.1] * 1000), 'B': ([1] * 1000), })
        self.df2 = concat([self.df2, (2 * self.df2), (3 * self.df2)])

    def time_frame_loc_dups(self):
        self.df2.loc[self.idx]


class frame_xs_mi_ix(object):
    goal_time = 0.2

    def setup(self):
        self.mi = MultiIndex.from_tuples([(x, y) for x in range(1000) for y in range(1000)])
        self.s = Series(np.random.randn(1000000), index=self.mi)
        self.df = DataFrame(self.s)

    def time_frame_xs_mi_ix(self):
        self.df.ix[999]


class indexing_dataframe_boolean(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(50000, 100))
        self.df2 = DataFrame(np.random.randn(50000, 100))

    def time_indexing_dataframe_boolean(self):
        (self.df > self.df2)


class indexing_dataframe_boolean_no_ne(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(50000, 100))
        self.df2 = DataFrame(np.random.randn(50000, 100))
        expr.set_use_numexpr(False)

    def time_indexing_dataframe_boolean_no_ne(self):
        (self.df > self.df2)

    def teardown(self):
        expr.set_use_numexpr(True)


class indexing_dataframe_boolean_rows(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(10000, 4), columns=['A', 'B', 'C', 'D'])
        self.indexer = (self.df['B'] > 0)
        self.obj_indexer = self.indexer.astype('O')

    def time_indexing_dataframe_boolean_rows(self):
        self.df[self.indexer]


class indexing_dataframe_boolean_rows_object(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(10000, 4), columns=['A', 'B', 'C', 'D'])
        self.indexer = (self.df['B'] > 0)
        self.obj_indexer = self.indexer.astype('O')

    def time_indexing_dataframe_boolean_rows_object(self):
        self.df[self.obj_indexer]


class indexing_dataframe_boolean_st(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(50000, 100))
        self.df2 = DataFrame(np.random.randn(50000, 100))
        expr.set_numexpr_threads(1)

    def time_indexing_dataframe_boolean_st(self):
        (self.df > self.df2)

    def teardown(self):
        expr.set_numexpr_threads()


class indexing_frame_get_value(object):
    goal_time = 0.2

    def setup(self):
        self.index = tm.makeStringIndex(1000)
        self.columns = tm.makeStringIndex(30)
        self.df = DataFrame(np.random.randn(1000, 30), index=self.index, columns=self.columns)
        self.idx = self.index[100]
        self.col = self.columns[10]

    def time_indexing_frame_get_value(self):
        self.df.get_value(self.idx, self.col)


class indexing_frame_get_value_ix(object):
    goal_time = 0.2

    def setup(self):
        self.index = tm.makeStringIndex(1000)
        self.columns = tm.makeStringIndex(30)
        self.df = DataFrame(np.random.randn(1000, 30), index=self.index, columns=self.columns)
        self.idx = self.index[100]
        self.col = self.columns[10]

    def time_indexing_frame_get_value_ix(self):
        self.df.ix[(self.idx, self.col)]


class indexing_panel_subset(object):
    goal_time = 0.2

    def setup(self):
        self.p = Panel(np.random.randn(100, 100, 100))
        self.inds = range(0, 100, 10)

    def time_indexing_panel_subset(self):
        self.p.ix[(self.inds, self.inds, self.inds)]


class multiindex_slicers(object):
    goal_time = 0.2

    def setup(self):
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

    def time_multiindex_slicers(self):
        self.mdt2.loc[self.idx[(self.test_A - self.eps_A):(self.test_A + self.eps_A), (self.test_B - self.eps_B):(self.test_B + self.eps_B), (self.test_C - self.eps_C):(self.test_C + self.eps_C), (self.test_D - self.eps_D):(self.test_D + self.eps_D)], :]


class series_getitem_array(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_getitem_array(self):
        self.s[np.arange(10000)]


class series_getitem_label_slice(object):
    goal_time = 0.2

    def setup(self):
        self.index = tm.makeStringIndex(1000000)
        self.s = Series(np.random.rand(1000000), index=self.index)
        self.lbl = self.s.index[800000]

    def time_series_getitem_label_slice(self):
        self.s[:self.lbl]


class series_getitem_list_like(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_getitem_list_like(self):
        self.s[[800000]]


class series_getitem_pos_slice(object):
    goal_time = 0.2

    def setup(self):
        self.index = tm.makeStringIndex(1000000)
        self.s = Series(np.random.rand(1000000), index=self.index)

    def time_series_getitem_pos_slice(self):
        self.s[:800000]


class series_getitem_scalar(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_getitem_scalar(self):
        self.s[800000]


class series_getitem_slice(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_getitem_slice(self):
        self.s[:800000]


class series_iloc_array(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_iloc_array(self):
        self.s.iloc[np.arange(10000)]


class series_iloc_list_like(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_iloc_list_like(self):
        self.s.iloc[[800000]]


class series_iloc_scalar(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_iloc_scalar(self):
        self.s.iloc[800000]


class series_iloc_slice(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_iloc_slice(self):
        self.s.iloc[:800000]


class series_ix_array(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_ix_array(self):
        self.s.ix[np.arange(10000)]


class series_ix_list_like(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_ix_list_like(self):
        self.s.ix[[800000]]


class series_ix_scalar(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_ix_scalar(self):
        self.s.ix[800000]


class series_ix_slice(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_ix_slice(self):
        self.s.ix[:800000]


class series_loc_array(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_loc_array(self):
        self.s.loc[np.arange(10000)]


class series_loc_list_like(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_loc_list_like(self):
        self.s.loc[[800000]]


class series_loc_scalar(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_loc_scalar(self):
        self.s.loc[800000]


class series_loc_slice(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.rand(1000000))

    def time_series_loc_slice(self):
        self.s.loc[:800000]


class series_xs_mi_ix(object):
    goal_time = 0.2

    def setup(self):
        self.mi = MultiIndex.from_tuples([(x, y) for x in range(1000) for y in range(1000)])
        self.s = Series(np.random.randn(1000000), index=self.mi)

    def time_series_xs_mi_ix(self):
        self.s.ix[999]


class sort_level_one(object):
    goal_time = 0.2

    def setup(self):
        self.a = np.repeat(np.arange(100), 1000)
        self.b = np.tile(np.arange(1000), 100)
        self.midx = MultiIndex.from_arrays([self.a, self.b])
        self.midx = self.midx.take(np.random.permutation(np.arange(100000)))

    def time_sort_level_one(self):
        self.midx.sortlevel(1)


class sort_level_zero(object):
    goal_time = 0.2

    def setup(self):
        self.a = np.repeat(np.arange(100), 1000)
        self.b = np.tile(np.arange(1000), 100)
        self.midx = MultiIndex.from_arrays([self.a, self.b])
        self.midx = self.midx.take(np.random.permutation(np.arange(100000)))

    def time_sort_level_zero(self):
        self.midx.sortlevel(0)