from .pandas_vb_common import *


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

    def time_getitem_lists(self):
        self.s[np.arange(10000).tolist()]

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

        self.mi_large = MultiIndex.from_product(
            [np.arange(1000), np.arange(20), list(string.ascii_letters)],
            names=['one', 'two', 'three'])
        self.mi_med = MultiIndex.from_product(
            [np.arange(1000), np.arange(10), list('A')],
            names=['one', 'two', 'three'])
        self.mi_small = MultiIndex.from_product(
            [np.arange(100), list('A'), list('A')],
            names=['one', 'two', 'three'])

        rng = np.random.RandomState(4)
        size = 1 << 16
        self.mi_unused_levels = pd.MultiIndex.from_arrays([
            rng.randint(0, 1 << 13, size),
            rng.randint(0, 1 << 10, size)])[rng.rand(size) < 0.1]

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

    def time_multiindex_large_get_loc(self):
        self.mi_large.get_loc((999, 19, 'Z'))

    def time_multiindex_large_get_loc_warm(self):
        for _ in range(1000):
            self.mi_large.get_loc((999, 19, 'Z'))

    def time_multiindex_med_get_loc(self):
        self.mi_med.get_loc((999, 9, 'A'))

    def time_multiindex_med_get_loc_warm(self):
        for _ in range(1000):
            self.mi_med.get_loc((999, 9, 'A'))

    def time_multiindex_string_get_loc(self):
        self.mi_small.get_loc((99, 'A', 'A'))

    def time_multiindex_small_get_loc_warm(self):
        for _ in range(1000):
            self.mi_small.get_loc((99, 'A', 'A'))

    def time_is_monotonic(self):
        self.miint.is_monotonic

    def time_remove_unused_levels(self):
        self.mi_unused_levels.remove_unused_levels()


class IntervalIndexing(object):
    goal_time = 0.2

    def setup(self):
        self.monotonic = Series(np.arange(1000000),
                                index=IntervalIndex.from_breaks(np.arange(1000001)))

    def time_getitem_scalar(self):
        self.monotonic[80000]

    def time_loc_scalar(self):
        self.monotonic.loc[80000]

    def time_getitem_list(self):
        self.monotonic[80000:]

    def time_loc_list(self):
        self.monotonic.loc[80000:]


class PanelIndexing(object):
    goal_time = 0.2

    def setup(self):
        self.p = Panel(np.random.randn(100, 100, 100))
        self.inds = range(0, 100, 10)

    def time_subset(self):
        self.p.ix[(self.inds, self.inds, self.inds)]


class IndexerLookup(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(range(10))

    def time_lookup_iloc(self):
        self.s.iloc

    def time_lookup_ix(self):
        self.s.ix

    def time_lookup_loc(self):
        self.s.loc


class BooleanRowSelect(object):

    goal_time = 0.2

    def setup(self):
        N = 10000
        np.random.seed(1234)
        self.df = DataFrame(np.random.randn(N, 100))
        self.bool_arr = np.zeros(N, dtype=bool)
        self.bool_arr[:1000] = True

    def time_frame_boolean_row_select(self):
        self.df[self.bool_arr]


class GetItemSingleColumn(object):

    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.df2 = DataFrame(np.random.randn(3000, 1), columns=['A'])
        self.df3 = DataFrame(np.random.randn(3000, 1))

    def time_frame_getitem_single_column_label(self):
        self.df2['A']

    def time_frame_getitem_single_column_int(self):
        self.df3[0]


class AssignTimeseriesIndex(object):

    goal_time = 0.2

    def setup(self):
        N = 100000
        np.random.seed(1234)
        dx = date_range('1/1/2000', periods=N, freq='H')
        self.df = DataFrame(np.random.randn(N, 1), columns=['A'], index=idx)

    def time_frame_assign_timeseries_index(self):
        self.df['date'] = self.df.index


class InsertColumns(object):

    goal_time = 0.2

    def setup(self):
        self.N = 10**3
        self.df = DataFrame(index=range(N))

    def time_insert(self):
        np.random.seed(1234)
        for i in range(100):
            self.df.insert(0, i, np.random.randn(self.N))

    def time_assign_with_setitem(self):
        np.random.seed(1234)
        for i in range(100):
            self.df[i] = np.random.randn(self.N)


