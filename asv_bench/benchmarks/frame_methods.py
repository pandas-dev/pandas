from pandas_vb_common import *


class frame_apply_axis_1(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 100))

    def time_frame_apply_axis_1(self):
        self.df.apply((lambda x: (x + 1)), axis=1)


class frame_apply_lambda_mean(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 100))

    def time_frame_apply_lambda_mean(self):
        self.df.apply((lambda x: x.sum()))


class frame_apply_np_mean(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 100))

    def time_frame_apply_np_mean(self):
        self.df.apply(np.mean)


class frame_apply_pass_thru(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 100))

    def time_frame_apply_pass_thru(self):
        self.df.apply((lambda x: x))


class frame_apply_ref_by_name(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 3), columns=list('ABC'))

    def time_frame_apply_ref_by_name(self):
        self.df.apply((lambda x: (x['A'] + x['B'])), axis=1)


class frame_apply_user_func(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.arange(1028.0))
        self.df = DataFrame({i: self.s for i in range(1028)})

    def time_frame_apply_user_func(self):
        self.df.apply((lambda x: np.corrcoef(x, self.s)[(0, 1)]))


class frame_assign_timeseries_index(object):
    goal_time = 0.2

    def setup(self):
        self.idx = date_range('1/1/2000', periods=100000, freq='D')
        self.df = DataFrame(randn(100000, 1), columns=['A'], index=self.idx)

        def f(x):
            self.x = self.x.copy()
            self.x['date'] = self.x.index

    def time_frame_assign_timeseries_index(self):
        f(self.df)


class frame_boolean_row_select(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 100))
        self.bool_arr = np.zeros(10000, dtype=bool)
        self.bool_arr[:1000] = True

    def time_frame_boolean_row_select(self):
        self.df[self.bool_arr]


class frame_count_level_axis0_mixed_dtypes_multi(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(10000, 1000)
        self.df = DataFrame(self.data)
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan
        self.df['foo'] = 'bar'
        self.df.index = MultiIndex.from_tuples(self.df.index.map((lambda x: (x, x))))
        self.df.columns = MultiIndex.from_tuples(self.df.columns.map((lambda x: (x, x))))

    def time_frame_count_level_axis0_mixed_dtypes_multi(self):
        self.df.count(axis=0, level=1)


class frame_count_level_axis0_multi(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(10000, 1000)
        self.df = DataFrame(self.data)
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan
        self.df.index = MultiIndex.from_tuples(self.df.index.map((lambda x: (x, x))))
        self.df.columns = MultiIndex.from_tuples(self.df.columns.map((lambda x: (x, x))))

    def time_frame_count_level_axis0_multi(self):
        self.df.count(axis=0, level=1)


class frame_count_level_axis1_mixed_dtypes_multi(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(10000, 1000)
        self.df = DataFrame(self.data)
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan
        self.df['foo'] = 'bar'
        self.df.index = MultiIndex.from_tuples(self.df.index.map((lambda x: (x, x))))
        self.df.columns = MultiIndex.from_tuples(self.df.columns.map((lambda x: (x, x))))

    def time_frame_count_level_axis1_mixed_dtypes_multi(self):
        self.df.count(axis=1, level=1)


class frame_count_level_axis1_multi(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(10000, 1000)
        self.df = DataFrame(self.data)
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan
        self.df.index = MultiIndex.from_tuples(self.df.index.map((lambda x: (x, x))))
        self.df.columns = MultiIndex.from_tuples(self.df.columns.map((lambda x: (x, x))))

    def time_frame_count_level_axis1_multi(self):
        self.df.count(axis=1, level=1)


class frame_dropna_axis0_all(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(10000, 1000)
        self.df = DataFrame(self.data)
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan

    def time_frame_dropna_axis0_all(self):
        self.df.dropna(how='all', axis=0)


class frame_dropna_axis0_all_mixed_dtypes(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(10000, 1000)
        self.df = DataFrame(self.data)
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan
        self.df['foo'] = 'bar'

    def time_frame_dropna_axis0_all_mixed_dtypes(self):
        self.df.dropna(how='all', axis=0)


class frame_dropna_axis0_any(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(10000, 1000)
        self.df = DataFrame(self.data)
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan

    def time_frame_dropna_axis0_any(self):
        self.df.dropna(how='any', axis=0)


class frame_dropna_axis0_any_mixed_dtypes(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(10000, 1000)
        self.df = DataFrame(self.data)
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan
        self.df['foo'] = 'bar'

    def time_frame_dropna_axis0_any_mixed_dtypes(self):
        self.df.dropna(how='any', axis=0)


class frame_dropna_axis1_all(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(10000, 1000)
        self.df = DataFrame(self.data)
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan

    def time_frame_dropna_axis1_all(self):
        self.df.dropna(how='all', axis=1)


class frame_dropna_axis1_all_mixed_dtypes(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(10000, 1000)
        self.df = DataFrame(self.data)
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan
        self.df['foo'] = 'bar'

    def time_frame_dropna_axis1_all_mixed_dtypes(self):
        self.df.dropna(how='all', axis=1)


class frame_dropna_axis1_any(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(10000, 1000)
        self.df = DataFrame(self.data)
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan

    def time_frame_dropna_axis1_any(self):
        self.df.dropna(how='any', axis=1)


class frame_dropna_axis1_any_mixed_dtypes(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(10000, 1000)
        self.df = DataFrame(self.data)
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan
        self.df['foo'] = 'bar'

    def time_frame_dropna_axis1_any_mixed_dtypes(self):
        self.df.dropna(how='any', axis=1)


class frame_dtypes(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 1000))

    def time_frame_dtypes(self):
        self.df.dtypes


class frame_duplicated(object):
    goal_time = 0.2

    def setup(self):
        self.n = (1 << 20)
        self.t = date_range('2015-01-01', freq='S', periods=(self.n // 64))
        self.xs = np.random.randn((self.n // 64)).round(2)
        self.df = DataFrame({'a': np.random.randint(((-1) << 8), (1 << 8), self.n), 'b': np.random.choice(self.t, self.n), 'c': np.random.choice(self.xs, self.n), })

    def time_frame_duplicated(self):
        self.df.duplicated()


class frame_fancy_lookup(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(10000, 8), columns=list('abcdefgh'))
        self.df['foo'] = 'bar'
        self.row_labels = list(self.df.index[::10])[:900]
        self.col_labels = (list(self.df.columns) * 100)
        self.row_labels_all = np.array((list(self.df.index) * len(self.df.columns)), dtype='object')
        self.col_labels_all = np.array((list(self.df.columns) * len(self.df.index)), dtype='object')

    def time_frame_fancy_lookup(self):
        self.df.lookup(self.row_labels, self.col_labels)


class frame_fancy_lookup_all(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(10000, 8), columns=list('abcdefgh'))
        self.df['foo'] = 'bar'
        self.row_labels = list(self.df.index[::10])[:900]
        self.col_labels = (list(self.df.columns) * 100)
        self.row_labels_all = np.array((list(self.df.index) * len(self.df.columns)), dtype='object')
        self.col_labels_all = np.array((list(self.df.columns) * len(self.df.index)), dtype='object')

    def time_frame_fancy_lookup_all(self):
        self.df.lookup(self.row_labels_all, self.col_labels_all)


class frame_fillna_inplace(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 100))
        self.df.values[::2] = np.nan

    def time_frame_fillna_inplace(self):
        self.df.fillna(0, inplace=True)


class frame_float_equal(object):
    goal_time = 0.2

    def setup(self):

        def make_pair(frame):
            self.df = frame
            self.df2 = self.df.copy()
            self.df2.ix[((-1), (-1))] = np.nan
            return (self.df, self.df2)

        def test_equal(name):
            (self.df, self.df2) = pairs[name]
            return self.df.equals(self.df)

        def test_unequal(name):
            (self.df, self.df2) = pairs[name]
            return self.df.equals(self.df2)
        self.float_df = DataFrame(np.random.randn(1000, 1000))
        self.object_df = DataFrame(([(['foo'] * 1000)] * 1000))
        self.nonunique_cols = self.object_df.copy()
        self.nonunique_cols.columns = (['A'] * len(self.nonunique_cols.columns))
        self.pairs = dict([(name, make_pair(frame)) for (name, frame) in (('float_df', self.float_df), ('object_df', self.object_df), ('nonunique_cols', self.nonunique_cols))])

    def time_frame_float_equal(self):
        test_equal('float_df')


class frame_float_unequal(object):
    goal_time = 0.2

    def setup(self):

        def make_pair(frame):
            self.df = frame
            self.df2 = self.df.copy()
            self.df2.ix[((-1), (-1))] = np.nan
            return (self.df, self.df2)

        def test_equal(name):
            (self.df, self.df2) = pairs[name]
            return self.df.equals(self.df)

        def test_unequal(name):
            (self.df, self.df2) = pairs[name]
            return self.df.equals(self.df2)
        self.float_df = DataFrame(np.random.randn(1000, 1000))
        self.object_df = DataFrame(([(['foo'] * 1000)] * 1000))
        self.nonunique_cols = self.object_df.copy()
        self.nonunique_cols.columns = (['A'] * len(self.nonunique_cols.columns))
        self.pairs = dict([(name, make_pair(frame)) for (name, frame) in (('float_df', self.float_df), ('object_df', self.object_df), ('nonunique_cols', self.nonunique_cols))])

    def time_frame_float_unequal(self):
        test_unequal('float_df')


class frame_from_records_generator(object):
    goal_time = 0.2

    def setup(self):

        def get_data(n=100000):
            return ((x, (x * 20), (x * 100)) for x in xrange(n))

    def time_frame_from_records_generator(self):
        self.df = DataFrame.from_records(get_data())


class frame_from_records_generator_nrows(object):
    goal_time = 0.2

    def setup(self):

        def get_data(n=100000):
            return ((x, (x * 20), (x * 100)) for x in xrange(n))

    def time_frame_from_records_generator_nrows(self):
        self.df = DataFrame.from_records(get_data(), nrows=1000)


class frame_get_dtype_counts(object):
    goal_time = 0.2

    def setup(self):
        self.df = pandas.DataFrame(np.random.randn(10, 10000))

    def time_frame_get_dtype_counts(self):
        self.df.get_dtype_counts()


class frame_getitem_single_column(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 1000))
        self.df2 = DataFrame(randn(3000, 1), columns=['A'])
        self.df3 = DataFrame(randn(3000, 1))

        def f():
            if hasattr(self.df, '_item_cache'):
                self.df._item_cache.clear()
            for (name, col) in self.df.iteritems():
                pass

        def g():
            for (name, col) in self.df.iteritems():
                pass

        def h():
            for i in xrange(10000):
                self.df2['A']

        def j():
            for i in xrange(10000):
                self.df3[0]

    def time_frame_getitem_single_column(self):
        h()


class frame_getitem_single_column2(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 1000))
        self.df2 = DataFrame(randn(3000, 1), columns=['A'])
        self.df3 = DataFrame(randn(3000, 1))

        def f():
            if hasattr(self.df, '_item_cache'):
                self.df._item_cache.clear()
            for (name, col) in self.df.iteritems():
                pass

        def g():
            for (name, col) in self.df.iteritems():
                pass

        def h():
            for i in xrange(10000):
                self.df2['A']

        def j():
            for i in xrange(10000):
                self.df3[0]

    def time_frame_getitem_single_column2(self):
        j()


class frame_html_repr_trunc_mi(object):
    goal_time = 0.2

    def setup(self):
        self.nrows = 10000
        self.data = randn(self.nrows, 10)
        self.idx = MultiIndex.from_arrays(np.tile(randn(3, (self.nrows / 100)), 100))
        self.df = DataFrame(self.data, index=self.idx)

    def time_frame_html_repr_trunc_mi(self):
        self.df._repr_html_()


class frame_html_repr_trunc_si(object):
    goal_time = 0.2

    def setup(self):
        self.nrows = 10000
        self.data = randn(self.nrows, 10)
        self.idx = randn(self.nrows)
        self.df = DataFrame(self.data, index=self.idx)

    def time_frame_html_repr_trunc_si(self):
        self.df._repr_html_()


class frame_insert_100_columns_begin(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000

        def f(K=100):
            self.df = DataFrame(index=range(self.N))
            self.new_col = np.random.randn(self.N)
            for i in range(K):
                self.df.insert(0, i, self.new_col)

    def time_frame_insert_100_columns_begin(self):
        f()


class frame_insert_500_columns_end(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000

        def f(K=500):
            self.df = DataFrame(index=range(self.N))
            self.new_col = np.random.randn(self.N)
            for i in range(K):
                self.df[i] = self.new_col

    def time_frame_insert_500_columns_end(self):
        f()


class frame_interpolate(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 100))
        self.df.values[::2] = np.nan

    def time_frame_interpolate(self):
        self.df.interpolate()


class frame_interpolate_some_good(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'A': np.arange(0, 10000), 'B': np.random.randint(0, 100, 10000), 'C': randn(10000), 'D': randn(10000), })
        self.df.loc[1::5, 'A'] = np.nan
        self.df.loc[1::5, 'C'] = np.nan

    def time_frame_interpolate_some_good(self):
        self.df.interpolate()


class frame_interpolate_some_good_infer(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'A': np.arange(0, 10000), 'B': np.random.randint(0, 100, 10000), 'C': randn(10000), 'D': randn(10000), })
        self.df.loc[1::5, 'A'] = np.nan
        self.df.loc[1::5, 'C'] = np.nan

    def time_frame_interpolate_some_good_infer(self):
        self.df.interpolate(downcast='infer')


class frame_isnull(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(1000, 1000)
        self.df = DataFrame(self.data)

    def time_frame_isnull(self):
        isnull(self.df)


class frame_iteritems(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 1000))
        self.df2 = DataFrame(randn(3000, 1), columns=['A'])
        self.df3 = DataFrame(randn(3000, 1))

        def f():
            if hasattr(self.df, '_item_cache'):
                self.df._item_cache.clear()
            for (name, col) in self.df.iteritems():
                pass

        def g():
            for (name, col) in self.df.iteritems():
                pass

        def h():
            for i in xrange(10000):
                self.df2['A']

        def j():
            for i in xrange(10000):
                self.df3[0]

    def time_frame_iteritems(self):
        f()


class frame_iteritems_cached(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 1000))
        self.df2 = DataFrame(randn(3000, 1), columns=['A'])
        self.df3 = DataFrame(randn(3000, 1))

        def f():
            if hasattr(self.df, '_item_cache'):
                self.df._item_cache.clear()
            for (name, col) in self.df.iteritems():
                pass

        def g():
            for (name, col) in self.df.iteritems():
                pass

        def h():
            for i in xrange(10000):
                self.df2['A']

        def j():
            for i in xrange(10000):
                self.df3[0]

    def time_frame_iteritems_cached(self):
        g()


class frame_mask_bools(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(1000, 500)
        self.df = DataFrame(self.data)
        self.df = self.df.where((self.df > 0))
        self.bools = (self.df > 0)
        self.mask = isnull(self.df)

    def time_frame_mask_bools(self):
        self.bools.mask(self.mask)


class frame_mask_floats(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(1000, 500)
        self.df = DataFrame(self.data)
        self.df = self.df.where((self.df > 0))
        self.bools = (self.df > 0)
        self.mask = isnull(self.df)

    def time_frame_mask_floats(self):
        self.bools.astype(float).mask(self.mask)


class frame_nonunique_equal(object):
    goal_time = 0.2

    def setup(self):

        def make_pair(frame):
            self.df = frame
            self.df2 = self.df.copy()
            self.df2.ix[((-1), (-1))] = np.nan
            return (self.df, self.df2)

        def test_equal(name):
            (self.df, self.df2) = pairs[name]
            return self.df.equals(self.df)

        def test_unequal(name):
            (self.df, self.df2) = pairs[name]
            return self.df.equals(self.df2)
        self.float_df = DataFrame(np.random.randn(1000, 1000))
        self.object_df = DataFrame(([(['foo'] * 1000)] * 1000))
        self.nonunique_cols = self.object_df.copy()
        self.nonunique_cols.columns = (['A'] * len(self.nonunique_cols.columns))
        self.pairs = dict([(name, make_pair(frame)) for (name, frame) in (('float_df', self.float_df), ('object_df', self.object_df), ('nonunique_cols', self.nonunique_cols))])

    def time_frame_nonunique_equal(self):
        test_equal('nonunique_cols')


class frame_nonunique_unequal(object):
    goal_time = 0.2

    def setup(self):

        def make_pair(frame):
            self.df = frame
            self.df2 = self.df.copy()
            self.df2.ix[((-1), (-1))] = np.nan
            return (self.df, self.df2)

        def test_equal(name):
            (self.df, self.df2) = pairs[name]
            return self.df.equals(self.df)

        def test_unequal(name):
            (self.df, self.df2) = pairs[name]
            return self.df.equals(self.df2)
        self.float_df = DataFrame(np.random.randn(1000, 1000))
        self.object_df = DataFrame(([(['foo'] * 1000)] * 1000))
        self.nonunique_cols = self.object_df.copy()
        self.nonunique_cols.columns = (['A'] * len(self.nonunique_cols.columns))
        self.pairs = dict([(name, make_pair(frame)) for (name, frame) in (('float_df', self.float_df), ('object_df', self.object_df), ('nonunique_cols', self.nonunique_cols))])

    def time_frame_nonunique_unequal(self):
        test_unequal('nonunique_cols')


class frame_object_equal(object):
    goal_time = 0.2

    def setup(self):

        def make_pair(frame):
            self.df = frame
            self.df2 = self.df.copy()
            self.df2.ix[((-1), (-1))] = np.nan
            return (self.df, self.df2)

        def test_equal(name):
            (self.df, self.df2) = pairs[name]
            return self.df.equals(self.df)

        def test_unequal(name):
            (self.df, self.df2) = pairs[name]
            return self.df.equals(self.df2)
        self.float_df = DataFrame(np.random.randn(1000, 1000))
        self.object_df = DataFrame(([(['foo'] * 1000)] * 1000))
        self.nonunique_cols = self.object_df.copy()
        self.nonunique_cols.columns = (['A'] * len(self.nonunique_cols.columns))
        self.pairs = dict([(name, make_pair(frame)) for (name, frame) in (('float_df', self.float_df), ('object_df', self.object_df), ('nonunique_cols', self.nonunique_cols))])

    def time_frame_object_equal(self):
        test_equal('object_df')


class frame_object_unequal(object):
    goal_time = 0.2

    def setup(self):

        def make_pair(frame):
            self.df = frame
            self.df2 = self.df.copy()
            self.df2.ix[((-1), (-1))] = np.nan
            return (self.df, self.df2)

        def test_equal(name):
            (self.df, self.df2) = pairs[name]
            return self.df.equals(self.df)

        def test_unequal(name):
            (self.df, self.df2) = pairs[name]
            return self.df.equals(self.df2)
        self.float_df = DataFrame(np.random.randn(1000, 1000))
        self.object_df = DataFrame(([(['foo'] * 1000)] * 1000))
        self.nonunique_cols = self.object_df.copy()
        self.nonunique_cols.columns = (['A'] * len(self.nonunique_cols.columns))
        self.pairs = dict([(name, make_pair(frame)) for (name, frame) in (('float_df', self.float_df), ('object_df', self.object_df), ('nonunique_cols', self.nonunique_cols))])

    def time_frame_object_unequal(self):
        test_unequal('object_df')


class frame_reindex_axis0(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 10000))
        self.idx = np.arange(4000, 7000)

    def time_frame_reindex_axis0(self):
        self.df.reindex(self.idx)


class frame_reindex_axis1(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 10000))
        self.idx = np.arange(4000, 7000)

    def time_frame_reindex_axis1(self):
        self.df.reindex(columns=self.idx)


class frame_reindex_both_axes(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 10000))
        self.idx = np.arange(4000, 7000)

    def time_frame_reindex_both_axes(self):
        self.df.reindex(index=self.idx, columns=self.idx)


class frame_reindex_both_axes_ix(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 10000))
        self.idx = np.arange(4000, 7000)

    def time_frame_reindex_both_axes_ix(self):
        self.df.ix[(self.idx, self.idx)]


class frame_reindex_upcast(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(dict([(c, {0: randint(0, 2, 1000).astype(np.bool_), 1: randint(0, 1000, 1000).astype(np.int16), 2: randint(0, 1000, 1000).astype(np.int32), 3: randint(0, 1000, 1000).astype(np.int64), }[randint(0, 4)]) for c in range(1000)]))

    def time_frame_reindex_upcast(self):
        self.df.reindex(permutation(range(1200)))


class frame_repr_tall(object):
    goal_time = 0.2

    def setup(self):
        self.df = pandas.DataFrame(np.random.randn(10000, 10))

    def time_frame_repr_tall(self):
        repr(self.df)


class frame_repr_wide(object):
    goal_time = 0.2

    def setup(self):
        self.df = pandas.DataFrame(np.random.randn(10, 10000))

    def time_frame_repr_wide(self):
        repr(self.df)


class frame_shift_axis0(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.rand(10000, 500))

    def time_frame_shift_axis0(self):
        self.df.shift(1, axis=0)


class frame_shift_axis_1(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.rand(10000, 500))

    def time_frame_shift_axis_1(self):
        self.df.shift(1, axis=1)


class frame_to_html_mixed(object):
    goal_time = 0.2

    def setup(self):
        self.nrows = 500
        self.df = DataFrame(randn(self.nrows, 10))
        self.df[0] = period_range('2000', '2010', self.nrows)
        self.df[1] = range(self.nrows)

    def time_frame_to_html_mixed(self):
        self.df.to_html()


class frame_to_string_floats(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(100, 10))

    def time_frame_to_string_floats(self):
        self.df.to_string()


class frame_xs_col(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(1, 100000))

    def time_frame_xs_col(self):
        self.df.xs(50000, axis=1)


class frame_xs_row(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(100000, 1))

    def time_frame_xs_row(self):
        self.df.xs(50000)


class series_string_vector_slice(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series((['abcdefg', np.nan] * 500000))

    def time_series_string_vector_slice(self):
        self.s.str[:5]