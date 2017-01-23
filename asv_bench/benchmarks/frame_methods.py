from .pandas_vb_common import *
import string


#----------------------------------------------------------------------
# lookup

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

    def time_frame_fancy_lookup_all(self):
        self.df.lookup(self.row_labels_all, self.col_labels_all)


#----------------------------------------------------------------------
# reindex

class Reindex(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 1000))
        self.idx = np.arange(4000, 7000)

        self.df2 = DataFrame(
            dict([(c, {0: randint(0, 2, 1000).astype(np.bool_),
                       1: randint(0, 1000, 1000).astype(
                           np.int16),
                       2: randint(0, 1000, 1000).astype(
                           np.int32),
                       3: randint(0, 1000, 1000).astype(
                           np.int64),}[randint(0, 4)]) for c in
                  range(1000)]))

    def time_reindex_axis0(self):
        self.df.reindex(self.idx)

    def time_reindex_axis1(self):
        self.df.reindex(columns=self.idx)

    def time_reindex_both_axes(self):
        self.df.reindex(index=self.idx, columns=self.idx)

    def time_reindex_both_axes_ix(self):
        self.df.ix[(self.idx, self.idx)]

    def time_reindex_upcast(self):
        self.df2.reindex(permutation(range(1200)))


#----------------------------------------------------------------------
# iteritems (monitor no-copying behaviour)

class Iteration(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 1000))
        self.df2 = DataFrame(np.random.randn(50000, 10))
        self.df3 = pd.DataFrame(np.random.randn(1000,5000),
                                columns=['C'+str(c) for c in range(5000)])

    def f(self):
        if hasattr(self.df, '_item_cache'):
            self.df._item_cache.clear()
        for (name, col) in self.df.iteritems():
            pass

    def g(self):
        for (name, col) in self.df.iteritems():
            pass

    def time_iteritems(self):
        self.f()

    def time_iteritems_cached(self):
        self.g()

    def time_iteritems_indexing(self):
        df = self.df3
        for col in df:
            df[col]

    def time_itertuples(self):
        for row in self.df2.itertuples():
            pass


#----------------------------------------------------------------------
# to_string, to_html, repr

class Formatting(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(100, 10))

        self.nrows = 500
        self.df2 = DataFrame(randn(self.nrows, 10))
        self.df2[0] = period_range('2000', '2010', self.nrows)
        self.df2[1] = range(self.nrows)

        self.nrows = 10000
        self.data = randn(self.nrows, 10)
        self.idx = MultiIndex.from_arrays(np.tile(randn(3, int(self.nrows / 100)), 100))
        self.df3 = DataFrame(self.data, index=self.idx)
        self.idx = randn(self.nrows)
        self.df4 = DataFrame(self.data, index=self.idx)

        self.df_tall = pandas.DataFrame(np.random.randn(10000, 10))

        self.df_wide = pandas.DataFrame(np.random.randn(10, 10000))

    def time_to_string_floats(self):
        self.df.to_string()

    def time_to_html_mixed(self):
        self.df2.to_html()

    def time_html_repr_trunc_mi(self):
        self.df3._repr_html_()

    def time_html_repr_trunc_si(self):
        self.df4._repr_html_()

    def time_repr_tall(self):
        repr(self.df_tall)

    def time_frame_repr_wide(self):
        repr(self.df_wide)


#----------------------------------------------------------------------
# nulls/masking


## masking

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

    def time_frame_mask_floats(self):
        self.bools.astype(float).mask(self.mask)


## isnull

class FrameIsnull(object):
    goal_time = 0.2

    def setup(self):
        self.df_no_null = DataFrame(np.random.randn(1000, 1000))

        np.random.seed(1234)
        self.sample = np.array([np.nan, 1.0])
        self.data = np.random.choice(self.sample, (1000, 1000))
        self.df = DataFrame(self.data)

        np.random.seed(1234)
        self.sample = np.array(list(string.ascii_lowercase) +
                               list(string.ascii_uppercase) +
                               list(string.whitespace))
        self.data = np.random.choice(self.sample, (1000, 1000))
        self.df_strings= DataFrame(self.data)

        np.random.seed(1234)
        self.sample = np.array([NaT, np.nan, None, np.datetime64('NaT'),
                                np.timedelta64('NaT'), 0, 1, 2.0, '', 'abcd'])
        self.data = np.random.choice(self.sample, (1000, 1000))
        self.df_obj = DataFrame(self.data)

    def time_isnull_floats_no_null(self):
        isnull(self.df_no_null)

    def time_isnull(self):
        isnull(self.df)

    def time_isnull_strngs(self):
        isnull(self.df_strings)

    def time_isnull_obj(self):
        isnull(self.df_obj)


# ----------------------------------------------------------------------
# fillna in place

class frame_fillna_inplace(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 100))
        self.df.values[::2] = np.nan

    def time_frame_fillna_inplace(self):
        self.df.fillna(0, inplace=True)



class frame_fillna_many_columns_pad(object):
    goal_time = 0.2

    def setup(self):
        self.values = np.random.randn(1000, 1000)
        self.values[::2] = np.nan
        self.df = DataFrame(self.values)

    def time_frame_fillna_many_columns_pad(self):
        self.df.fillna(method='pad')



class Dropna(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(10000, 1000)
        self.df = DataFrame(self.data)
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan
        self.df_mixed = self.df.copy()
        self.df_mixed['foo'] = 'bar'

        self.df_mi = self.df.copy()
        self.df_mi.index = MultiIndex.from_tuples(self.df_mi.index.map((lambda x: (x, x))))
        self.df_mi.columns = MultiIndex.from_tuples(self.df_mi.columns.map((lambda x: (x, x))))

        self.df_mixed_mi = self.df_mixed.copy()
        self.df_mixed_mi.index = MultiIndex.from_tuples(self.df_mixed_mi.index.map((lambda x: (x, x))))
        self.df_mixed_mi.columns = MultiIndex.from_tuples(self.df_mixed_mi.columns.map((lambda x: (x, x))))

    def time_dropna_axis0_all(self):
        self.df.dropna(how='all', axis=0)

    def time_dropna_axis0_any(self):
        self.df.dropna(how='any', axis=0)

    def time_dropna_axis1_all(self):
        self.df.dropna(how='all', axis=1)

    def time_dropna_axis1_any(self):
        self.df.dropna(how='any', axis=1)

    def time_dropna_axis0_all_mixed_dtypes(self):
        self.df_mixed.dropna(how='all', axis=0)

    def time_dropna_axis0_any_mixed_dtypes(self):
        self.df_mixed.dropna(how='any', axis=0)

    def time_dropna_axis1_all_mixed_dtypes(self):
        self.df_mixed.dropna(how='all', axis=1)

    def time_dropna_axis1_any_mixed_dtypes(self):
        self.df_mixed.dropna(how='any', axis=1)

    def time_count_level_axis0_multi(self):
        self.df_mi.count(axis=0, level=1)

    def time_count_level_axis1_multi(self):
        self.df_mi.count(axis=1, level=1)

    def time_count_level_axis0_mixed_dtypes_multi(self):
        self.df_mixed_mi.count(axis=0, level=1)

    def time_count_level_axis1_mixed_dtypes_multi(self):
        self.df_mixed_mi.count(axis=1, level=1)


class Apply(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 100))

        self.s = Series(np.arange(1028.0))
        self.df2 = DataFrame({i: self.s for i in range(1028)})

        self.df3 = DataFrame(np.random.randn(1000, 3), columns=list('ABC'))

    def time_apply_user_func(self):
        self.df2.apply((lambda x: np.corrcoef(x, self.s)[(0, 1)]))

    def time_apply_axis_1(self):
        self.df.apply((lambda x: (x + 1)), axis=1)

    def time_apply_lambda_mean(self):
        self.df.apply((lambda x: x.mean()))

    def time_apply_np_mean(self):
        self.df.apply(np.mean)

    def time_apply_pass_thru(self):
        self.df.apply((lambda x: x))

    def time_apply_ref_by_name(self):
        self.df3.apply((lambda x: (x['A'] + x['B'])), axis=1)


#----------------------------------------------------------------------
# dtypes

class frame_dtypes(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 1000))

    def time_frame_dtypes(self):
        self.df.dtypes

#----------------------------------------------------------------------
# equals

class Equals(object):
    goal_time = 0.2

    def setup(self):
        self.float_df = DataFrame(np.random.randn(1000, 1000))
        self.object_df = DataFrame(([(['foo'] * 1000)] * 1000))
        self.nonunique_cols = self.object_df.copy()
        self.nonunique_cols.columns = (['A'] * len(self.nonunique_cols.columns))
        self.pairs = dict([(name, self.make_pair(frame)) for (name, frame) in (
            ('float_df', self.float_df), ('object_df', self.object_df),
            ('nonunique_cols', self.nonunique_cols))])

    def make_pair(self, frame):
        self.df = frame
        self.df2 = self.df.copy()
        self.df2.ix[((-1), (-1))] = np.nan
        return (self.df, self.df2)

    def test_equal(self, name):
        (self.df, self.df2) = self.pairs[name]
        return self.df.equals(self.df)

    def test_unequal(self, name):
        (self.df, self.df2) = self.pairs[name]
        return self.df.equals(self.df2)

    def time_frame_float_equal(self):
        self.test_equal('float_df')

    def time_frame_float_unequal(self):
        self.test_unequal('float_df')

    def time_frame_nonunique_equal(self):
        self.test_equal('nonunique_cols')

    def time_frame_nonunique_unequal(self):
        self.test_unequal('nonunique_cols')

    def time_frame_object_equal(self):
        self.test_equal('object_df')

    def time_frame_object_unequal(self):
        self.test_unequal('object_df')


class Interpolate(object):
    goal_time = 0.2

    def setup(self):
        # this is the worst case, where every column has NaNs.
        self.df = DataFrame(randn(10000, 100))
        self.df.values[::2] = np.nan

        self.df2 = DataFrame(
            {'A': np.arange(0, 10000), 'B': np.random.randint(0, 100, 10000),
             'C': randn(10000), 'D': randn(10000),})
        self.df2.loc[1::5, 'A'] = np.nan
        self.df2.loc[1::5, 'C'] = np.nan

    def time_interpolate(self):
        self.df.interpolate()

    def time_interpolate_some_good(self):
        self.df2.interpolate()

    def time_interpolate_some_good_infer(self):
        self.df2.interpolate(downcast='infer')


class Shift(object):
    # frame shift speedup issue-5609
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.rand(10000, 500))

    def time_shift_axis0(self):
        self.df.shift(1, axis=0)

    def time_shift_axis_1(self):
        self.df.shift(1, axis=1)


#-----------------------------------------------------------------------------
# from_records issue-6700

class frame_from_records_generator(object):
    goal_time = 0.2

    def get_data(self, n=100000):
        return ((x, (x * 20), (x * 100)) for x in range(n))

    def time_frame_from_records_generator(self):
        self.df = DataFrame.from_records(self.get_data())

    def time_frame_from_records_generator_nrows(self):
        self.df = DataFrame.from_records(self.get_data(), nrows=1000)



#-----------------------------------------------------------------------------
# nunique

class frame_nunique(object):

    def setup(self):
        self.data = np.random.randn(10000, 1000)
        self.df = DataFrame(self.data)

    def time_frame_nunique(self):
        self.df.nunique()



#-----------------------------------------------------------------------------
# duplicated

class frame_duplicated(object):
    goal_time = 0.2

    def setup(self):
        self.n = (1 << 20)
        self.t = date_range('2015-01-01', freq='S', periods=(self.n // 64))
        self.xs = np.random.randn((self.n // 64)).round(2)
        self.df = DataFrame({'a': np.random.randint(((-1) << 8), (1 << 8), self.n), 'b': np.random.choice(self.t, self.n), 'c': np.random.choice(self.xs, self.n), })

        self.df2 = DataFrame(np.random.randn(1000, 100).astype(str))

    def time_frame_duplicated(self):
        self.df.duplicated()

    def time_frame_duplicated_wide(self):
        self.df2.T.duplicated()

















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


class frame_sort_index(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(1000000, 2), columns=list('AB'))

    def time_frame_sort_index(self):
        self.df.sort_index()


class frame_sort_index_by_columns(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.K = 10
        self.key1 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.key2 = tm.makeStringIndex(self.N).values.repeat(self.K)
        self.df = DataFrame({'key1': self.key1, 'key2': self.key2, 'value': np.random.randn((self.N * self.K)), })
        self.col_array_list = list(self.df.values.T)

    def time_frame_sort_index_by_columns(self):
        self.df.sort_index(by=['key1', 'key2'])


class frame_quantile_axis1(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 3),
                            columns=list('ABC'))

    def time_frame_quantile_axis1(self):
        self.df.quantile([0.1, 0.5], axis=1)


#----------------------------------------------------------------------
# boolean indexing

class frame_boolean_row_select(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 100))
        self.bool_arr = np.zeros(10000, dtype=bool)
        self.bool_arr[:1000] = True

    def time_frame_boolean_row_select(self):
        self.df[self.bool_arr]

class frame_getitem_single_column(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(randn(10000, 1000))
        self.df2 = DataFrame(randn(3000, 1), columns=['A'])
        self.df3 = DataFrame(randn(3000, 1))

    def h(self):
        for i in range(10000):
            self.df2['A']

    def j(self):
        for i in range(10000):
            self.df3[0]

    def time_frame_getitem_single_column(self):
        self.h()

    def time_frame_getitem_single_column2(self):
        self.j()


#----------------------------------------------------------------------
# assignment

class frame_assign_timeseries_index(object):
    goal_time = 0.2

    def setup(self):
        self.idx = date_range('1/1/2000', periods=100000, freq='D')
        self.df = DataFrame(randn(100000, 1), columns=['A'], index=self.idx)

    def time_frame_assign_timeseries_index(self):
        self.f(self.df)

    def f(self, df):
        self.x = self.df.copy()
        self.x['date'] = self.x.index



# insert many columns

class frame_insert_100_columns_begin(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000

    def f(self, K=100):
        self.df = DataFrame(index=range(self.N))
        self.new_col = np.random.randn(self.N)
        for i in range(K):
            self.df.insert(0, i, self.new_col)

    def g(self, K=500):
        self.df = DataFrame(index=range(self.N))
        self.new_col = np.random.randn(self.N)
        for i in range(K):
            self.df[i] = self.new_col

    def time_frame_insert_100_columns_begin(self):
        self.f()

    def time_frame_insert_500_columns_end(self):
        self.g()



#----------------------------------------------------------------------
# strings methods, #2602

class series_string_vector_slice(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series((['abcdefg', np.nan] * 500000))

    def time_series_string_vector_slice(self):
        self.s.str[:5]


#----------------------------------------------------------------------
# df.info() and get_dtype_counts() # 2807

class frame_get_dtype_counts(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(10, 10000))

    def time_frame_get_dtype_counts(self):
        self.df.get_dtype_counts()


class frame_nlargest(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 3),
                            columns=list('ABC'))

    def time_frame_nlargest(self):
        self.df.nlargest(100, 'A')
