import string
import warnings

import numpy as np
import pandas.util.testing as tm
from pandas import (DataFrame, Series, MultiIndex, date_range, period_range,
                    isnull, NaT)

from .pandas_vb_common import setup  # noqa


class GetNumericData(object):

    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(10000, 25))
        self.df['foo'] = 'bar'
        self.df['bar'] = 'baz'
        with warnings.catch_warnings(record=True):
            self.df = self.df.consolidate()

    def time_frame_get_numeric_data(self):
        self.df._get_numeric_data()


class Lookup(object):

    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(10000, 8),
                            columns=list('abcdefgh'))
        self.df['foo'] = 'bar'
        self.row_labels = list(self.df.index[::10])[:900]
        self.col_labels = list(self.df.columns) * 100
        self.row_labels_all = np.array(
            list(self.df.index) * len(self.df.columns), dtype='object')
        self.col_labels_all = np.array(
            list(self.df.columns) * len(self.df.index), dtype='object')

    def time_frame_fancy_lookup(self):
        self.df.lookup(self.row_labels, self.col_labels)

    def time_frame_fancy_lookup_all(self):
        self.df.lookup(self.row_labels_all, self.col_labels_all)


class Reindex(object):

    goal_time = 0.2

    def setup(self):
        N = 10**3
        self.df = DataFrame(np.random.randn(N * 10, N))
        self.idx = np.arange(4 * N, 7 * N)
        self.df2 = DataFrame(
            {c: {0: np.random.randint(0, 2, N).astype(np.bool_),
                 1: np.random.randint(0, N, N).astype(np.int16),
                 2: np.random.randint(0, N, N).astype(np.int32),
                 3: np.random.randint(0, N, N).astype(np.int64)}
                [np.random.randint(0, 4)] for c in range(N)})

    def time_reindex_axis0(self):
        self.df.reindex(self.idx)

    def time_reindex_axis1(self):
        self.df.reindex(columns=self.idx)

    def time_reindex_both_axes(self):
        self.df.reindex(index=self.idx, columns=self.idx)

    def time_reindex_both_axes_ix(self):
        self.df.ix[self.idx, self.idx]

    def time_reindex_upcast(self):
        self.df2.reindex(np.random.permutation(range(1200)))


class Iteration(object):

    goal_time = 0.2

    def setup(self):
        N = 1000
        self.df = DataFrame(np.random.randn(N * 10, N))
        self.df2 = DataFrame(np.random.randn(N * 50, 10))
        self.df3 = DataFrame(np.random.randn(N, 5 * N),
                             columns=['C' + str(c) for c in range(N * 5)])

    def time_iteritems(self):
        # (monitor no-copying behaviour)
        if hasattr(self.df, '_item_cache'):
            self.df._item_cache.clear()
        for name, col in self.df.iteritems():
            pass

    def time_iteritems_cached(self):
        for name, col in self.df.iteritems():
            pass

    def time_iteritems_indexing(self):
        for col in self.df3:
            self.df3[col]

    def time_itertuples(self):
        for row in self.df2.itertuples():
            pass

    def time_iterrows(self):
        for row in self.df.iterrows():
            pass


class ToString(object):

    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(100, 10))

    def time_to_string_floats(self):
        self.df.to_string()


class ToHTML(object):

    goal_time = 0.2

    def setup(self):
        nrows = 500
        self.df2 = DataFrame(np.random.randn(nrows, 10))
        self.df2[0] = period_range('2000', periods=nrows)
        self.df2[1] = range(nrows)

    def time_to_html_mixed(self):
        self.df2.to_html()


class Repr(object):

    goal_time = 0.2

    def setup(self):
        nrows = 10000
        data = np.random.randn(nrows, 10)
        arrays = np.tile(np.random.randn(3, int(nrows / 100)), 100)
        idx = MultiIndex.from_arrays(arrays)
        self.df3 = DataFrame(data, index=idx)
        self.df4 = DataFrame(data, index=np.random.randn(nrows))
        self.df_tall = DataFrame(np.random.randn(nrows, 10))
        self.df_wide = DataFrame(np.random.randn(10, nrows))

    def time_html_repr_trunc_mi(self):
        self.df3._repr_html_()

    def time_html_repr_trunc_si(self):
        self.df4._repr_html_()

    def time_repr_tall(self):
        repr(self.df_tall)

    def time_frame_repr_wide(self):
        repr(self.df_wide)


class MaskBool(object):

    goal_time = 0.2

    def setup(self):
        data = np.random.randn(1000, 500)
        df = DataFrame(data)
        df = df.where(df > 0)
        self.bools = df > 0
        self.mask = isnull(df)

    def time_frame_mask_bools(self):
        self.bools.mask(self.mask)

    def time_frame_mask_floats(self):
        self.bools.astype(float).mask(self.mask)


class Isnull(object):

    goal_time = 0.2

    def setup(self):
        N = 10**3
        self.df_no_null = DataFrame(np.random.randn(N, N))

        sample = np.array([np.nan, 1.0])
        data = np.random.choice(sample, (N, N))
        self.df = DataFrame(data)

        sample = np.array(list(string.ascii_letters + string.whitespace))
        data = np.random.choice(sample, (N, N))
        self.df_strings = DataFrame(data)

        sample = np.array([NaT, np.nan, None, np.datetime64('NaT'),
                           np.timedelta64('NaT'), 0, 1, 2.0, '', 'abcd'])
        data = np.random.choice(sample, (N, N))
        self.df_obj = DataFrame(data)

    def time_isnull_floats_no_null(self):
        isnull(self.df_no_null)

    def time_isnull(self):
        isnull(self.df)

    def time_isnull_strngs(self):
        isnull(self.df_strings)

    def time_isnull_obj(self):
        isnull(self.df_obj)


class Fillna(object):

    goal_time = 0.2
    params = ([True, False], ['pad', 'bfill'])
    param_names = ['inplace', 'method']

    def setup(self, inplace, method):
        values = np.random.randn(10000, 100)
        values[::2] = np.nan
        self.df = DataFrame(values)

    def time_frame_fillna(self, inplace, method):
        self.df.fillna(inplace=inplace, method=method)


class Dropna(object):

    goal_time = 0.2
    params = (['all', 'any'], [0, 1])
    param_names = ['how', 'axis']

    def setup(self, how, axis):
        self.df = DataFrame(np.random.randn(10000, 1000))
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan
        self.df_mixed = self.df.copy()
        self.df_mixed['foo'] = 'bar'

    def time_dropna(self, how, axis):
        self.df.dropna(how=how, axis=axis)

    def time_dropna_axis_mixed_dtypes(self, how, axis):
        self.df_mixed.dropna(how=how, axis=axis)


class Count(object):

    goal_time = 0.2

    params = [0, 1]
    param_names = ['axis']

    def setup(self, axis):
        self.df = DataFrame(np.random.randn(10000, 1000))
        self.df.ix[50:1000, 20:50] = np.nan
        self.df.ix[2000:3000] = np.nan
        self.df.ix[:, 60:70] = np.nan
        self.df_mixed = self.df.copy()
        self.df_mixed['foo'] = 'bar'

        self.df.index = MultiIndex.from_arrays([self.df.index, self.df.index])
        self.df.columns = MultiIndex.from_arrays([self.df.columns,
                                                  self.df.columns])
        self.df_mixed.index = MultiIndex.from_arrays([self.df_mixed.index,
                                                      self.df_mixed.index])
        self.df_mixed.columns = MultiIndex.from_arrays([self.df_mixed.columns,
                                                        self.df_mixed.columns])

    def time_count_level_multi(self, axis):
        self.df.count(axis=axis, level=1)

    def time_count_level_mixed_dtypes_multi(self, axis):
        self.df_mixed.count(axis=axis, level=1)


class Apply(object):

    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 100))

        self.s = Series(np.arange(1028.0))
        self.df2 = DataFrame({i: self.s for i in range(1028)})
        self.df3 = DataFrame(np.random.randn(1000, 3), columns=list('ABC'))

    def time_apply_user_func(self):
        self.df2.apply(lambda x: np.corrcoef(x, self.s)[(0, 1)])

    def time_apply_axis_1(self):
        self.df.apply(lambda x: x + 1, axis=1)

    def time_apply_lambda_mean(self):
        self.df.apply(lambda x: x.mean())

    def time_apply_np_mean(self):
        self.df.apply(np.mean)

    def time_apply_pass_thru(self):
        self.df.apply(lambda x: x)

    def time_apply_ref_by_name(self):
        self.df3.apply(lambda x: x['A'] + x['B'], axis=1)


class Dtypes(object):

    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 1000))

    def time_frame_dtypes(self):
        self.df.dtypes


class Equals(object):

    goal_time = 0.2

    def setup(self):
        N = 10**3
        self.float_df = DataFrame(np.random.randn(N, N))
        self.float_df_nan = self.float_df.copy()
        self.float_df_nan.iloc[-1, -1] = np.nan

        self.object_df = DataFrame('foo', index=range(N), columns=range(N))
        self.object_df_nan = self.object_df.copy()
        self.object_df_nan.iloc[-1, -1] = np.nan

        self.nonunique_cols = self.object_df.copy()
        self.nonunique_cols.columns = ['A'] * len(self.nonunique_cols.columns)
        self.nonunique_cols_nan = self.nonunique_cols.copy()
        self.nonunique_cols_nan.iloc[-1, -1] = np.nan

    def time_frame_float_equal(self):
        self.float_df.equals(self.float_df)

    def time_frame_float_unequal(self):
        self.float_df.equals(self.float_df_nan)

    def time_frame_nonunique_equal(self):
        self.nonunique_cols.equals(self.nonunique_cols)

    def time_frame_nonunique_unequal(self):
        self.nonunique_cols.equals(self.nonunique_cols_nan)

    def time_frame_object_equal(self):
        self.object_df.equals(self.object_df)

    def time_frame_object_unequal(self):
        self.object_df.equals(self.object_df_nan)


class Interpolate(object):

    goal_time = 0.2
    params = [None, 'infer']
    param_names = ['downcast']

    def setup(self, downcast):
        N = 10000
        # this is the worst case, where every column has NaNs.
        self.df = DataFrame(np.random.randn(N, 100))
        self.df.values[::2] = np.nan

        self.df2 = DataFrame({'A': np.arange(0, N),
                              'B': np.random.randint(0, 100, N),
                              'C': np.random.randn(N),
                              'D': np.random.randn(N)})
        self.df2.loc[1::5, 'A'] = np.nan
        self.df2.loc[1::5, 'C'] = np.nan

    def time_interpolate(self, downcast):
        self.df.interpolate(downcast=downcast)

    def time_interpolate_some_good(self, downcast):
        self.df2.interpolate(downcast=downcast)


class Shift(object):
    # frame shift speedup issue-5609
    goal_time = 0.2
    params = [0, 1]
    param_names = ['axis']

    def setup(self, axis):
        self.df = DataFrame(np.random.rand(10000, 500))

    def time_shift(self, axis):
        self.df.shift(1, axis=axis)


class Nunique(object):

    def setup(self):
        self.df = DataFrame(np.random.randn(10000, 1000))

    def time_frame_nunique(self):
        self.df.nunique()


class Duplicated(object):

    goal_time = 0.2

    def setup(self):
        n = (1 << 20)
        t = date_range('2015-01-01', freq='S', periods=(n // 64))
        xs = np.random.randn(n // 64).round(2)
        self.df = DataFrame({'a': np.random.randint(-1 << 8, 1 << 8, n),
                             'b': np.random.choice(t, n),
                             'c': np.random.choice(xs, n)})
        self.df2 = DataFrame(np.random.randn(1000, 100).astype(str)).T

    def time_frame_duplicated(self):
        self.df.duplicated()

    def time_frame_duplicated_wide(self):
        self.df2.duplicated()


class XS(object):

    goal_time = 0.2
    params = [0, 1]
    param_names = ['axis']

    def setup(self, axis):
        self.N = 10**4
        self.df = DataFrame(np.random.randn(self.N, self.N))

    def time_frame_xs(self, axis):
        self.df.xs(self.N / 2, axis=axis)


class SortValues(object):

    goal_time = 0.2
    params = [True, False]
    param_names = ['ascending']

    def setup(self, ascending):
        self.df = DataFrame(np.random.randn(1000000, 2), columns=list('AB'))

    def time_frame_sort_values(self, ascending):
        self.df.sort_values(by='A', ascending=ascending)


class SortIndexByColumns(object):

    goal_time = 0.2

    def setup(self):
        N = 10000
        K = 10
        self.df = DataFrame({'key1': tm.makeStringIndex(N).values.repeat(K),
                             'key2': tm.makeStringIndex(N).values.repeat(K),
                             'value': np.random.randn(N * K)})

    def time_frame_sort_values_by_columns(self):
        self.df.sort_values(by=['key1', 'key2'])


class Quantile(object):

    goal_time = 0.2
    params = [0, 1]
    param_names = ['axis']

    def setup(self, axis):
        self.df = DataFrame(np.random.randn(1000, 3), columns=list('ABC'))

    def time_frame_quantile(self, axis):
        self.df.quantile([0.1, 0.5], axis=axis)


class GetDtypeCounts(object):
    # 2807
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(10, 10000))

    def time_frame_get_dtype_counts(self):
        self.df.get_dtype_counts()

    def time_info(self):
        self.df.info()


class NSort(object):

    goal_time = 0.2
    params = ['first', 'last', 'all']
    param_names = ['keep']

    def setup(self, keep):
        self.df = DataFrame(np.random.randn(1000, 3), columns=list('ABC'))

    def time_nlargest(self, keep):
        self.df.nlargest(100, 'A', keep=keep)

    def time_nsmallest(self, keep):
        self.df.nsmallest(100, 'A', keep=keep)


class Describe(object):

    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({
            'a': np.random.randint(0, 100, int(1e6)),
            'b': np.random.randint(0, 100, int(1e6)),
            'c': np.random.randint(0, 100, int(1e6))
        })

    def time_series_describe(self):
        self.df['a'].describe()

    def time_dataframe_describe(self):
        self.df.describe()
