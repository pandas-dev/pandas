# -*- coding: utf-8 -*-
# pylint: disable-msg=W0612,E1101
import nose

from pandas import DataFrame, Series
from pandas.core.sparse import SparseDataFrame
import pandas as pd

from numpy import nan
import numpy as np

from pandas.util.testing import assert_frame_equal

from pandas.core.reshape import (melt, lreshape, get_dummies, wide_to_long)
import pandas.util.testing as tm
from pandas.compat import range, u

_multiprocess_can_split_ = True


class TestMelt(tm.TestCase):

    def setUp(self):
        self.df = tm.makeTimeDataFrame()[:10]
        self.df['id1'] = (self.df['A'] > 0).astype(np.int64)
        self.df['id2'] = (self.df['B'] > 0).astype(np.int64)

        self.var_name = 'var'
        self.value_name = 'val'

        self.df1 = pd.DataFrame([[1.067683, -1.110463, 0.20867
                                  ], [-1.321405, 0.368915, -1.055342],
                                 [-0.807333, 0.08298, -0.873361]])
        self.df1.columns = [list('ABC'), list('abc')]
        self.df1.columns.names = ['CAP', 'low']

    def test_default_col_names(self):
        result = melt(self.df)
        self.assertEqual(result.columns.tolist(), ['variable', 'value'])

        result1 = melt(self.df, id_vars=['id1'])
        self.assertEqual(result1.columns.tolist(), ['id1', 'variable', 'value'
                                                    ])

        result2 = melt(self.df, id_vars=['id1', 'id2'])
        self.assertEqual(result2.columns.tolist(), ['id1', 'id2', 'variable',
                                                    'value'])

    def test_value_vars(self):
        result3 = melt(self.df, id_vars=['id1', 'id2'], value_vars='A')
        self.assertEqual(len(result3), 10)

        result4 = melt(self.df, id_vars=['id1', 'id2'], value_vars=['A', 'B'])
        expected4 = DataFrame({'id1': self.df['id1'].tolist() * 2,
                               'id2': self.df['id2'].tolist() * 2,
                               'variable': ['A'] * 10 + ['B'] * 10,
                               'value': (self.df['A'].tolist() +
                                         self.df['B'].tolist())},
                              columns=['id1', 'id2', 'variable', 'value'])
        tm.assert_frame_equal(result4, expected4)

    def test_custom_var_name(self):
        result5 = melt(self.df, var_name=self.var_name)
        self.assertEqual(result5.columns.tolist(), ['var', 'value'])

        result6 = melt(self.df, id_vars=['id1'], var_name=self.var_name)
        self.assertEqual(result6.columns.tolist(), ['id1', 'var', 'value'])

        result7 = melt(self.df, id_vars=['id1', 'id2'], var_name=self.var_name)
        self.assertEqual(result7.columns.tolist(), ['id1', 'id2', 'var',
                                                    'value'])

        result8 = melt(self.df, id_vars=['id1', 'id2'], value_vars='A',
                       var_name=self.var_name)
        self.assertEqual(result8.columns.tolist(), ['id1', 'id2', 'var',
                                                    'value'])

        result9 = melt(self.df, id_vars=['id1', 'id2'], value_vars=['A', 'B'],
                       var_name=self.var_name)
        expected9 = DataFrame({'id1': self.df['id1'].tolist() * 2,
                               'id2': self.df['id2'].tolist() * 2,
                               self.var_name: ['A'] * 10 + ['B'] * 10,
                               'value': (self.df['A'].tolist() +
                                         self.df['B'].tolist())},
                              columns=['id1', 'id2', self.var_name, 'value'])
        tm.assert_frame_equal(result9, expected9)

    def test_custom_value_name(self):
        result10 = melt(self.df, value_name=self.value_name)
        self.assertEqual(result10.columns.tolist(), ['variable', 'val'])

        result11 = melt(self.df, id_vars=['id1'], value_name=self.value_name)
        self.assertEqual(result11.columns.tolist(), ['id1', 'variable', 'val'])

        result12 = melt(self.df, id_vars=['id1', 'id2'],
                        value_name=self.value_name)
        self.assertEqual(result12.columns.tolist(), ['id1', 'id2', 'variable',
                                                     'val'])

        result13 = melt(self.df, id_vars=['id1', 'id2'], value_vars='A',
                        value_name=self.value_name)
        self.assertEqual(result13.columns.tolist(), ['id1', 'id2', 'variable',
                                                     'val'])

        result14 = melt(self.df, id_vars=['id1', 'id2'], value_vars=['A', 'B'],
                        value_name=self.value_name)
        expected14 = DataFrame({'id1': self.df['id1'].tolist() * 2,
                                'id2': self.df['id2'].tolist() * 2,
                                'variable': ['A'] * 10 + ['B'] * 10,
                                self.value_name: (self.df['A'].tolist() +
                                                  self.df['B'].tolist())},
                               columns=['id1', 'id2', 'variable',
                                        self.value_name])
        tm.assert_frame_equal(result14, expected14)

    def test_custom_var_and_value_name(self):

        result15 = melt(self.df, var_name=self.var_name,
                        value_name=self.value_name)
        self.assertEqual(result15.columns.tolist(), ['var', 'val'])

        result16 = melt(self.df, id_vars=['id1'], var_name=self.var_name,
                        value_name=self.value_name)
        self.assertEqual(result16.columns.tolist(), ['id1', 'var', 'val'])

        result17 = melt(self.df, id_vars=['id1', 'id2'],
                        var_name=self.var_name, value_name=self.value_name)
        self.assertEqual(result17.columns.tolist(), ['id1', 'id2', 'var', 'val'
                                                     ])

        result18 = melt(self.df, id_vars=['id1', 'id2'], value_vars='A',
                        var_name=self.var_name, value_name=self.value_name)
        self.assertEqual(result18.columns.tolist(), ['id1', 'id2', 'var', 'val'
                                                     ])

        result19 = melt(self.df, id_vars=['id1', 'id2'], value_vars=['A', 'B'],
                        var_name=self.var_name, value_name=self.value_name)
        expected19 = DataFrame({'id1': self.df['id1'].tolist() * 2,
                                'id2': self.df['id2'].tolist() * 2,
                                self.var_name: ['A'] * 10 + ['B'] * 10,
                                self.value_name: (self.df['A'].tolist() +
                                                  self.df['B'].tolist())},
                               columns=['id1', 'id2', self.var_name,
                                        self.value_name])
        tm.assert_frame_equal(result19, expected19)

        df20 = self.df.copy()
        df20.columns.name = 'foo'
        result20 = melt(df20)
        self.assertEqual(result20.columns.tolist(), ['foo', 'value'])

    def test_col_level(self):
        res1 = melt(self.df1, col_level=0)
        res2 = melt(self.df1, col_level='CAP')
        self.assertEqual(res1.columns.tolist(), ['CAP', 'value'])
        self.assertEqual(res2.columns.tolist(), ['CAP', 'value'])

    def test_multiindex(self):
        res = pd.melt(self.df1)
        self.assertEqual(res.columns.tolist(), ['CAP', 'low', 'value'])


class TestGetDummies(tm.TestCase):

    sparse = False

    def setUp(self):
        self.df = DataFrame({'A': ['a', 'b', 'a'],
                             'B': ['b', 'b', 'c'],
                             'C': [1, 2, 3]})

    def test_basic(self):
        s_list = list('abc')
        s_series = Series(s_list)
        s_series_index = Series(s_list, list('ABC'))

        expected = DataFrame({'a': {0: 1,
                                    1: 0,
                                    2: 0},
                              'b': {0: 0,
                                    1: 1,
                                    2: 0},
                              'c': {0: 0,
                                    1: 0,
                                    2: 1}}, dtype=np.uint8)
        assert_frame_equal(get_dummies(s_list, sparse=self.sparse), expected)
        assert_frame_equal(get_dummies(s_series, sparse=self.sparse), expected)

        expected.index = list('ABC')
        assert_frame_equal(
            get_dummies(s_series_index, sparse=self.sparse), expected)

    def test_basic_types(self):
        # GH 10531
        s_list = list('abc')
        s_series = Series(s_list)
        s_df = DataFrame({'a': [0, 1, 0, 1, 2],
                          'b': ['A', 'A', 'B', 'C', 'C'],
                          'c': [2, 3, 3, 3, 2]})

        if not self.sparse:
            exp_df_type = DataFrame
            exp_blk_type = pd.core.internals.IntBlock
        else:
            exp_df_type = SparseDataFrame
            exp_blk_type = pd.core.internals.SparseBlock

        self.assertEqual(
            type(get_dummies(s_list, sparse=self.sparse)), exp_df_type)
        self.assertEqual(
            type(get_dummies(s_series, sparse=self.sparse)), exp_df_type)

        r = get_dummies(s_df, sparse=self.sparse, columns=s_df.columns)
        self.assertEqual(type(r), exp_df_type)

        r = get_dummies(s_df, sparse=self.sparse, columns=['a'])
        self.assertEqual(type(r[['a_0']]._data.blocks[0]), exp_blk_type)
        self.assertEqual(type(r[['a_1']]._data.blocks[0]), exp_blk_type)
        self.assertEqual(type(r[['a_2']]._data.blocks[0]), exp_blk_type)

    def test_just_na(self):
        just_na_list = [np.nan]
        just_na_series = Series(just_na_list)
        just_na_series_index = Series(just_na_list, index=['A'])

        res_list = get_dummies(just_na_list, sparse=self.sparse)
        res_series = get_dummies(just_na_series, sparse=self.sparse)
        res_series_index = get_dummies(just_na_series_index,
                                       sparse=self.sparse)

        self.assertEqual(res_list.empty, True)
        self.assertEqual(res_series.empty, True)
        self.assertEqual(res_series_index.empty, True)

        self.assertEqual(res_list.index.tolist(), [0])
        self.assertEqual(res_series.index.tolist(), [0])
        self.assertEqual(res_series_index.index.tolist(), ['A'])

    def test_include_na(self):
        s = ['a', 'b', np.nan]
        res = get_dummies(s, sparse=self.sparse)
        exp = DataFrame({'a': {0: 1, 1: 0, 2: 0},
                         'b': {0: 0, 1: 1, 2: 0}}, dtype=np.uint8)
        assert_frame_equal(res, exp)

        # Sparse dataframes do not allow nan labelled columns, see #GH8822
        res_na = get_dummies(s, dummy_na=True, sparse=self.sparse)
        exp_na = DataFrame({nan: {0: 0, 1: 0, 2: 1},
                            'a': {0: 1, 1: 0, 2: 0},
                            'b': {0: 0, 1: 1, 2: 0}},
                           dtype=np.uint8)
        exp_na = exp_na.reindex_axis(['a', 'b', nan], 1)
        # hack (NaN handling in assert_index_equal)
        exp_na.columns = res_na.columns
        assert_frame_equal(res_na, exp_na)

        res_just_na = get_dummies([nan], dummy_na=True, sparse=self.sparse)
        exp_just_na = DataFrame(Series(1, index=[0]), columns=[nan],
                                dtype=np.uint8)
        tm.assert_numpy_array_equal(res_just_na.values, exp_just_na.values)

    def test_unicode(self
                     ):  # See GH 6885 - get_dummies chokes on unicode values
        import unicodedata
        e = 'e'
        eacute = unicodedata.lookup('LATIN SMALL LETTER E WITH ACUTE')
        s = [e, eacute, eacute]
        res = get_dummies(s, prefix='letter', sparse=self.sparse)
        exp = DataFrame({'letter_e': {0: 1,
                                      1: 0,
                                      2: 0},
                         u('letter_%s') % eacute: {0: 0,
                                                   1: 1,
                                                   2: 1}},
                        dtype=np.uint8)
        assert_frame_equal(res, exp)

    def test_dataframe_dummies_all_obj(self):
        df = self.df[['A', 'B']]
        result = get_dummies(df, sparse=self.sparse)
        expected = DataFrame({'A_a': [1, 0, 1],
                              'A_b': [0, 1, 0],
                              'B_b': [1, 1, 0],
                              'B_c': [0, 0, 1]}, dtype=np.uint8)
        assert_frame_equal(result, expected)

    def test_dataframe_dummies_mix_default(self):
        df = self.df
        result = get_dummies(df, sparse=self.sparse)
        expected = DataFrame({'C': [1, 2, 3],
                              'A_a': [1, 0, 1],
                              'A_b': [0, 1, 0],
                              'B_b': [1, 1, 0],
                              'B_c': [0, 0, 1]})
        cols = ['A_a', 'A_b', 'B_b', 'B_c']
        expected[cols] = expected[cols].astype(np.uint8)
        expected = expected[['C', 'A_a', 'A_b', 'B_b', 'B_c']]
        assert_frame_equal(result, expected)

    def test_dataframe_dummies_prefix_list(self):
        prefixes = ['from_A', 'from_B']
        df = DataFrame({'A': ['a', 'b', 'a'],
                        'B': ['b', 'b', 'c'],
                        'C': [1, 2, 3]})
        result = get_dummies(df, prefix=prefixes, sparse=self.sparse)
        expected = DataFrame({'C': [1, 2, 3],
                              'from_A_a': [1, 0, 1],
                              'from_A_b': [0, 1, 0],
                              'from_B_b': [1, 1, 0],
                              'from_B_c': [0, 0, 1]})
        cols = expected.columns[1:]
        expected[cols] = expected[cols].astype(np.uint8)
        expected = expected[['C', 'from_A_a', 'from_A_b', 'from_B_b',
                             'from_B_c']]
        assert_frame_equal(result, expected)

    def test_dataframe_dummies_prefix_str(self):
        # not that you should do this...
        df = self.df
        result = get_dummies(df, prefix='bad', sparse=self.sparse)
        expected = DataFrame([[1, 1, 0, 1, 0],
                              [2, 0, 1, 1, 0],
                              [3, 1, 0, 0, 1]],
                             columns=['C', 'bad_a', 'bad_b', 'bad_b', 'bad_c'],
                             dtype=np.uint8)
        expected = expected.astype({"C": np.int64})
        assert_frame_equal(result, expected)

    def test_dataframe_dummies_subset(self):
        df = self.df
        result = get_dummies(df, prefix=['from_A'], columns=['A'],
                             sparse=self.sparse)
        expected = DataFrame({'from_A_a': [1, 0, 1],
                              'from_A_b': [0, 1, 0],
                              'B': ['b', 'b', 'c'],
                              'C': [1, 2, 3]})
        cols = ['from_A_a', 'from_A_b']
        expected[cols] = expected[cols].astype(np.uint8)
        assert_frame_equal(result, expected)

    def test_dataframe_dummies_prefix_sep(self):
        df = self.df
        result = get_dummies(df, prefix_sep='..', sparse=self.sparse)
        expected = DataFrame({'C': [1, 2, 3],
                              'A..a': [1, 0, 1],
                              'A..b': [0, 1, 0],
                              'B..b': [1, 1, 0],
                              'B..c': [0, 0, 1]})
        expected = expected[['C', 'A..a', 'A..b', 'B..b', 'B..c']]
        cols = expected.columns[1:]
        expected[cols] = expected[cols].astype(np.uint8)
        assert_frame_equal(result, expected)

        result = get_dummies(df, prefix_sep=['..', '__'], sparse=self.sparse)
        expected = expected.rename(columns={'B..b': 'B__b', 'B..c': 'B__c'})
        assert_frame_equal(result, expected)

        result = get_dummies(df, prefix_sep={'A': '..',
                                             'B': '__'}, sparse=self.sparse)
        assert_frame_equal(result, expected)

    def test_dataframe_dummies_prefix_bad_length(self):
        with tm.assertRaises(ValueError):
            get_dummies(self.df, prefix=['too few'], sparse=self.sparse)

    def test_dataframe_dummies_prefix_sep_bad_length(self):
        with tm.assertRaises(ValueError):
            get_dummies(self.df, prefix_sep=['bad'], sparse=self.sparse)

    def test_dataframe_dummies_prefix_dict(self):
        prefixes = {'A': 'from_A', 'B': 'from_B'}
        df = DataFrame({'A': ['a', 'b', 'a'],
                        'B': ['b', 'b', 'c'],
                        'C': [1, 2, 3]})
        result = get_dummies(df, prefix=prefixes, sparse=self.sparse)
        expected = DataFrame({'from_A_a': [1, 0, 1],
                              'from_A_b': [0, 1, 0],
                              'from_B_b': [1, 1, 0],
                              'from_B_c': [0, 0, 1],
                              'C': [1, 2, 3]})
        cols = ['from_A_a', 'from_A_b', 'from_B_b', 'from_B_c']
        expected[cols] = expected[cols].astype(np.uint8)
        assert_frame_equal(result, expected)

    def test_dataframe_dummies_with_na(self):
        df = self.df
        df.loc[3, :] = [np.nan, np.nan, np.nan]
        result = get_dummies(df, dummy_na=True, sparse=self.sparse)
        expected = DataFrame({'C': [1, 2, 3, np.nan],
                              'A_a': [1, 0, 1, 0],
                              'A_b': [0, 1, 0, 0],
                              'A_nan': [0, 0, 0, 1],
                              'B_b': [1, 1, 0, 0],
                              'B_c': [0, 0, 1, 0],
                              'B_nan': [0, 0, 0, 1]})
        cols = ['A_a', 'A_b', 'A_nan', 'B_b', 'B_c', 'B_nan']
        expected[cols] = expected[cols].astype(np.uint8)
        expected = expected[['C', 'A_a', 'A_b', 'A_nan',
                             'B_b', 'B_c', 'B_nan']]
        assert_frame_equal(result, expected)

        result = get_dummies(df, dummy_na=False, sparse=self.sparse)
        expected = expected[['C', 'A_a', 'A_b', 'B_b', 'B_c']]
        assert_frame_equal(result, expected)

    def test_dataframe_dummies_with_categorical(self):
        df = self.df
        df['cat'] = pd.Categorical(['x', 'y', 'y'])
        result = get_dummies(df, sparse=self.sparse)
        expected = DataFrame({'C': [1, 2, 3],
                              'A_a': [1, 0, 1],
                              'A_b': [0, 1, 0],
                              'B_b': [1, 1, 0],
                              'B_c': [0, 0, 1],
                              'cat_x': [1, 0, 0],
                              'cat_y': [0, 1, 1]})
        cols = ['A_a', 'A_b', 'B_b', 'B_c', 'cat_x', 'cat_y']
        expected[cols] = expected[cols].astype(np.uint8)
        expected = expected[['C', 'A_a', 'A_b', 'B_b', 'B_c',
                             'cat_x', 'cat_y']]
        assert_frame_equal(result, expected)

    # GH12402 Add a new parameter `drop_first` to avoid collinearity
    def test_basic_drop_first(self):
        # Basic case
        s_list = list('abc')
        s_series = Series(s_list)
        s_series_index = Series(s_list, list('ABC'))

        expected = DataFrame({'b': {0: 0,
                                    1: 1,
                                    2: 0},
                              'c': {0: 0,
                                    1: 0,
                                    2: 1}}, dtype=np.uint8)

        result = get_dummies(s_list, sparse=self.sparse, drop_first=True)
        assert_frame_equal(result, expected)

        result = get_dummies(s_series, sparse=self.sparse, drop_first=True)
        assert_frame_equal(result, expected)

        expected.index = list('ABC')
        result = get_dummies(s_series_index, sparse=self.sparse,
                             drop_first=True)
        assert_frame_equal(result, expected)

    def test_basic_drop_first_one_level(self):
        # Test the case that categorical variable only has one level.
        s_list = list('aaa')
        s_series = Series(s_list)
        s_series_index = Series(s_list, list('ABC'))

        expected = DataFrame(index=np.arange(3))

        result = get_dummies(s_list, sparse=self.sparse, drop_first=True)
        assert_frame_equal(result, expected)

        result = get_dummies(s_series, sparse=self.sparse, drop_first=True)
        assert_frame_equal(result, expected)

        expected = DataFrame(index=list('ABC'))
        result = get_dummies(s_series_index, sparse=self.sparse,
                             drop_first=True)
        assert_frame_equal(result, expected)

    def test_basic_drop_first_NA(self):
        # Test NA hadling together with drop_first
        s_NA = ['a', 'b', np.nan]
        res = get_dummies(s_NA, sparse=self.sparse, drop_first=True)
        exp = DataFrame({'b': {0: 0,
                               1: 1,
                               2: 0}}, dtype=np.uint8)
        assert_frame_equal(res, exp)

        res_na = get_dummies(s_NA, dummy_na=True, sparse=self.sparse,
                             drop_first=True)
        exp_na = DataFrame({'b': {0: 0,
                                  1: 1,
                                  2: 0},
                            nan: {0: 0,
                                  1: 0,
                                  2: 1}}, dtype=np.uint8).reindex_axis(
                                      ['b', nan], 1)
        assert_frame_equal(res_na, exp_na)

        res_just_na = get_dummies([nan], dummy_na=True, sparse=self.sparse,
                                  drop_first=True)
        exp_just_na = DataFrame(index=np.arange(1))
        assert_frame_equal(res_just_na, exp_just_na)

    def test_dataframe_dummies_drop_first(self):
        df = self.df[['A', 'B']]
        result = get_dummies(df, sparse=self.sparse, drop_first=True)
        expected = DataFrame({'A_b': [0, 1, 0],
                              'B_c': [0, 0, 1]}, dtype=np.uint8)
        assert_frame_equal(result, expected)

    def test_dataframe_dummies_drop_first_with_categorical(self):
        df = self.df
        df['cat'] = pd.Categorical(['x', 'y', 'y'])
        result = get_dummies(df, sparse=self.sparse, drop_first=True)
        expected = DataFrame({'C': [1, 2, 3],
                              'A_b': [0, 1, 0],
                              'B_c': [0, 0, 1],
                              'cat_y': [0, 1, 1]})
        cols = ['A_b', 'B_c', 'cat_y']
        expected[cols] = expected[cols].astype(np.uint8)
        expected = expected[['C', 'A_b', 'B_c', 'cat_y']]
        assert_frame_equal(result, expected)

    def test_dataframe_dummies_drop_first_with_na(self):
        df = self.df
        df.loc[3, :] = [np.nan, np.nan, np.nan]
        result = get_dummies(df, dummy_na=True, sparse=self.sparse,
                             drop_first=True)
        expected = DataFrame({'C': [1, 2, 3, np.nan],
                              'A_b': [0, 1, 0, 0],
                              'A_nan': [0, 0, 0, 1],
                              'B_c': [0, 0, 1, 0],
                              'B_nan': [0, 0, 0, 1]})
        cols = ['A_b', 'A_nan', 'B_c', 'B_nan']
        expected[cols] = expected[cols].astype(np.uint8)

        expected = expected[['C', 'A_b', 'A_nan', 'B_c', 'B_nan']]
        assert_frame_equal(result, expected)

        result = get_dummies(df, dummy_na=False, sparse=self.sparse,
                             drop_first=True)
        expected = expected[['C', 'A_b', 'B_c']]
        assert_frame_equal(result, expected)

    def test_int_int(self):
        data = Series([1, 2, 1])
        result = pd.get_dummies(data)
        expected = DataFrame([[1, 0], [0, 1], [1, 0]], columns=[1, 2],
                             dtype=np.uint8)
        tm.assert_frame_equal(result, expected)

        data = Series(pd.Categorical(['a', 'b', 'a']))
        result = pd.get_dummies(data)
        expected = DataFrame([[1, 0], [0, 1], [1, 0]],
                             columns=pd.Categorical(['a', 'b']),
                             dtype=np.uint8)
        tm.assert_frame_equal(result, expected)

    def test_int_df(self):
        data = DataFrame(
            {'A': [1, 2, 1],
             'B': pd.Categorical(['a', 'b', 'a']),
             'C': [1, 2, 1],
             'D': [1., 2., 1.]
             }
        )
        columns = ['C', 'D', 'A_1', 'A_2', 'B_a', 'B_b']
        expected = DataFrame([
            [1, 1., 1, 0, 1, 0],
            [2, 2., 0, 1, 0, 1],
            [1, 1., 1, 0, 1, 0]
        ], columns=columns)
        expected[columns[2:]] = expected[columns[2:]].astype(np.uint8)
        result = pd.get_dummies(data, columns=['A', 'B'])
        tm.assert_frame_equal(result, expected)

    def test_dataframe_dummies_preserve_categorical_dtype(self):
        # GH13854
        for ordered in [False, True]:
            cat = pd.Categorical(list("xy"), categories=list("xyz"),
                                 ordered=ordered)
            result = get_dummies(cat)

            data = np.array([[1, 0, 0], [0, 1, 0]], dtype=np.uint8)
            cols = pd.CategoricalIndex(cat.categories,
                                       categories=cat.categories,
                                       ordered=ordered)
            expected = DataFrame(data, columns=cols)

            tm.assert_frame_equal(result, expected)


class TestGetDummiesSparse(TestGetDummies):
    sparse = True


class TestMakeAxisDummies(tm.TestCase):

    def test_preserve_categorical_dtype(self):
        # GH13854
        for ordered in [False, True]:
            cidx = pd.CategoricalIndex(list("xyz"), ordered=ordered)
            midx = pd.MultiIndex(levels=[['a'], cidx],
                                 labels=[[0, 0], [0, 1]])
            df = DataFrame([[10, 11]], index=midx)

            expected = DataFrame([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                                 index=midx, columns=cidx)

            from pandas.core.reshape import make_axis_dummies
            result = make_axis_dummies(df)
            tm.assert_frame_equal(result, expected)

            result = make_axis_dummies(df, transform=lambda x: x)
            tm.assert_frame_equal(result, expected)


class TestLreshape(tm.TestCase):

    def test_pairs(self):
        data = {'birthdt': ['08jan2009', '20dec2008', '30dec2008', '21dec2008',
                            '11jan2009'],
                'birthwt': [1766, 3301, 1454, 3139, 4133],
                'id': [101, 102, 103, 104, 105],
                'sex': ['Male', 'Female', 'Female', 'Female', 'Female'],
                'visitdt1': ['11jan2009', '22dec2008', '04jan2009',
                             '29dec2008', '20jan2009'],
                'visitdt2':
                ['21jan2009', nan, '22jan2009', '31dec2008', '03feb2009'],
                'visitdt3': ['05feb2009', nan, nan, '02jan2009', '15feb2009'],
                'wt1': [1823, 3338, 1549, 3298, 4306],
                'wt2': [2011.0, nan, 1892.0, 3338.0, 4575.0],
                'wt3': [2293.0, nan, nan, 3377.0, 4805.0]}

        df = DataFrame(data)

        spec = {'visitdt': ['visitdt%d' % i for i in range(1, 4)],
                'wt': ['wt%d' % i for i in range(1, 4)]}
        result = lreshape(df, spec)

        exp_data = {'birthdt':
                    ['08jan2009', '20dec2008', '30dec2008', '21dec2008',
                     '11jan2009', '08jan2009', '30dec2008', '21dec2008',
                     '11jan2009', '08jan2009', '21dec2008', '11jan2009'],
                    'birthwt': [1766, 3301, 1454, 3139, 4133, 1766, 1454, 3139,
                                4133, 1766, 3139, 4133],
                    'id': [101, 102, 103, 104, 105, 101, 103, 104, 105, 101,
                           104, 105],
                    'sex': ['Male', 'Female', 'Female', 'Female', 'Female',
                            'Male', 'Female', 'Female', 'Female', 'Male',
                            'Female', 'Female'],
                    'visitdt': ['11jan2009', '22dec2008', '04jan2009',
                                '29dec2008', '20jan2009', '21jan2009',
                                '22jan2009', '31dec2008', '03feb2009',
                                '05feb2009', '02jan2009', '15feb2009'],
                    'wt': [1823.0, 3338.0, 1549.0, 3298.0, 4306.0, 2011.0,
                           1892.0, 3338.0, 4575.0, 2293.0, 3377.0, 4805.0]}
        exp = DataFrame(exp_data, columns=result.columns)
        tm.assert_frame_equal(result, exp)

        result = lreshape(df, spec, dropna=False)
        exp_data = {'birthdt':
                    ['08jan2009', '20dec2008', '30dec2008', '21dec2008',
                     '11jan2009', '08jan2009', '20dec2008', '30dec2008',
                     '21dec2008', '11jan2009', '08jan2009', '20dec2008',
                     '30dec2008', '21dec2008', '11jan2009'],
                    'birthwt': [1766, 3301, 1454, 3139, 4133, 1766, 3301, 1454,
                                3139, 4133, 1766, 3301, 1454, 3139, 4133],
                    'id': [101, 102, 103, 104, 105, 101, 102, 103, 104, 105,
                           101, 102, 103, 104, 105],
                    'sex': ['Male', 'Female', 'Female', 'Female', 'Female',
                            'Male', 'Female', 'Female', 'Female', 'Female',
                            'Male', 'Female', 'Female', 'Female', 'Female'],
                    'visitdt': ['11jan2009', '22dec2008', '04jan2009',
                                '29dec2008', '20jan2009', '21jan2009', nan,
                                '22jan2009', '31dec2008', '03feb2009',
                                '05feb2009', nan, nan, '02jan2009',
                                '15feb2009'],
                    'wt': [1823.0, 3338.0, 1549.0, 3298.0, 4306.0, 2011.0, nan,
                           1892.0, 3338.0, 4575.0, 2293.0, nan, nan, 3377.0,
                           4805.0]}
        exp = DataFrame(exp_data, columns=result.columns)
        tm.assert_frame_equal(result, exp)

        spec = {'visitdt': ['visitdt%d' % i for i in range(1, 3)],
                'wt': ['wt%d' % i for i in range(1, 4)]}
        self.assertRaises(ValueError, lreshape, df, spec)


class TestWideToLong(tm.TestCase):

    def test_simple(self):
        np.random.seed(123)
        x = np.random.randn(3)
        df = pd.DataFrame({"A1970": {0: "a",
                                     1: "b",
                                     2: "c"},
                           "A1980": {0: "d",
                                     1: "e",
                                     2: "f"},
                           "B1970": {0: 2.5,
                                     1: 1.2,
                                     2: .7},
                           "B1980": {0: 3.2,
                                     1: 1.3,
                                     2: .1},
                           "X": dict(zip(
                               range(3), x))})
        df["id"] = df.index
        exp_data = {"X": x.tolist() + x.tolist(),
                    "A": ['a', 'b', 'c', 'd', 'e', 'f'],
                    "B": [2.5, 1.2, 0.7, 3.2, 1.3, 0.1],
                    "year": [1970, 1970, 1970, 1980, 1980, 1980],
                    "id": [0, 1, 2, 0, 1, 2]}
        exp_frame = DataFrame(exp_data)
        exp_frame = exp_frame.set_index(['id', 'year'])[["X", "A", "B"]]
        long_frame = wide_to_long(df, ["A", "B"], i="id", j="year")
        tm.assert_frame_equal(long_frame, exp_frame)

    def test_stubs(self):
        # GH9204
        df = pd.DataFrame([[0, 1, 2, 3, 8], [4, 5, 6, 7, 9]])
        df.columns = ['id', 'inc1', 'inc2', 'edu1', 'edu2']
        stubs = ['inc', 'edu']

        # TODO: unused?
        df_long = pd.wide_to_long(df, stubs, i='id', j='age')  # noqa

        self.assertEqual(stubs, ['inc', 'edu'])


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
