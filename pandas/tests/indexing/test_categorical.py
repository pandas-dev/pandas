# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from pandas import Series, DataFrame
from pandas.util.testing import assert_series_equal, assert_frame_equal
from pandas.util import testing as tm


class TestCategoricalIndex(tm.TestCase):

    def setUp(self):

        self.df = DataFrame({'A': np.arange(6, dtype='int64'),
                             'B': Series(list('aabbca')).astype(
                                 'category', categories=list(
                                     'cab'))}).set_index('B')
        self.df2 = DataFrame({'A': np.arange(6, dtype='int64'),
                              'B': Series(list('aabbca')).astype(
                                  'category', categories=list(
                                      'cabe'))}).set_index('B')
        self.df3 = DataFrame({'A': np.arange(6, dtype='int64'),
                              'B': (Series([1, 1, 2, 1, 3, 2])
                                    .astype('category', categories=[3, 2, 1],
                                            ordered=True))}).set_index('B')
        self.df4 = DataFrame({'A': np.arange(6, dtype='int64'),
                              'B': (Series([1, 1, 2, 1, 3, 2])
                                    .astype('category', categories=[3, 2, 1],
                                            ordered=False))}).set_index('B')

    def test_loc_scalar(self):
        result = self.df.loc['a']
        expected = (DataFrame({'A': [0, 1, 5],
                               'B': (Series(list('aaa'))
                                     .astype('category',
                                             categories=list('cab')))})
                    .set_index('B'))
        assert_frame_equal(result, expected)

        df = self.df.copy()
        df.loc['a'] = 20
        expected = (DataFrame({'A': [20, 20, 2, 3, 4, 20],
                               'B': (Series(list('aabbca'))
                                     .astype('category',
                                             categories=list('cab')))})
                    .set_index('B'))
        assert_frame_equal(df, expected)

        # value not in the categories
        self.assertRaises(KeyError, lambda: df.loc['d'])

        def f():
            df.loc['d'] = 10

        self.assertRaises(TypeError, f)

        def f():
            df.loc['d', 'A'] = 10

        self.assertRaises(TypeError, f)

        def f():
            df.loc['d', 'C'] = 10

        self.assertRaises(TypeError, f)

    def test_loc_listlike(self):

        # list of labels
        result = self.df.loc[['c', 'a']]
        expected = self.df.iloc[[4, 0, 1, 5]]
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.loc[['a', 'b', 'e']]
        exp_index = pd.CategoricalIndex(
            list('aaabbe'), categories=list('cabe'), name='B')
        expected = DataFrame({'A': [0, 1, 5, 2, 3, np.nan]}, index=exp_index)
        assert_frame_equal(result, expected, check_index_type=True)

        # element in the categories but not in the values
        self.assertRaises(KeyError, lambda: self.df2.loc['e'])

        # assign is ok
        df = self.df2.copy()
        df.loc['e'] = 20
        result = df.loc[['a', 'b', 'e']]
        exp_index = pd.CategoricalIndex(
            list('aaabbe'), categories=list('cabe'), name='B')
        expected = DataFrame({'A': [0, 1, 5, 2, 3, 20]}, index=exp_index)
        assert_frame_equal(result, expected)

        df = self.df2.copy()
        result = df.loc[['a', 'b', 'e']]
        exp_index = pd.CategoricalIndex(
            list('aaabbe'), categories=list('cabe'), name='B')
        expected = DataFrame({'A': [0, 1, 5, 2, 3, np.nan]}, index=exp_index)
        assert_frame_equal(result, expected, check_index_type=True)

        # not all labels in the categories
        self.assertRaises(KeyError, lambda: self.df2.loc[['a', 'd']])

    def test_loc_listlike_dtypes(self):
        # GH 11586

        # unique categories and codes
        index = pd.CategoricalIndex(['a', 'b', 'c'])
        df = DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]}, index=index)

        # unique slice
        res = df.loc[['a', 'b']]
        exp_index = pd.CategoricalIndex(['a', 'b'],
                                        categories=index.categories)
        exp = DataFrame({'A': [1, 2], 'B': [4, 5]}, index=exp_index)
        tm.assert_frame_equal(res, exp, check_index_type=True)

        # duplicated slice
        res = df.loc[['a', 'a', 'b']]

        exp_index = pd.CategoricalIndex(['a', 'a', 'b'],
                                        categories=index.categories)
        exp = DataFrame({'A': [1, 1, 2], 'B': [4, 4, 5]}, index=exp_index)
        tm.assert_frame_equal(res, exp, check_index_type=True)

        with tm.assertRaisesRegexp(
                KeyError,
                'a list-indexer must only include values that are '
                'in the categories'):
            df.loc[['a', 'x']]

        # duplicated categories and codes
        index = pd.CategoricalIndex(['a', 'b', 'a'])
        df = DataFrame({'A': [1, 2, 3], 'B': [4, 5, 6]}, index=index)

        # unique slice
        res = df.loc[['a', 'b']]
        exp = DataFrame({'A': [1, 3, 2],
                         'B': [4, 6, 5]},
                        index=pd.CategoricalIndex(['a', 'a', 'b']))
        tm.assert_frame_equal(res, exp, check_index_type=True)

        # duplicated slice
        res = df.loc[['a', 'a', 'b']]
        exp = DataFrame(
            {'A': [1, 3, 1, 3, 2],
             'B': [4, 6, 4, 6, 5
                   ]}, index=pd.CategoricalIndex(['a', 'a', 'a', 'a', 'b']))
        tm.assert_frame_equal(res, exp, check_index_type=True)

        with tm.assertRaisesRegexp(
                KeyError,
                'a list-indexer must only include values '
                'that are in the categories'):
            df.loc[['a', 'x']]

        # contains unused category
        index = pd.CategoricalIndex(
            ['a', 'b', 'a', 'c'], categories=list('abcde'))
        df = DataFrame({'A': [1, 2, 3, 4], 'B': [5, 6, 7, 8]}, index=index)

        res = df.loc[['a', 'b']]
        exp = DataFrame({'A': [1, 3, 2],
                         'B': [5, 7, 6]}, index=pd.CategoricalIndex(
                             ['a', 'a', 'b'], categories=list('abcde')))
        tm.assert_frame_equal(res, exp, check_index_type=True)

        res = df.loc[['a', 'e']]
        exp = DataFrame({'A': [1, 3, np.nan], 'B': [5, 7, np.nan]},
                        index=pd.CategoricalIndex(['a', 'a', 'e'],
                                                  categories=list('abcde')))
        tm.assert_frame_equal(res, exp, check_index_type=True)

        # duplicated slice
        res = df.loc[['a', 'a', 'b']]
        exp = DataFrame({'A': [1, 3, 1, 3, 2], 'B': [5, 7, 5, 7, 6]},
                        index=pd.CategoricalIndex(['a', 'a', 'a', 'a', 'b'],
                                                  categories=list('abcde')))
        tm.assert_frame_equal(res, exp, check_index_type=True)

        with tm.assertRaisesRegexp(
                KeyError,
                'a list-indexer must only include values '
                'that are in the categories'):
            df.loc[['a', 'x']]

    def test_ix_categorical_index(self):
        # GH 12531
        df = pd.DataFrame(np.random.randn(3, 3),
                          index=list('ABC'), columns=list('XYZ'))
        cdf = df.copy()
        cdf.index = pd.CategoricalIndex(df.index)
        cdf.columns = pd.CategoricalIndex(df.columns)

        expect = pd.Series(df.ix['A', :], index=cdf.columns, name='A')
        assert_series_equal(cdf.ix['A', :], expect)

        expect = pd.Series(df.ix[:, 'X'], index=cdf.index, name='X')
        assert_series_equal(cdf.ix[:, 'X'], expect)

        exp_index = pd.CategoricalIndex(list('AB'), categories=['A', 'B', 'C'])
        expect = pd.DataFrame(df.ix[['A', 'B'], :], columns=cdf.columns,
                              index=exp_index)
        assert_frame_equal(cdf.ix[['A', 'B'], :], expect)

        exp_columns = pd.CategoricalIndex(list('XY'),
                                          categories=['X', 'Y', 'Z'])
        expect = pd.DataFrame(df.ix[:, ['X', 'Y']], index=cdf.index,
                              columns=exp_columns)
        assert_frame_equal(cdf.ix[:, ['X', 'Y']], expect)

        # non-unique
        df = pd.DataFrame(np.random.randn(3, 3),
                          index=list('ABA'), columns=list('XYX'))
        cdf = df.copy()
        cdf.index = pd.CategoricalIndex(df.index)
        cdf.columns = pd.CategoricalIndex(df.columns)

        exp_index = pd.CategoricalIndex(list('AA'), categories=['A', 'B'])
        expect = pd.DataFrame(df.ix['A', :], columns=cdf.columns,
                              index=exp_index)
        assert_frame_equal(cdf.ix['A', :], expect)

        exp_columns = pd.CategoricalIndex(list('XX'), categories=['X', 'Y'])
        expect = pd.DataFrame(df.ix[:, 'X'], index=cdf.index,
                              columns=exp_columns)
        assert_frame_equal(cdf.ix[:, 'X'], expect)

        expect = pd.DataFrame(df.ix[['A', 'B'], :], columns=cdf.columns,
                              index=pd.CategoricalIndex(list('AAB')))
        assert_frame_equal(cdf.ix[['A', 'B'], :], expect)

        expect = pd.DataFrame(df.ix[:, ['X', 'Y']], index=cdf.index,
                              columns=pd.CategoricalIndex(list('XXY')))
        assert_frame_equal(cdf.ix[:, ['X', 'Y']], expect)

    def test_read_only_source(self):
        # GH 10043
        rw_array = np.eye(10)
        rw_df = DataFrame(rw_array)

        ro_array = np.eye(10)
        ro_array.setflags(write=False)
        ro_df = DataFrame(ro_array)

        assert_frame_equal(rw_df.iloc[[1, 2, 3]], ro_df.iloc[[1, 2, 3]])
        assert_frame_equal(rw_df.iloc[[1]], ro_df.iloc[[1]])
        assert_series_equal(rw_df.iloc[1], ro_df.iloc[1])
        assert_frame_equal(rw_df.iloc[1:3], ro_df.iloc[1:3])

        assert_frame_equal(rw_df.loc[[1, 2, 3]], ro_df.loc[[1, 2, 3]])
        assert_frame_equal(rw_df.loc[[1]], ro_df.loc[[1]])
        assert_series_equal(rw_df.loc[1], ro_df.loc[1])
        assert_frame_equal(rw_df.loc[1:3], ro_df.loc[1:3])

    def test_reindexing(self):

        # reindexing
        # convert to a regular index
        result = self.df2.reindex(['a', 'b', 'e'])
        expected = DataFrame({'A': [0, 1, 5, 2, 3, np.nan],
                              'B': Series(list('aaabbe'))}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(['a', 'b'])
        expected = DataFrame({'A': [0, 1, 5, 2, 3],
                              'B': Series(list('aaabb'))}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(['e'])
        expected = DataFrame({'A': [np.nan],
                              'B': Series(['e'])}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(['d'])
        expected = DataFrame({'A': [np.nan],
                              'B': Series(['d'])}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        # since we are actually reindexing with a Categorical
        # then return a Categorical
        cats = list('cabe')

        result = self.df2.reindex(pd.Categorical(['a', 'd'], categories=cats))
        expected = DataFrame({'A': [0, 1, 5, np.nan],
                              'B': Series(list('aaad')).astype(
                                  'category', categories=cats)}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(pd.Categorical(['a'], categories=cats))
        expected = DataFrame({'A': [0, 1, 5],
                              'B': Series(list('aaa')).astype(
                                  'category', categories=cats)}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(['a', 'b', 'e'])
        expected = DataFrame({'A': [0, 1, 5, 2, 3, np.nan],
                              'B': Series(list('aaabbe'))}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(['a', 'b'])
        expected = DataFrame({'A': [0, 1, 5, 2, 3],
                              'B': Series(list('aaabb'))}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(['e'])
        expected = DataFrame({'A': [np.nan],
                              'B': Series(['e'])}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        # give back the type of categorical that we received
        result = self.df2.reindex(pd.Categorical(
            ['a', 'd'], categories=cats, ordered=True))
        expected = DataFrame(
            {'A': [0, 1, 5, np.nan],
             'B': Series(list('aaad')).astype('category', categories=cats,
                                              ordered=True)}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        result = self.df2.reindex(pd.Categorical(
            ['a', 'd'], categories=['a', 'd']))
        expected = DataFrame({'A': [0, 1, 5, np.nan],
                              'B': Series(list('aaad')).astype(
                                  'category', categories=['a', 'd'
                                                          ])}).set_index('B')
        assert_frame_equal(result, expected, check_index_type=True)

        # passed duplicate indexers are not allowed
        self.assertRaises(ValueError, lambda: self.df2.reindex(['a', 'a']))

        # args NotImplemented ATM
        self.assertRaises(NotImplementedError,
                          lambda: self.df2.reindex(['a'], method='ffill'))
        self.assertRaises(NotImplementedError,
                          lambda: self.df2.reindex(['a'], level=1))
        self.assertRaises(NotImplementedError,
                          lambda: self.df2.reindex(['a'], limit=2))

    def test_loc_slice(self):
        # slicing
        # not implemented ATM
        # GH9748

        self.assertRaises(TypeError, lambda: self.df.loc[1:5])

        # result = df.loc[1:5]
        # expected = df.iloc[[1,2,3,4]]
        # assert_frame_equal(result, expected)

    def test_boolean_selection(self):

        df3 = self.df3
        df4 = self.df4

        result = df3[df3.index == 'a']
        expected = df3.iloc[[]]
        assert_frame_equal(result, expected)

        result = df4[df4.index == 'a']
        expected = df4.iloc[[]]
        assert_frame_equal(result, expected)

        result = df3[df3.index == 1]
        expected = df3.iloc[[0, 1, 3]]
        assert_frame_equal(result, expected)

        result = df4[df4.index == 1]
        expected = df4.iloc[[0, 1, 3]]
        assert_frame_equal(result, expected)

        # since we have an ordered categorical

        # CategoricalIndex([1, 1, 2, 1, 3, 2],
        #         categories=[3, 2, 1],
        #         ordered=True,
        #         name=u'B')
        result = df3[df3.index < 2]
        expected = df3.iloc[[4]]
        assert_frame_equal(result, expected)

        result = df3[df3.index > 1]
        expected = df3.iloc[[]]
        assert_frame_equal(result, expected)

        # unordered
        # cannot be compared

        # CategoricalIndex([1, 1, 2, 1, 3, 2],
        #         categories=[3, 2, 1],
        #         ordered=False,
        #         name=u'B')
        self.assertRaises(TypeError, lambda: df4[df4.index < 2])
        self.assertRaises(TypeError, lambda: df4[df4.index > 1])

    def test_indexing_with_category(self):

        # https://github.com/pandas-dev/pandas/issues/12564
        # consistent result if comparing as Dataframe

        cat = DataFrame({'A': ['foo', 'bar', 'baz']})
        exp = DataFrame({'A': [True, False, False]})

        res = (cat[['A']] == 'foo')
        tm.assert_frame_equal(res, exp)

        cat['A'] = cat['A'].astype('category')

        res = (cat[['A']] == 'foo')
        tm.assert_frame_equal(res, exp)
