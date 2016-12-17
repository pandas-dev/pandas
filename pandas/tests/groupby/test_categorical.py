# -*- coding: utf-8 -*-
from __future__ import print_function
import nose
from numpy import nan


from pandas.core.index import Index, MultiIndex, CategoricalIndex
from pandas.core.api import DataFrame, Categorical

from pandas.core.series import Series

from pandas.util.testing import (assert_frame_equal, assert_series_equal
                                 )

from pandas.compat import (lmap)

from pandas import compat

import pandas.core.common as com
import numpy as np

import pandas.util.testing as tm
import pandas as pd


class TestGroupByCategorical(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.ts = tm.makeTimeSeries()

        self.seriesd = tm.getSeriesData()
        self.tsd = tm.getTimeSeriesData()
        self.frame = DataFrame(self.seriesd)
        self.tsframe = DataFrame(self.tsd)

        self.df = DataFrame(
            {'A': ['foo', 'bar', 'foo', 'bar', 'foo', 'bar', 'foo', 'foo'],
             'B': ['one', 'one', 'two', 'three', 'two', 'two', 'one', 'three'],
             'C': np.random.randn(8),
             'D': np.random.randn(8)})

        self.df_mixed_floats = DataFrame(
            {'A': ['foo', 'bar', 'foo', 'bar', 'foo', 'bar', 'foo', 'foo'],
             'B': ['one', 'one', 'two', 'three', 'two', 'two', 'one', 'three'],
             'C': np.random.randn(8),
             'D': np.array(
                 np.random.randn(8), dtype='float32')})

        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'], ['one', 'two',
                                                                  'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        self.mframe = DataFrame(np.random.randn(10, 3), index=index,
                                columns=['A', 'B', 'C'])

        self.three_group = DataFrame(
            {'A': ['foo', 'foo', 'foo', 'foo', 'bar', 'bar', 'bar', 'bar',
                   'foo', 'foo', 'foo'],
             'B': ['one', 'one', 'one', 'two', 'one', 'one', 'one', 'two',
                   'two', 'two', 'one'],
             'C': ['dull', 'dull', 'shiny', 'dull', 'dull', 'shiny', 'shiny',
                   'dull', 'shiny', 'shiny', 'shiny'],
             'D': np.random.randn(11),
             'E': np.random.randn(11),
             'F': np.random.randn(11)})

    def test_apply_use_categorical_name(self):
        from pandas import qcut
        cats = qcut(self.df.C, 4)

        def get_stats(group):
            return {'min': group.min(),
                    'max': group.max(),
                    'count': group.count(),
                    'mean': group.mean()}

        result = self.df.groupby(cats).D.apply(get_stats)
        self.assertEqual(result.index.names[0], 'C')

    def test_apply_categorical_data(self):
        # GH 10138
        for ordered in [True, False]:
            dense = Categorical(list('abc'), ordered=ordered)
            # 'b' is in the categories but not in the list
            missing = Categorical(
                list('aaa'), categories=['a', 'b'], ordered=ordered)
            values = np.arange(len(dense))
            df = DataFrame({'missing': missing,
                            'dense': dense,
                            'values': values})
            grouped = df.groupby(['missing', 'dense'])

            # missing category 'b' should still exist in the output index
            idx = MultiIndex.from_product(
                [Categorical(['a', 'b'], ordered=ordered),
                 Categorical(['a', 'b', 'c'], ordered=ordered)],
                names=['missing', 'dense'])
            expected = DataFrame([0, 1, 2, np.nan, np.nan, np.nan],
                                 index=idx,
                                 columns=['values'])

            assert_frame_equal(grouped.apply(lambda x: np.mean(x)), expected)
            assert_frame_equal(grouped.mean(), expected)
            assert_frame_equal(grouped.agg(np.mean), expected)

            # but for transform we should still get back the original index
            idx = MultiIndex.from_product([['a'], ['a', 'b', 'c']],
                                          names=['missing', 'dense'])
            expected = Series(1, index=idx)
            assert_series_equal(grouped.apply(lambda x: 1), expected)

    def test_groupby_categorical(self):
        levels = ['foo', 'bar', 'baz', 'qux']
        codes = np.random.randint(0, 4, size=100)

        cats = Categorical.from_codes(codes, levels, ordered=True)

        data = DataFrame(np.random.randn(100, 4))

        result = data.groupby(cats).mean()

        expected = data.groupby(np.asarray(cats)).mean()
        exp_idx = CategoricalIndex(levels, categories=cats.categories,
                                   ordered=True)
        expected = expected.reindex(exp_idx)

        assert_frame_equal(result, expected)

        grouped = data.groupby(cats)
        desc_result = grouped.describe()

        idx = cats.codes.argsort()
        ord_labels = np.asarray(cats).take(idx)
        ord_data = data.take(idx)

        exp_cats = Categorical(ord_labels, ordered=True,
                               categories=['foo', 'bar', 'baz', 'qux'])
        expected = ord_data.groupby(exp_cats, sort=False).describe()
        expected.index.names = [None, None]
        assert_frame_equal(desc_result, expected)

        # GH 10460
        expc = Categorical.from_codes(np.arange(4).repeat(8),
                                      levels, ordered=True)
        exp = CategoricalIndex(expc)
        self.assert_index_equal(desc_result.index.get_level_values(0), exp)
        exp = Index(['count', 'mean', 'std', 'min', '25%', '50%',
                     '75%', 'max'] * 4)
        self.assert_index_equal(desc_result.index.get_level_values(1), exp)

    def test_groupby_datetime_categorical(self):
        # GH9049: ensure backward compatibility
        levels = pd.date_range('2014-01-01', periods=4)
        codes = np.random.randint(0, 4, size=100)

        cats = Categorical.from_codes(codes, levels, ordered=True)

        data = DataFrame(np.random.randn(100, 4))
        result = data.groupby(cats).mean()

        expected = data.groupby(np.asarray(cats)).mean()
        expected = expected.reindex(levels)
        expected.index = CategoricalIndex(expected.index,
                                          categories=expected.index,
                                          ordered=True)

        assert_frame_equal(result, expected)

        grouped = data.groupby(cats)
        desc_result = grouped.describe()

        idx = cats.codes.argsort()
        ord_labels = cats.take_nd(idx)
        ord_data = data.take(idx)
        expected = ord_data.groupby(ord_labels).describe()
        expected.index.names = [None, None]
        assert_frame_equal(desc_result, expected)
        tm.assert_index_equal(desc_result.index, expected.index)
        tm.assert_index_equal(
            desc_result.index.get_level_values(0),
            expected.index.get_level_values(0))

        # GH 10460
        expc = Categorical.from_codes(
            np.arange(4).repeat(8), levels, ordered=True)
        exp = CategoricalIndex(expc)
        self.assert_index_equal(desc_result.index.get_level_values(0), exp)
        exp = Index(['count', 'mean', 'std', 'min', '25%', '50%',
                     '75%', 'max'] * 4)
        self.assert_index_equal(desc_result.index.get_level_values(1), exp)

    def test_groupby_categorical_index(self):

        levels = ['foo', 'bar', 'baz', 'qux']
        codes = np.random.randint(0, 4, size=20)
        cats = Categorical.from_codes(codes, levels, ordered=True)
        df = DataFrame(
            np.repeat(
                np.arange(20), 4).reshape(-1, 4), columns=list('abcd'))
        df['cats'] = cats

        # with a cat index
        result = df.set_index('cats').groupby(level=0).sum()
        expected = df[list('abcd')].groupby(cats.codes).sum()
        expected.index = CategoricalIndex(
            Categorical.from_codes(
                [0, 1, 2, 3], levels, ordered=True), name='cats')
        assert_frame_equal(result, expected)

        # with a cat column, should produce a cat index
        result = df.groupby('cats').sum()
        expected = df[list('abcd')].groupby(cats.codes).sum()
        expected.index = CategoricalIndex(
            Categorical.from_codes(
                [0, 1, 2, 3], levels, ordered=True), name='cats')
        assert_frame_equal(result, expected)

    def test_groupby_describe_categorical_columns(self):
        # GH 11558
        cats = pd.CategoricalIndex(['qux', 'foo', 'baz', 'bar'],
                                   categories=['foo', 'bar', 'baz', 'qux'],
                                   ordered=True)
        df = DataFrame(np.random.randn(20, 4), columns=cats)
        result = df.groupby([1, 2, 3, 4] * 5).describe()

        tm.assert_index_equal(result.columns, cats)
        tm.assert_categorical_equal(result.columns.values, cats.values)

    def test_groupby_unstack_categorical(self):
        # GH11558 (example is taken from the original issue)
        df = pd.DataFrame({'a': range(10),
                           'medium': ['A', 'B'] * 5,
                           'artist': list('XYXXY') * 2})
        df['medium'] = df['medium'].astype('category')

        gcat = df.groupby(['artist', 'medium'])['a'].count().unstack()
        result = gcat.describe()

        exp_columns = pd.CategoricalIndex(['A', 'B'], ordered=False,
                                          name='medium')
        tm.assert_index_equal(result.columns, exp_columns)
        tm.assert_categorical_equal(result.columns.values, exp_columns.values)

        result = gcat['A'] + gcat['B']
        expected = pd.Series([6, 4], index=pd.Index(['X', 'Y'], name='artist'))
        tm.assert_series_equal(result, expected)

    def test_groupby_categorical_unequal_len(self):
        # GH3011
        series = Series([np.nan, np.nan, 1, 1, 2, 2, 3, 3, 4, 4])
        # The raises only happens with categorical, not with series of types
        # category
        bins = pd.cut(series.dropna().values, 4)

        # len(bins) != len(series) here
        self.assertRaises(ValueError, lambda: series.groupby(bins).mean())

    def test_groupby_categorical_two_columns(self):

        # https://github.com/pandas-dev/pandas/issues/8138
        d = {'cat':
             pd.Categorical(["a", "b", "a", "b"], categories=["a", "b", "c"],
                            ordered=True),
             'ints': [1, 1, 2, 2],
             'val': [10, 20, 30, 40]}
        test = pd.DataFrame(d)

        # Grouping on a single column
        groups_single_key = test.groupby("cat")
        res = groups_single_key.agg('mean')

        exp_index = pd.CategoricalIndex(["a", "b", "c"], name="cat",
                                        ordered=True)
        exp = DataFrame({"ints": [1.5, 1.5, np.nan], "val": [20, 30, np.nan]},
                        index=exp_index)
        tm.assert_frame_equal(res, exp)

        # Grouping on two columns
        groups_double_key = test.groupby(["cat", "ints"])
        res = groups_double_key.agg('mean')
        exp = DataFrame({"val": [10, 30, 20, 40, np.nan, np.nan],
                         "cat": pd.Categorical(["a", "a", "b", "b", "c", "c"],
                                               ordered=True),
                         "ints": [1, 2, 1, 2, 1, 2]}).set_index(["cat", "ints"
                                                                 ])
        tm.assert_frame_equal(res, exp)

        # GH 10132
        for key in [('a', 1), ('b', 2), ('b', 1), ('a', 2)]:
            c, i = key
            result = groups_double_key.get_group(key)
            expected = test[(test.cat == c) & (test.ints == i)]
            assert_frame_equal(result, expected)

        d = {'C1': [3, 3, 4, 5], 'C2': [1, 2, 3, 4], 'C3': [10, 100, 200, 34]}
        test = pd.DataFrame(d)
        values = pd.cut(test['C1'], [1, 2, 3, 6])
        values.name = "cat"
        groups_double_key = test.groupby([values, 'C2'])

        res = groups_double_key.agg('mean')
        nan = np.nan
        idx = MultiIndex.from_product(
            [Categorical(["(1, 2]", "(2, 3]", "(3, 6]"], ordered=True),
             [1, 2, 3, 4]],
            names=["cat", "C2"])
        exp = DataFrame({"C1": [nan, nan, nan, nan, 3, 3,
                                nan, nan, nan, nan, 4, 5],
                         "C3": [nan, nan, nan, nan, 10, 100,
                                nan, nan, nan, nan, 200, 34]}, index=idx)
        tm.assert_frame_equal(res, exp)

    def test_groupby_multi_categorical_as_index(self):
        # GH13204
        df = DataFrame({'cat': Categorical([1, 2, 2], [1, 2, 3]),
                        'A': [10, 11, 11],
                        'B': [101, 102, 103]})
        result = df.groupby(['cat', 'A'], as_index=False).sum()
        expected = DataFrame({'cat': Categorical([1, 1, 2, 2, 3, 3]),
                              'A': [10, 11, 10, 11, 10, 11],
                              'B': [101.0, nan, nan, 205.0, nan, nan]},
                             columns=['cat', 'A', 'B'])
        tm.assert_frame_equal(result, expected)

        # function grouper
        f = lambda r: df.loc[r, 'A']
        result = df.groupby(['cat', f], as_index=False).sum()
        expected = DataFrame({'cat': Categorical([1, 1, 2, 2, 3, 3]),
                              'A': [10.0, nan, nan, 22.0, nan, nan],
                              'B': [101.0, nan, nan, 205.0, nan, nan]},
                             columns=['cat', 'A', 'B'])
        tm.assert_frame_equal(result, expected)

        # another not in-axis grouper (conflicting names in index)
        s = Series(['a', 'b', 'b'], name='cat')
        result = df.groupby(['cat', s], as_index=False).sum()
        expected = DataFrame({'cat': Categorical([1, 1, 2, 2, 3, 3]),
                              'A': [10.0, nan, nan, 22.0, nan, nan],
                              'B': [101.0, nan, nan, 205.0, nan, nan]},
                             columns=['cat', 'A', 'B'])
        tm.assert_frame_equal(result, expected)

        # is original index dropped?
        expected = DataFrame({'cat': Categorical([1, 1, 2, 2, 3, 3]),
                              'A': [10, 11, 10, 11, 10, 11],
                              'B': [101.0, nan, nan, 205.0, nan, nan]},
                             columns=['cat', 'A', 'B'])

        group_columns = ['cat', 'A']

        for name in [None, 'X', 'B', 'cat']:
            df.index = Index(list("abc"), name=name)

            if name in group_columns and name in df.index.names:
                with tm.assert_produces_warning(FutureWarning,
                                                check_stacklevel=False):
                    result = df.groupby(group_columns, as_index=False).sum()

            else:
                result = df.groupby(group_columns, as_index=False).sum()

            tm.assert_frame_equal(result, expected, check_index_type=True)

    def test_groupby_preserve_categorical_dtype(self):
        # GH13743, GH13854
        df = DataFrame({'A': [1, 2, 1, 1, 2],
                        'B': [10, 16, 22, 28, 34],
                        'C1': Categorical(list("abaab"),
                                          categories=list("bac"),
                                          ordered=False),
                        'C2': Categorical(list("abaab"),
                                          categories=list("bac"),
                                          ordered=True)})
        # single grouper
        exp_full = DataFrame({'A': [2.0, 1.0, np.nan],
                              'B': [25.0, 20.0, np.nan],
                              'C1': Categorical(list("bac"),
                                                categories=list("bac"),
                                                ordered=False),
                              'C2': Categorical(list("bac"),
                                                categories=list("bac"),
                                                ordered=True)})
        for col in ['C1', 'C2']:
            result1 = df.groupby(by=col, as_index=False).mean()
            result2 = df.groupby(by=col, as_index=True).mean().reset_index()
            expected = exp_full.reindex(columns=result1.columns)
            tm.assert_frame_equal(result1, expected)
            tm.assert_frame_equal(result2, expected)

        # multiple grouper
        exp_full = DataFrame({'A': [1, 1, 1, 2, 2, 2],
                              'B': [np.nan, 20.0, np.nan, 25.0, np.nan,
                                    np.nan],
                              'C1': Categorical(list("bacbac"),
                                                categories=list("bac"),
                                                ordered=False),
                              'C2': Categorical(list("bacbac"),
                                                categories=list("bac"),
                                                ordered=True)})
        for cols in [['A', 'C1'], ['A', 'C2']]:
            result1 = df.groupby(by=cols, as_index=False).mean()
            result2 = df.groupby(by=cols, as_index=True).mean().reset_index()
            expected = exp_full.reindex(columns=result1.columns)
            tm.assert_frame_equal(result1, expected)
            tm.assert_frame_equal(result2, expected)

    def test_groupby_categorical_no_compress(self):
        data = Series(np.random.randn(9))

        codes = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
        cats = Categorical.from_codes(codes, [0, 1, 2], ordered=True)

        result = data.groupby(cats).mean()
        exp = data.groupby(codes).mean()

        exp.index = CategoricalIndex(exp.index, categories=cats.categories,
                                     ordered=cats.ordered)
        assert_series_equal(result, exp)

        codes = np.array([0, 0, 0, 1, 1, 1, 3, 3, 3])
        cats = Categorical.from_codes(codes, [0, 1, 2, 3], ordered=True)

        result = data.groupby(cats).mean()
        exp = data.groupby(codes).mean().reindex(cats.categories)
        exp.index = CategoricalIndex(exp.index, categories=cats.categories,
                                     ordered=cats.ordered)
        assert_series_equal(result, exp)

        cats = Categorical(["a", "a", "a", "b", "b", "b", "c", "c", "c"],
                           categories=["a", "b", "c", "d"], ordered=True)
        data = DataFrame({"a": [1, 1, 1, 2, 2, 2, 3, 4, 5], "b": cats})

        result = data.groupby("b").mean()
        result = result["a"].values
        exp = np.array([1, 2, 4, np.nan])
        self.assert_numpy_array_equal(result, exp)


def assert_fp_equal(a, b):
    assert (np.abs(a - b) < 1e-12).all()


def _check_groupby(df, result, keys, field, f=lambda x: x.sum()):
    tups = lmap(tuple, df[keys].values)
    tups = com._asarray_tuplesafe(tups)
    expected = f(df.groupby(tups)[field])
    for k, v in compat.iteritems(expected):
        assert (result[k] == v)


def test_decons():
    from pandas.core.groupby import decons_group_index, get_group_index

    def testit(label_list, shape):
        group_index = get_group_index(label_list, shape, sort=True, xnull=True)
        label_list2 = decons_group_index(group_index, shape)

        for a, b in zip(label_list, label_list2):
            assert (np.array_equal(a, b))

    shape = (4, 5, 6)
    label_list = [np.tile([0, 1, 2, 3, 0, 1, 2, 3], 100), np.tile(
        [0, 2, 4, 3, 0, 1, 2, 3], 100), np.tile(
            [5, 1, 0, 2, 3, 0, 5, 4], 100)]
    testit(label_list, shape)

    shape = (10000, 10000)
    label_list = [np.tile(np.arange(10000), 5), np.tile(np.arange(10000), 5)]
    testit(label_list, shape)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure', '-s'
                         ], exit=False)
