# -*- coding: utf-8 -*-

import pytest
from distutils.version import LooseVersion

import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas import (Categorical, Index, Series, DataFrame, CategoricalIndex)
from pandas.core.dtypes.dtypes import CategoricalDtype


class TestCategoricalGeneric(object):

    def setup_method(self, method):
        self.factor = Categorical(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c'],
                                  ordered=True)

    def test_categories_none(self):
        factor = Categorical(['a', 'b', 'b', 'a',
                              'a', 'c', 'c', 'c'], ordered=True)
        tm.assert_categorical_equal(factor, self.factor)

    def test_describe(self):
        # string type
        desc = self.factor.describe()
        assert self.factor.ordered
        exp_index = CategoricalIndex(['a', 'b', 'c'], name='categories',
                                     ordered=self.factor.ordered)
        expected = DataFrame({'counts': [3, 2, 3],
                              'freqs': [3 / 8., 2 / 8., 3 / 8.]},
                             index=exp_index)
        tm.assert_frame_equal(desc, expected)

        # check unused categories
        cat = self.factor.copy()
        cat.set_categories(["a", "b", "c", "d"], inplace=True)
        desc = cat.describe()

        exp_index = CategoricalIndex(
            list('abcd'), ordered=self.factor.ordered, name='categories')
        expected = DataFrame({'counts': [3, 2, 3, 0],
                              'freqs': [3 / 8., 2 / 8., 3 / 8., 0]},
                             index=exp_index)
        tm.assert_frame_equal(desc, expected)

        # check an integer one
        cat = Categorical([1, 2, 3, 1, 2, 3, 3, 2, 1, 1, 1])
        desc = cat.describe()
        exp_index = CategoricalIndex([1, 2, 3], ordered=cat.ordered,
                                     name='categories')
        expected = DataFrame({'counts': [5, 3, 3],
                              'freqs': [5 / 11., 3 / 11., 3 / 11.]},
                             index=exp_index)
        tm.assert_frame_equal(desc, expected)

        # https://github.com/pandas-dev/pandas/issues/3678
        # describe should work with NaN
        cat = Categorical([np.nan, 1, 2, 2])
        desc = cat.describe()
        expected = DataFrame({'counts': [1, 2, 1],
                              'freqs': [1 / 4., 2 / 4., 1 / 4.]},
                             index=CategoricalIndex([1, 2, np.nan],
                                                    categories=[1, 2],
                                                    name='categories'))
        tm.assert_frame_equal(desc, expected)

    def test_getitem(self):
        assert self.factor[0] == 'a'
        assert self.factor[-1] == 'c'

        subf = self.factor[[0, 1, 2]]
        tm.assert_numpy_array_equal(subf._codes,
                                    np.array([0, 1, 1], dtype=np.int8))

        subf = self.factor[np.asarray(self.factor) == 'c']
        tm.assert_numpy_array_equal(subf._codes,
                                    np.array([2, 2, 2], dtype=np.int8))

    def test_setitem(self):

        # int/positional
        c = self.factor.copy()
        c[0] = 'b'
        assert c[0] == 'b'
        c[-1] = 'a'
        assert c[-1] == 'a'

        # boolean
        c = self.factor.copy()
        indexer = np.zeros(len(c), dtype='bool')
        indexer[0] = True
        indexer[-1] = True
        c[indexer] = 'c'
        expected = Categorical(['c', 'b', 'b', 'a', 'a', 'c', 'c', 'c'],
                               ordered=True)

        tm.assert_categorical_equal(c, expected)

    def test_set_categories_inplace(self):
        cat = self.factor.copy()
        cat.set_categories(['a', 'b', 'c', 'd'], inplace=True)
        tm.assert_index_equal(cat.categories, Index(['a', 'b', 'c', 'd']))

    def test_comparisons(self):

        result = self.factor[self.factor == 'a']
        expected = self.factor[np.asarray(self.factor) == 'a']
        tm.assert_categorical_equal(result, expected)

        result = self.factor[self.factor != 'a']
        expected = self.factor[np.asarray(self.factor) != 'a']
        tm.assert_categorical_equal(result, expected)

        result = self.factor[self.factor < 'c']
        expected = self.factor[np.asarray(self.factor) < 'c']
        tm.assert_categorical_equal(result, expected)

        result = self.factor[self.factor > 'a']
        expected = self.factor[np.asarray(self.factor) > 'a']
        tm.assert_categorical_equal(result, expected)

        result = self.factor[self.factor >= 'b']
        expected = self.factor[np.asarray(self.factor) >= 'b']
        tm.assert_categorical_equal(result, expected)

        result = self.factor[self.factor <= 'b']
        expected = self.factor[np.asarray(self.factor) <= 'b']
        tm.assert_categorical_equal(result, expected)

        n = len(self.factor)

        other = self.factor[np.random.permutation(n)]
        result = self.factor == other
        expected = np.asarray(self.factor) == np.asarray(other)
        tm.assert_numpy_array_equal(result, expected)

        result = self.factor == 'd'
        expected = np.repeat(False, len(self.factor))
        tm.assert_numpy_array_equal(result, expected)

        # comparisons with categoricals
        cat_rev = Categorical(
            ["a", "b", "c"], categories=["c", "b", "a"], ordered=True)
        cat_rev_base = Categorical(
            ["b", "b", "b"], categories=["c", "b", "a"], ordered=True)
        cat = Categorical(["a", "b", "c"], ordered=True)
        cat_base = Categorical(
            ["b", "b", "b"], categories=cat.categories, ordered=True)

        # comparisons need to take categories ordering into account
        res_rev = cat_rev > cat_rev_base
        exp_rev = np.array([True, False, False])
        tm.assert_numpy_array_equal(res_rev, exp_rev)

        res_rev = cat_rev < cat_rev_base
        exp_rev = np.array([False, False, True])
        tm.assert_numpy_array_equal(res_rev, exp_rev)

        res = cat > cat_base
        exp = np.array([False, False, True])
        tm.assert_numpy_array_equal(res, exp)

        # Only categories with same categories can be compared
        def f():
            cat > cat_rev

        pytest.raises(TypeError, f)

        cat_rev_base2 = Categorical(
            ["b", "b", "b"], categories=["c", "b", "a", "d"])

        def f():
            cat_rev > cat_rev_base2

        pytest.raises(TypeError, f)

        # Only categories with same ordering information can be compared
        cat_unorderd = cat.set_ordered(False)
        assert not (cat > cat).any()

        def f():
            cat > cat_unorderd

        pytest.raises(TypeError, f)

        # comparison (in both directions) with Series will raise
        s = Series(["b", "b", "b"])
        pytest.raises(TypeError, lambda: cat > s)
        pytest.raises(TypeError, lambda: cat_rev > s)
        pytest.raises(TypeError, lambda: s < cat)
        pytest.raises(TypeError, lambda: s < cat_rev)

        # comparison with numpy.array will raise in both direction, but only on
        # newer numpy versions
        a = np.array(["b", "b", "b"])
        pytest.raises(TypeError, lambda: cat > a)
        pytest.raises(TypeError, lambda: cat_rev > a)

        # The following work via '__array_priority__ = 1000'
        # works only on numpy >= 1.7.1
        if LooseVersion(np.__version__) > "1.7.1":
            pytest.raises(TypeError, lambda: a < cat)
            pytest.raises(TypeError, lambda: a < cat_rev)

        # Make sure that unequal comparison take the categories order in
        # account
        cat_rev = Categorical(
            list("abc"), categories=list("cba"), ordered=True)
        exp = np.array([True, False, False])
        res = cat_rev > "b"
        tm.assert_numpy_array_equal(res, exp)

    def test_print(self):
        expected = ["[a, b, b, a, a, c, c, c]",
                    "Categories (3, object): [a < b < c]"]
        expected = "\n".join(expected)
        actual = repr(self.factor)
        assert actual == expected


class TestCategoricalGenericBlock(object):

    def setup_method(self, method):
        self.factor = Categorical(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c'])

        df = DataFrame({'value': np.random.randint(0, 10000, 100)})
        labels = ["{0} - {1}".format(i, i + 499) for i in range(0, 10000, 500)]
        cat_labels = Categorical(labels, labels)

        df = df.sort_values(by=['value'], ascending=True)
        df['value_group'] = pd.cut(df.value, range(0, 10500, 500),
                                   right=False, labels=cat_labels)
        self.cat = df

    def test_basic(self):

        # test basic creation / coercion of categoricals
        s = Series(self.factor, name='A')
        assert s.dtype == 'category'
        assert len(s) == len(self.factor)
        str(s.values)
        str(s)

        # in a frame
        df = DataFrame({'A': self.factor})
        result = df['A']
        tm.assert_series_equal(result, s)
        result = df.iloc[:, 0]
        tm.assert_series_equal(result, s)
        assert len(df) == len(self.factor)
        str(df.values)
        str(df)

        df = DataFrame({'A': s})
        result = df['A']
        tm.assert_series_equal(result, s)
        assert len(df) == len(self.factor)
        str(df.values)
        str(df)

        # multiples
        df = DataFrame({'A': s, 'B': s, 'C': 1})
        result1 = df['A']
        result2 = df['B']
        tm.assert_series_equal(result1, s)
        tm.assert_series_equal(result2, s, check_names=False)
        assert result2.name == 'B'
        assert len(df) == len(self.factor)
        str(df.values)
        str(df)

        # GH8623
        x = DataFrame([[1, 'John P. Doe'], [2, 'Jane Dove'],
                       [1, 'John P. Doe']],
                      columns=['person_id', 'person_name'])
        x['person_name'] = Categorical(x.person_name
                                       )  # doing this breaks transform

        expected = x.iloc[0].person_name
        result = x.person_name.iloc[0]
        assert result == expected

        result = x.person_name[0]
        assert result == expected

        result = x.person_name.loc[0]
        assert result == expected

    def test_describe(self):

        # Categoricals should not show up together with numerical columns
        result = self.cat.describe()
        assert len(result.columns) == 1

        # In a frame, describe() for the cat should be the same as for string
        # arrays (count, unique, top, freq)

        cat = Categorical(["a", "b", "b", "b"], categories=['a', 'b', 'c'],
                          ordered=True)
        s = Series(cat)
        result = s.describe()
        expected = Series([4, 2, "b", 3],
                          index=['count', 'unique', 'top', 'freq'])
        tm.assert_series_equal(result, expected)

        cat = Series(Categorical(["a", "b", "c", "c"]))
        df3 = DataFrame({"cat": cat, "s": ["a", "b", "c", "c"]})
        res = df3.describe()
        tm.assert_numpy_array_equal(res["cat"].values, res["s"].values)

    def test_groupby_sort(self):

        # http://stackoverflow.com/questions/23814368/sorting-pandas-categorical-labels-after-groupby
        # This should result in a properly sorted Series so that the plot
        # has a sorted x axis
        # self.cat.groupby(['value_group'])['value_group'].count().plot(kind='bar')

        res = self.cat.groupby(['value_group'])['value_group'].count()
        exp = res[sorted(res.index, key=lambda x: float(x.split()[0]))]
        exp.index = CategoricalIndex(exp.index, name=exp.index.name)
        tm.assert_series_equal(res, exp)

    def test_astype_to_other(self):

        s = self.cat['value_group']
        expected = s
        tm.assert_series_equal(s.astype('category'), expected)
        tm.assert_series_equal(s.astype(CategoricalDtype()), expected)
        pytest.raises(ValueError, lambda: s.astype('float64'))

        cat = Series(Categorical(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c']))
        exp = Series(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c'])
        tm.assert_series_equal(cat.astype('str'), exp)
        s2 = Series(Categorical(['1', '2', '3', '4']))
        exp2 = Series([1, 2, 3, 4]).astype(int)
        tm.assert_series_equal(s2.astype('int'), exp2)

        # object don't sort correctly, so just compare that we have the same
        # values
        def cmp(a, b):
            tm.assert_almost_equal(
                np.sort(np.unique(a)), np.sort(np.unique(b)))

        expected = Series(np.array(s.values), name='value_group')
        cmp(s.astype('object'), expected)
        cmp(s.astype(np.object_), expected)

        # array conversion
        tm.assert_almost_equal(np.array(s), np.array(s.values))

        # valid conversion
        for valid in [lambda x: x.astype('category'),
                      lambda x: x.astype(CategoricalDtype()),
                      lambda x: x.astype('object').astype('category'),
                      lambda x: x.astype('object').astype(
                          CategoricalDtype())
                      ]:

            result = valid(s)
            # compare series values
            # internal .categories can't be compared because it is sorted
            tm.assert_series_equal(result, s, check_categorical=False)

        # invalid conversion (these are NOT a dtype)
        for invalid in [lambda x: x.astype(Categorical),
                        lambda x: x.astype('object').astype(Categorical)]:
            pytest.raises(TypeError, lambda: invalid(s))

    def test_numeric_like_ops(self):

        # numeric ops should not succeed
        for op in ['__add__', '__sub__', '__mul__', '__truediv__']:
            pytest.raises(TypeError,
                          lambda: getattr(self.cat, op)(self.cat))

        # reduction ops should not succeed (unless specifically defined, e.g.
        # min/max)
        s = self.cat['value_group']
        for op in ['kurt', 'skew', 'var', 'std', 'mean', 'sum', 'median']:
            pytest.raises(TypeError,
                          lambda: getattr(s, op)(numeric_only=False))

        # mad technically works because it takes always the numeric data

        # numpy ops
        s = Series(Categorical([1, 2, 3, 4]))
        pytest.raises(TypeError, lambda: np.sum(s))

        # numeric ops on a Series
        for op in ['__add__', '__sub__', '__mul__', '__truediv__']:
            pytest.raises(TypeError, lambda: getattr(s, op)(2))

        # invalid ufunc
        pytest.raises(TypeError, lambda: np.log(s))
