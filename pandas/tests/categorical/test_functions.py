# -*- coding: utf-8 -*-

import pytest

import numpy as np

import pandas.util.testing as tm
from pandas import (Categorical, Index, Series, DataFrame, CategoricalIndex)
from pandas.tests.categorical.common import TestCategorical


class TestCategoricalGeneric(TestCategorical):

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

        # Make sure that unequal comparison take the categories order in
        # account
        cat_rev = Categorical(
            list("abc"), categories=list("cba"), ordered=True)
        exp = np.array([True, False, False])
        res = cat_rev > "b"
        tm.assert_numpy_array_equal(res, exp)
