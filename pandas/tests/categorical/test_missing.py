# -*- coding: utf-8 -*-

import numpy as np

import pandas.util.testing as tm
from pandas import (Categorical, Index, Series, CategoricalIndex, isna)
from pandas.compat import lrange
from pandas.core.dtypes.dtypes import CategoricalDtype


class TestCategoricalMissing(object):

    def test_na_flags_int_categories(self):
        # #1457

        categories = lrange(10)
        labels = np.random.randint(0, 10, 20)
        labels[::5] = -1

        cat = Categorical(labels, categories, fastpath=True)
        repr(cat)

        tm.assert_numpy_array_equal(isna(cat), labels == -1)

    def test_nan_handling(self):

        # Nans are represented as -1 in codes
        c = Categorical(["a", "b", np.nan, "a"])
        tm.assert_index_equal(c.categories, Index(["a", "b"]))
        tm.assert_numpy_array_equal(c._codes, np.array([0, 1, -1, 0],
                                                       dtype=np.int8))
        c[1] = np.nan
        tm.assert_index_equal(c.categories, Index(["a", "b"]))
        tm.assert_numpy_array_equal(c._codes, np.array([0, -1, -1, 0],
                                                       dtype=np.int8))

        # Adding nan to categories should make assigned nan point to the
        # category!
        c = Categorical(["a", "b", np.nan, "a"])
        tm.assert_index_equal(c.categories, Index(["a", "b"]))
        tm.assert_numpy_array_equal(c._codes, np.array([0, 1, -1, 0],
                                                       dtype=np.int8))

    def test_set_dtype_nans(self):
        c = Categorical(['a', 'b', np.nan])
        result = c._set_dtype(CategoricalDtype(['a', 'c']))
        tm.assert_numpy_array_equal(result.codes, np.array([0, -1, -1],
                                                           dtype='int8'))

    def test_isna(self):
        exp = np.array([False, False, True])
        c = Categorical(["a", "b", np.nan])
        res = c.isna()

        tm.assert_numpy_array_equal(res, exp)

    def test_set_item_nan(self):
        cat = Categorical([1, 2, 3])
        cat[1] = np.nan

        exp = Categorical([1, np.nan, 3], categories=[1, 2, 3])
        tm.assert_categorical_equal(cat, exp)


class TestCategoricalBlockMissing(object):

    def test_value_counts_with_nan(self):
        # see gh-9443

        # sanity check
        s = Series(["a", "b", "a"], dtype="category")
        exp = Series([2, 1], index=CategoricalIndex(["a", "b"]))

        res = s.value_counts(dropna=True)
        tm.assert_series_equal(res, exp)

        res = s.value_counts(dropna=True)
        tm.assert_series_equal(res, exp)

        # same Series via two different constructions --> same behaviour
        series = [
            Series(["a", "b", None, "a", None, None], dtype="category"),
            Series(Categorical(["a", "b", None, "a", None, None],
                               categories=["a", "b"]))
        ]

        for s in series:
            # None is a NaN value, so we exclude its count here
            exp = Series([2, 1], index=CategoricalIndex(["a", "b"]))
            res = s.value_counts(dropna=True)
            tm.assert_series_equal(res, exp)

            # we don't exclude the count of None and sort by counts
            exp = Series([3, 2, 1], index=CategoricalIndex([np.nan, "a", "b"]))
            res = s.value_counts(dropna=False)
            tm.assert_series_equal(res, exp)

            # When we aren't sorting by counts, and np.nan isn't a
            # category, it should be last.
            exp = Series([2, 1, 3], index=CategoricalIndex(["a", "b", np.nan]))
            res = s.value_counts(dropna=False, sort=False)
            tm.assert_series_equal(res, exp)

    def test_nan_handling(self):

        # NaNs are represented as -1 in labels
        s = Series(Categorical(["a", "b", np.nan, "a"]))
        tm.assert_index_equal(s.cat.categories, Index(["a", "b"]))
        tm.assert_numpy_array_equal(s.values.codes,
                                    np.array([0, 1, -1, 0], dtype=np.int8))
