# -*- coding: utf-8 -*-

import numpy as np

import pandas.util.testing as tm
from pandas import (Categorical, Index, Series, DataFrame)


class TestCategoricalSort(object):

    def test_argsort(self):
        c = Categorical([5, 3, 1, 4, 2], ordered=True)

        expected = np.array([2, 4, 1, 3, 0])
        tm.assert_numpy_array_equal(c.argsort(ascending=True), expected,
                                    check_dtype=False)

        expected = expected[::-1]
        tm.assert_numpy_array_equal(c.argsort(ascending=False), expected,
                                    check_dtype=False)

    def test_numpy_argsort(self):
        c = Categorical([5, 3, 1, 4, 2], ordered=True)

        expected = np.array([2, 4, 1, 3, 0])
        tm.assert_numpy_array_equal(np.argsort(c), expected,
                                    check_dtype=False)

        tm.assert_numpy_array_equal(np.argsort(c, kind='mergesort'), expected,
                                    check_dtype=False)

        msg = "the 'axis' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, np.argsort,
                               c, axis=0)

        msg = "the 'order' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, np.argsort,
                               c, order='C')

    def test_sort_values(self):

        # unordered cats are sortable
        cat = Categorical(["a", "b", "b", "a"], ordered=False)
        cat.sort_values()

        cat = Categorical(["a", "c", "b", "d"], ordered=True)

        # sort_values
        res = cat.sort_values()
        exp = np.array(["a", "b", "c", "d"], dtype=object)
        tm.assert_numpy_array_equal(res.__array__(), exp)
        tm.assert_index_equal(res.categories, cat.categories)

        cat = Categorical(["a", "c", "b", "d"],
                          categories=["a", "b", "c", "d"], ordered=True)
        res = cat.sort_values()
        exp = np.array(["a", "b", "c", "d"], dtype=object)
        tm.assert_numpy_array_equal(res.__array__(), exp)
        tm.assert_index_equal(res.categories, cat.categories)

        res = cat.sort_values(ascending=False)
        exp = np.array(["d", "c", "b", "a"], dtype=object)
        tm.assert_numpy_array_equal(res.__array__(), exp)
        tm.assert_index_equal(res.categories, cat.categories)

        # sort (inplace order)
        cat1 = cat.copy()
        cat1.sort_values(inplace=True)
        exp = np.array(["a", "b", "c", "d"], dtype=object)
        tm.assert_numpy_array_equal(cat1.__array__(), exp)
        tm.assert_index_equal(res.categories, cat.categories)

        # reverse
        cat = Categorical(["a", "c", "c", "b", "d"], ordered=True)
        res = cat.sort_values(ascending=False)
        exp_val = np.array(["d", "c", "c", "b", "a"], dtype=object)
        exp_categories = Index(["a", "b", "c", "d"])
        tm.assert_numpy_array_equal(res.__array__(), exp_val)
        tm.assert_index_equal(res.categories, exp_categories)

    def test_sort_values_na_position(self):
        # see gh-12882
        cat = Categorical([5, 2, np.nan, 2, np.nan], ordered=True)
        exp_categories = Index([2, 5])

        exp = np.array([2.0, 2.0, 5.0, np.nan, np.nan])
        res = cat.sort_values()  # default arguments
        tm.assert_numpy_array_equal(res.__array__(), exp)
        tm.assert_index_equal(res.categories, exp_categories)

        exp = np.array([np.nan, np.nan, 2.0, 2.0, 5.0])
        res = cat.sort_values(ascending=True, na_position='first')
        tm.assert_numpy_array_equal(res.__array__(), exp)
        tm.assert_index_equal(res.categories, exp_categories)

        exp = np.array([np.nan, np.nan, 5.0, 2.0, 2.0])
        res = cat.sort_values(ascending=False, na_position='first')
        tm.assert_numpy_array_equal(res.__array__(), exp)
        tm.assert_index_equal(res.categories, exp_categories)

        exp = np.array([2.0, 2.0, 5.0, np.nan, np.nan])
        res = cat.sort_values(ascending=True, na_position='last')
        tm.assert_numpy_array_equal(res.__array__(), exp)
        tm.assert_index_equal(res.categories, exp_categories)

        exp = np.array([5.0, 2.0, 2.0, np.nan, np.nan])
        res = cat.sort_values(ascending=False, na_position='last')
        tm.assert_numpy_array_equal(res.__array__(), exp)
        tm.assert_index_equal(res.categories, exp_categories)

        cat = Categorical(["a", "c", "b", "d", np.nan], ordered=True)
        res = cat.sort_values(ascending=False, na_position='last')
        exp_val = np.array(["d", "c", "b", "a", np.nan], dtype=object)
        exp_categories = Index(["a", "b", "c", "d"])
        tm.assert_numpy_array_equal(res.__array__(), exp_val)
        tm.assert_index_equal(res.categories, exp_categories)

        cat = Categorical(["a", "c", "b", "d", np.nan], ordered=True)
        res = cat.sort_values(ascending=False, na_position='first')
        exp_val = np.array([np.nan, "d", "c", "b", "a"], dtype=object)
        exp_categories = Index(["a", "b", "c", "d"])
        tm.assert_numpy_array_equal(res.__array__(), exp_val)
        tm.assert_index_equal(res.categories, exp_categories)


class TestCategoricalBlockSort(object):

    def test_sort_values(self):

        c = Categorical(["a", "b", "b", "a"], ordered=False)
        cat = Series(c.copy())

        # sort in the categories order
        expected = Series(
            Categorical(["a", "a", "b", "b"],
                        ordered=False), index=[0, 3, 1, 2])
        result = cat.sort_values()
        tm.assert_series_equal(result, expected)

        cat = Series(Categorical(["a", "c", "b", "d"], ordered=True))
        res = cat.sort_values()
        exp = np.array(["a", "b", "c", "d"], dtype=np.object_)
        tm.assert_numpy_array_equal(res.__array__(), exp)

        cat = Series(Categorical(["a", "c", "b", "d"], categories=[
                     "a", "b", "c", "d"], ordered=True))
        res = cat.sort_values()
        exp = np.array(["a", "b", "c", "d"], dtype=np.object_)
        tm.assert_numpy_array_equal(res.__array__(), exp)

        res = cat.sort_values(ascending=False)
        exp = np.array(["d", "c", "b", "a"], dtype=np.object_)
        tm.assert_numpy_array_equal(res.__array__(), exp)

        raw_cat1 = Categorical(["a", "b", "c", "d"],
                               categories=["a", "b", "c", "d"], ordered=False)
        raw_cat2 = Categorical(["a", "b", "c", "d"],
                               categories=["d", "c", "b", "a"], ordered=True)
        s = ["a", "b", "c", "d"]
        df = DataFrame({"unsort": raw_cat1,
                        "sort": raw_cat2,
                        "string": s,
                        "values": [1, 2, 3, 4]})

        # Cats must be sorted in a dataframe
        res = df.sort_values(by=["string"], ascending=False)
        exp = np.array(["d", "c", "b", "a"], dtype=np.object_)
        tm.assert_numpy_array_equal(res["sort"].values.__array__(), exp)
        assert res["sort"].dtype == "category"

        res = df.sort_values(by=["sort"], ascending=False)
        exp = df.sort_values(by=["string"], ascending=True)
        tm.assert_series_equal(res["values"], exp["values"])
        assert res["sort"].dtype == "category"
        assert res["unsort"].dtype == "category"

        # unordered cat, but we allow this
        df.sort_values(by=["unsort"], ascending=False)

        # multi-columns sort
        # GH 7848
        df = DataFrame({"id": [6, 5, 4, 3, 2, 1],
                        "raw_grade": ['a', 'b', 'b', 'a', 'a', 'e']})
        df["grade"] = Categorical(df["raw_grade"], ordered=True)
        df['grade'] = df['grade'].cat.set_categories(['b', 'e', 'a'])

        # sorts 'grade' according to the order of the categories
        result = df.sort_values(by=['grade'])
        expected = df.iloc[[1, 2, 5, 0, 3, 4]]
        tm.assert_frame_equal(result, expected)

        # multi
        result = df.sort_values(by=['grade', 'id'])
        expected = df.iloc[[2, 1, 5, 4, 3, 0]]
        tm.assert_frame_equal(result, expected)
