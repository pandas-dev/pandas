# -*- coding: utf-8 -*-

import pytest
import sys

import numpy as np

import pandas.util.testing as tm
from pandas import Categorical, CategoricalIndex, Index, Series, DataFrame

from pandas.compat import PYPY
from pandas.core.categorical import _recode_for_categories
from pandas.tests.categorical.common import TestCategorical


class TestCategoricalAPI(object):

    def test_searchsorted(self):
        # https://github.com/pandas-dev/pandas/issues/8420
        # https://github.com/pandas-dev/pandas/issues/14522

        c1 = Categorical(['cheese', 'milk', 'apple', 'bread', 'bread'],
                         categories=['cheese', 'milk', 'apple', 'bread'],
                         ordered=True)
        s1 = Series(c1)
        c2 = Categorical(['cheese', 'milk', 'apple', 'bread', 'bread'],
                         categories=['cheese', 'milk', 'apple', 'bread'],
                         ordered=False)
        s2 = Series(c2)

        # Searching for single item argument, side='left' (default)
        res_cat = c1.searchsorted('apple')
        res_ser = s1.searchsorted('apple')
        exp = np.array([2], dtype=np.intp)
        tm.assert_numpy_array_equal(res_cat, exp)
        tm.assert_numpy_array_equal(res_ser, exp)

        # Searching for single item array, side='left' (default)
        res_cat = c1.searchsorted(['bread'])
        res_ser = s1.searchsorted(['bread'])
        exp = np.array([3], dtype=np.intp)
        tm.assert_numpy_array_equal(res_cat, exp)
        tm.assert_numpy_array_equal(res_ser, exp)

        # Searching for several items array, side='right'
        res_cat = c1.searchsorted(['apple', 'bread'], side='right')
        res_ser = s1.searchsorted(['apple', 'bread'], side='right')
        exp = np.array([3, 5], dtype=np.intp)
        tm.assert_numpy_array_equal(res_cat, exp)
        tm.assert_numpy_array_equal(res_ser, exp)

        # Searching for a single value that is not from the Categorical
        pytest.raises(ValueError, lambda: c1.searchsorted('cucumber'))
        pytest.raises(ValueError, lambda: s1.searchsorted('cucumber'))

        # Searching for multiple values one of each is not from the Categorical
        pytest.raises(ValueError,
                      lambda: c1.searchsorted(['bread', 'cucumber']))
        pytest.raises(ValueError,
                      lambda: s1.searchsorted(['bread', 'cucumber']))

        # searchsorted call for unordered Categorical
        pytest.raises(ValueError, lambda: c2.searchsorted('apple'))
        pytest.raises(ValueError, lambda: s2.searchsorted('apple'))

        with tm.assert_produces_warning(FutureWarning):
            res = c1.searchsorted(v=['bread'])
            exp = np.array([3], dtype=np.intp)
            tm.assert_numpy_array_equal(res, exp)

    def test_ordered_api(self):
        # GH 9347
        cat1 = Categorical(list('acb'), ordered=False)
        tm.assert_index_equal(cat1.categories, Index(['a', 'b', 'c']))
        assert not cat1.ordered

        cat2 = Categorical(list('acb'), categories=list('bca'), ordered=False)
        tm.assert_index_equal(cat2.categories, Index(['b', 'c', 'a']))
        assert not cat2.ordered

        cat3 = Categorical(list('acb'), ordered=True)
        tm.assert_index_equal(cat3.categories, Index(['a', 'b', 'c']))
        assert cat3.ordered

        cat4 = Categorical(list('acb'), categories=list('bca'), ordered=True)
        tm.assert_index_equal(cat4.categories, Index(['b', 'c', 'a']))
        assert cat4.ordered

    def test_set_ordered(self):

        cat = Categorical(["a", "b", "c", "a"], ordered=True)
        cat2 = cat.as_unordered()
        assert not cat2.ordered
        cat2 = cat.as_ordered()
        assert cat2.ordered
        cat2.as_unordered(inplace=True)
        assert not cat2.ordered
        cat2.as_ordered(inplace=True)
        assert cat2.ordered

        assert cat2.set_ordered(True).ordered
        assert not cat2.set_ordered(False).ordered
        cat2.set_ordered(True, inplace=True)
        assert cat2.ordered
        cat2.set_ordered(False, inplace=True)
        assert not cat2.ordered

        # removed in 0.19.0
        msg = "can\'t set attribute"
        with tm.assert_raises_regex(AttributeError, msg):
            cat.ordered = True
        with tm.assert_raises_regex(AttributeError, msg):
            cat.ordered = False

    def test_rename_categories(self):
        cat = Categorical(["a", "b", "c", "a"])

        # inplace=False: the old one must not be changed
        res = cat.rename_categories([1, 2, 3])
        tm.assert_numpy_array_equal(res.__array__(), np.array([1, 2, 3, 1],
                                                              dtype=np.int64))
        tm.assert_index_equal(res.categories, Index([1, 2, 3]))

        exp_cat = np.array(["a", "b", "c", "a"], dtype=np.object_)
        tm.assert_numpy_array_equal(cat.__array__(), exp_cat)

        exp_cat = Index(["a", "b", "c"])
        tm.assert_index_equal(cat.categories, exp_cat)
        res = cat.rename_categories([1, 2, 3], inplace=True)

        # and now inplace
        assert res is None
        tm.assert_numpy_array_equal(cat.__array__(), np.array([1, 2, 3, 1],
                                                              dtype=np.int64))
        tm.assert_index_equal(cat.categories, Index([1, 2, 3]))

        # Lengthen
        with pytest.raises(ValueError):
            cat.rename_categories([1, 2, 3, 4])

        # Shorten
        with pytest.raises(ValueError):
            cat.rename_categories([1, 2])

    def test_rename_categories_series(self):
        # https://github.com/pandas-dev/pandas/issues/17981
        c = Categorical(['a', 'b'])
        xpr = "Treating Series 'new_categories' as a list-like "
        with tm.assert_produces_warning(FutureWarning) as rec:
            result = c.rename_categories(Series([0, 1]))

        assert len(rec) == 1
        assert xpr in str(rec[0].message)
        expected = Categorical([0, 1])
        tm.assert_categorical_equal(result, expected)

    def test_rename_categories_dict(self):
        # GH 17336
        cat = Categorical(['a', 'b', 'c', 'd'])
        res = cat.rename_categories({'a': 4, 'b': 3, 'c': 2, 'd': 1})
        expected = Index([4, 3, 2, 1])
        tm.assert_index_equal(res.categories, expected)

        # Test for inplace
        res = cat.rename_categories({'a': 4, 'b': 3, 'c': 2, 'd': 1},
                                    inplace=True)
        assert res is None
        tm.assert_index_equal(cat.categories, expected)

        # Test for dicts of smaller length
        cat = Categorical(['a', 'b', 'c', 'd'])
        res = cat.rename_categories({'a': 1, 'c': 3})

        expected = Index([1, 'b', 3, 'd'])
        tm.assert_index_equal(res.categories, expected)

        # Test for dicts with bigger length
        cat = Categorical(['a', 'b', 'c', 'd'])
        res = cat.rename_categories({'a': 1, 'b': 2, 'c': 3,
                                     'd': 4, 'e': 5, 'f': 6})
        expected = Index([1, 2, 3, 4])
        tm.assert_index_equal(res.categories, expected)

        # Test for dicts with no items from old categories
        cat = Categorical(['a', 'b', 'c', 'd'])
        res = cat.rename_categories({'f': 1, 'g': 3})

        expected = Index(['a', 'b', 'c', 'd'])
        tm.assert_index_equal(res.categories, expected)

    @pytest.mark.parametrize('codes, old, new, expected', [
        ([0, 1], ['a', 'b'], ['a', 'b'], [0, 1]),
        ([0, 1], ['b', 'a'], ['b', 'a'], [0, 1]),
        ([0, 1], ['a', 'b'], ['b', 'a'], [1, 0]),
        ([0, 1], ['b', 'a'], ['a', 'b'], [1, 0]),
        ([0, 1, 0, 1], ['a', 'b'], ['a', 'b', 'c'], [0, 1, 0, 1]),
        ([0, 1, 2, 2], ['a', 'b', 'c'], ['a', 'b'], [0, 1, -1, -1]),
        ([0, 1, -1], ['a', 'b', 'c'], ['a', 'b', 'c'], [0, 1, -1]),
        ([0, 1, -1], ['a', 'b', 'c'], ['b'], [-1, 0, -1]),
        ([0, 1, -1], ['a', 'b', 'c'], ['d'], [-1, -1, -1]),
        ([0, 1, -1], ['a', 'b', 'c'], [], [-1, -1, -1]),
        ([-1, -1], [], ['a', 'b'], [-1, -1]),
        ([1, 0], ['b', 'a'], ['a', 'b'], [0, 1]),
    ])
    def test_recode_to_categories(self, codes, old, new, expected):
        codes = np.asanyarray(codes, dtype=np.int8)
        expected = np.asanyarray(expected, dtype=np.int8)
        old = Index(old)
        new = Index(new)
        result = _recode_for_categories(codes, old, new)
        tm.assert_numpy_array_equal(result, expected)

    def test_recode_to_categories_large(self):
        N = 1000
        codes = np.arange(N)
        old = Index(codes)
        expected = np.arange(N - 1, -1, -1, dtype=np.int16)
        new = Index(expected)
        result = _recode_for_categories(codes, old, new)
        tm.assert_numpy_array_equal(result, expected)

    def test_reorder_categories(self):
        cat = Categorical(["a", "b", "c", "a"], ordered=True)
        old = cat.copy()
        new = Categorical(["a", "b", "c", "a"], categories=["c", "b", "a"],
                          ordered=True)

        # first inplace == False
        res = cat.reorder_categories(["c", "b", "a"])
        # cat must be the same as before
        tm.assert_categorical_equal(cat, old)
        # only res is changed
        tm.assert_categorical_equal(res, new)

        # inplace == True
        res = cat.reorder_categories(["c", "b", "a"], inplace=True)
        assert res is None
        tm.assert_categorical_equal(cat, new)

        # not all "old" included in "new"
        cat = Categorical(["a", "b", "c", "a"], ordered=True)

        def f():
            cat.reorder_categories(["a"])

        pytest.raises(ValueError, f)

        # still not all "old" in "new"
        def f():
            cat.reorder_categories(["a", "b", "d"])

        pytest.raises(ValueError, f)

        # all "old" included in "new", but too long
        def f():
            cat.reorder_categories(["a", "b", "c", "d"])

        pytest.raises(ValueError, f)

    def test_add_categories(self):
        cat = Categorical(["a", "b", "c", "a"], ordered=True)
        old = cat.copy()
        new = Categorical(["a", "b", "c", "a"],
                          categories=["a", "b", "c", "d"], ordered=True)

        # first inplace == False
        res = cat.add_categories("d")
        tm.assert_categorical_equal(cat, old)
        tm.assert_categorical_equal(res, new)

        res = cat.add_categories(["d"])
        tm.assert_categorical_equal(cat, old)
        tm.assert_categorical_equal(res, new)

        # inplace == True
        res = cat.add_categories("d", inplace=True)
        tm.assert_categorical_equal(cat, new)
        assert res is None

        # new is in old categories
        def f():
            cat.add_categories(["d"])

        pytest.raises(ValueError, f)

        # GH 9927
        cat = Categorical(list("abc"), ordered=True)
        expected = Categorical(
            list("abc"), categories=list("abcde"), ordered=True)
        # test with Series, np.array, index, list
        res = cat.add_categories(Series(["d", "e"]))
        tm.assert_categorical_equal(res, expected)
        res = cat.add_categories(np.array(["d", "e"]))
        tm.assert_categorical_equal(res, expected)
        res = cat.add_categories(Index(["d", "e"]))
        tm.assert_categorical_equal(res, expected)
        res = cat.add_categories(["d", "e"])
        tm.assert_categorical_equal(res, expected)

    def test_set_categories(self):
        cat = Categorical(["a", "b", "c", "a"], ordered=True)
        exp_categories = Index(["c", "b", "a"])
        exp_values = np.array(["a", "b", "c", "a"], dtype=np.object_)

        res = cat.set_categories(["c", "b", "a"], inplace=True)
        tm.assert_index_equal(cat.categories, exp_categories)
        tm.assert_numpy_array_equal(cat.__array__(), exp_values)
        assert res is None

        res = cat.set_categories(["a", "b", "c"])
        # cat must be the same as before
        tm.assert_index_equal(cat.categories, exp_categories)
        tm.assert_numpy_array_equal(cat.__array__(), exp_values)
        # only res is changed
        exp_categories_back = Index(["a", "b", "c"])
        tm.assert_index_equal(res.categories, exp_categories_back)
        tm.assert_numpy_array_equal(res.__array__(), exp_values)

        # not all "old" included in "new" -> all not included ones are now
        # np.nan
        cat = Categorical(["a", "b", "c", "a"], ordered=True)
        res = cat.set_categories(["a"])
        tm.assert_numpy_array_equal(res.codes, np.array([0, -1, -1, 0],
                                                        dtype=np.int8))

        # still not all "old" in "new"
        res = cat.set_categories(["a", "b", "d"])
        tm.assert_numpy_array_equal(res.codes, np.array([0, 1, -1, 0],
                                                        dtype=np.int8))
        tm.assert_index_equal(res.categories, Index(["a", "b", "d"]))

        # all "old" included in "new"
        cat = cat.set_categories(["a", "b", "c", "d"])
        exp_categories = Index(["a", "b", "c", "d"])
        tm.assert_index_equal(cat.categories, exp_categories)

        # internals...
        c = Categorical([1, 2, 3, 4, 1], categories=[1, 2, 3, 4], ordered=True)
        tm.assert_numpy_array_equal(c._codes, np.array([0, 1, 2, 3, 0],
                                                       dtype=np.int8))
        tm.assert_index_equal(c.categories, Index([1, 2, 3, 4]))

        exp = np.array([1, 2, 3, 4, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(c.get_values(), exp)

        # all "pointers" to '4' must be changed from 3 to 0,...
        c = c.set_categories([4, 3, 2, 1])

        # positions are changed
        tm.assert_numpy_array_equal(c._codes, np.array([3, 2, 1, 0, 3],
                                                       dtype=np.int8))

        # categories are now in new order
        tm.assert_index_equal(c.categories, Index([4, 3, 2, 1]))

        # output is the same
        exp = np.array([1, 2, 3, 4, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(c.get_values(), exp)
        assert c.min() == 4
        assert c.max() == 1

        # set_categories should set the ordering if specified
        c2 = c.set_categories([4, 3, 2, 1], ordered=False)
        assert not c2.ordered

        tm.assert_numpy_array_equal(c.get_values(), c2.get_values())

        # set_categories should pass thru the ordering
        c2 = c.set_ordered(False).set_categories([4, 3, 2, 1])
        assert not c2.ordered

        tm.assert_numpy_array_equal(c.get_values(), c2.get_values())

    @pytest.mark.parametrize('values, categories, new_categories', [
        # No NaNs, same cats, same order
        (['a', 'b', 'a'], ['a', 'b'], ['a', 'b'],),
        # No NaNs, same cats, different order
        (['a', 'b', 'a'], ['a', 'b'], ['b', 'a'],),
        # Same, unsorted
        (['b', 'a', 'a'], ['a', 'b'], ['a', 'b'],),
        # No NaNs, same cats, different order
        (['b', 'a', 'a'], ['a', 'b'], ['b', 'a'],),
        # NaNs
        (['a', 'b', 'c'], ['a', 'b'], ['a', 'b']),
        (['a', 'b', 'c'], ['a', 'b'], ['b', 'a']),
        (['b', 'a', 'c'], ['a', 'b'], ['a', 'b']),
        (['b', 'a', 'c'], ['a', 'b'], ['a', 'b']),
        # Introduce NaNs
        (['a', 'b', 'c'], ['a', 'b'], ['a']),
        (['a', 'b', 'c'], ['a', 'b'], ['b']),
        (['b', 'a', 'c'], ['a', 'b'], ['a']),
        (['b', 'a', 'c'], ['a', 'b'], ['a']),
        # No overlap
        (['a', 'b', 'c'], ['a', 'b'], ['d', 'e']),
    ])
    @pytest.mark.parametrize('ordered', [True, False])
    def test_set_categories_many(self, values, categories, new_categories,
                                 ordered):
        c = Categorical(values, categories)
        expected = Categorical(values, new_categories, ordered)
        result = c.set_categories(new_categories, ordered=ordered)
        tm.assert_categorical_equal(result, expected)

    def test_set_categories_private(self):
        cat = Categorical(['a', 'b', 'c'], categories=['a', 'b', 'c', 'd'])
        cat._set_categories(['a', 'c', 'd', 'e'])
        expected = Categorical(['a', 'c', 'd'], categories=list('acde'))
        tm.assert_categorical_equal(cat, expected)

        # fastpath
        cat = Categorical(['a', 'b', 'c'], categories=['a', 'b', 'c', 'd'])
        cat._set_categories(['a', 'c', 'd', 'e'], fastpath=True)
        expected = Categorical(['a', 'c', 'd'], categories=list('acde'))
        tm.assert_categorical_equal(cat, expected)

    def test_remove_categories(self):
        cat = Categorical(["a", "b", "c", "a"], ordered=True)
        old = cat.copy()
        new = Categorical(["a", "b", np.nan, "a"], categories=["a", "b"],
                          ordered=True)

        # first inplace == False
        res = cat.remove_categories("c")
        tm.assert_categorical_equal(cat, old)
        tm.assert_categorical_equal(res, new)

        res = cat.remove_categories(["c"])
        tm.assert_categorical_equal(cat, old)
        tm.assert_categorical_equal(res, new)

        # inplace == True
        res = cat.remove_categories("c", inplace=True)
        tm.assert_categorical_equal(cat, new)
        assert res is None

        # removal is not in categories
        def f():
            cat.remove_categories(["c"])

        pytest.raises(ValueError, f)

    def test_remove_unused_categories(self):
        c = Categorical(["a", "b", "c", "d", "a"],
                        categories=["a", "b", "c", "d", "e"])
        exp_categories_all = Index(["a", "b", "c", "d", "e"])
        exp_categories_dropped = Index(["a", "b", "c", "d"])

        tm.assert_index_equal(c.categories, exp_categories_all)

        res = c.remove_unused_categories()
        tm.assert_index_equal(res.categories, exp_categories_dropped)
        tm.assert_index_equal(c.categories, exp_categories_all)

        res = c.remove_unused_categories(inplace=True)
        tm.assert_index_equal(c.categories, exp_categories_dropped)
        assert res is None

        # with NaN values (GH11599)
        c = Categorical(["a", "b", "c", np.nan],
                        categories=["a", "b", "c", "d", "e"])
        res = c.remove_unused_categories()
        tm.assert_index_equal(res.categories,
                              Index(np.array(["a", "b", "c"])))
        exp_codes = np.array([0, 1, 2, -1], dtype=np.int8)
        tm.assert_numpy_array_equal(res.codes, exp_codes)
        tm.assert_index_equal(c.categories, exp_categories_all)

        val = ['F', np.nan, 'D', 'B', 'D', 'F', np.nan]
        cat = Categorical(values=val, categories=list('ABCDEFG'))
        out = cat.remove_unused_categories()
        tm.assert_index_equal(out.categories, Index(['B', 'D', 'F']))
        exp_codes = np.array([2, -1, 1, 0, 1, 2, -1], dtype=np.int8)
        tm.assert_numpy_array_equal(out.codes, exp_codes)
        assert out.get_values().tolist() == val

        alpha = list('abcdefghijklmnopqrstuvwxyz')
        val = np.random.choice(alpha[::2], 10000).astype('object')
        val[np.random.choice(len(val), 100)] = np.nan

        cat = Categorical(values=val, categories=alpha)
        out = cat.remove_unused_categories()
        assert out.get_values().tolist() == val.tolist()

    def test_codes_immutable(self):

        # Codes should be read only
        c = Categorical(["a", "b", "c", "a", np.nan])
        exp = np.array([0, 1, 2, 0, -1], dtype='int8')
        tm.assert_numpy_array_equal(c.codes, exp)

        # Assignments to codes should raise
        def f():
            c.codes = np.array([0, 1, 2, 0, 1], dtype='int8')

        pytest.raises(ValueError, f)

        # changes in the codes array should raise
        # np 1.6.1 raises RuntimeError rather than ValueError
        codes = c.codes

        def f():
            codes[4] = 1

        pytest.raises(ValueError, f)

        # But even after getting the codes, the original array should still be
        # writeable!
        c[4] = "a"
        exp = np.array([0, 1, 2, 0, 0], dtype='int8')
        tm.assert_numpy_array_equal(c.codes, exp)
        c._codes[4] = 2
        exp = np.array([0, 1, 2, 0, 2], dtype='int8')
        tm.assert_numpy_array_equal(c.codes, exp)

    def test_unique(self):
        # categories are reordered based on value when ordered=False
        cat = Categorical(["a", "b"])
        exp = Index(["a", "b"])
        res = cat.unique()
        tm.assert_index_equal(res.categories, exp)
        tm.assert_categorical_equal(res, cat)

        cat = Categorical(["a", "b", "a", "a"], categories=["a", "b", "c"])
        res = cat.unique()
        tm.assert_index_equal(res.categories, exp)
        tm.assert_categorical_equal(res, Categorical(exp))

        cat = Categorical(["c", "a", "b", "a", "a"],
                          categories=["a", "b", "c"])
        exp = Index(["c", "a", "b"])
        res = cat.unique()
        tm.assert_index_equal(res.categories, exp)
        exp_cat = Categorical(exp, categories=['c', 'a', 'b'])
        tm.assert_categorical_equal(res, exp_cat)

        # nan must be removed
        cat = Categorical(["b", np.nan, "b", np.nan, "a"],
                          categories=["a", "b", "c"])
        res = cat.unique()
        exp = Index(["b", "a"])
        tm.assert_index_equal(res.categories, exp)
        exp_cat = Categorical(["b", np.nan, "a"], categories=["b", "a"])
        tm.assert_categorical_equal(res, exp_cat)

    def test_unique_ordered(self):
        # keep categories order when ordered=True
        cat = Categorical(['b', 'a', 'b'], categories=['a', 'b'], ordered=True)
        res = cat.unique()
        exp_cat = Categorical(['b', 'a'], categories=['a', 'b'], ordered=True)
        tm.assert_categorical_equal(res, exp_cat)

        cat = Categorical(['c', 'b', 'a', 'a'], categories=['a', 'b', 'c'],
                          ordered=True)
        res = cat.unique()
        exp_cat = Categorical(['c', 'b', 'a'], categories=['a', 'b', 'c'],
                              ordered=True)
        tm.assert_categorical_equal(res, exp_cat)

        cat = Categorical(['b', 'a', 'a'], categories=['a', 'b', 'c'],
                          ordered=True)
        res = cat.unique()
        exp_cat = Categorical(['b', 'a'], categories=['a', 'b'], ordered=True)
        tm.assert_categorical_equal(res, exp_cat)

        cat = Categorical(['b', 'b', np.nan, 'a'], categories=['a', 'b', 'c'],
                          ordered=True)
        res = cat.unique()
        exp_cat = Categorical(['b', np.nan, 'a'], categories=['a', 'b'],
                              ordered=True)
        tm.assert_categorical_equal(res, exp_cat)

    def test_unique_index_series(self):
        c = Categorical([3, 1, 2, 2, 1], categories=[3, 2, 1])
        # Categorical.unique sorts categories by appearance order
        # if ordered=False
        exp = Categorical([3, 1, 2], categories=[3, 1, 2])
        tm.assert_categorical_equal(c.unique(), exp)

        tm.assert_index_equal(Index(c).unique(), Index(exp))
        tm.assert_categorical_equal(Series(c).unique(), exp)

        c = Categorical([1, 1, 2, 2], categories=[3, 2, 1])
        exp = Categorical([1, 2], categories=[1, 2])
        tm.assert_categorical_equal(c.unique(), exp)
        tm.assert_index_equal(Index(c).unique(), Index(exp))
        tm.assert_categorical_equal(Series(c).unique(), exp)

        c = Categorical([3, 1, 2, 2, 1], categories=[3, 2, 1], ordered=True)
        # Categorical.unique keeps categories order if ordered=True
        exp = Categorical([3, 1, 2], categories=[3, 2, 1], ordered=True)
        tm.assert_categorical_equal(c.unique(), exp)

        tm.assert_index_equal(Index(c).unique(), Index(exp))
        tm.assert_categorical_equal(Series(c).unique(), exp)

    def test_shift(self):
        # GH 9416
        cat = Categorical(['a', 'b', 'c', 'd', 'a'])

        # shift forward
        sp1 = cat.shift(1)
        xp1 = Categorical([np.nan, 'a', 'b', 'c', 'd'])
        tm.assert_categorical_equal(sp1, xp1)
        tm.assert_categorical_equal(cat[:-1], sp1[1:])

        # shift back
        sn2 = cat.shift(-2)
        xp2 = Categorical(['c', 'd', 'a', np.nan, np.nan],
                          categories=['a', 'b', 'c', 'd'])
        tm.assert_categorical_equal(sn2, xp2)
        tm.assert_categorical_equal(cat[2:], sn2[:-2])

        # shift by zero
        tm.assert_categorical_equal(cat, cat.shift(0))

    def test_nbytes(self):
        cat = Categorical([1, 2, 3])
        exp = 3 + 3 * 8  # 3 int8s for values + 3 int64s for categories
        assert cat.nbytes == exp

    def test_memory_usage(self):
        cat = Categorical([1, 2, 3])

        # .categories is an index, so we include the hashtable
        assert 0 < cat.nbytes <= cat.memory_usage()
        assert 0 < cat.nbytes <= cat.memory_usage(deep=True)

        cat = Categorical(['foo', 'foo', 'bar'])
        assert cat.memory_usage(deep=True) > cat.nbytes

        if not PYPY:
            # sys.getsizeof will call the .memory_usage with
            # deep=True, and add on some GC overhead
            diff = cat.memory_usage(deep=True) - sys.getsizeof(cat)
            assert abs(diff) < 100

    def test_deprecated_labels(self):
        # TODO: labels is deprecated and should be removed in 0.18 or 2017,
        # whatever is earlier
        cat = Categorical([1, 2, 3, np.nan], categories=[1, 2, 3])
        exp = cat.codes
        with tm.assert_produces_warning(FutureWarning):
            res = cat.labels
        tm.assert_numpy_array_equal(res, exp)

    def test_deprecated_from_array(self):
        # GH13854, `.from_array` is deprecated
        with tm.assert_produces_warning(FutureWarning):
            Categorical.from_array([0, 1])

    def test_map(self):
        c = Categorical(list('ABABC'), categories=list('CBA'), ordered=True)
        result = c.map(lambda x: x.lower())
        exp = Categorical(list('ababc'), categories=list('cba'), ordered=True)
        tm.assert_categorical_equal(result, exp)

        c = Categorical(list('ABABC'), categories=list('ABC'), ordered=False)
        result = c.map(lambda x: x.lower())
        exp = Categorical(list('ababc'), categories=list('abc'), ordered=False)
        tm.assert_categorical_equal(result, exp)

        result = c.map(lambda x: 1)
        # GH 12766: Return an index not an array
        tm.assert_index_equal(result, Index(np.array([1] * 5, dtype=np.int64)))

    def test_validate_inplace(self):
        cat = Categorical(['A', 'B', 'B', 'C', 'A'])
        invalid_values = [1, "True", [1, 2, 3], 5.0]

        for value in invalid_values:
            with pytest.raises(ValueError):
                cat.set_ordered(value=True, inplace=value)

            with pytest.raises(ValueError):
                cat.as_ordered(inplace=value)

            with pytest.raises(ValueError):
                cat.as_unordered(inplace=value)

            with pytest.raises(ValueError):
                cat.set_categories(['X', 'Y', 'Z'], rename=True, inplace=value)

            with pytest.raises(ValueError):
                cat.rename_categories(['X', 'Y', 'Z'], inplace=value)

            with pytest.raises(ValueError):
                cat.reorder_categories(
                    ['X', 'Y', 'Z'], ordered=True, inplace=value)

            with pytest.raises(ValueError):
                cat.add_categories(
                    new_categories=['D', 'E', 'F'], inplace=value)

            with pytest.raises(ValueError):
                cat.remove_categories(removals=['D', 'E', 'F'], inplace=value)

            with pytest.raises(ValueError):
                cat.remove_unused_categories(inplace=value)

            with pytest.raises(ValueError):
                cat.sort_values(inplace=value)

    def test_repeat(self):
        # GH10183
        cat = Categorical(["a", "b"], categories=["a", "b"])
        exp = Categorical(["a", "a", "b", "b"], categories=["a", "b"])
        res = cat.repeat(2)
        tm.assert_categorical_equal(res, exp)

    def test_numpy_repeat(self):
        cat = Categorical(["a", "b"], categories=["a", "b"])
        exp = Categorical(["a", "a", "b", "b"], categories=["a", "b"])
        tm.assert_categorical_equal(np.repeat(cat, 2), exp)

        msg = "the 'axis' parameter is not supported"
        tm.assert_raises_regex(ValueError, msg, np.repeat, cat, 2, axis=1)

    def test_astype_categorical(self):

        cat = Categorical(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c'])
        tm.assert_categorical_equal(cat, cat.astype('category'))
        tm.assert_almost_equal(np.array(cat), cat.astype('object'))

        pytest.raises(ValueError, lambda: cat.astype(float))

    def test_isna(self):
        exp = np.array([False, False, True])
        c = Categorical(["a", "b", np.nan])
        res = c.isna()

        tm.assert_numpy_array_equal(res, exp)


class TestCategoricalAPIWithFactor(TestCategorical):

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
