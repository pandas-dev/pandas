# -*- coding: utf-8 -*-

from warnings import catch_warnings
import pytest
import sys

import numpy as np

import pandas as pd
import pandas.compat as compat
import pandas.util.testing as tm
from pandas import (Categorical, Index, Series, DataFrame, Timestamp,
                    CategoricalIndex, DatetimeIndex, date_range,
                    period_range, PeriodIndex, timedelta_range,
                    TimedeltaIndex)
from pandas.core.dtypes.dtypes import CategoricalDtype

from pandas.compat import PYPY
from pandas.core.categorical import _recode_for_categories


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

    def test_min_max(self):

        # unordered cats have no min/max
        cat = Categorical(["a", "b", "c", "d"], ordered=False)
        pytest.raises(TypeError, lambda: cat.min())
        pytest.raises(TypeError, lambda: cat.max())
        cat = Categorical(["a", "b", "c", "d"], ordered=True)
        _min = cat.min()
        _max = cat.max()
        assert _min == "a"
        assert _max == "d"
        cat = Categorical(["a", "b", "c", "d"],
                          categories=['d', 'c', 'b', 'a'], ordered=True)
        _min = cat.min()
        _max = cat.max()
        assert _min == "d"
        assert _max == "a"
        cat = Categorical([np.nan, "b", "c", np.nan],
                          categories=['d', 'c', 'b', 'a'], ordered=True)
        _min = cat.min()
        _max = cat.max()
        assert np.isnan(_min)
        assert _max == "b"

        _min = cat.min(numeric_only=True)
        assert _min == "c"
        _max = cat.max(numeric_only=True)
        assert _max == "b"

        cat = Categorical([np.nan, 1, 2, np.nan], categories=[5, 4, 3, 2, 1],
                          ordered=True)
        _min = cat.min()
        _max = cat.max()
        assert np.isnan(_min)
        assert _max == 1

        _min = cat.min(numeric_only=True)
        assert _min == 2
        _max = cat.max(numeric_only=True)
        assert _max == 1

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

    def test_mode(self):
        s = Categorical([1, 1, 2, 4, 5, 5, 5], categories=[5, 4, 3, 2, 1],
                        ordered=True)
        res = s.mode()
        exp = Categorical([5], categories=[5, 4, 3, 2, 1], ordered=True)
        tm.assert_categorical_equal(res, exp)
        s = Categorical([1, 1, 1, 4, 5, 5, 5], categories=[5, 4, 3, 2, 1],
                        ordered=True)
        res = s.mode()
        exp = Categorical([5, 1], categories=[5, 4, 3, 2, 1], ordered=True)
        tm.assert_categorical_equal(res, exp)
        s = Categorical([1, 2, 3, 4, 5], categories=[5, 4, 3, 2, 1],
                        ordered=True)
        res = s.mode()
        exp = Categorical([5, 4, 3, 2, 1],
                          categories=[5, 4, 3, 2, 1], ordered=True)
        tm.assert_categorical_equal(res, exp)
        # NaN should not become the mode!
        s = Categorical([np.nan, np.nan, np.nan, 4, 5],
                        categories=[5, 4, 3, 2, 1], ordered=True)
        res = s.mode()
        exp = Categorical([5, 4], categories=[5, 4, 3, 2, 1], ordered=True)
        tm.assert_categorical_equal(res, exp)
        s = Categorical([np.nan, np.nan, np.nan, 4, 5, 4],
                        categories=[5, 4, 3, 2, 1], ordered=True)
        res = s.mode()
        exp = Categorical([4], categories=[5, 4, 3, 2, 1], ordered=True)
        tm.assert_categorical_equal(res, exp)
        s = Categorical([np.nan, np.nan, 4, 5, 4], categories=[5, 4, 3, 2, 1],
                        ordered=True)
        res = s.mode()
        exp = Categorical([4], categories=[5, 4, 3, 2, 1], ordered=True)
        tm.assert_categorical_equal(res, exp)

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

    @pytest.mark.xfail(reason="Imaginary values not supported in Categorical")
    def test_imaginary(self):
        values = [1, 2, 3 + 1j]
        c1 = Categorical(values)
        tm.assert_index_equal(c1.categories, Index(values))
        tm.assert_numpy_array_equal(np.array(c1), np.array(values))


class TestBlockCategoricalAPI(object):

    def test_reshaping(self):

        with catch_warnings(record=True):
            p = tm.makePanel()
            p['str'] = 'foo'
            df = p.to_frame()

        df['category'] = df['str'].astype('category')
        result = df['category'].unstack()

        c = Categorical(['foo'] * len(p.major_axis))
        expected = DataFrame({'A': c.copy(),
                              'B': c.copy(),
                              'C': c.copy(),
                              'D': c.copy()},
                             columns=Index(list('ABCD'), name='minor'),
                             index=p.major_axis.set_names('major'))
        tm.assert_frame_equal(result, expected)

    def test_sideeffects_free(self):
        # Passing a categorical to a Series and then changing values in either
        # the series or the categorical should not change the values in the
        # other one, IF you specify copy!
        cat = Categorical(["a", "b", "c", "a"])
        s = Series(cat, copy=True)
        assert s.cat is not cat
        s.cat.categories = [1, 2, 3]
        exp_s = np.array([1, 2, 3, 1], dtype=np.int64)
        exp_cat = np.array(["a", "b", "c", "a"], dtype=np.object_)
        tm.assert_numpy_array_equal(s.__array__(), exp_s)
        tm.assert_numpy_array_equal(cat.__array__(), exp_cat)

        # setting
        s[0] = 2
        exp_s2 = np.array([2, 2, 3, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(s.__array__(), exp_s2)
        tm.assert_numpy_array_equal(cat.__array__(), exp_cat)

        # however, copy is False by default
        # so this WILL change values
        cat = Categorical(["a", "b", "c", "a"])
        s = Series(cat)
        assert s.values is cat
        s.cat.categories = [1, 2, 3]
        exp_s = np.array([1, 2, 3, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(s.__array__(), exp_s)
        tm.assert_numpy_array_equal(cat.__array__(), exp_s)

        s[0] = 2
        exp_s2 = np.array([2, 2, 3, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(s.__array__(), exp_s2)
        tm.assert_numpy_array_equal(cat.__array__(), exp_s2)

    def test_cat_accessor(self):
        s = Series(Categorical(["a", "b", np.nan, "a"]))
        tm.assert_index_equal(s.cat.categories, Index(["a", "b"]))
        assert not s.cat.ordered, False

        exp = Categorical(["a", "b", np.nan, "a"], categories=["b", "a"])
        s.cat.set_categories(["b", "a"], inplace=True)
        tm.assert_categorical_equal(s.values, exp)

        res = s.cat.set_categories(["b", "a"])
        tm.assert_categorical_equal(res.values, exp)

        s[:] = "a"
        s = s.cat.remove_unused_categories()
        tm.assert_index_equal(s.cat.categories, Index(["a"]))

    def test_sequence_like(self):

        # GH 7839
        # make sure can iterate
        df = DataFrame({"id": [1, 2, 3, 4, 5, 6],
                        "raw_grade": ['a', 'b', 'b', 'a', 'a', 'e']})
        df['grade'] = Categorical(df['raw_grade'])

        # basic sequencing testing
        result = list(df.grade.values)
        expected = np.array(df.grade.values).tolist()
        tm.assert_almost_equal(result, expected)

        # iteration
        for t in df.itertuples(index=False):
            str(t)

        for row, s in df.iterrows():
            str(s)

        for c, col in df.iteritems():
            str(s)

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

    def test_reshape(self):
        cat = Categorical([], categories=["a", "b"])
        tm.assert_produces_warning(FutureWarning, cat.reshape, 0)

        with tm.assert_produces_warning(FutureWarning):
            cat = Categorical([], categories=["a", "b"])
            tm.assert_categorical_equal(cat.reshape(0), cat)

        with tm.assert_produces_warning(FutureWarning):
            cat = Categorical([], categories=["a", "b"])
            tm.assert_categorical_equal(cat.reshape((5, -1)), cat)

        with tm.assert_produces_warning(FutureWarning):
            cat = Categorical(["a", "b"], categories=["a", "b"])
            tm.assert_categorical_equal(cat.reshape(cat.shape), cat)

        with tm.assert_produces_warning(FutureWarning):
            cat = Categorical(["a", "b"], categories=["a", "b"])
            tm.assert_categorical_equal(cat.reshape(cat.size), cat)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            msg = "can only specify one unknown dimension"
            cat = Categorical(["a", "b"], categories=["a", "b"])
            tm.assert_raises_regex(ValueError, msg, cat.reshape, (-2, -1))

    def test_numpy_reshape(self):
        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            cat = Categorical(["a", "b"], categories=["a", "b"])
            tm.assert_categorical_equal(np.reshape(cat, cat.shape), cat)

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            msg = "the 'order' parameter is not supported"
            tm.assert_raises_regex(ValueError, msg, np.reshape,
                                   cat, cat.shape, order='F')

    def test_astype_categorical(self):

        cat = Categorical(['a', 'b', 'b', 'a', 'a', 'c', 'c', 'c'])
        tm.assert_categorical_equal(cat, cat.astype('category'))
        tm.assert_almost_equal(np.array(cat), cat.astype('object'))

        pytest.raises(ValueError, lambda: cat.astype(float))

    def test_cat_tab_completition(self):
        # test the tab completion display
        ok_for_cat = ['categories', 'codes', 'ordered', 'set_categories',
                      'add_categories', 'remove_categories',
                      'rename_categories', 'reorder_categories',
                      'remove_unused_categories', 'as_ordered', 'as_unordered']

        def get_dir(s):
            results = [r for r in s.cat.__dir__() if not r.startswith('_')]
            return list(sorted(set(results)))

        s = Series(list('aabbcde')).astype('category')
        results = get_dir(s)
        tm.assert_almost_equal(results, list(sorted(set(ok_for_cat))))

    @pytest.mark.parametrize(
        "dtype",
        ["int_", "uint", "float_", "unicode_", "timedelta64[h]",
         pytest.param("datetime64[D]",
                      marks=pytest.mark.xfail(reason="issue7996"))]
    )
    @pytest.mark.parametrize("is_ordered", [True, False])
    def test_drop_duplicates_non_bool(self, dtype, is_ordered):
        cat_array = np.array([1, 2, 3, 4, 5], dtype=np.dtype(dtype))

        # Test case 1
        input1 = np.array([1, 2, 3, 3], dtype=np.dtype(dtype))
        tc1 = Series(Categorical(input1, categories=cat_array,
                                 ordered=is_ordered))

        expected = Series([False, False, False, True])
        tm.assert_series_equal(tc1.duplicated(), expected)
        tm.assert_series_equal(tc1.drop_duplicates(), tc1[~expected])
        sc = tc1.copy()
        sc.drop_duplicates(inplace=True)
        tm.assert_series_equal(sc, tc1[~expected])

        expected = Series([False, False, True, False])
        tm.assert_series_equal(tc1.duplicated(keep='last'), expected)
        tm.assert_series_equal(tc1.drop_duplicates(keep='last'),
                               tc1[~expected])
        sc = tc1.copy()
        sc.drop_duplicates(keep='last', inplace=True)
        tm.assert_series_equal(sc, tc1[~expected])

        expected = Series([False, False, True, True])
        tm.assert_series_equal(tc1.duplicated(keep=False), expected)
        tm.assert_series_equal(tc1.drop_duplicates(keep=False), tc1[~expected])
        sc = tc1.copy()
        sc.drop_duplicates(keep=False, inplace=True)
        tm.assert_series_equal(sc, tc1[~expected])

        # Test case 2
        input2 = np.array([1, 2, 3, 5, 3, 2, 4], dtype=np.dtype(dtype))
        tc2 = Series(Categorical(
            input2, categories=cat_array, ordered=is_ordered)
        )

        expected = Series([False, False, False, False, True, True, False])
        tm.assert_series_equal(tc2.duplicated(), expected)
        tm.assert_series_equal(tc2.drop_duplicates(), tc2[~expected])
        sc = tc2.copy()
        sc.drop_duplicates(inplace=True)
        tm.assert_series_equal(sc, tc2[~expected])

        expected = Series([False, True, True, False, False, False, False])
        tm.assert_series_equal(tc2.duplicated(keep='last'), expected)
        tm.assert_series_equal(tc2.drop_duplicates(keep='last'),
                               tc2[~expected])
        sc = tc2.copy()
        sc.drop_duplicates(keep='last', inplace=True)
        tm.assert_series_equal(sc, tc2[~expected])

        expected = Series([False, True, True, False, True, True, False])
        tm.assert_series_equal(tc2.duplicated(keep=False), expected)
        tm.assert_series_equal(tc2.drop_duplicates(keep=False), tc2[~expected])
        sc = tc2.copy()
        sc.drop_duplicates(keep=False, inplace=True)
        tm.assert_series_equal(sc, tc2[~expected])

    @pytest.mark.parametrize("is_ordered", [True, False])
    def test_drop_duplicates_bool(self, is_ordered):
        tc = Series(Categorical([True, False, True, False],
                                categories=[True, False], ordered=is_ordered))

        expected = Series([False, False, True, True])
        tm.assert_series_equal(tc.duplicated(), expected)
        tm.assert_series_equal(tc.drop_duplicates(), tc[~expected])
        sc = tc.copy()
        sc.drop_duplicates(inplace=True)
        tm.assert_series_equal(sc, tc[~expected])

        expected = Series([True, True, False, False])
        tm.assert_series_equal(tc.duplicated(keep='last'), expected)
        tm.assert_series_equal(tc.drop_duplicates(keep='last'), tc[~expected])
        sc = tc.copy()
        sc.drop_duplicates(keep='last', inplace=True)
        tm.assert_series_equal(sc, tc[~expected])

        expected = Series([True, True, True, True])
        tm.assert_series_equal(tc.duplicated(keep=False), expected)
        tm.assert_series_equal(tc.drop_duplicates(keep=False), tc[~expected])
        sc = tc.copy()
        sc.drop_duplicates(keep=False, inplace=True)
        tm.assert_series_equal(sc, tc[~expected])

    def test_reindex(self):

        index = date_range('20000101', periods=3)

        # reindexing to an invalid Categorical
        s = Series(['a', 'b', 'c'], dtype='category')
        result = s.reindex(index)
        expected = Series(Categorical(values=[np.nan, np.nan, np.nan],
                                      categories=['a', 'b', 'c']))
        expected.index = index
        tm.assert_series_equal(result, expected)

        # partial reindexing
        expected = Series(Categorical(values=['b', 'c'], categories=['a', 'b',
                                                                     'c']))
        expected.index = [1, 2]
        result = s.reindex([1, 2])
        tm.assert_series_equal(result, expected)

        expected = Series(Categorical(
            values=['c', np.nan], categories=['a', 'b', 'c']))
        expected.index = [2, 3]
        result = s.reindex([2, 3])
        tm.assert_series_equal(result, expected)

    def test_series_delegations(self):

        # invalid accessor
        pytest.raises(AttributeError, lambda: Series([1, 2, 3]).cat)
        tm.assert_raises_regex(
            AttributeError,
            r"Can only use .cat accessor with a 'category' dtype",
            lambda: Series([1, 2, 3]).cat)
        pytest.raises(AttributeError, lambda: Series(['a', 'b', 'c']).cat)
        pytest.raises(AttributeError, lambda: Series(np.arange(5.)).cat)
        pytest.raises(AttributeError,
                      lambda: Series([Timestamp('20130101')]).cat)

        # Series should delegate calls to '.categories', '.codes', '.ordered'
        # and the methods '.set_categories()' 'drop_unused_categories()' to the
        # categorical
        s = Series(Categorical(["a", "b", "c", "a"], ordered=True))
        exp_categories = Index(["a", "b", "c"])
        tm.assert_index_equal(s.cat.categories, exp_categories)
        s.cat.categories = [1, 2, 3]
        exp_categories = Index([1, 2, 3])
        tm.assert_index_equal(s.cat.categories, exp_categories)

        exp_codes = Series([0, 1, 2, 0], dtype='int8')
        tm.assert_series_equal(s.cat.codes, exp_codes)

        assert s.cat.ordered
        s = s.cat.as_unordered()
        assert not s.cat.ordered
        s.cat.as_ordered(inplace=True)
        assert s.cat.ordered

        # reorder
        s = Series(Categorical(["a", "b", "c", "a"], ordered=True))
        exp_categories = Index(["c", "b", "a"])
        exp_values = np.array(["a", "b", "c", "a"], dtype=np.object_)
        s = s.cat.set_categories(["c", "b", "a"])
        tm.assert_index_equal(s.cat.categories, exp_categories)
        tm.assert_numpy_array_equal(s.values.__array__(), exp_values)
        tm.assert_numpy_array_equal(s.__array__(), exp_values)

        # remove unused categories
        s = Series(Categorical(["a", "b", "b", "a"], categories=["a", "b", "c"
                                                                 ]))
        exp_categories = Index(["a", "b"])
        exp_values = np.array(["a", "b", "b", "a"], dtype=np.object_)
        s = s.cat.remove_unused_categories()
        tm.assert_index_equal(s.cat.categories, exp_categories)
        tm.assert_numpy_array_equal(s.values.__array__(), exp_values)
        tm.assert_numpy_array_equal(s.__array__(), exp_values)

        # This method is likely to be confused, so test that it raises an error
        # on wrong inputs:
        def f():
            s.set_categories([4, 3, 2, 1])

        pytest.raises(Exception, f)
        # right: s.cat.set_categories([4,3,2,1])

    def test_series_functions_no_warnings(self):
        df = DataFrame({'value': np.random.randint(0, 100, 20)})
        labels = ["{0} - {1}".format(i, i + 9) for i in range(0, 100, 10)]
        with tm.assert_produces_warning(False):
            df['group'] = pd.cut(df.value, range(0, 105, 10), right=False,
                                 labels=labels)

    def test_assignment_to_dataframe(self):
        # assignment
        df = DataFrame({'value': np.array(
            np.random.randint(0, 10000, 100), dtype='int32')})
        labels = Categorical(["{0} - {1}".format(i, i + 499)
                              for i in range(0, 10000, 500)])

        df = df.sort_values(by=['value'], ascending=True)
        s = pd.cut(df.value, range(0, 10500, 500), right=False, labels=labels)
        d = s.values
        df['D'] = d
        str(df)

        result = df.dtypes
        expected = Series(
            [np.dtype('int32'), CategoricalDtype(categories=labels,
                                                 ordered=False)],
            index=['value', 'D'])
        tm.assert_series_equal(result, expected)

        df['E'] = s
        str(df)

        result = df.dtypes
        expected = Series([np.dtype('int32'),
                           CategoricalDtype(categories=labels, ordered=False),
                           CategoricalDtype(categories=labels, ordered=False)],
                          index=['value', 'D', 'E'])
        tm.assert_series_equal(result, expected)

        result1 = df['D']
        result2 = df['E']
        tm.assert_categorical_equal(result1._data._block.values, d)

        # sorting
        s.name = 'E'
        tm.assert_series_equal(result2.sort_index(), s.sort_index())

        cat = Categorical([1, 2, 3, 10], categories=[1, 2, 3, 4, 10])
        df = DataFrame(Series(cat))

    def test_categorical_frame(self):
        # normal DataFrame
        dt = date_range('2011-01-01 09:00', freq='H', periods=5,
                        tz='US/Eastern')
        p = period_range('2011-01', freq='M', periods=5)
        df = DataFrame({'dt': dt, 'p': p})
        exp = """                         dt       p
0 2011-01-01 09:00:00-05:00 2011-01
1 2011-01-01 10:00:00-05:00 2011-02
2 2011-01-01 11:00:00-05:00 2011-03
3 2011-01-01 12:00:00-05:00 2011-04
4 2011-01-01 13:00:00-05:00 2011-05"""

        df = DataFrame({'dt': Categorical(dt), 'p': Categorical(p)})
        assert repr(df) == exp

    def test_info(self):

        # make sure it works
        n = 2500
        df = DataFrame({'int64': np.random.randint(100, size=n)})
        df['category'] = Series(np.array(list('abcdefghij')).take(
            np.random.randint(0, 10, size=n))).astype('category')
        df.isna()
        buf = compat.StringIO()
        df.info(buf=buf)

        df2 = df[df['category'] == 'd']
        buf = compat.StringIO()
        df2.info(buf=buf)

    def test_min_max(self):
        # unordered cats have no min/max
        cat = Series(Categorical(["a", "b", "c", "d"], ordered=False))
        pytest.raises(TypeError, lambda: cat.min())
        pytest.raises(TypeError, lambda: cat.max())

        cat = Series(Categorical(["a", "b", "c", "d"], ordered=True))
        _min = cat.min()
        _max = cat.max()
        assert _min == "a"
        assert _max == "d"

        cat = Series(Categorical(["a", "b", "c", "d"], categories=[
                     'd', 'c', 'b', 'a'], ordered=True))
        _min = cat.min()
        _max = cat.max()
        assert _min == "d"
        assert _max == "a"

        cat = Series(Categorical(
            [np.nan, "b", "c", np.nan], categories=['d', 'c', 'b', 'a'
                                                    ], ordered=True))
        _min = cat.min()
        _max = cat.max()
        assert np.isnan(_min)
        assert _max == "b"

        cat = Series(Categorical(
            [np.nan, 1, 2, np.nan], categories=[5, 4, 3, 2, 1], ordered=True))
        _min = cat.min()
        _max = cat.max()
        assert np.isnan(_min)
        assert _max == 1

    def test_mode(self):
        s = Series(Categorical([1, 1, 2, 4, 5, 5, 5],
                               categories=[5, 4, 3, 2, 1], ordered=True))
        res = s.mode()
        exp = Series(Categorical([5], categories=[
                     5, 4, 3, 2, 1], ordered=True))
        tm.assert_series_equal(res, exp)
        s = Series(Categorical([1, 1, 1, 4, 5, 5, 5],
                               categories=[5, 4, 3, 2, 1], ordered=True))
        res = s.mode()
        exp = Series(Categorical([5, 1], categories=[
                     5, 4, 3, 2, 1], ordered=True))
        tm.assert_series_equal(res, exp)
        s = Series(Categorical([1, 2, 3, 4, 5], categories=[5, 4, 3, 2, 1],
                               ordered=True))
        res = s.mode()
        exp = Series(Categorical([5, 4, 3, 2, 1], categories=[5, 4, 3, 2, 1],
                                 ordered=True))
        tm.assert_series_equal(res, exp)

    def test_value_counts(self):
        # GH 12835
        cats = Categorical(list('abcccb'), categories=list('cabd'))
        s = Series(cats, name='xxx')
        res = s.value_counts(sort=False)

        exp_index = CategoricalIndex(list('cabd'), categories=cats.categories)
        exp = Series([3, 1, 2, 0], name='xxx', index=exp_index)
        tm.assert_series_equal(res, exp)

        res = s.value_counts(sort=True)

        exp_index = CategoricalIndex(list('cbad'), categories=cats.categories)
        exp = Series([3, 2, 1, 0], name='xxx', index=exp_index)
        tm.assert_series_equal(res, exp)

        # check object dtype handles the Series.name as the same
        # (tested in test_base.py)
        s = Series(["a", "b", "c", "c", "c", "b"], name='xxx')
        res = s.value_counts()
        exp = Series([3, 2, 1], name='xxx', index=["c", "b", "a"])
        tm.assert_series_equal(res, exp)

    def test_groupby(self):

        cats = Categorical(["a", "a", "a", "b", "b", "b", "c", "c", "c"],
                           categories=["a", "b", "c", "d"], ordered=True)
        data = DataFrame({"a": [1, 1, 1, 2, 2, 2, 3, 4, 5], "b": cats})

        exp_index = CategoricalIndex(list('abcd'), name='b', ordered=True)
        expected = DataFrame({'a': [1, 2, 4, np.nan]}, index=exp_index)
        result = data.groupby("b").mean()
        tm.assert_frame_equal(result, expected)

        raw_cat1 = Categorical(["a", "a", "b", "b"],
                               categories=["a", "b", "z"], ordered=True)
        raw_cat2 = Categorical(["c", "d", "c", "d"],
                               categories=["c", "d", "y"], ordered=True)
        df = DataFrame({"A": raw_cat1, "B": raw_cat2, "values": [1, 2, 3, 4]})

        # single grouper
        gb = df.groupby("A")
        exp_idx = CategoricalIndex(['a', 'b', 'z'], name='A', ordered=True)
        expected = DataFrame({'values': Series([3, 7, np.nan], index=exp_idx)})
        result = gb.sum()
        tm.assert_frame_equal(result, expected)

        # multiple groupers
        gb = df.groupby(['A', 'B'])
        exp_index = pd.MultiIndex.from_product(
            [Categorical(["a", "b", "z"], ordered=True),
             Categorical(["c", "d", "y"], ordered=True)],
            names=['A', 'B'])
        expected = DataFrame({'values': [1, 2, np.nan, 3, 4, np.nan,
                                         np.nan, np.nan, np.nan]},
                             index=exp_index)
        result = gb.sum()
        tm.assert_frame_equal(result, expected)

        # multiple groupers with a non-cat
        df = df.copy()
        df['C'] = ['foo', 'bar'] * 2
        gb = df.groupby(['A', 'B', 'C'])
        exp_index = pd.MultiIndex.from_product(
            [Categorical(["a", "b", "z"], ordered=True),
             Categorical(["c", "d", "y"], ordered=True),
             ['foo', 'bar']],
            names=['A', 'B', 'C'])
        expected = DataFrame({'values': Series(
            np.nan, index=exp_index)}).sort_index()
        expected.iloc[[1, 2, 7, 8], 0] = [1, 2, 3, 4]
        result = gb.sum()
        tm.assert_frame_equal(result, expected)

        # GH 8623
        x = DataFrame([[1, 'John P. Doe'], [2, 'Jane Dove'],
                       [1, 'John P. Doe']],
                      columns=['person_id', 'person_name'])
        x['person_name'] = Categorical(x.person_name)

        g = x.groupby(['person_id'])
        result = g.transform(lambda x: x)
        tm.assert_frame_equal(result, x[['person_name']])

        result = x.drop_duplicates('person_name')
        expected = x.iloc[[0, 1]]
        tm.assert_frame_equal(result, expected)

        def f(x):
            return x.drop_duplicates('person_name').iloc[0]

        result = g.apply(f)
        expected = x.iloc[[0, 1]].copy()
        expected.index = Index([1, 2], name='person_id')
        expected['person_name'] = expected['person_name'].astype('object')
        tm.assert_frame_equal(result, expected)

        # GH 9921
        # Monotonic
        df = DataFrame({"a": [5, 15, 25]})
        c = pd.cut(df.a, bins=[0, 10, 20, 30, 40])

        result = df.a.groupby(c).transform(sum)
        tm.assert_series_equal(result, df['a'])

        tm.assert_series_equal(
            df.a.groupby(c).transform(lambda xs: np.sum(xs)), df['a'])
        tm.assert_frame_equal(df.groupby(c).transform(sum), df[['a']])
        tm.assert_frame_equal(
            df.groupby(c).transform(lambda xs: np.max(xs)), df[['a']])

        # Filter
        tm.assert_series_equal(df.a.groupby(c).filter(np.all), df['a'])
        tm.assert_frame_equal(df.groupby(c).filter(np.all), df)

        # Non-monotonic
        df = DataFrame({"a": [5, 15, 25, -5]})
        c = pd.cut(df.a, bins=[-10, 0, 10, 20, 30, 40])

        result = df.a.groupby(c).transform(sum)
        tm.assert_series_equal(result, df['a'])

        tm.assert_series_equal(
            df.a.groupby(c).transform(lambda xs: np.sum(xs)), df['a'])
        tm.assert_frame_equal(df.groupby(c).transform(sum), df[['a']])
        tm.assert_frame_equal(
            df.groupby(c).transform(lambda xs: np.sum(xs)), df[['a']])

        # GH 9603
        df = DataFrame({'a': [1, 0, 0, 0]})
        c = pd.cut(df.a, [0, 1, 2, 3, 4], labels=Categorical(list('abcd')))
        result = df.groupby(c).apply(len)

        exp_index = CategoricalIndex(
            c.values.categories, ordered=c.values.ordered)
        expected = Series([1, 0, 0, 0], index=exp_index)
        expected.index.name = 'a'
        tm.assert_series_equal(result, expected)

    def test_pivot_table(self):

        raw_cat1 = Categorical(["a", "a", "b", "b"],
                               categories=["a", "b", "z"], ordered=True)
        raw_cat2 = Categorical(["c", "d", "c", "d"],
                               categories=["c", "d", "y"], ordered=True)
        df = DataFrame({"A": raw_cat1, "B": raw_cat2, "values": [1, 2, 3, 4]})
        result = pd.pivot_table(df, values='values', index=['A', 'B'])

        exp_index = pd.MultiIndex.from_product(
            [Categorical(["a", "b", "z"], ordered=True),
             Categorical(["c", "d", "y"], ordered=True)],
            names=['A', 'B'])
        expected = DataFrame(
            {'values': [1, 2, np.nan, 3, 4, np.nan, np.nan, np.nan, np.nan]},
            index=exp_index)
        tm.assert_frame_equal(result, expected)

    def test_count(self):

        s = Series(Categorical([np.nan, 1, 2, np.nan],
                               categories=[5, 4, 3, 2, 1], ordered=True))
        result = s.count()
        assert result == 2

    def test_concat_append(self):
        cat = Categorical(["a", "b"], categories=["a", "b"])
        vals = [1, 2]
        df = DataFrame({"cats": cat, "vals": vals})
        cat2 = Categorical(["a", "b", "a", "b"], categories=["a", "b"])
        vals2 = [1, 2, 1, 2]
        exp = DataFrame({"cats": cat2, "vals": vals2},
                        index=Index([0, 1, 0, 1]))

        tm.assert_frame_equal(pd.concat([df, df]), exp)
        tm.assert_frame_equal(df.append(df), exp)

        # GH 13524 can concat different categories
        cat3 = Categorical(["a", "b"], categories=["a", "b", "c"])
        vals3 = [1, 2]
        df_different_categories = DataFrame({"cats": cat3, "vals": vals3})

        res = pd.concat([df, df_different_categories], ignore_index=True)
        exp = DataFrame({"cats": list('abab'), "vals": [1, 2, 1, 2]})
        tm.assert_frame_equal(res, exp)

        res = df.append(df_different_categories, ignore_index=True)
        tm.assert_frame_equal(res, exp)

    def test_concat_append_gh7864(self):
        # GH 7864
        # make sure ordering is preserverd
        df = DataFrame({"id": [1, 2, 3, 4, 5, 6], "raw_grade": list('abbaae')})
        df["grade"] = Categorical(df["raw_grade"])
        df['grade'].cat.set_categories(['e', 'a', 'b'])

        df1 = df[0:3]
        df2 = df[3:]

        tm.assert_index_equal(df['grade'].cat.categories,
                              df1['grade'].cat.categories)
        tm.assert_index_equal(df['grade'].cat.categories,
                              df2['grade'].cat.categories)

        dfx = pd.concat([df1, df2])
        tm.assert_index_equal(df['grade'].cat.categories,
                              dfx['grade'].cat.categories)

        dfa = df1.append(df2)
        tm.assert_index_equal(df['grade'].cat.categories,
                              dfa['grade'].cat.categories)

    def test_concat_preserve(self):

        # GH 8641  series concat not preserving category dtype
        # GH 13524 can concat different categories
        s = Series(list('abc'), dtype='category')
        s2 = Series(list('abd'), dtype='category')

        exp = Series(list('abcabd'))
        res = pd.concat([s, s2], ignore_index=True)
        tm.assert_series_equal(res, exp)

        exp = Series(list('abcabc'), dtype='category')
        res = pd.concat([s, s], ignore_index=True)
        tm.assert_series_equal(res, exp)

        exp = Series(list('abcabc'), index=[0, 1, 2, 0, 1, 2],
                     dtype='category')
        res = pd.concat([s, s])
        tm.assert_series_equal(res, exp)

        a = Series(np.arange(6, dtype='int64'))
        b = Series(list('aabbca'))

        df2 = DataFrame({'A': a,
                         'B': b.astype(CategoricalDtype(list('cab')))})
        res = pd.concat([df2, df2])
        exp = DataFrame(
            {'A': pd.concat([a, a]),
             'B': pd.concat([b, b]).astype(CategoricalDtype(list('cab')))})
        tm.assert_frame_equal(res, exp)

    def test_categorical_index_preserver(self):

        a = Series(np.arange(6, dtype='int64'))
        b = Series(list('aabbca'))

        df2 = DataFrame({'A': a,
                         'B': b.astype(CategoricalDtype(list('cab')))
                         }).set_index('B')
        result = pd.concat([df2, df2])
        expected = DataFrame(
            {'A': pd.concat([a, a]),
             'B': pd.concat([b, b]).astype(CategoricalDtype(list('cab')))
             }).set_index('B')
        tm.assert_frame_equal(result, expected)

        # wrong catgories
        df3 = DataFrame({'A': a, 'B': Categorical(b, categories=list('abe'))
                         }).set_index('B')
        pytest.raises(TypeError, lambda: pd.concat([df2, df3]))

    def test_merge(self):
        # GH 9426

        right = DataFrame({'c': {0: 'a',
                                 1: 'b',
                                 2: 'c',
                                 3: 'd',
                                 4: 'e'},
                           'd': {0: 'null',
                                 1: 'null',
                                 2: 'null',
                                 3: 'null',
                                 4: 'null'}})
        left = DataFrame({'a': {0: 'f',
                                1: 'f',
                                2: 'f',
                                3: 'f',
                                4: 'f'},
                          'b': {0: 'g',
                                1: 'g',
                                2: 'g',
                                3: 'g',
                                4: 'g'}})
        df = pd.merge(left, right, how='left', left_on='b', right_on='c')

        # object-object
        expected = df.copy()

        # object-cat
        # note that we propagate the category
        # because we don't have any matching rows
        cright = right.copy()
        cright['d'] = cright['d'].astype('category')
        result = pd.merge(left, cright, how='left', left_on='b', right_on='c')
        expected['d'] = expected['d'].astype(CategoricalDtype(['null']))
        tm.assert_frame_equal(result, expected)

        # cat-object
        cleft = left.copy()
        cleft['b'] = cleft['b'].astype('category')
        result = pd.merge(cleft, cright, how='left', left_on='b', right_on='c')
        tm.assert_frame_equal(result, expected)

        # cat-cat
        cright = right.copy()
        cright['d'] = cright['d'].astype('category')
        cleft = left.copy()
        cleft['b'] = cleft['b'].astype('category')
        result = pd.merge(cleft, cright, how='left', left_on='b', right_on='c')
        tm.assert_frame_equal(result, expected)

    def test_to_records(self):

        # GH8626

        # dict creation
        df = DataFrame({'A': list('abc')}, dtype='category')
        expected = Series(list('abc'), dtype='category', name='A')
        tm.assert_series_equal(df['A'], expected)

        # list-like creation
        df = DataFrame(list('abc'), dtype='category')
        expected = Series(list('abc'), dtype='category', name=0)
        tm.assert_series_equal(df[0], expected)

        # to record array
        # this coerces
        result = df.to_records()
        expected = np.rec.array([(0, 'a'), (1, 'b'), (2, 'c')],
                                dtype=[('index', '=i8'), ('0', 'O')])
        tm.assert_almost_equal(result, expected)

    def test_cat_accessor_api(self):
        # GH 9322
        from pandas.core.categorical import CategoricalAccessor
        assert Series.cat is CategoricalAccessor
        s = Series(list('aabbcde')).astype('category')
        assert isinstance(s.cat, CategoricalAccessor)

        invalid = Series([1])
        with tm.assert_raises_regex(AttributeError,
                                    "only use .cat accessor"):
            invalid.cat
        assert not hasattr(invalid, 'cat')

    def test_cat_accessor_no_new_attributes(self):
        # https://github.com/pandas-dev/pandas/issues/10673
        c = Series(list('aabbcde')).astype('category')
        with tm.assert_raises_regex(AttributeError,
                                    "You cannot add any new attribute"):
            c.cat.xlabel = "a"

    def test_str_accessor_api_for_categorical(self):
        # https://github.com/pandas-dev/pandas/issues/10661
        from pandas.core.strings import StringMethods
        s = Series(list('aabb'))
        s = s + " " + s
        c = s.astype('category')
        assert isinstance(c.str, StringMethods)

        # str functions, which need special arguments
        special_func_defs = [
            ('cat', (list("zyxw"),), {"sep": ","}),
            ('center', (10,), {}),
            ('contains', ("a",), {}),
            ('count', ("a",), {}),
            ('decode', ("UTF-8",), {}),
            ('encode', ("UTF-8",), {}),
            ('endswith', ("a",), {}),
            ('extract', ("([a-z]*) ",), {"expand": False}),
            ('extract', ("([a-z]*) ",), {"expand": True}),
            ('extractall', ("([a-z]*) ",), {}),
            ('find', ("a",), {}),
            ('findall', ("a",), {}),
            ('index', (" ",), {}),
            ('ljust', (10,), {}),
            ('match', ("a"), {}),  # deprecated...
            ('normalize', ("NFC",), {}),
            ('pad', (10,), {}),
            ('partition', (" ",), {"expand": False}),  # not default
            ('partition', (" ",), {"expand": True}),  # default
            ('repeat', (3,), {}),
            ('replace', ("a", "z"), {}),
            ('rfind', ("a",), {}),
            ('rindex', (" ",), {}),
            ('rjust', (10,), {}),
            ('rpartition', (" ",), {"expand": False}),  # not default
            ('rpartition', (" ",), {"expand": True}),  # default
            ('slice', (0, 1), {}),
            ('slice_replace', (0, 1, "z"), {}),
            ('split', (" ",), {"expand": False}),  # default
            ('split', (" ",), {"expand": True}),  # not default
            ('startswith', ("a",), {}),
            ('wrap', (2,), {}),
            ('zfill', (10,), {})
        ]
        _special_func_names = [f[0] for f in special_func_defs]

        # * get, join: they need a individual elements of type lists, but
        #   we can't make a categorical with lists as individual categories.
        #   -> `s.str.split(" ").astype("category")` will error!
        # * `translate` has different interfaces for py2 vs. py3
        _ignore_names = ["get", "join", "translate"]

        str_func_names = [f for f in dir(s.str) if not (
            f.startswith("_") or
            f in _special_func_names or
            f in _ignore_names)]

        func_defs = [(f, (), {}) for f in str_func_names]
        func_defs.extend(special_func_defs)

        for func, args, kwargs in func_defs:
            res = getattr(c.str, func)(*args, **kwargs)
            exp = getattr(s.str, func)(*args, **kwargs)

            if isinstance(res, DataFrame):
                tm.assert_frame_equal(res, exp)
            else:
                tm.assert_series_equal(res, exp)

        invalid = Series([1, 2, 3]).astype('category')
        with tm.assert_raises_regex(AttributeError,
                                    "Can only use .str "
                                    "accessor with string"):
            invalid.str
        assert not hasattr(invalid, 'str')

    def test_dt_accessor_api_for_categorical(self):
        # https://github.com/pandas-dev/pandas/issues/10661
        from pandas.core.indexes.accessors import Properties

        s_dr = Series(date_range('1/1/2015', periods=5, tz="MET"))
        c_dr = s_dr.astype("category")

        s_pr = Series(period_range('1/1/2015', freq='D', periods=5))
        c_pr = s_pr.astype("category")

        s_tdr = Series(timedelta_range('1 days', '10 days'))
        c_tdr = s_tdr.astype("category")

        # only testing field (like .day)
        # and bool (is_month_start)
        get_ops = lambda x: x._datetimelike_ops

        test_data = [
            ("Datetime", get_ops(DatetimeIndex), s_dr, c_dr),
            ("Period", get_ops(PeriodIndex), s_pr, c_pr),
            ("Timedelta", get_ops(TimedeltaIndex), s_tdr, c_tdr)]

        assert isinstance(c_dr.dt, Properties)

        special_func_defs = [
            ('strftime', ("%Y-%m-%d",), {}),
            ('tz_convert', ("EST",), {}),
            ('round', ("D",), {}),
            ('floor', ("D",), {}),
            ('ceil', ("D",), {}),
            ('asfreq', ("D",), {}),
            # ('tz_localize', ("UTC",), {}),
        ]
        _special_func_names = [f[0] for f in special_func_defs]

        # the series is already localized
        _ignore_names = ['tz_localize', 'components']

        for name, attr_names, s, c in test_data:
            func_names = [f
                          for f in dir(s.dt)
                          if not (f.startswith("_") or f in attr_names or f in
                                  _special_func_names or f in _ignore_names)]

            func_defs = [(f, (), {}) for f in func_names]
            for f_def in special_func_defs:
                if f_def[0] in dir(s.dt):
                    func_defs.append(f_def)

            for func, args, kwargs in func_defs:
                res = getattr(c.dt, func)(*args, **kwargs)
                exp = getattr(s.dt, func)(*args, **kwargs)

                if isinstance(res, DataFrame):
                    tm.assert_frame_equal(res, exp)
                elif isinstance(res, Series):
                    tm.assert_series_equal(res, exp)
                else:
                    tm.assert_almost_equal(res, exp)

            for attr in attr_names:
                try:
                    res = getattr(c.dt, attr)
                    exp = getattr(s.dt, attr)
                except Exception as e:
                    print(name, attr)
                    raise e

            if isinstance(res, DataFrame):
                tm.assert_frame_equal(res, exp)
            elif isinstance(res, Series):
                tm.assert_series_equal(res, exp)
            else:
                tm.assert_almost_equal(res, exp)

        invalid = Series([1, 2, 3]).astype('category')
        with tm.assert_raises_regex(
                AttributeError, "Can only use .dt accessor with datetimelike"):
            invalid.dt
        assert not hasattr(invalid, 'str')

    def test_concat_categorical(self):
        # See GH 10177
        df1 = DataFrame(np.arange(18, dtype='int64').reshape(6, 3),
                        columns=["a", "b", "c"])

        df2 = DataFrame(np.arange(14, dtype='int64').reshape(7, 2),
                        columns=["a", "c"])

        cat_values = ["one", "one", "two", "one", "two", "two", "one"]
        df2['h'] = Series(Categorical(cat_values))

        res = pd.concat((df1, df2), axis=0, ignore_index=True)
        exp = DataFrame({'a': [0, 3, 6, 9, 12, 15, 0, 2, 4, 6, 8, 10, 12],
                         'b': [1, 4, 7, 10, 13, 16, np.nan, np.nan, np.nan,
                               np.nan, np.nan, np.nan, np.nan],
                         'c': [2, 5, 8, 11, 14, 17, 1, 3, 5, 7, 9, 11, 13],
                         'h': [None] * 6 + cat_values})
        tm.assert_frame_equal(res, exp)
