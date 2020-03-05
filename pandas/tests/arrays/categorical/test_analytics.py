import sys

import numpy as np
import pytest

from pandas.compat import PYPY

from pandas import Categorical, Index, NaT, Series, date_range
import pandas._testing as tm
from pandas.api.types import is_scalar


class TestCategoricalAnalytics:
    @pytest.mark.parametrize("aggregation", ["min", "max"])
    def test_min_max_not_ordered_raises(self, aggregation):
        # unordered cats have no min/max
        cat = Categorical(["a", "b", "c", "d"], ordered=False)
        msg = f"Categorical is not ordered for operation {aggregation}"
        agg_func = getattr(cat, aggregation)

        with pytest.raises(TypeError, match=msg):
            agg_func()

    def test_min_max_ordered(self):
        cat = Categorical(["a", "b", "c", "d"], ordered=True)
        _min = cat.min()
        _max = cat.max()
        assert _min == "a"
        assert _max == "d"

        cat = Categorical(
            ["a", "b", "c", "d"], categories=["d", "c", "b", "a"], ordered=True
        )
        _min = cat.min()
        _max = cat.max()
        assert _min == "d"
        assert _max == "a"

    @pytest.mark.parametrize(
        "categories,expected",
        [
            (list("ABC"), np.NaN),
            ([1, 2, 3], np.NaN),
            pytest.param(
                Series(date_range("2020-01-01", periods=3), dtype="category"),
                NaT,
                marks=pytest.mark.xfail(
                    reason="https://github.com/pandas-dev/pandas/issues/29962"
                ),
            ),
        ],
    )
    @pytest.mark.parametrize("aggregation", ["min", "max"])
    def test_min_max_ordered_empty(self, categories, expected, aggregation):
        # GH 30227
        cat = Categorical([], categories=list("ABC"), ordered=True)

        agg_func = getattr(cat, aggregation)
        result = agg_func()
        assert result is expected

    @pytest.mark.parametrize("skipna", [True, False])
    def test_min_max_with_nan(self, skipna):
        # GH 25303
        cat = Categorical(
            [np.nan, "b", "c", np.nan], categories=["d", "c", "b", "a"], ordered=True
        )
        _min = cat.min(skipna=skipna)
        _max = cat.max(skipna=skipna)

        if skipna is False:
            assert np.isnan(_min)
            assert np.isnan(_max)
        else:
            assert _min == "c"
            assert _max == "b"

        cat = Categorical(
            [np.nan, 1, 2, np.nan], categories=[5, 4, 3, 2, 1], ordered=True
        )
        _min = cat.min(skipna=skipna)
        _max = cat.max(skipna=skipna)

        if skipna is False:
            assert np.isnan(_min)
            assert np.isnan(_max)
        else:
            assert _min == 2
            assert _max == 1

    @pytest.mark.parametrize("method", ["min", "max"])
    def test_deprecate_numeric_only_min_max(self, method):
        # GH 25303
        cat = Categorical(
            [np.nan, 1, 2, np.nan], categories=[5, 4, 3, 2, 1], ordered=True
        )
        with tm.assert_produces_warning(expected_warning=FutureWarning):
            getattr(cat, method)(numeric_only=True)

    @pytest.mark.parametrize(
        "values,categories,exp_mode",
        [
            ([1, 1, 2, 4, 5, 5, 5], [5, 4, 3, 2, 1], [5]),
            ([1, 1, 1, 4, 5, 5, 5], [5, 4, 3, 2, 1], [5, 1]),
            ([1, 2, 3, 4, 5], [5, 4, 3, 2, 1], [5, 4, 3, 2, 1]),
            ([np.nan, np.nan, np.nan, 4, 5], [5, 4, 3, 2, 1], [5, 4]),
            ([np.nan, np.nan, np.nan, 4, 5, 4], [5, 4, 3, 2, 1], [4]),
            ([np.nan, np.nan, 4, 5, 4], [5, 4, 3, 2, 1], [4]),
        ],
    )
    def test_mode(self, values, categories, exp_mode):
        s = Categorical(values, categories=categories, ordered=True)
        res = s.mode()
        exp = Categorical(exp_mode, categories=categories, ordered=True)
        tm.assert_categorical_equal(res, exp)

    def test_searchsorted(self, ordered_fixture):
        # https://github.com/pandas-dev/pandas/issues/8420
        # https://github.com/pandas-dev/pandas/issues/14522

        cat = Categorical(
            ["cheese", "milk", "apple", "bread", "bread"],
            categories=["cheese", "milk", "apple", "bread"],
            ordered=ordered_fixture,
        )
        ser = Series(cat)

        # Searching for single item argument, side='left' (default)
        res_cat = cat.searchsorted("apple")
        assert res_cat == 2
        assert is_scalar(res_cat)

        res_ser = ser.searchsorted("apple")
        assert res_ser == 2
        assert is_scalar(res_ser)

        # Searching for single item array, side='left' (default)
        res_cat = cat.searchsorted(["bread"])
        res_ser = ser.searchsorted(["bread"])
        exp = np.array([3], dtype=np.intp)
        tm.assert_numpy_array_equal(res_cat, exp)
        tm.assert_numpy_array_equal(res_ser, exp)

        # Searching for several items array, side='right'
        res_cat = cat.searchsorted(["apple", "bread"], side="right")
        res_ser = ser.searchsorted(["apple", "bread"], side="right")
        exp = np.array([3, 5], dtype=np.intp)
        tm.assert_numpy_array_equal(res_cat, exp)
        tm.assert_numpy_array_equal(res_ser, exp)

        # Searching for a single value that is not from the Categorical
        with pytest.raises(KeyError, match="cucumber"):
            cat.searchsorted("cucumber")
        with pytest.raises(KeyError, match="cucumber"):
            ser.searchsorted("cucumber")

        # Searching for multiple values one of each is not from the Categorical
        with pytest.raises(KeyError, match="cucumber"):
            cat.searchsorted(["bread", "cucumber"])
        with pytest.raises(KeyError, match="cucumber"):
            ser.searchsorted(["bread", "cucumber"])

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

        cat = Categorical(["c", "a", "b", "a", "a"], categories=["a", "b", "c"])
        exp = Index(["c", "a", "b"])
        res = cat.unique()
        tm.assert_index_equal(res.categories, exp)
        exp_cat = Categorical(exp, categories=["c", "a", "b"])
        tm.assert_categorical_equal(res, exp_cat)

        # nan must be removed
        cat = Categorical(["b", np.nan, "b", np.nan, "a"], categories=["a", "b", "c"])
        res = cat.unique()
        exp = Index(["b", "a"])
        tm.assert_index_equal(res.categories, exp)
        exp_cat = Categorical(["b", np.nan, "a"], categories=["b", "a"])
        tm.assert_categorical_equal(res, exp_cat)

    def test_unique_ordered(self):
        # keep categories order when ordered=True
        cat = Categorical(["b", "a", "b"], categories=["a", "b"], ordered=True)
        res = cat.unique()
        exp_cat = Categorical(["b", "a"], categories=["a", "b"], ordered=True)
        tm.assert_categorical_equal(res, exp_cat)

        cat = Categorical(
            ["c", "b", "a", "a"], categories=["a", "b", "c"], ordered=True
        )
        res = cat.unique()
        exp_cat = Categorical(["c", "b", "a"], categories=["a", "b", "c"], ordered=True)
        tm.assert_categorical_equal(res, exp_cat)

        cat = Categorical(["b", "a", "a"], categories=["a", "b", "c"], ordered=True)
        res = cat.unique()
        exp_cat = Categorical(["b", "a"], categories=["a", "b"], ordered=True)
        tm.assert_categorical_equal(res, exp_cat)

        cat = Categorical(
            ["b", "b", np.nan, "a"], categories=["a", "b", "c"], ordered=True
        )
        res = cat.unique()
        exp_cat = Categorical(["b", np.nan, "a"], categories=["a", "b"], ordered=True)
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
        cat = Categorical(["a", "b", "c", "d", "a"])

        # shift forward
        sp1 = cat.shift(1)
        xp1 = Categorical([np.nan, "a", "b", "c", "d"])
        tm.assert_categorical_equal(sp1, xp1)
        tm.assert_categorical_equal(cat[:-1], sp1[1:])

        # shift back
        sn2 = cat.shift(-2)
        xp2 = Categorical(
            ["c", "d", "a", np.nan, np.nan], categories=["a", "b", "c", "d"]
        )
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

        cat = Categorical(["foo", "foo", "bar"])
        assert cat.memory_usage(deep=True) > cat.nbytes

        if not PYPY:
            # sys.getsizeof will call the .memory_usage with
            # deep=True, and add on some GC overhead
            diff = cat.memory_usage(deep=True) - sys.getsizeof(cat)
            assert abs(diff) < 100

    def test_map(self):
        c = Categorical(list("ABABC"), categories=list("CBA"), ordered=True)
        result = c.map(lambda x: x.lower())
        exp = Categorical(list("ababc"), categories=list("cba"), ordered=True)
        tm.assert_categorical_equal(result, exp)

        c = Categorical(list("ABABC"), categories=list("ABC"), ordered=False)
        result = c.map(lambda x: x.lower())
        exp = Categorical(list("ababc"), categories=list("abc"), ordered=False)
        tm.assert_categorical_equal(result, exp)

        result = c.map(lambda x: 1)
        # GH 12766: Return an index not an array
        tm.assert_index_equal(result, Index(np.array([1] * 5, dtype=np.int64)))

    @pytest.mark.parametrize("value", [1, "True", [1, 2, 3], 5.0])
    def test_validate_inplace_raises(self, value):
        cat = Categorical(["A", "B", "B", "C", "A"])
        msg = (
            'For argument "inplace" expected type bool, '
            f"received type {type(value).__name__}"
        )
        with pytest.raises(ValueError, match=msg):
            cat.set_ordered(value=True, inplace=value)

        with pytest.raises(ValueError, match=msg):
            cat.as_ordered(inplace=value)

        with pytest.raises(ValueError, match=msg):
            cat.as_unordered(inplace=value)

        with pytest.raises(ValueError, match=msg):
            cat.set_categories(["X", "Y", "Z"], rename=True, inplace=value)

        with pytest.raises(ValueError, match=msg):
            cat.rename_categories(["X", "Y", "Z"], inplace=value)

        with pytest.raises(ValueError, match=msg):
            cat.reorder_categories(["X", "Y", "Z"], ordered=True, inplace=value)

        with pytest.raises(ValueError, match=msg):
            cat.add_categories(new_categories=["D", "E", "F"], inplace=value)

        with pytest.raises(ValueError, match=msg):
            cat.remove_categories(removals=["D", "E", "F"], inplace=value)

        with pytest.raises(ValueError, match=msg):
            cat.remove_unused_categories(inplace=value)

        with pytest.raises(ValueError, match=msg):
            cat.sort_values(inplace=value)

    def test_isna(self):
        exp = np.array([False, False, True])
        c = Categorical(["a", "b", np.nan])
        res = c.isna()

        tm.assert_numpy_array_equal(res, exp)
