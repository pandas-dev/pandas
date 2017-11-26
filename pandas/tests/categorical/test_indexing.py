# -*- coding: utf-8 -*-

import pytest

import numpy as np

import pandas as pd
import pandas.compat as compat
import pandas.util.testing as tm
from pandas.core.dtypes.dtypes import CategoricalDtype
from pandas import (Categorical, Index, Series, DataFrame, PeriodIndex,
                    Interval)
from pandas.core.dtypes.common import (
    is_categorical_dtype)


class TestCategoricalIndexing(object):

    def test_getitem_listlike(self):

        # GH 9469
        # properly coerce the input indexers
        np.random.seed(1)
        c = Categorical(np.random.randint(0, 5, size=150000).astype(np.int8))
        result = c.codes[np.array([100000]).astype(np.int64)]
        expected = c[np.array([100000]).astype(np.int64)].codes
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize(
        "method",
        [
            lambda x: x.cat.set_categories([1, 2, 3]),
            lambda x: x.cat.reorder_categories([2, 3, 1], ordered=True),
            lambda x: x.cat.rename_categories([1, 2, 3]),
            lambda x: x.cat.remove_unused_categories(),
            lambda x: x.cat.remove_categories([2]),
            lambda x: x.cat.add_categories([4]),
            lambda x: x.cat.as_ordered(),
            lambda x: x.cat.as_unordered(),
        ])
    def test_getname_categorical_accessor(self, method):
        # GH 17509
        s = Series([1, 2, 3], name='A').astype('category')
        expected = 'A'
        result = method(s).name
        assert result == expected

    def test_getitem_category_type(self):
        # GH 14580
        # test iloc() on Series with Categorical data

        s = Series([1, 2, 3]).astype('category')

        # get slice
        result = s.iloc[0:2]
        expected = Series([1, 2]).astype(CategoricalDtype([1, 2, 3]))
        tm.assert_series_equal(result, expected)

        # get list of indexes
        result = s.iloc[[0, 1]]
        expected = Series([1, 2]).astype(CategoricalDtype([1, 2, 3]))
        tm.assert_series_equal(result, expected)

        # get boolean array
        result = s.iloc[[True, False, False]]
        expected = Series([1]).astype(CategoricalDtype([1, 2, 3]))
        tm.assert_series_equal(result, expected)

    def test_setitem_listlike(self):

        # GH 9469
        # properly coerce the input indexers
        np.random.seed(1)
        c = Categorical(np.random.randint(0, 5, size=150000).astype(
            np.int8)).add_categories([-1000])
        indexer = np.array([100000]).astype(np.int64)
        c[indexer] = -1000

        # we are asserting the code result here
        # which maps to the -1000 category
        result = c.codes[np.array([100000]).astype(np.int64)]
        tm.assert_numpy_array_equal(result, np.array([5], dtype='int8'))

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

    def test_slicing_directly(self):
        cat = Categorical(["a", "b", "c", "d", "a", "b", "c"])
        sliced = cat[3]
        assert sliced == "d"
        sliced = cat[3:5]
        expected = Categorical(["d", "a"], categories=['a', 'b', 'c', 'd'])
        tm.assert_numpy_array_equal(sliced._codes, expected._codes)
        tm.assert_index_equal(sliced.categories, expected.categories)

    def test_periodindex(self):
        idx1 = PeriodIndex(['2014-01', '2014-01', '2014-02', '2014-02',
                            '2014-03', '2014-03'], freq='M')

        cat1 = Categorical(idx1)
        str(cat1)
        exp_arr = np.array([0, 0, 1, 1, 2, 2], dtype=np.int8)
        exp_idx = PeriodIndex(['2014-01', '2014-02', '2014-03'], freq='M')
        tm.assert_numpy_array_equal(cat1._codes, exp_arr)
        tm.assert_index_equal(cat1.categories, exp_idx)

        idx2 = PeriodIndex(['2014-03', '2014-03', '2014-02', '2014-01',
                            '2014-03', '2014-01'], freq='M')
        cat2 = Categorical(idx2, ordered=True)
        str(cat2)
        exp_arr = np.array([2, 2, 1, 0, 2, 0], dtype=np.int8)
        exp_idx2 = PeriodIndex(['2014-01', '2014-02', '2014-03'], freq='M')
        tm.assert_numpy_array_equal(cat2._codes, exp_arr)
        tm.assert_index_equal(cat2.categories, exp_idx2)

        idx3 = PeriodIndex(['2013-12', '2013-11', '2013-10', '2013-09',
                            '2013-08', '2013-07', '2013-05'], freq='M')
        cat3 = Categorical(idx3, ordered=True)
        exp_arr = np.array([6, 5, 4, 3, 2, 1, 0], dtype=np.int8)
        exp_idx = PeriodIndex(['2013-05', '2013-07', '2013-08', '2013-09',
                               '2013-10', '2013-11', '2013-12'], freq='M')
        tm.assert_numpy_array_equal(cat3._codes, exp_arr)
        tm.assert_index_equal(cat3.categories, exp_idx)

    def test_categories_assigments(self):
        s = Categorical(["a", "b", "c", "a"])
        exp = np.array([1, 2, 3, 1], dtype=np.int64)
        s.categories = [1, 2, 3]
        tm.assert_numpy_array_equal(s.__array__(), exp)
        tm.assert_index_equal(s.categories, Index([1, 2, 3]))

        # lengthen
        def f():
            s.categories = [1, 2, 3, 4]

        pytest.raises(ValueError, f)

        # shorten
        def f():
            s.categories = [1, 2]

        pytest.raises(ValueError, f)

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


class TestCategoricalBlockIndexing(object):

    def test_slicing(self):
        cat = Series(Categorical([1, 2, 3, 4]))
        reversed = cat[::-1]
        exp = np.array([4, 3, 2, 1], dtype=np.int64)
        tm.assert_numpy_array_equal(reversed.__array__(), exp)

        df = DataFrame({'value': (np.arange(100) + 1).astype('int64')})
        df['D'] = pd.cut(df.value, bins=[0, 25, 50, 75, 100])

        expected = Series([11, Interval(0, 25)], index=['value', 'D'], name=10)
        result = df.iloc[10]
        tm.assert_series_equal(result, expected)

        expected = DataFrame({'value': np.arange(11, 21).astype('int64')},
                             index=np.arange(10, 20).astype('int64'))
        expected['D'] = pd.cut(expected.value, bins=[0, 25, 50, 75, 100])
        result = df.iloc[10:20]
        tm.assert_frame_equal(result, expected)

        expected = Series([9, Interval(0, 25)], index=['value', 'D'], name=8)
        result = df.loc[8]
        tm.assert_series_equal(result, expected)

    def test_slicing_and_getting_ops(self):

        # systematically test the slicing operations:
        #  for all slicing ops:
        #   - returning a dataframe
        #   - returning a column
        #   - returning a row
        #   - returning a single value

        cats = Categorical(
            ["a", "c", "b", "c", "c", "c", "c"], categories=["a", "b", "c"])
        idx = Index(["h", "i", "j", "k", "l", "m", "n"])
        values = [1, 2, 3, 4, 5, 6, 7]
        df = DataFrame({"cats": cats, "values": values}, index=idx)

        # the expected values
        cats2 = Categorical(["b", "c"], categories=["a", "b", "c"])
        idx2 = Index(["j", "k"])
        values2 = [3, 4]

        # 2:4,: | "j":"k",:
        exp_df = DataFrame({"cats": cats2, "values": values2}, index=idx2)

        # :,"cats" | :,0
        exp_col = Series(cats, index=idx, name='cats')

        # "j",: | 2,:
        exp_row = Series(["b", 3], index=["cats", "values"], dtype="object",
                         name="j")

        # "j","cats | 2,0
        exp_val = "b"

        # iloc
        # frame
        res_df = df.iloc[2:4, :]
        tm.assert_frame_equal(res_df, exp_df)
        assert is_categorical_dtype(res_df["cats"])

        # row
        res_row = df.iloc[2, :]
        tm.assert_series_equal(res_row, exp_row)
        assert isinstance(res_row["cats"], compat.string_types)

        # col
        res_col = df.iloc[:, 0]
        tm.assert_series_equal(res_col, exp_col)
        assert is_categorical_dtype(res_col)

        # single value
        res_val = df.iloc[2, 0]
        assert res_val == exp_val

        # loc
        # frame
        res_df = df.loc["j":"k", :]
        tm.assert_frame_equal(res_df, exp_df)
        assert is_categorical_dtype(res_df["cats"])

        # row
        res_row = df.loc["j", :]
        tm.assert_series_equal(res_row, exp_row)
        assert isinstance(res_row["cats"], compat.string_types)

        # col
        res_col = df.loc[:, "cats"]
        tm.assert_series_equal(res_col, exp_col)
        assert is_categorical_dtype(res_col)

        # single value
        res_val = df.loc["j", "cats"]
        assert res_val == exp_val

        # ix
        # frame
        # res_df = df.loc["j":"k",[0,1]] # doesn't work?
        res_df = df.loc["j":"k", :]
        tm.assert_frame_equal(res_df, exp_df)
        assert is_categorical_dtype(res_df["cats"])

        # row
        res_row = df.loc["j", :]
        tm.assert_series_equal(res_row, exp_row)
        assert isinstance(res_row["cats"], compat.string_types)

        # col
        res_col = df.loc[:, "cats"]
        tm.assert_series_equal(res_col, exp_col)
        assert is_categorical_dtype(res_col)

        # single value
        res_val = df.loc["j", df.columns[0]]
        assert res_val == exp_val

        # iat
        res_val = df.iat[2, 0]
        assert res_val == exp_val

        # at
        res_val = df.at["j", "cats"]
        assert res_val == exp_val

        # fancy indexing
        exp_fancy = df.iloc[[2]]

        res_fancy = df[df["cats"] == "b"]
        tm.assert_frame_equal(res_fancy, exp_fancy)
        res_fancy = df[df["values"] == 3]
        tm.assert_frame_equal(res_fancy, exp_fancy)

        # get_value
        res_val = df.at["j", "cats"]
        assert res_val == exp_val

        # i : int, slice, or sequence of integers
        res_row = df.iloc[2]
        tm.assert_series_equal(res_row, exp_row)
        assert isinstance(res_row["cats"], compat.string_types)

        res_df = df.iloc[slice(2, 4)]
        tm.assert_frame_equal(res_df, exp_df)
        assert is_categorical_dtype(res_df["cats"])

        res_df = df.iloc[[2, 3]]
        tm.assert_frame_equal(res_df, exp_df)
        assert is_categorical_dtype(res_df["cats"])

        res_col = df.iloc[:, 0]
        tm.assert_series_equal(res_col, exp_col)
        assert is_categorical_dtype(res_col)

        res_df = df.iloc[:, slice(0, 2)]
        tm.assert_frame_equal(res_df, df)
        assert is_categorical_dtype(res_df["cats"])

        res_df = df.iloc[:, [0, 1]]
        tm.assert_frame_equal(res_df, df)
        assert is_categorical_dtype(res_df["cats"])

    def test_slicing_doc_examples(self):

        # GH 7918
        cats = Categorical(["a", "b", "b", "b", "c", "c", "c"],
                           categories=["a", "b", "c"])
        idx = Index(["h", "i", "j", "k", "l", "m", "n", ])
        values = [1, 2, 2, 2, 3, 4, 5]
        df = DataFrame({"cats": cats, "values": values}, index=idx)

        result = df.iloc[2:4, :]
        expected = DataFrame(
            {"cats": Categorical(['b', 'b'], categories=['a', 'b', 'c']),
             "values": [2, 2]}, index=['j', 'k'])
        tm.assert_frame_equal(result, expected)

        result = df.iloc[2:4, :].dtypes
        expected = Series(['category', 'int64'], ['cats', 'values'])
        tm.assert_series_equal(result, expected)

        result = df.loc["h":"j", "cats"]
        expected = Series(Categorical(['a', 'b', 'b'],
                                      categories=['a', 'b', 'c']),
                          index=['h', 'i', 'j'], name='cats')
        tm.assert_series_equal(result, expected)

        result = df.loc["h":"j", df.columns[0:1]]
        expected = DataFrame({'cats': Categorical(['a', 'b', 'b'],
                                                  categories=['a', 'b', 'c'])},
                             index=['h', 'i', 'j'])
        tm.assert_frame_equal(result, expected)

    def test_assigning_ops(self):
        # systematically test the assigning operations:
        # for all slicing ops:
        #  for value in categories and value not in categories:

        #   - assign a single value -> exp_single_cats_value

        #   - assign a complete row (mixed values) -> exp_single_row

        # assign multiple rows (mixed values) (-> array) -> exp_multi_row

        # assign a part of a column with dtype == categorical ->
        # exp_parts_cats_col

        # assign a part of a column with dtype != categorical ->
        # exp_parts_cats_col

        cats = Categorical(["a", "a", "a", "a", "a", "a", "a"],
                           categories=["a", "b"])
        idx = Index(["h", "i", "j", "k", "l", "m", "n"])
        values = [1, 1, 1, 1, 1, 1, 1]
        orig = DataFrame({"cats": cats, "values": values}, index=idx)

        # the expected values
        # changed single row
        cats1 = Categorical(["a", "a", "b", "a", "a", "a", "a"],
                            categories=["a", "b"])
        idx1 = Index(["h", "i", "j", "k", "l", "m", "n"])
        values1 = [1, 1, 2, 1, 1, 1, 1]
        exp_single_row = DataFrame({"cats": cats1,
                                    "values": values1}, index=idx1)

        # changed multiple rows
        cats2 = Categorical(["a", "a", "b", "b", "a", "a", "a"],
                            categories=["a", "b"])
        idx2 = Index(["h", "i", "j", "k", "l", "m", "n"])
        values2 = [1, 1, 2, 2, 1, 1, 1]
        exp_multi_row = DataFrame({"cats": cats2,
                                   "values": values2}, index=idx2)

        # changed part of the cats column
        cats3 = Categorical(
            ["a", "a", "b", "b", "a", "a", "a"], categories=["a", "b"])
        idx3 = Index(["h", "i", "j", "k", "l", "m", "n"])
        values3 = [1, 1, 1, 1, 1, 1, 1]
        exp_parts_cats_col = DataFrame({"cats": cats3,
                                        "values": values3}, index=idx3)

        # changed single value in cats col
        cats4 = Categorical(
            ["a", "a", "b", "a", "a", "a", "a"], categories=["a", "b"])
        idx4 = Index(["h", "i", "j", "k", "l", "m", "n"])
        values4 = [1, 1, 1, 1, 1, 1, 1]
        exp_single_cats_value = DataFrame({"cats": cats4,
                                           "values": values4}, index=idx4)

        #  iloc
        # ###############
        #   - assign a single value -> exp_single_cats_value
        df = orig.copy()
        df.iloc[2, 0] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        df = orig.copy()
        df.iloc[df.index == "j", 0] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        #   - assign a single value not in the current categories set
        def f():
            df = orig.copy()
            df.iloc[2, 0] = "c"

        pytest.raises(ValueError, f)

        #   - assign a complete row (mixed values) -> exp_single_row
        df = orig.copy()
        df.iloc[2, :] = ["b", 2]
        tm.assert_frame_equal(df, exp_single_row)

        #   - assign a complete row (mixed values) not in categories set
        def f():
            df = orig.copy()
            df.iloc[2, :] = ["c", 2]

        pytest.raises(ValueError, f)

        #   - assign multiple rows (mixed values) -> exp_multi_row
        df = orig.copy()
        df.iloc[2:4, :] = [["b", 2], ["b", 2]]
        tm.assert_frame_equal(df, exp_multi_row)

        def f():
            df = orig.copy()
            df.iloc[2:4, :] = [["c", 2], ["c", 2]]

        pytest.raises(ValueError, f)

        # assign a part of a column with dtype == categorical ->
        # exp_parts_cats_col
        df = orig.copy()
        df.iloc[2:4, 0] = Categorical(["b", "b"], categories=["a", "b"])
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with pytest.raises(ValueError):
            # different categories -> not sure if this should fail or pass
            df = orig.copy()
            df.iloc[2:4, 0] = Categorical(list('bb'), categories=list('abc'))

        with pytest.raises(ValueError):
            # different values
            df = orig.copy()
            df.iloc[2:4, 0] = Categorical(list('cc'), categories=list('abc'))

        # assign a part of a column with dtype != categorical ->
        # exp_parts_cats_col
        df = orig.copy()
        df.iloc[2:4, 0] = ["b", "b"]
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with pytest.raises(ValueError):
            df.iloc[2:4, 0] = ["c", "c"]

        #  loc
        # ##############
        #   - assign a single value -> exp_single_cats_value
        df = orig.copy()
        df.loc["j", "cats"] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        df = orig.copy()
        df.loc[df.index == "j", "cats"] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        #   - assign a single value not in the current categories set
        def f():
            df = orig.copy()
            df.loc["j", "cats"] = "c"

        pytest.raises(ValueError, f)

        #   - assign a complete row (mixed values) -> exp_single_row
        df = orig.copy()
        df.loc["j", :] = ["b", 2]
        tm.assert_frame_equal(df, exp_single_row)

        #   - assign a complete row (mixed values) not in categories set
        def f():
            df = orig.copy()
            df.loc["j", :] = ["c", 2]

        pytest.raises(ValueError, f)

        #   - assign multiple rows (mixed values) -> exp_multi_row
        df = orig.copy()
        df.loc["j":"k", :] = [["b", 2], ["b", 2]]
        tm.assert_frame_equal(df, exp_multi_row)

        def f():
            df = orig.copy()
            df.loc["j":"k", :] = [["c", 2], ["c", 2]]

        pytest.raises(ValueError, f)

        # assign a part of a column with dtype == categorical ->
        # exp_parts_cats_col
        df = orig.copy()
        df.loc["j":"k", "cats"] = Categorical(
            ["b", "b"], categories=["a", "b"])
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with pytest.raises(ValueError):
            # different categories -> not sure if this should fail or pass
            df = orig.copy()
            df.loc["j":"k", "cats"] = Categorical(
                ["b", "b"], categories=["a", "b", "c"])

        with pytest.raises(ValueError):
            # different values
            df = orig.copy()
            df.loc["j":"k", "cats"] = Categorical(
                ["c", "c"], categories=["a", "b", "c"])

        # assign a part of a column with dtype != categorical ->
        # exp_parts_cats_col
        df = orig.copy()
        df.loc["j":"k", "cats"] = ["b", "b"]
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with pytest.raises(ValueError):
            df.loc["j":"k", "cats"] = ["c", "c"]

        #  loc
        # ##############
        #   - assign a single value -> exp_single_cats_value
        df = orig.copy()
        df.loc["j", df.columns[0]] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        df = orig.copy()
        df.loc[df.index == "j", df.columns[0]] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        #   - assign a single value not in the current categories set
        def f():
            df = orig.copy()
            df.loc["j", df.columns[0]] = "c"

        pytest.raises(ValueError, f)

        #   - assign a complete row (mixed values) -> exp_single_row
        df = orig.copy()
        df.loc["j", :] = ["b", 2]
        tm.assert_frame_equal(df, exp_single_row)

        #   - assign a complete row (mixed values) not in categories set
        def f():
            df = orig.copy()
            df.loc["j", :] = ["c", 2]

        pytest.raises(ValueError, f)

        #   - assign multiple rows (mixed values) -> exp_multi_row
        df = orig.copy()
        df.loc["j":"k", :] = [["b", 2], ["b", 2]]
        tm.assert_frame_equal(df, exp_multi_row)

        def f():
            df = orig.copy()
            df.loc["j":"k", :] = [["c", 2], ["c", 2]]

        pytest.raises(ValueError, f)

        # assign a part of a column with dtype == categorical ->
        # exp_parts_cats_col
        df = orig.copy()
        df.loc["j":"k", df.columns[0]] = Categorical(
            ["b", "b"], categories=["a", "b"])
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with pytest.raises(ValueError):
            # different categories -> not sure if this should fail or pass
            df = orig.copy()
            df.loc["j":"k", df.columns[0]] = Categorical(
                ["b", "b"], categories=["a", "b", "c"])

        with pytest.raises(ValueError):
            # different values
            df = orig.copy()
            df.loc["j":"k", df.columns[0]] = Categorical(
                ["c", "c"], categories=["a", "b", "c"])

        # assign a part of a column with dtype != categorical ->
        # exp_parts_cats_col
        df = orig.copy()
        df.loc["j":"k", df.columns[0]] = ["b", "b"]
        tm.assert_frame_equal(df, exp_parts_cats_col)

        with pytest.raises(ValueError):
            df.loc["j":"k", df.columns[0]] = ["c", "c"]

        # iat
        df = orig.copy()
        df.iat[2, 0] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        #   - assign a single value not in the current categories set
        def f():
            df = orig.copy()
            df.iat[2, 0] = "c"

        pytest.raises(ValueError, f)

        # at
        #   - assign a single value -> exp_single_cats_value
        df = orig.copy()
        df.at["j", "cats"] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        #   - assign a single value not in the current categories set
        def f():
            df = orig.copy()
            df.at["j", "cats"] = "c"

        pytest.raises(ValueError, f)

        # fancy indexing
        catsf = Categorical(["a", "a", "c", "c", "a", "a", "a"],
                            categories=["a", "b", "c"])
        idxf = Index(["h", "i", "j", "k", "l", "m", "n"])
        valuesf = [1, 1, 3, 3, 1, 1, 1]
        df = DataFrame({"cats": catsf, "values": valuesf}, index=idxf)

        exp_fancy = exp_multi_row.copy()
        exp_fancy["cats"].cat.set_categories(["a", "b", "c"], inplace=True)

        df[df["cats"] == "c"] = ["b", 2]
        # category c is kept in .categories
        tm.assert_frame_equal(df, exp_fancy)

        # set_value
        df = orig.copy()
        df.at["j", "cats"] = "b"
        tm.assert_frame_equal(df, exp_single_cats_value)

        def f():
            df = orig.copy()
            df.at["j", "cats"] = "c"

        pytest.raises(ValueError, f)

        # Assigning a Category to parts of a int/... column uses the values of
        # the Catgorical
        df = DataFrame({"a": [1, 1, 1, 1, 1], "b": list("aaaaa")})
        exp = DataFrame({"a": [1, "b", "b", 1, 1], "b": list("aabba")})
        df.loc[1:2, "a"] = Categorical(["b", "b"], categories=["a", "b"])
        df.loc[2:3, "b"] = Categorical(["b", "b"], categories=["a", "b"])
        tm.assert_frame_equal(df, exp)

        # Series
        orig = Series(Categorical(["b", "b"], categories=["a", "b"]))
        s = orig.copy()
        s[:] = "a"
        exp = Series(Categorical(["a", "a"], categories=["a", "b"]))
        tm.assert_series_equal(s, exp)

        s = orig.copy()
        s[1] = "a"
        exp = Series(Categorical(["b", "a"], categories=["a", "b"]))
        tm.assert_series_equal(s, exp)

        s = orig.copy()
        s[s.index > 0] = "a"
        exp = Series(Categorical(["b", "a"], categories=["a", "b"]))
        tm.assert_series_equal(s, exp)

        s = orig.copy()
        s[[False, True]] = "a"
        exp = Series(Categorical(["b", "a"], categories=["a", "b"]))
        tm.assert_series_equal(s, exp)

        s = orig.copy()
        s.index = ["x", "y"]
        s["y"] = "a"
        exp = Series(Categorical(["b", "a"], categories=["a", "b"]),
                     index=["x", "y"])
        tm.assert_series_equal(s, exp)

        # ensure that one can set something to np.nan
        s = Series(Categorical([1, 2, 3]))
        exp = Series(Categorical([1, np.nan, 3], categories=[1, 2, 3]))
        s[1] = np.nan
        tm.assert_series_equal(s, exp)
