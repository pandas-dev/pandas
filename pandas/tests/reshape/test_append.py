from itertools import product

import numpy as np
import pytest

import pandas as pd
from pandas.core.indexes.base import InvalidIndexError
from pandas.util.testing import assert_frame_equal


indexes = [
    # indexes listed here must be sorted

    # base
    pd.Index(['A', 'B', 'C']),

    # numeric
    pd.RangeIndex(3),
    pd.Int64Index([3, 4, 5]),
    pd.UInt64Index([6, 7, 8]),
    pd.Float64Index([3.5, 4.5, 5.5]),

    # datetime
    pd.to_datetime(['2013-01-01', '2013-01-10', '2013-01-15']),
    pd.to_timedelta(['1 day', '2 days', '3 days']),
    pd.PeriodIndex(start='2000', periods=3),

    # interval
    pd.interval_range(start=0, end=3),

    # categorical
    pd.CategoricalIndex('A B C'.split()),
    pd.CategoricalIndex('D E F'.split(), ordered=True),

    # multi-index
    pd.MultiIndex.from_arrays(['A B C'.split(), 'D E F'.split()]),
]


index_sort_groups = [
    # When indexes from the same group are joined, the result is sortable.
    # When indexes from different groups are joined, the result is not
    # sortable.

    [  # joining produces a string index
     pd.Index(['A', 'B', 'C']),
     pd.CategoricalIndex('A B C'.split()),
     pd.CategoricalIndex('D E F'.split(), ordered=True)],

    [  # numeric indexes
     pd.RangeIndex(3),
     pd.Int64Index([3, 4, 5]),
     pd.UInt64Index([6, 7, 8]),
     pd.Float64Index([3.5, 4.5, 5.5])],

    [pd.to_datetime(['2013-01-01', '2013-01-10', '2013-01-15'])],
    [pd.to_timedelta(['1 day', '2 days', '3 days'])],
    [pd.PeriodIndex(start='2000', periods=3)],
    [pd.interval_range(start=0, end=3)],
    [pd.MultiIndex.from_arrays(['A B C'.split(), 'D E F'.split()])],
]


def cls_name(obj):
    return obj.__class__.__name__


@pytest.fixture(params=[True, False])
def sort(request):
    """Boolean sort keyword for DataFrame.append
    """
    return request.param


class TestAppendBasic(object):
    def test_different_types_of_input(self, sort):
        # There are 7 types of accepted input by append:
        #
        # dict
        # Series
        # DataFrame
        # empty list
        # list of dicts
        # list of Series
        # list of DataFrames
        #
        # Using one or another should always be interchangeable.

        # append to dict
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        map = {
            0: 7,
            1: 8,
            2: 9
        }
        result = df.append(map, ignore_index=True, sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        assert_frame_equal(result, expected)

        # append to Series
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        ser = pd.Series([7, 8, 9])
        result = df.append(ser, ignore_index=True, sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        assert_frame_equal(result, expected)

        # append to DataFrame
        df1 = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        df2 = pd.DataFrame([[7, 8, 9]])
        result = df1.append(df2, ignore_index=True, sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        assert_frame_equal(result, expected)

        # append to empty list
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        result = df1.append([], sort=sort)
        expected = df
        assert_frame_equal(result, expected)

        # append to list of dicts
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        map = {
            0: 7,
            1: 8,
            2: 9
        }
        result = df.append([map], ignore_index=True, sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        assert_frame_equal(result, expected)

        # append to list of Series
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        ser = pd.Series([7, 8, 9])
        result = df.append([ser], ignore_index=True, sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        assert_frame_equal(result, expected)

        # append to list of DataFrames
        df1 = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        df2 = pd.DataFrame([[7, 8, 9]])
        result = df1.append([df2], ignore_index=True, sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        assert_frame_equal(result, expected)

        # append to list of dicts (2 dicts)
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        map = {
            0: 7,
            1: 8,
            2: 9
        }
        result = df.append([map, map], ignore_index=True, sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9], [7, 8, 9]])
        assert_frame_equal(result, expected)

        # append to list of Series (2 series)
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        ser = pd.Series([7, 8, 9])
        result = df.append([ser, ser], ignore_index=True, sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9], [7, 8, 9]])
        assert_frame_equal(result, expected)

        # append to list of DataFrames (2 dframes)
        df1 = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        df2 = pd.DataFrame([[7, 8, 9]])
        result = df1.append([df2, df2], ignore_index=True, sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9], [7, 8, 9]])
        assert_frame_equal(result, expected)

    def test_bad_input_type(self, sort):
        # When appending a bad input type, the function
        # should raise an exception.

        bad_input_msg = r'The value of other must be .*'
        mixed_list_msg = r'When other is a list, its .*'

        # integer input
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError, match=bad_input_msg):
            df.append(1, ignore_index=True, sort=sort)

        # string input
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError, match=bad_input_msg):
            df.append("1 2 3", ignore_index=True, sort=sort)

        # tuple input
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError, match=bad_input_msg):
            df.append((df, ), ignore_index=True, sort=sort)

        # list of integers
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError, match=bad_input_msg):
            df.append([1], ignore_index=True, sort=sort)

        # list of strings
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError, match=bad_input_msg):
            df.append(["1 2 3"], ignore_index=True, sort=sort)

        # list of lists
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError, match=bad_input_msg):
            df.append([[df]], ignore_index=True, sort=sort)

        # list of tuples
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        with pytest.raises(TypeError, match=bad_input_msg):
            df.append([(df, )], ignore_index=True, sort=sort)

        # mixed list
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        ser = pd.Series([7, 8, 9])
        dict = {
            0: 10,
            1: 11,
            2: 12
        }
        with pytest.raises(TypeError, match=mixed_list_msg):
            df.append([ser, dict], ignore_index=True, sort=sort)
        with pytest.raises(TypeError, match=mixed_list_msg):
            df.append([dict, ser], ignore_index=True, sort=sort)

        # mixed list with bad first element
        # (when the first element is bad, display the
        #  bad input msg instead of the mixed list one)
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        ser = pd.Series([7, 8, 9])
        with pytest.raises(TypeError, match=bad_input_msg):
            df.append([1, ser, ser], ignore_index=True, sort=sort)

        # mixed list with bad second element
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        ser = pd.Series([7, 8, 9])
        with pytest.raises(TypeError, match=mixed_list_msg):
            df.append([ser, 1, ser], ignore_index=True, sort=sort)

        # mixed list with bad third element
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        ser = pd.Series([7, 8, 9])
        with pytest.raises(TypeError, match=mixed_list_msg):
            df.append([ser, ser, 1], ignore_index=True, sort=sort)

    def test_no_unecessary_upcast(self, sort):
        # GH: 22621
        # When appending, the resulting columns should
        # not be float64 without necessity.

        # basic
        df1 = pd.DataFrame([[1, 2, 3]])
        df2 = pd.DataFrame([[4, 5, 6]], index=[1])
        result = df1.append(df2, sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        assert_frame_equal(result, expected)

        # 0 rows 0 columns
        df1 = pd.DataFrame([[1, 2, 3]])
        df2 = pd.DataFrame()
        result = df1.append(df2, sort=sort)
        expected = df1.copy()
        assert_frame_equal(result, expected)

        df1 = pd.DataFrame()
        df2 = pd.DataFrame([[1, 2, 3]])
        result = df1.append(df2, sort=sort)
        expected = df2.copy()
        assert_frame_equal(result, expected)

        # 0 rows 2 columns
        # (the original dtype (object) of the empty columns
        #  must be preserved)
        df1 = pd.DataFrame([[1, 2, 3]], columns=[0, 1, 2])
        df2 = pd.DataFrame(columns=[3, 4])
        result = df1.append(df2, sort=sort)
        expected = pd.DataFrame([[1, 2, 3, np.nan, np.nan]])
        expected[[3, 4]] = expected[[3, 4]].astype(object)
        assert_frame_equal(result, expected)

        df1 = pd.DataFrame(columns=[0, 1])
        df2 = pd.DataFrame([[1, 2, 3]], columns=[2, 3, 4])
        result = df1.append(df2, sort=sort)
        expected = pd.DataFrame([[np.nan, np.nan, 1, 2, 3]])
        expected[[0, 1]] = expected[[0, 1]].astype(object)
        assert_frame_equal(result, expected)

        # big.append(small)
        big = pd.DataFrame([[1, 2, 3]])
        small = pd.DataFrame([[4, 5]], index=[1])
        result = big.append(small, sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, np.nan]])
        assert_frame_equal(result, expected)

        # small.append(big)
        small = pd.DataFrame([[1, 2]])
        big = pd.DataFrame([[3, 4, 5]], index=[1])
        result = small.append(big, sort=sort)
        expected = pd.DataFrame([[1, 2, np.nan], [3, 4, 5]])
        assert_frame_equal(result, expected)


class TestAppendColumnsIndex(object):
    @pytest.mark.parametrize('idx_name3', [None, 'foo', 'bar', 'baz'])
    @pytest.mark.parametrize('idx_name2', [None, 'foo', 'bar', 'baz'])
    @pytest.mark.parametrize('idx_name1', [None, 'foo', 'bar', 'baz'])
    def test_preserve_index_name(self, sort, idx_name1, idx_name2, idx_name3):
        # When appending, the name of the indexes
        # of the base DataFrame must always be
        # preserved in the result.

        df1 = pd.DataFrame([[1, 2, 3]])
        df2 = pd.DataFrame([[4, 5, 6]], index=[1])
        df3 = pd.DataFrame([[7, 8, 9]], index=[2])

        df1.columns.name = idx_name1
        df2.columns.name = idx_name2
        df3.columns.name = idx_name3

        # append []
        result = df1.append([], sort=sort)
        expected = df1.copy()
        assert_frame_equal(result, expected)

        # append [df]
        result = df1.append([df2], sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        expected.columns.name = idx_name1
        assert_frame_equal(result, expected)

        # append [df, df]
        result = df1.append([df2, df3], sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        expected.columns.name = idx_name1
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize('index', indexes, ids=cls_name)
    def test_preserve_index_type(self, sort, index):
        # when there's only one index type in the inputs,
        # it must be preserved in the output.

        # basic
        df1 = pd.DataFrame([[1, 2, 3]], columns=index)
        df2 = pd.DataFrame([[4, 5, 6]], index=[1], columns=index)
        result = df1.append(df2, sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6]], columns=index)
        assert_frame_equal(result, expected)

        # big.append(small)
        big = pd.DataFrame([[1, 2, 3]], columns=index)
        small = pd.DataFrame([[4, 5]], index=[1], columns=index[:2])
        result = big.append(small, sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, np.nan]], columns=index)
        assert_frame_equal(result, expected)

        # small.append(big)
        small = pd.DataFrame([[1, 2]], columns=index[:2])
        big = pd.DataFrame([[3, 4, 5]], index=[1], columns=index)
        result = small.append(big, sort=sort)
        expected = pd.DataFrame([[1, 2, np.nan], [3, 4, 5]], columns=index)
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize('index2', indexes, ids=cls_name)
    @pytest.mark.parametrize('index1', indexes, ids=cls_name)
    def test_preserve_index_values_without_sort(self, index1, index2):
        # When appending indexes of different types, we want
        # the resulting index to preserve the exact indexes
        # values.

        # Related to GH13626
        from pandas.core.dtypes.generic import (
            ABCDatetimeIndex, ABCMultiIndex, ABCTimedeltaIndex
        )
        if isinstance(index1, ABCMultiIndex):
            if isinstance(index2, ABCDatetimeIndex):
                pytest.xfail("MultiIndex + DatetimeIndex produces bad value")
            if isinstance(index2, ABCTimedeltaIndex):
                pytest.xfail("MultiIndex + TimedeltaIndex produces bad value")

        df1 = pd.DataFrame([[1, 2, 3]], columns=index1)
        df2 = pd.DataFrame([[4, 5, 6]], columns=index2, index=[1])
        result = df1.append(df2, sort=False)
        for value in index1:
            assert value in result.columns
        for value in index2:
            assert value in result.columns

    @pytest.mark.parametrize(
        'index1, index2',
        [(i1, i2)
            for group in index_sort_groups
            for i1, i2 in product(group, repeat=2)],
        ids=cls_name
    )
    def test_preserve_index_values_with_sort(self, index1, index2):
        # When appending indexes of different types, we want
        # the resulting index to preserve the exact indexes
        # values.

        df1 = pd.DataFrame([[1, 2, 3]], columns=index1)
        df2 = pd.DataFrame([[4, 5, 6]], columns=index2, index=[1])
        result = df1.append(df2, sort=True)
        for value in index1:
            assert value in result.columns
        for value in index2:
            assert value in result.columns

    def test_raise_on_duplicates(self, sort):
        # Append should not allow DataFrames with repeated
        # column names (or series with repeated row names).

        # dupe on base
        df1 = pd.DataFrame([[1, 2, 3]], columns=['A', 'B', 'B'])
        df2 = pd.DataFrame([[1, 2, 3]], columns=['A', 'B', 'C'])
        with pytest.raises(InvalidIndexError):
            df1.append([], sort=sort)
        with pytest.raises(InvalidIndexError):
            df1.append([df2], sort=sort)
        with pytest.raises(InvalidIndexError):
            df1.append([df2, df2], sort=sort)

        # dupe on other
        df1 = pd.DataFrame([[1, 2, 3]], columns=['A', 'B', 'C'])
        df2 = pd.DataFrame([[1, 2, 3]], columns=['A', 'B', 'B'])
        with pytest.raises(InvalidIndexError):
            df1.append([df2], sort=sort)
        with pytest.raises(InvalidIndexError):
            df1.append([df2, df2], sort=sort)

        # dupe on both
        # (we could avoid raising errors here, but, to keep the api
        #  consistent, we don't)
        df1 = pd.DataFrame([[1, 2, 3]], columns=['A', 'B', 'B'])
        df2 = pd.DataFrame([[1, 2, 3]], columns=['A', 'B', 'B'])
        with pytest.raises(InvalidIndexError):
            df1.append([], sort=sort)
        with pytest.raises(InvalidIndexError):
            df1.append([df2], sort=sort)
        with pytest.raises(InvalidIndexError):
            df1.append([df2, df2], sort=sort)

    def test_nosort_basic(self):
        # When sort=False, the resulting columns come
        # in the order that they appear in the inputs.

        nan = np.nan

        # NUMERIC INDEX TESTS

        # append []
        df = pd.DataFrame([[1, 2, 3]], columns=[0, 1, 2])
        result = df.append([], sort=False)
        expected = df[[0, 1, 2]]
        assert_frame_equal(result, expected)

        df = pd.DataFrame([[1, 2, 3]], columns=[2, 1, 0])
        result = df.append([], sort=False)
        expected = df[[2, 1, 0]]
        assert_frame_equal(result, expected)

        # append [df]
        df1 = pd.DataFrame([[1, 2]], columns=[0.0, 1.0])
        df2 = pd.DataFrame([[1, 2]], columns=[0.5, 1.5], index=[1])
        result = df1.append(df2, sort=False)
        expected = pd.DataFrame([[1, 2, nan, nan],
                                 [nan, nan, 1, 2]],
                                columns=[0.0, 1.0, 0.5, 1.5])
        assert_frame_equal(result, expected)

        # append [df, df]
        df1 = pd.DataFrame([[1, 2]], columns=[0.0, 1.0])
        df2 = pd.DataFrame([[1, 2]], columns=[0.3, 1.3], index=[1])
        df3 = pd.DataFrame([[1, 2]], columns=[0.6, 1.6], index=[2])
        result = df1.append([df2, df3], sort=False)
        expected = pd.DataFrame([[1, 2, nan, nan, nan, nan],
                                 [nan, nan, 1, 2, nan, nan],
                                 [nan, nan, nan, nan, 1, 2]],
                                columns=[0.0, 1.0, 0.3, 1.3, 0.6, 1.6])
        assert_frame_equal(result, expected)

        # STRING INDEX TESTS

        # append []
        df = pd.DataFrame([[1, 2, 3]], columns=['a', 'b', 'c'])
        result = df.append([], sort=False)
        expected = df[['a', 'b', 'c']]
        assert_frame_equal(result, expected)

        df = pd.DataFrame([[1, 2, 3]], columns=['c', 'b', 'a'])
        result = df.append([], sort=False)
        expected = df[['c', 'b', 'a']]
        assert_frame_equal(result, expected)

        # append [df]
        df1 = pd.DataFrame([[1, 2]], columns=['a', 'c'])
        df2 = pd.DataFrame([[1, 2]], columns=['b', 'd'], index=[1])
        result = df1.append(df2, sort=False)
        expected = pd.DataFrame([[1, 2, nan, nan],
                                 [nan, nan, 1, 2]],
                                columns=['a', 'c', 'b', 'd'])
        assert_frame_equal(result, expected)

        # append [df, df]
        df1 = pd.DataFrame([[1, 2]], columns=['a', 'd'])
        df2 = pd.DataFrame([[1, 2]], columns=['b', 'e'], index=[1])
        df3 = pd.DataFrame([[1, 2]], columns=['c', 'f'], index=[2])
        result = df1.append([df2, df3], sort=False)
        expected = pd.DataFrame([[1, 2, nan, nan, nan, nan],
                                 [nan, nan, 1, 2, nan, nan],
                                 [nan, nan, nan, nan, 1, 2]],
                                columns=['a', 'd', 'b', 'e', 'c', 'f'])
        assert_frame_equal(result, expected)

    def test_sort_basic(self):
        # When sort=True, the resulting columns must come
        # out sorted.

        nan = np.nan

        # NUMERIC INDEX TESTS

        # append []
        df = pd.DataFrame([[1, 2, 3]], columns=[0, 1, 2])
        result = df.append([], sort=True)
        expected = df[[0, 1, 2]]
        assert_frame_equal(result, expected)

        df = pd.DataFrame([[1, 2, 3]], columns=[2, 1, 0])
        result = df.append([], sort=True)
        expected = df[[0, 1, 2]]
        assert_frame_equal(result, expected)

        # append [df]
        df1 = pd.DataFrame([[1, 2]], columns=[0.0, 1.0])
        df2 = pd.DataFrame([[1, 2]], columns=[0.5, 1.5], index=[1])
        result = df1.append(df2, sort=True)
        expected = pd.DataFrame([[1, nan, 2, nan],
                                 [nan, 1, nan, 2]],
                                columns=[0.0, 0.5, 1.0, 1.5])
        assert_frame_equal(result, expected)

        # append [df, df]
        df1 = pd.DataFrame([[1, 2]], columns=[0.0, 1.0])
        df2 = pd.DataFrame([[1, 2]], columns=[0.3, 1.3], index=[1])
        df3 = pd.DataFrame([[1, 2]], columns=[0.6, 1.6], index=[2])
        result = df1.append([df2, df3], sort=True)
        expected = pd.DataFrame([[1, nan, nan, 2, nan, nan],
                                 [nan, 1, nan, nan, 2, nan],
                                 [nan, nan, 1, nan, nan, 2]],
                                columns=[0.0, 0.3, 0.6, 1.0, 1.3, 1.6])
        assert_frame_equal(result, expected)

        # STRING INDEX TESTS

        # append []
        df = pd.DataFrame([[1, 2, 3]], columns=['a', 'b', 'c'])
        result = df.append([], sort=True)
        expected = df[['a', 'b', 'c']]
        assert_frame_equal(result, expected)

        df = pd.DataFrame([[1, 2, 3]], columns=['c', 'b', 'a'])
        result = df.append([], sort=True)
        expected = df[['a', 'b', 'c']]
        assert_frame_equal(result, expected)

        # append [df]
        df1 = pd.DataFrame([[1, 2]], columns=['a', 'c'])
        df2 = pd.DataFrame([[1, 2]], columns=['b', 'd'], index=[1])
        result = df1.append(df2, sort=True)
        expected = pd.DataFrame([[1, nan, 2, nan],
                                 [nan, 1, nan, 2]],
                                columns=['a', 'b', 'c', 'd'])
        assert_frame_equal(result, expected)

        # append [df, df]
        df1 = pd.DataFrame([[1, 2]], columns=['a', 'd'])
        df2 = pd.DataFrame([[1, 2]], columns=['b', 'e'], index=[1])
        df3 = pd.DataFrame([[1, 2]], columns=['c', 'f'], index=[2])
        result = df1.append([df2, df3], sort=True)
        expected = pd.DataFrame([[1, nan, nan, 2, nan, nan],
                                 [nan, 1, nan, nan, 2, nan],
                                 [nan, nan, 1, nan, nan, 2]],
                                columns=['a', 'b', 'c', 'd', 'e', 'f'])
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize('index2', indexes, ids=cls_name)
    @pytest.mark.parametrize('index1', indexes, ids=cls_name)
    def test_index_types_without_sort(self, index1, index2):
        # We should be able to append to a DataFrame
        # regardless of the type of its index.

        # TODO: check end of append and create tests (empty / IntervalIndex)
        # TODO: implement different way for df.append([])
        from pandas.core.dtypes.generic import ABCIntervalIndex
        if isinstance(index1, ABCIntervalIndex):
            pytest.xfail("Cannot do df[interval] for IntervalIndex")

        # the code below should not raise any exceptions
        df1 = pd.DataFrame([[1, 2, 3]], columns=index1)
        df2 = pd.DataFrame([[4, 5, 6]], columns=index2, index=[1])
        df1.append([], sort=False)
        df1.append([df2], sort=False)
        df1.append([df2, df2], sort=False)

    @pytest.mark.parametrize(
        'index1, index2',
        [(i1, i2)
            for group in index_sort_groups
            for i1, i2 in product(group, repeat=2)],
        ids=cls_name
    )
    def test_index_types_with_possible_sort(self, index1, index2):
        # When the result of joining two indexes is sortable,
        # we should not raise any exceptions.

        # TODO: check end of append and create tests (empty / IntervalIndex)
        # TODO: implement different way for df.append([])
        from pandas.core.dtypes.generic import ABCIntervalIndex
        if isinstance(index1, ABCIntervalIndex):
            pytest.xfail("Cannot do df[interval] for IntervalIndex")

        df1 = pd.DataFrame([[1, 2, 3]], columns=index1)
        df2 = pd.DataFrame([[4, 5, 6]], columns=index2, index=[1])
        df1.append([], sort=True)  # sorts the original frame
        df1.append([df2], sort=True)
        df1.append([df2, df2], sort=True)

    @pytest.mark.parametrize(
        'index1, index2',
        [(i1, i2)
            for g1, g2 in product(index_sort_groups, repeat=2)
                # different sort groups
                if type(g1[0]) != type(g2[0])
            for i1, i2 in product(g1, g2)],
        ids=cls_name
    )
    def test_index_types_with_impossible_sort(self, index1, index2):
        # When the result of joining two indexes is not sortable,
        # we should raise an exception.

        # TODO: check end of append and create tests (empty / IntervalIndex)
        # TODO: implement different way for df.append([])
        from pandas.core.dtypes.generic import ABCIntervalIndex
        if isinstance(index1, ABCIntervalIndex):
            pytest.xfail("Cannot do df[interval] for IntervalIndex")

        err_msg = r'The resulting columns could not be sorted\..*'

        df1 = pd.DataFrame([[1, 2, 3]], columns=index1)
        df2 = pd.DataFrame([[4, 5, 6]], columns=index2, index=[1])

        with pytest.raises(TypeError, match=err_msg):
            df1.append([df2], sort=True)
        with pytest.raises(TypeError, match=err_msg):
            df1.append([df2, df2], sort=True)


class TestAppendRowsIndex(object):
    @pytest.mark.parametrize('idx_name3', [None, 'foo', 'bar', 'baz'])
    @pytest.mark.parametrize('idx_name2', [None, 'foo', 'bar', 'baz'])
    @pytest.mark.parametrize('idx_name1', [None, 'foo', 'bar', 'baz'])
    def test_preserve_index_name(self, sort, idx_name1, idx_name2, idx_name3):
        # When appending, the name of the indexes
        # of the base DataFrame must always be
        # preserved in the result.

        df1 = pd.DataFrame([[1, 2, 3]])
        df2 = pd.DataFrame([[4, 5, 6]], index=[1])
        df3 = pd.DataFrame([[7, 8, 9]], index=[2])

        df1.index.name = idx_name1
        df2.index.name = idx_name2
        df3.index.name = idx_name3

        # append []
        result = df1.append([], sort=sort)
        expected = df1.copy()
        assert_frame_equal(result, expected)

        # append [df]
        result = df1.append([df2], sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        expected.index.name = idx_name1
        assert_frame_equal(result, expected)

        # append [df, df]
        result = df1.append([df2, df3], sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        expected.index.name = idx_name1
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize('index', indexes, ids=cls_name)
    def test_preserve_index_type(self, sort, index):
        # when there's only one index type in the inputs,
        # it must be preserved in the output.

        index1 = index[:1]
        index2 = index[1:2]
        index_comb = index1.append(index2)

        df1 = pd.DataFrame([[1, 2, 3]], index=index1)
        df2 = pd.DataFrame([[4, 5, 6]], index=index2)
        result = df1.append(df2, sort=sort)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=index_comb)
        assert_frame_equal(result, expected)

    @pytest.mark.parametrize('index2', indexes, ids=cls_name)
    @pytest.mark.parametrize('index1', indexes, ids=cls_name)
    def test_preserve_index_values(self, sort, index1, index2):
        # When appending indexes of different types, we want
        # the resulting index to preserve the exact indexes
        # values.

        # Related to GH13626
        from pandas.core.dtypes.generic import (
            ABCDatetimeIndex, ABCMultiIndex, ABCTimedeltaIndex
        )
        if isinstance(index1, ABCMultiIndex):
            if isinstance(index2, ABCDatetimeIndex):
                pytest.xfail("MultiIndex + DatetimeIndex produces bad value")
            if isinstance(index2, ABCTimedeltaIndex):
                pytest.xfail("MultiIndex + TimedeltaIndex produces bad value")

        # Concat raises a TypeError when appending a CategoricalIndex
        # with another type
        from pandas.core.dtypes.generic import ABCCategoricalIndex
        if isinstance(index1, ABCCategoricalIndex):
            pytest.xfail("Cannot have a CategoricalIndex append to another typ")

        df1 = pd.DataFrame([[1, 2, 3]], index=index1[:1])
        df2 = pd.DataFrame([[4, 5, 6]], index=index2[:1])
        result = df1.append(df2, sort=sort)
        assert index1[0] in result.index
        assert index2[0] in result.index

    def test_duplicates_without_verify_integrity(self):
        # When verify_integrity=False, the function should
        # allow duplicate values in the rows index.

        raise NotImplementedError

    def test_duplicates_with_verify_integrity(self):
        # When verify_integrity=True, the function should
        # not allow duplicate values in the rows index (whether
        # in the input or output).

        raise NotImplementedError

    def test_ignore_index(self):
        # When ignore_index=True, the function should completely
        # ignore the input indexes and generate one that is brand
        # new (RangeIndex).

        raise NotImplementedError

    def test_warning_ignore_index_and_verify_integrity(self):
        # It makes no sense to set verify_integrity=True when
        # ignore_index=True. To warn of a possible user
        # misunderstanding, append should raise a warning in
        # this situation.

        raise NotImplementedError


class TestAppendDangling(object):
    """Tests that have not been concretized yet
    """

    def test_append_unnamed_series_raises(self, sort):
        dict_msg = 'Can only append a dict if ignore_index=True'
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        dict = {
            0: 7,
            1: 8,
            2: 9
        }
        with pytest.raises(TypeError, match=dict_msg):
            df.append(dict, sort=sort)
        with pytest.raises(TypeError, match=dict_msg):
            df.append([dict], sort=sort)
        with pytest.raises(TypeError, match=dict_msg):
            df.append([dict, dict], sort=sort)

        series_msg = 'Can only append a Series if ignore_index=True or .*'
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6]])
        series = pd.Series([7, 8, 9])
        with pytest.raises(TypeError, match=series_msg):
            df.append(series, sort=sort)
        with pytest.raises(TypeError, match=series_msg):
            df.append([series], sort=sort)
        with pytest.raises(TypeError, match=series_msg):
            df.append([series, series], sort=sort)

    indexes = [
        None,
        pd.Index([0, 1]),
        pd.Index(['a', 'b']),
        pd.Index(['a', 'b'], name='foo')
    ]

    @pytest.mark.parametrize('index1', indexes, ids=lambda x: repr(x))
    @pytest.mark.parametrize('index2', indexes, ids=lambda x: repr(x))
    def test_append_ignore_index(self, sort, index1, index2):
        # when appending with ignore_index=True,
        # all index content must be forgotten
        df1 = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=index1)
        df2 = pd.DataFrame([[1, 2, 3], [4, 5, 6]], index=index2)

        result = df1.append(df2, ignore_index=True)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6],
                                 [1, 2, 3], [4, 5, 6]])
        assert_frame_equal(result, expected)

        result = df1.append([df2], ignore_index=True)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6],
                                 [1, 2, 3], [4, 5, 6]])
        assert_frame_equal(result, expected)

        result = df1.append([df2, df2], ignore_index=True)
        expected = pd.DataFrame([[1, 2, 3], [4, 5, 6],
                                 [1, 2, 3], [4, 5, 6],
                                 [1, 2, 3], [4, 5, 6]])
        assert_frame_equal(result, expected)
