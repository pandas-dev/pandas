import numpy as np
import pytest

import pandas as pd
from pandas import CategoricalDtype, DataFrame, IntervalIndex, MultiIndex, Series
import pandas.util.testing as tm


class TestDataFrameSortIndex:
    def test_sort_index_nan(self):
        # GH#3917

        # Test DataFrame with nan label
        df = DataFrame(
            {"A": [1, 2, np.nan, 1, 6, 8, 4], "B": [9, np.nan, 5, 2, 5, 4, 5]},
            index=[1, 2, 3, 4, 5, 6, np.nan],
        )

        # NaN label, ascending=True, na_position='last'
        sorted_df = df.sort_index(kind="quicksort", ascending=True, na_position="last")
        expected = DataFrame(
            {"A": [1, 2, np.nan, 1, 6, 8, 4], "B": [9, np.nan, 5, 2, 5, 4, 5]},
            index=[1, 2, 3, 4, 5, 6, np.nan],
        )
        tm.assert_frame_equal(sorted_df, expected)

        # NaN label, ascending=True, na_position='first'
        sorted_df = df.sort_index(na_position="first")
        expected = DataFrame(
            {"A": [4, 1, 2, np.nan, 1, 6, 8], "B": [5, 9, np.nan, 5, 2, 5, 4]},
            index=[np.nan, 1, 2, 3, 4, 5, 6],
        )
        tm.assert_frame_equal(sorted_df, expected)

        # NaN label, ascending=False, na_position='last'
        sorted_df = df.sort_index(kind="quicksort", ascending=False)
        expected = DataFrame(
            {"A": [8, 6, 1, np.nan, 2, 1, 4], "B": [4, 5, 2, 5, np.nan, 9, 5]},
            index=[6, 5, 4, 3, 2, 1, np.nan],
        )
        tm.assert_frame_equal(sorted_df, expected)

        # NaN label, ascending=False, na_position='first'
        sorted_df = df.sort_index(
            kind="quicksort", ascending=False, na_position="first"
        )
        expected = DataFrame(
            {"A": [4, 8, 6, 1, np.nan, 2, 1], "B": [5, 4, 5, 2, 5, np.nan, 9]},
            index=[np.nan, 6, 5, 4, 3, 2, 1],
        )
        tm.assert_frame_equal(sorted_df, expected)

    def test_sort_index_multi_index(self):
        # GH#25775, testing that sorting by index works with a multi-index.
        df = DataFrame(
            {"a": [3, 1, 2], "b": [0, 0, 0], "c": [0, 1, 2], "d": list("abc")}
        )
        result = df.set_index(list("abc")).sort_index(level=list("ba"))

        expected = DataFrame(
            {"a": [1, 2, 3], "b": [0, 0, 0], "c": [1, 2, 0], "d": list("bca")}
        )
        expected = expected.set_index(list("abc"))

        tm.assert_frame_equal(result, expected)

    def test_sort_index_inplace(self):
        frame = DataFrame(
            np.random.randn(4, 4), index=[1, 2, 3, 4], columns=["A", "B", "C", "D"]
        )

        # axis=0
        unordered = frame.loc[[3, 2, 4, 1]]
        a_id = id(unordered["A"])
        df = unordered.copy()
        df.sort_index(inplace=True)
        expected = frame
        tm.assert_frame_equal(df, expected)
        assert a_id != id(df["A"])

        df = unordered.copy()
        df.sort_index(ascending=False, inplace=True)
        expected = frame[::-1]
        tm.assert_frame_equal(df, expected)

        # axis=1
        unordered = frame.loc[:, ["D", "B", "C", "A"]]
        df = unordered.copy()
        df.sort_index(axis=1, inplace=True)
        expected = frame
        tm.assert_frame_equal(df, expected)

        df = unordered.copy()
        df.sort_index(axis=1, ascending=False, inplace=True)
        expected = frame.iloc[:, ::-1]
        tm.assert_frame_equal(df, expected)

    def test_sort_index_different_sortorder(self):
        A = np.arange(20).repeat(5)
        B = np.tile(np.arange(5), 20)

        indexer = np.random.permutation(100)
        A = A.take(indexer)
        B = B.take(indexer)

        df = DataFrame({"A": A, "B": B, "C": np.random.randn(100)})

        ex_indexer = np.lexsort((df.B.max() - df.B, df.A))
        expected = df.take(ex_indexer)

        # test with multiindex, too
        idf = df.set_index(["A", "B"])

        result = idf.sort_index(ascending=[1, 0])
        expected = idf.take(ex_indexer)
        tm.assert_frame_equal(result, expected)

        # also, Series!
        result = idf["C"].sort_index(ascending=[1, 0])
        tm.assert_series_equal(result, expected["C"])

    def test_sort_index_level(self):
        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list("ABC"))
        df = DataFrame([[1, 2], [3, 4]], mi)

        result = df.sort_index(level="A", sort_remaining=False)
        expected = df
        tm.assert_frame_equal(result, expected)

        result = df.sort_index(level=["A", "B"], sort_remaining=False)
        expected = df
        tm.assert_frame_equal(result, expected)

        # Error thrown by sort_index when
        # first index is sorted last (GH#26053)
        result = df.sort_index(level=["C", "B", "A"])
        expected = df.iloc[[1, 0]]
        tm.assert_frame_equal(result, expected)

        result = df.sort_index(level=["B", "C", "A"])
        expected = df.iloc[[1, 0]]
        tm.assert_frame_equal(result, expected)

        result = df.sort_index(level=["C", "A"])
        expected = df.iloc[[1, 0]]
        tm.assert_frame_equal(result, expected)

    def test_sort_index_categorical_index(self):

        df = DataFrame(
            {
                "A": np.arange(6, dtype="int64"),
                "B": Series(list("aabbca")).astype(CategoricalDtype(list("cab"))),
            }
        ).set_index("B")

        result = df.sort_index()
        expected = df.iloc[[4, 0, 1, 5, 2, 3]]
        tm.assert_frame_equal(result, expected)

        result = df.sort_index(ascending=False)
        expected = df.iloc[[2, 3, 0, 1, 5, 4]]
        tm.assert_frame_equal(result, expected)

    def test_sort_index(self):
        # GH#13496

        frame = DataFrame(
            np.arange(16).reshape(4, 4),
            index=[1, 2, 3, 4],
            columns=["A", "B", "C", "D"],
        )

        # axis=0 : sort rows by index labels
        unordered = frame.loc[[3, 2, 4, 1]]
        result = unordered.sort_index(axis=0)
        expected = frame
        tm.assert_frame_equal(result, expected)

        result = unordered.sort_index(ascending=False)
        expected = frame[::-1]
        tm.assert_frame_equal(result, expected)

        # axis=1 : sort columns by column names
        unordered = frame.iloc[:, [2, 1, 3, 0]]
        result = unordered.sort_index(axis=1)
        tm.assert_frame_equal(result, frame)

        result = unordered.sort_index(axis=1, ascending=False)
        expected = frame.iloc[:, ::-1]
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("level", ["A", 0])  # GH#21052
    def test_sort_index_multiindex(self, level):
        # GH#13496

        # sort rows by specified level of multi-index
        mi = MultiIndex.from_tuples(
            [[2, 1, 3], [2, 1, 2], [1, 1, 1]], names=list("ABC")
        )
        df = DataFrame([[1, 2], [3, 4], [5, 6]], index=mi)

        expected_mi = MultiIndex.from_tuples(
            [[1, 1, 1], [2, 1, 2], [2, 1, 3]], names=list("ABC")
        )
        expected = pd.DataFrame([[5, 6], [3, 4], [1, 2]], index=expected_mi)
        result = df.sort_index(level=level)
        tm.assert_frame_equal(result, expected)

        # sort_remaining=False
        expected_mi = MultiIndex.from_tuples(
            [[1, 1, 1], [2, 1, 3], [2, 1, 2]], names=list("ABC")
        )
        expected = pd.DataFrame([[5, 6], [1, 2], [3, 4]], index=expected_mi)
        result = df.sort_index(level=level, sort_remaining=False)
        tm.assert_frame_equal(result, expected)

    def test_sort_index_intervalindex(self):
        # this is a de-facto sort via unstack
        # confirming that we sort in the order of the bins
        y = Series(np.random.randn(100))
        x1 = Series(np.sign(np.random.randn(100)))
        x2 = pd.cut(Series(np.random.randn(100)), bins=[-3, -0.5, 0, 0.5, 3])
        model = pd.concat([y, x1, x2], axis=1, keys=["Y", "X1", "X2"])

        result = model.groupby(["X1", "X2"], observed=True).mean().unstack()
        expected = IntervalIndex.from_tuples(
            [(-3.0, -0.5), (-0.5, 0.0), (0.0, 0.5), (0.5, 3.0)], closed="right"
        )
        result = result.columns.levels[1].categories
        tm.assert_index_equal(result, expected)
