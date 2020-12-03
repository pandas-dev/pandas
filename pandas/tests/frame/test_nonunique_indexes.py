import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, MultiIndex, Series, date_range
import pandas._testing as tm


def check(result, expected=None):
    if expected is not None:
        tm.assert_frame_equal(result, expected)
    result.dtypes
    str(result)


class TestDataFrameNonuniqueIndexes:
    def test_column_dups_operations(self):

        # assignment
        # GH 3687
        arr = np.random.randn(3, 2)
        idx = list(range(2))
        df = DataFrame(arr, columns=["A", "A"])
        df.columns = idx
        expected = DataFrame(arr, columns=idx)
        check(df, expected)

        idx = date_range("20130101", periods=4, freq="Q-NOV")
        df = DataFrame(
            [[1, 1, 1, 5], [1, 1, 2, 5], [2, 1, 3, 5]], columns=["a", "a", "a", "a"]
        )
        df.columns = idx
        expected = DataFrame([[1, 1, 1, 5], [1, 1, 2, 5], [2, 1, 3, 5]], columns=idx)
        check(df, expected)

        # insert
        df = DataFrame(
            [[1, 1, 1, 5], [1, 1, 2, 5], [2, 1, 3, 5]],
            columns=["foo", "bar", "foo", "hello"],
        )
        df["string"] = "bah"
        expected = DataFrame(
            [[1, 1, 1, 5, "bah"], [1, 1, 2, 5, "bah"], [2, 1, 3, 5, "bah"]],
            columns=["foo", "bar", "foo", "hello", "string"],
        )
        check(df, expected)
        with pytest.raises(ValueError, match="Length of value"):
            df.insert(0, "AnotherColumn", range(len(df.index) - 1))

        # insert same dtype
        df["foo2"] = 3
        expected = DataFrame(
            [[1, 1, 1, 5, "bah", 3], [1, 1, 2, 5, "bah", 3], [2, 1, 3, 5, "bah", 3]],
            columns=["foo", "bar", "foo", "hello", "string", "foo2"],
        )
        check(df, expected)

        # set (non-dup)
        df["foo2"] = 4
        expected = DataFrame(
            [[1, 1, 1, 5, "bah", 4], [1, 1, 2, 5, "bah", 4], [2, 1, 3, 5, "bah", 4]],
            columns=["foo", "bar", "foo", "hello", "string", "foo2"],
        )
        check(df, expected)
        df["foo2"] = 3

        # delete (non dup)
        del df["bar"]
        expected = DataFrame(
            [[1, 1, 5, "bah", 3], [1, 2, 5, "bah", 3], [2, 3, 5, "bah", 3]],
            columns=["foo", "foo", "hello", "string", "foo2"],
        )
        check(df, expected)

        # try to delete again (its not consolidated)
        del df["hello"]
        expected = DataFrame(
            [[1, 1, "bah", 3], [1, 2, "bah", 3], [2, 3, "bah", 3]],
            columns=["foo", "foo", "string", "foo2"],
        )
        check(df, expected)

        # consolidate
        df = df._consolidate()
        expected = DataFrame(
            [[1, 1, "bah", 3], [1, 2, "bah", 3], [2, 3, "bah", 3]],
            columns=["foo", "foo", "string", "foo2"],
        )
        check(df, expected)

        # insert
        df.insert(2, "new_col", 5.0)
        expected = DataFrame(
            [[1, 1, 5.0, "bah", 3], [1, 2, 5.0, "bah", 3], [2, 3, 5.0, "bah", 3]],
            columns=["foo", "foo", "new_col", "string", "foo2"],
        )
        check(df, expected)

        # insert a dup
        with pytest.raises(ValueError, match="cannot insert"):
            df.insert(2, "new_col", 4.0)

        df.insert(2, "new_col", 4.0, allow_duplicates=True)
        expected = DataFrame(
            [
                [1, 1, 4.0, 5.0, "bah", 3],
                [1, 2, 4.0, 5.0, "bah", 3],
                [2, 3, 4.0, 5.0, "bah", 3],
            ],
            columns=["foo", "foo", "new_col", "new_col", "string", "foo2"],
        )
        check(df, expected)

        # delete (dup)
        del df["foo"]
        expected = DataFrame(
            [[4.0, 5.0, "bah", 3], [4.0, 5.0, "bah", 3], [4.0, 5.0, "bah", 3]],
            columns=["new_col", "new_col", "string", "foo2"],
        )
        tm.assert_frame_equal(df, expected)

        # dup across dtypes
        df = DataFrame(
            [[1, 1, 1.0, 5], [1, 1, 2.0, 5], [2, 1, 3.0, 5]],
            columns=["foo", "bar", "foo", "hello"],
        )
        check(df)

        df["foo2"] = 7.0
        expected = DataFrame(
            [[1, 1, 1.0, 5, 7.0], [1, 1, 2.0, 5, 7.0], [2, 1, 3.0, 5, 7.0]],
            columns=["foo", "bar", "foo", "hello", "foo2"],
        )
        check(df, expected)

        result = df["foo"]
        expected = DataFrame([[1, 1.0], [1, 2.0], [2, 3.0]], columns=["foo", "foo"])
        check(result, expected)

        # multiple replacements
        df["foo"] = "string"
        expected = DataFrame(
            [
                ["string", 1, "string", 5, 7.0],
                ["string", 1, "string", 5, 7.0],
                ["string", 1, "string", 5, 7.0],
            ],
            columns=["foo", "bar", "foo", "hello", "foo2"],
        )
        check(df, expected)

        del df["foo"]
        expected = DataFrame(
            [[1, 5, 7.0], [1, 5, 7.0], [1, 5, 7.0]], columns=["bar", "hello", "foo2"]
        )
        check(df, expected)

        # values
        df = DataFrame([[1, 2.5], [3, 4.5]], index=[1, 2], columns=["x", "x"])
        result = df.values
        expected = np.array([[1, 2.5], [3, 4.5]])
        assert (result == expected).all().all()

        # rename, GH 4403
        df4 = DataFrame(
            {"RT": [0.0454], "TClose": [22.02], "TExg": [0.0422]},
            index=MultiIndex.from_tuples(
                [(600809, 20130331)], names=["STK_ID", "RPT_Date"]
            ),
        )

        df5 = DataFrame(
            {
                "RPT_Date": [20120930, 20121231, 20130331],
                "STK_ID": [600809] * 3,
                "STK_Name": ["饡驦", "饡驦", "饡驦"],
                "TClose": [38.05, 41.66, 30.01],
            },
            index=MultiIndex.from_tuples(
                [(600809, 20120930), (600809, 20121231), (600809, 20130331)],
                names=["STK_ID", "RPT_Date"],
            ),
        )

        k = pd.merge(df4, df5, how="inner", left_index=True, right_index=True)
        result = k.rename(columns={"TClose_x": "TClose", "TClose_y": "QT_Close"})
        str(result)
        result.dtypes

        expected = DataFrame(
            [[0.0454, 22.02, 0.0422, 20130331, 600809, "饡驦", 30.01]],
            columns=[
                "RT",
                "TClose",
                "TExg",
                "RPT_Date",
                "STK_ID",
                "STK_Name",
                "QT_Close",
            ],
        ).set_index(["STK_ID", "RPT_Date"], drop=False)
        tm.assert_frame_equal(result, expected)

        # reindex is invalid!
        df = DataFrame(
            [[1, 5, 7.0], [1, 5, 7.0], [1, 5, 7.0]], columns=["bar", "a", "a"]
        )
        msg = "cannot reindex from a duplicate axis"
        with pytest.raises(ValueError, match=msg):
            df.reindex(columns=["bar"])
        with pytest.raises(ValueError, match=msg):
            df.reindex(columns=["bar", "foo"])

        # drop
        df = DataFrame(
            [[1, 5, 7.0], [1, 5, 7.0], [1, 5, 7.0]], columns=["bar", "a", "a"]
        )
        result = df.drop(["a"], axis=1)
        expected = DataFrame([[1], [1], [1]], columns=["bar"])
        check(result, expected)
        result = df.drop("a", axis=1)
        check(result, expected)

        # describe
        df = DataFrame(
            [[1, 1, 1], [2, 2, 2], [3, 3, 3]],
            columns=["bar", "a", "a"],
            dtype="float64",
        )
        result = df.describe()
        s = df.iloc[:, 0].describe()
        expected = pd.concat([s, s, s], keys=df.columns, axis=1)
        check(result, expected)

        # check column dups with index equal and not equal to df's index
        df = DataFrame(
            np.random.randn(5, 3),
            index=["a", "b", "c", "d", "e"],
            columns=["A", "B", "A"],
        )
        for index in [df.index, pd.Index(list("edcba"))]:
            this_df = df.copy()
            expected_ser = Series(index.values, index=this_df.index)
            expected_df = DataFrame(
                {"A": expected_ser, "B": this_df["B"], "A": expected_ser},
                columns=["A", "B", "A"],
            )
            this_df["A"] = index
            check(this_df, expected_df)

        # operations
        for op in ["__add__", "__mul__", "__sub__", "__truediv__"]:
            df = DataFrame({"A": np.arange(10), "B": np.random.rand(10)})
            expected = getattr(df, op)(df)
            expected.columns = ["A", "A"]
            df.columns = ["A", "A"]
            result = getattr(df, op)(df)
            check(result, expected)

        # multiple assignments that change dtypes
        # the location indexer is a slice
        # GH 6120
        df = DataFrame(np.random.randn(5, 2), columns=["that", "that"])
        expected = DataFrame(1.0, index=range(5), columns=["that", "that"])

        df["that"] = 1.0
        check(df, expected)

        df = DataFrame(np.random.rand(5, 2), columns=["that", "that"])
        expected = DataFrame(1, index=range(5), columns=["that", "that"])

        df["that"] = 1
        check(df, expected)

    def test_column_dups2(self):

        # drop buggy GH 6240
        df = DataFrame(
            {
                "A": np.random.randn(5),
                "B": np.random.randn(5),
                "C": np.random.randn(5),
                "D": ["a", "b", "c", "d", "e"],
            }
        )

        expected = df.take([0, 1, 1], axis=1)
        df2 = df.take([2, 0, 1, 2, 1], axis=1)
        result = df2.drop("C", axis=1)
        tm.assert_frame_equal(result, expected)

        # dropna
        df = DataFrame(
            {
                "A": np.random.randn(5),
                "B": np.random.randn(5),
                "C": np.random.randn(5),
                "D": ["a", "b", "c", "d", "e"],
            }
        )
        df.iloc[2, [0, 1, 2]] = np.nan
        df.iloc[0, 0] = np.nan
        df.iloc[1, 1] = np.nan
        df.iloc[:, 3] = np.nan
        expected = df.dropna(subset=["A", "B", "C"], how="all")
        expected.columns = ["A", "A", "B", "C"]

        df.columns = ["A", "A", "B", "C"]

        result = df.dropna(subset=["A", "C"], how="all")
        tm.assert_frame_equal(result, expected)

    def test_getitem_boolean_series_with_duplicate_columns(self):
        # boolean indexing
        # GH 4879
        dups = ["A", "A", "C", "D"]
        df = DataFrame(
            np.arange(12).reshape(3, 4), columns=["A", "B", "C", "D"], dtype="float64"
        )
        expected = df[df.C > 6]
        expected.columns = dups
        df = DataFrame(np.arange(12).reshape(3, 4), columns=dups, dtype="float64")
        result = df[df.C > 6]
        check(result, expected)

    def test_getitem_boolean_frame_with_duplicate_columns(self):
        dups = ["A", "A", "C", "D"]

        # where
        df = DataFrame(
            np.arange(12).reshape(3, 4), columns=["A", "B", "C", "D"], dtype="float64"
        )
        # `df > 6` is a DataFrame with the same shape+alignment as df
        expected = df[df > 6]
        expected.columns = dups
        df = DataFrame(np.arange(12).reshape(3, 4), columns=dups, dtype="float64")
        result = df[df > 6]
        check(result, expected)

    def test_getitem_boolean_frame_unaligned_with_duplicate_columns(self):
        # `df.A > 6` is a DataFrame with a different shape from df
        dups = ["A", "A", "C", "D"]

        # boolean with the duplicate raises
        df = DataFrame(np.arange(12).reshape(3, 4), columns=dups, dtype="float64")
        msg = "cannot reindex from a duplicate axis"
        with pytest.raises(ValueError, match=msg):
            df[df.A > 6]

    def test_column_dups_indexing(self):

        # dup aligning operations should work
        # GH 5185
        df1 = DataFrame([1, 2, 3, 4, 5], index=[1, 2, 1, 2, 3])
        df2 = DataFrame([1, 2, 3], index=[1, 2, 3])
        expected = DataFrame([0, 2, 0, 2, 2], index=[1, 1, 2, 2, 3])
        result = df1.sub(df2)
        tm.assert_frame_equal(result, expected)

        # equality
        df1 = DataFrame([[1, 2], [2, np.nan], [3, 4], [4, 4]], columns=["A", "B"])
        df2 = DataFrame([[0, 1], [2, 4], [2, np.nan], [4, 5]], columns=["A", "A"])

        # not-comparing like-labelled
        msg = "Can only compare identically-labeled DataFrame objects"
        with pytest.raises(ValueError, match=msg):
            df1 == df2

        df1r = df1.reindex_like(df2)
        result = df1r == df2
        expected = DataFrame(
            [[False, True], [True, False], [False, False], [True, False]],
            columns=["A", "A"],
        )
        tm.assert_frame_equal(result, expected)

        # mixed column selection
        # GH 5639
        dfbool = DataFrame(
            {
                "one": Series([True, True, False], index=["a", "b", "c"]),
                "two": Series([False, False, True, False], index=["a", "b", "c", "d"]),
                "three": Series([False, True, True, True], index=["a", "b", "c", "d"]),
            }
        )
        expected = pd.concat([dfbool["one"], dfbool["three"], dfbool["one"]], axis=1)
        result = dfbool[["one", "three", "one"]]
        check(result, expected)

        # multi-axis dups
        # GH 6121
        df = DataFrame(
            np.arange(25.0).reshape(5, 5),
            index=["a", "b", "c", "d", "e"],
            columns=["A", "B", "C", "D", "E"],
        )
        z = df[["A", "C", "A"]].copy()
        expected = z.loc[["a", "c", "a"]]

        df = DataFrame(
            np.arange(25.0).reshape(5, 5),
            index=["a", "b", "c", "d", "e"],
            columns=["A", "B", "C", "D", "E"],
        )
        z = df[["A", "C", "A"]]
        result = z.loc[["a", "c", "a"]]
        check(result, expected)

    def test_columns_with_dups(self):
        # GH 3468 related

        # basic
        df = DataFrame([[1, 2]], columns=["a", "a"])
        df.columns = ["a", "a.1"]
        str(df)
        expected = DataFrame([[1, 2]], columns=["a", "a.1"])
        tm.assert_frame_equal(df, expected)

        df = DataFrame([[1, 2, 3]], columns=["b", "a", "a"])
        df.columns = ["b", "a", "a.1"]
        str(df)
        expected = DataFrame([[1, 2, 3]], columns=["b", "a", "a.1"])
        tm.assert_frame_equal(df, expected)

        # with a dup index
        df = DataFrame([[1, 2]], columns=["a", "a"])
        df.columns = ["b", "b"]
        str(df)
        expected = DataFrame([[1, 2]], columns=["b", "b"])
        tm.assert_frame_equal(df, expected)

        # multi-dtype
        df = DataFrame(
            [[1, 2, 1.0, 2.0, 3.0, "foo", "bar"]],
            columns=["a", "a", "b", "b", "d", "c", "c"],
        )
        df.columns = list("ABCDEFG")
        str(df)
        expected = DataFrame(
            [[1, 2, 1.0, 2.0, 3.0, "foo", "bar"]], columns=list("ABCDEFG")
        )
        tm.assert_frame_equal(df, expected)

        df = DataFrame([[1, 2, "foo", "bar"]], columns=["a", "a", "a", "a"])
        df.columns = ["a", "a.1", "a.2", "a.3"]
        str(df)
        expected = DataFrame([[1, 2, "foo", "bar"]], columns=["a", "a.1", "a.2", "a.3"])
        tm.assert_frame_equal(df, expected)

        # dups across blocks
        df_float = DataFrame(np.random.randn(10, 3), dtype="float64")
        df_int = DataFrame(np.random.randn(10, 3), dtype="int64")
        df_bool = DataFrame(True, index=df_float.index, columns=df_float.columns)
        df_object = DataFrame("foo", index=df_float.index, columns=df_float.columns)
        df_dt = DataFrame(
            pd.Timestamp("20010101"), index=df_float.index, columns=df_float.columns
        )
        df = pd.concat([df_float, df_int, df_bool, df_object, df_dt], axis=1)

        assert len(df._mgr.blknos) == len(df.columns)
        assert len(df._mgr.blklocs) == len(df.columns)

        # testing iloc
        for i in range(len(df.columns)):
            df.iloc[:, i]

        # dup columns across dtype GH 2079/2194
        vals = [[1, -1, 2.0], [2, -2, 3.0]]
        rs = DataFrame(vals, columns=["A", "A", "B"])
        xp = DataFrame(vals)
        xp.columns = ["A", "A", "B"]
        tm.assert_frame_equal(rs, xp)

    def test_set_value_by_index(self):
        # See gh-12344
        df = DataFrame(np.arange(9).reshape(3, 3).T)
        df.columns = list("AAA")
        expected = df.iloc[:, 2]

        df.iloc[:, 0] = 3
        tm.assert_series_equal(df.iloc[:, 2], expected)

        df = DataFrame(np.arange(9).reshape(3, 3).T)
        df.columns = [2, float(2), str(2)]
        expected = df.iloc[:, 1]

        df.iloc[:, 0] = 3
        tm.assert_series_equal(df.iloc[:, 1], expected)
