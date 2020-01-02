""" test positional based indexing with iloc """

from datetime import datetime
from warnings import catch_warnings, simplefilter

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series, concat, date_range, isna
from pandas.api.types import is_scalar
from pandas.core.indexing import IndexingError
from pandas.tests.indexing.common import Base
import pandas.util.testing as tm


class TestiLoc(Base):
    def test_iloc_exceeds_bounds(self):

        # GH6296
        # iloc should allow indexers that exceed the bounds
        df = DataFrame(np.random.random_sample((20, 5)), columns=list("ABCDE"))

        # lists of positions should raise IndexError!
        msg = "positional indexers are out-of-bounds"
        with pytest.raises(IndexError, match=msg):
            df.iloc[:, [0, 1, 2, 3, 4, 5]]
        with pytest.raises(IndexError, match=msg):
            df.iloc[[1, 30]]
        with pytest.raises(IndexError, match=msg):
            df.iloc[[1, -30]]
        with pytest.raises(IndexError, match=msg):
            df.iloc[[100]]

        s = df["A"]
        with pytest.raises(IndexError, match=msg):
            s.iloc[[100]]
        with pytest.raises(IndexError, match=msg):
            s.iloc[[-100]]

        # still raise on a single indexer
        msg = "single positional indexer is out-of-bounds"
        with pytest.raises(IndexError, match=msg):
            df.iloc[30]
        with pytest.raises(IndexError, match=msg):
            df.iloc[-30]

        # GH10779
        # single positive/negative indexer exceeding Series bounds should raise
        # an IndexError
        with pytest.raises(IndexError, match=msg):
            s.iloc[30]
        with pytest.raises(IndexError, match=msg):
            s.iloc[-30]

        # slices are ok
        result = df.iloc[:, 4:10]  # 0 < start < len < stop
        expected = df.iloc[:, 4:]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, -4:-10]  # stop < 0 < start < len
        expected = df.iloc[:, :0]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, 10:4:-1]  # 0 < stop < len < start (down)
        expected = df.iloc[:, :4:-1]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, 4:-10:-1]  # stop < 0 < start < len (down)
        expected = df.iloc[:, 4::-1]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, -10:4]  # start < 0 < stop < len
        expected = df.iloc[:, :4]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, 10:4]  # 0 < stop < len < start
        expected = df.iloc[:, :0]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, -10:-11:-1]  # stop < start < 0 < len (down)
        expected = df.iloc[:, :0]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, 10:11]  # 0 < len < start < stop
        expected = df.iloc[:, :0]
        tm.assert_frame_equal(result, expected)

        # slice bounds exceeding is ok
        result = s.iloc[18:30]
        expected = s.iloc[18:]
        tm.assert_series_equal(result, expected)

        result = s.iloc[30:]
        expected = s.iloc[:0]
        tm.assert_series_equal(result, expected)

        result = s.iloc[30::-1]
        expected = s.iloc[::-1]
        tm.assert_series_equal(result, expected)

        # doc example
        def check(result, expected):
            str(result)
            result.dtypes
            tm.assert_frame_equal(result, expected)

        dfl = DataFrame(np.random.randn(5, 2), columns=list("AB"))
        check(dfl.iloc[:, 2:3], DataFrame(index=dfl.index))
        check(dfl.iloc[:, 1:3], dfl.iloc[:, [1]])
        check(dfl.iloc[4:6], dfl.iloc[[4]])

        msg = "positional indexers are out-of-bounds"
        with pytest.raises(IndexError, match=msg):
            dfl.iloc[[4, 5, 6]]
        msg = "single positional indexer is out-of-bounds"
        with pytest.raises(IndexError, match=msg):
            dfl.iloc[:, 4]

    @pytest.mark.parametrize("index,columns", [(np.arange(20), list("ABCDE"))])
    @pytest.mark.parametrize(
        "index_vals,column_vals",
        [
            ([slice(None), ["A", "D"]]),
            (["1", "2"], slice(None)),
            ([datetime(2019, 1, 1)], slice(None)),
        ],
    )
    def test_iloc_non_integer_raises(self, index, columns, index_vals, column_vals):
        # GH 25753
        df = DataFrame(
            np.random.randn(len(index), len(columns)), index=index, columns=columns
        )
        msg = ".iloc requires numeric indexers, got"
        with pytest.raises(IndexError, match=msg):
            df.iloc[index_vals, column_vals]

    def test_iloc_getitem_int(self):
        # integer
        self.check_result(
            "iloc",
            2,
            "iloc",
            2,
            typs=["labels", "mixed", "ts", "floats", "empty"],
            fails=IndexError,
        )

    def test_iloc_getitem_neg_int(self):
        # neg integer
        self.check_result(
            "iloc",
            -1,
            "iloc",
            -1,
            typs=["labels", "mixed", "ts", "floats", "empty"],
            fails=IndexError,
        )

    @pytest.mark.parametrize("dims", [1, 2])
    def test_iloc_getitem_invalid_scalar(self, dims):
        # GH 21982

        if dims == 1:
            s = Series(np.arange(10))
        else:
            s = DataFrame(np.arange(100).reshape(10, 10))

        with pytest.raises(TypeError, match="Cannot index by location index"):
            s.iloc["a"]

    def test_iloc_array_not_mutating_negative_indices(self):

        # GH 21867
        array_with_neg_numbers = np.array([1, 2, -1])
        array_copy = array_with_neg_numbers.copy()
        df = pd.DataFrame(
            {"A": [100, 101, 102], "B": [103, 104, 105], "C": [106, 107, 108]},
            index=[1, 2, 3],
        )
        df.iloc[array_with_neg_numbers]
        tm.assert_numpy_array_equal(array_with_neg_numbers, array_copy)
        df.iloc[:, array_with_neg_numbers]
        tm.assert_numpy_array_equal(array_with_neg_numbers, array_copy)

    def test_iloc_getitem_list_int(self):
        self.check_result(
            "iloc",
            [0, 1, 2],
            "iloc",
            [0, 1, 2],
            typs=["labels", "mixed", "ts", "floats", "empty"],
            fails=IndexError,
        )

        # array of ints (GH5006), make sure that a single indexer is returning
        # the correct type

    def test_iloc_getitem_neg_int_can_reach_first_index(self):
        # GH10547 and GH10779
        # negative integers should be able to reach index 0
        df = DataFrame({"A": [2, 3, 5], "B": [7, 11, 13]})
        s = df["A"]

        expected = df.iloc[0]
        result = df.iloc[-3]
        tm.assert_series_equal(result, expected)

        expected = df.iloc[[0]]
        result = df.iloc[[-3]]
        tm.assert_frame_equal(result, expected)

        expected = s.iloc[0]
        result = s.iloc[-3]
        assert result == expected

        expected = s.iloc[[0]]
        result = s.iloc[[-3]]
        tm.assert_series_equal(result, expected)

        # check the length 1 Series case highlighted in GH10547
        expected = Series(["a"], index=["A"])
        result = expected.iloc[[-1]]
        tm.assert_series_equal(result, expected)

    def test_iloc_getitem_dups(self):
        # GH 6766
        df1 = DataFrame([{"A": None, "B": 1}, {"A": 2, "B": 2}])
        df2 = DataFrame([{"A": 3, "B": 3}, {"A": 4, "B": 4}])
        df = concat([df1, df2], axis=1)

        # cross-sectional indexing
        result = df.iloc[0, 0]
        assert isna(result)

        result = df.iloc[0, :]
        expected = Series([np.nan, 1, 3, 3], index=["A", "B", "A", "B"], name=0)
        tm.assert_series_equal(result, expected)

    def test_iloc_getitem_array(self):
        # TODO: test something here?
        pass

    def test_iloc_getitem_bool(self):
        # TODO: test something here?
        pass

    @pytest.mark.parametrize("index", [[True, False], [True, False, True, False]])
    def test_iloc_getitem_bool_diff_len(self, index):
        # GH26658
        s = Series([1, 2, 3])
        with pytest.raises(
            IndexError,
            match=("Item wrong length {} instead of {}.".format(len(index), len(s))),
        ):
            _ = s.iloc[index]

    def test_iloc_getitem_slice(self):
        # TODO: test something here?
        pass

    def test_iloc_getitem_slice_dups(self):

        df1 = DataFrame(np.random.randn(10, 4), columns=["A", "A", "B", "B"])
        df2 = DataFrame(
            np.random.randint(0, 10, size=20).reshape(10, 2), columns=["A", "C"]
        )

        # axis=1
        df = concat([df1, df2], axis=1)
        tm.assert_frame_equal(df.iloc[:, :4], df1)
        tm.assert_frame_equal(df.iloc[:, 4:], df2)

        df = concat([df2, df1], axis=1)
        tm.assert_frame_equal(df.iloc[:, :2], df2)
        tm.assert_frame_equal(df.iloc[:, 2:], df1)

        exp = concat([df2, df1.iloc[:, [0]]], axis=1)
        tm.assert_frame_equal(df.iloc[:, 0:3], exp)

        # axis=0
        df = concat([df, df], axis=0)
        tm.assert_frame_equal(df.iloc[0:10, :2], df2)
        tm.assert_frame_equal(df.iloc[0:10, 2:], df1)
        tm.assert_frame_equal(df.iloc[10:, :2], df2)
        tm.assert_frame_equal(df.iloc[10:, 2:], df1)

    def test_iloc_setitem(self):
        df = self.frame_ints

        df.iloc[1, 1] = 1
        result = df.iloc[1, 1]
        assert result == 1

        df.iloc[:, 2:3] = 0
        expected = df.iloc[:, 2:3]
        result = df.iloc[:, 2:3]
        tm.assert_frame_equal(result, expected)

        # GH5771
        s = Series(0, index=[4, 5, 6])
        s.iloc[1:2] += 1
        expected = Series([0, 1, 0], index=[4, 5, 6])
        tm.assert_series_equal(s, expected)

    def test_iloc_setitem_list(self):

        # setitem with an iloc list
        df = DataFrame(
            np.arange(9).reshape((3, 3)), index=["A", "B", "C"], columns=["A", "B", "C"]
        )
        df.iloc[[0, 1], [1, 2]]
        df.iloc[[0, 1], [1, 2]] += 100

        expected = DataFrame(
            np.array([0, 101, 102, 3, 104, 105, 6, 7, 8]).reshape((3, 3)),
            index=["A", "B", "C"],
            columns=["A", "B", "C"],
        )
        tm.assert_frame_equal(df, expected)

    def test_iloc_setitem_pandas_object(self):
        # GH 17193
        s_orig = Series([0, 1, 2, 3])
        expected = Series([0, -1, -2, 3])

        s = s_orig.copy()
        s.iloc[Series([1, 2])] = [-1, -2]
        tm.assert_series_equal(s, expected)

        s = s_orig.copy()
        s.iloc[pd.Index([1, 2])] = [-1, -2]
        tm.assert_series_equal(s, expected)

    def test_iloc_setitem_dups(self):

        # GH 6766
        # iloc with a mask aligning from another iloc
        df1 = DataFrame([{"A": None, "B": 1}, {"A": 2, "B": 2}])
        df2 = DataFrame([{"A": 3, "B": 3}, {"A": 4, "B": 4}])
        df = concat([df1, df2], axis=1)

        expected = df.fillna(3)
        expected["A"] = expected["A"].astype("float64")
        inds = np.isnan(df.iloc[:, 0])
        mask = inds[inds].index
        df.iloc[mask, 0] = df.iloc[mask, 2]
        tm.assert_frame_equal(df, expected)

        # del a dup column across blocks
        expected = DataFrame({0: [1, 2], 1: [3, 4]})
        expected.columns = ["B", "B"]
        del df["A"]
        tm.assert_frame_equal(df, expected)

        # assign back to self
        df.iloc[[0, 1], [0, 1]] = df.iloc[[0, 1], [0, 1]]
        tm.assert_frame_equal(df, expected)

        # reversed x 2
        df.iloc[[1, 0], [0, 1]] = df.iloc[[1, 0], [0, 1]].reset_index(drop=True)
        df.iloc[[1, 0], [0, 1]] = df.iloc[[1, 0], [0, 1]].reset_index(drop=True)
        tm.assert_frame_equal(df, expected)

    # TODO: GH#27620 this test used to compare iloc against ix; check if this
    #  is redundant with another test comparing iloc against loc
    def test_iloc_getitem_frame(self):
        df = DataFrame(
            np.random.randn(10, 4), index=range(0, 20, 2), columns=range(0, 8, 2)
        )

        result = df.iloc[2]
        exp = df.loc[4]
        tm.assert_series_equal(result, exp)

        result = df.iloc[2, 2]
        exp = df.loc[4, 4]
        assert result == exp

        # slice
        result = df.iloc[4:8]
        expected = df.loc[8:14]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[:, 2:3]
        expected = df.loc[:, 4:5]
        tm.assert_frame_equal(result, expected)

        # list of integers
        result = df.iloc[[0, 1, 3]]
        expected = df.loc[[0, 2, 6]]
        tm.assert_frame_equal(result, expected)

        result = df.iloc[[0, 1, 3], [0, 1]]
        expected = df.loc[[0, 2, 6], [0, 2]]
        tm.assert_frame_equal(result, expected)

        # neg indices
        result = df.iloc[[-1, 1, 3], [-1, 1]]
        expected = df.loc[[18, 2, 6], [6, 2]]
        tm.assert_frame_equal(result, expected)

        # dups indices
        result = df.iloc[[-1, -1, 1, 3], [-1, 1]]
        expected = df.loc[[18, 18, 2, 6], [6, 2]]
        tm.assert_frame_equal(result, expected)

        # with index-like
        s = Series(index=range(1, 5), dtype=object)
        result = df.iloc[s.index]
        expected = df.loc[[2, 4, 6, 8]]
        tm.assert_frame_equal(result, expected)

    def test_iloc_getitem_labelled_frame(self):
        # try with labelled frame
        df = DataFrame(
            np.random.randn(10, 4), index=list("abcdefghij"), columns=list("ABCD")
        )

        result = df.iloc[1, 1]
        exp = df.loc["b", "B"]
        assert result == exp

        result = df.iloc[:, 2:3]
        expected = df.loc[:, ["C"]]
        tm.assert_frame_equal(result, expected)

        # negative indexing
        result = df.iloc[-1, -1]
        exp = df.loc["j", "D"]
        assert result == exp

        # out-of-bounds exception
        msg = "single positional indexer is out-of-bounds"
        with pytest.raises(IndexError, match=msg):
            df.iloc[10, 5]

        # trying to use a label
        msg = (
            r"Location based indexing can only have \[integer, integer"
            r" slice \(START point is INCLUDED, END point is EXCLUDED\),"
            r" listlike of integers, boolean array\] types"
        )
        with pytest.raises(ValueError, match=msg):
            df.iloc["j", "D"]

    def test_iloc_getitem_doc_issue(self):

        # multi axis slicing issue with single block
        # surfaced in GH 6059

        arr = np.random.randn(6, 4)
        index = date_range("20130101", periods=6)
        columns = list("ABCD")
        df = DataFrame(arr, index=index, columns=columns)

        # defines ref_locs
        df.describe()

        result = df.iloc[3:5, 0:2]
        str(result)
        result.dtypes

        expected = DataFrame(arr[3:5, 0:2], index=index[3:5], columns=columns[0:2])
        tm.assert_frame_equal(result, expected)

        # for dups
        df.columns = list("aaaa")
        result = df.iloc[3:5, 0:2]
        str(result)
        result.dtypes

        expected = DataFrame(arr[3:5, 0:2], index=index[3:5], columns=list("aa"))
        tm.assert_frame_equal(result, expected)

        # related
        arr = np.random.randn(6, 4)
        index = list(range(0, 12, 2))
        columns = list(range(0, 8, 2))
        df = DataFrame(arr, index=index, columns=columns)

        df._data.blocks[0].mgr_locs
        result = df.iloc[1:5, 2:4]
        str(result)
        result.dtypes
        expected = DataFrame(arr[1:5, 2:4], index=index[1:5], columns=columns[2:4])
        tm.assert_frame_equal(result, expected)

    def test_iloc_setitem_series(self):
        df = DataFrame(
            np.random.randn(10, 4), index=list("abcdefghij"), columns=list("ABCD")
        )

        df.iloc[1, 1] = 1
        result = df.iloc[1, 1]
        assert result == 1

        df.iloc[:, 2:3] = 0
        expected = df.iloc[:, 2:3]
        result = df.iloc[:, 2:3]
        tm.assert_frame_equal(result, expected)

        s = Series(np.random.randn(10), index=range(0, 20, 2))

        s.iloc[1] = 1
        result = s.iloc[1]
        assert result == 1

        s.iloc[:4] = 0
        expected = s.iloc[:4]
        result = s.iloc[:4]
        tm.assert_series_equal(result, expected)

        s = Series([-1] * 6)
        s.iloc[0::2] = [0, 2, 4]
        s.iloc[1::2] = [1, 3, 5]
        result = s
        expected = Series([0, 1, 2, 3, 4, 5])
        tm.assert_series_equal(result, expected)

    def test_iloc_setitem_list_of_lists(self):

        # GH 7551
        # list-of-list is set incorrectly in mixed vs. single dtyped frames
        df = DataFrame(
            dict(A=np.arange(5, dtype="int64"), B=np.arange(5, 10, dtype="int64"))
        )
        df.iloc[2:4] = [[10, 11], [12, 13]]
        expected = DataFrame(dict(A=[0, 1, 10, 12, 4], B=[5, 6, 11, 13, 9]))
        tm.assert_frame_equal(df, expected)

        df = DataFrame(dict(A=list("abcde"), B=np.arange(5, 10, dtype="int64")))
        df.iloc[2:4] = [["x", 11], ["y", 13]]
        expected = DataFrame(dict(A=["a", "b", "x", "y", "e"], B=[5, 6, 11, 13, 9]))
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize("indexer", [[0], slice(None, 1, None), np.array([0])])
    @pytest.mark.parametrize("value", [["Z"], np.array(["Z"])])
    def test_iloc_setitem_with_scalar_index(self, indexer, value):
        # GH #19474
        # assigning like "df.iloc[0, [0]] = ['Z']" should be evaluated
        # elementwisely, not using "setter('A', ['Z'])".

        df = pd.DataFrame([[1, 2], [3, 4]], columns=["A", "B"])
        df.iloc[0, indexer] = value
        result = df.iloc[0, 0]

        assert is_scalar(result) and result == "Z"

    def test_iloc_mask(self):

        # GH 3631, iloc with a mask (of a series) should raise
        df = DataFrame(list(range(5)), index=list("ABCDE"), columns=["a"])
        mask = df.a % 2 == 0
        msg = "iLocation based boolean indexing cannot use an indexable as a mask"
        with pytest.raises(ValueError, match=msg):
            df.iloc[mask]
        mask.index = range(len(mask))
        msg = "iLocation based boolean indexing on an integer type is not available"
        with pytest.raises(NotImplementedError, match=msg):
            df.iloc[mask]

        # ndarray ok
        result = df.iloc[np.array([True] * len(mask), dtype=bool)]
        tm.assert_frame_equal(result, df)

        # the possibilities
        locs = np.arange(4)
        nums = 2 ** locs
        reps = [bin(num) for num in nums]
        df = DataFrame({"locs": locs, "nums": nums}, reps)

        expected = {
            (None, ""): "0b1100",
            (None, ".loc"): "0b1100",
            (None, ".iloc"): "0b1100",
            ("index", ""): "0b11",
            ("index", ".loc"): "0b11",
            ("index", ".iloc"): (
                "iLocation based boolean indexing cannot use an indexable as a mask"
            ),
            ("locs", ""): "Unalignable boolean Series provided as indexer "
            "(index of the boolean Series and of the indexed "
            "object do not match).",
            ("locs", ".loc"): "Unalignable boolean Series provided as indexer "
            "(index of the boolean Series and of the "
            "indexed object do not match).",
            ("locs", ".iloc"): (
                "iLocation based boolean indexing on an "
                "integer type is not available"
            ),
        }

        # UserWarnings from reindex of a boolean mask
        with catch_warnings(record=True):
            simplefilter("ignore", UserWarning)
            result = dict()
            for idx in [None, "index", "locs"]:
                mask = (df.nums > 2).values
                if idx:
                    mask = Series(mask, list(reversed(getattr(df, idx))))
                for method in ["", ".loc", ".iloc"]:
                    try:
                        if method:
                            accessor = getattr(df, method[1:])
                        else:
                            accessor = df
                        ans = str(bin(accessor[mask]["nums"].sum()))
                    except (ValueError, IndexingError, NotImplementedError) as e:
                        ans = str(e)

                    key = tuple([idx, method])
                    r = expected.get(key)
                    if r != ans:
                        raise AssertionError(
                            "[{key}] does not match [{ans}], received [{r}]".format(
                                key=key, ans=ans, r=r
                            )
                        )

    def test_iloc_non_unique_indexing(self):

        # GH 4017, non-unique indexing (on the axis)
        df = DataFrame({"A": [0.1] * 3000, "B": [1] * 3000})
        idx = np.arange(30) * 99
        expected = df.iloc[idx]

        df3 = concat([df, 2 * df, 3 * df])
        result = df3.iloc[idx]

        tm.assert_frame_equal(result, expected)

        df2 = DataFrame({"A": [0.1] * 1000, "B": [1] * 1000})
        df2 = concat([df2, 2 * df2, 3 * df2])

        with pytest.raises(KeyError, match="with any missing labels"):
            df2.loc[idx]

    def test_iloc_empty_list_indexer_is_ok(self):

        df = tm.makeCustomDataframe(5, 2)
        # vertical empty
        tm.assert_frame_equal(
            df.iloc[:, []],
            df.iloc[:, :0],
            check_index_type=True,
            check_column_type=True,
        )
        # horizontal empty
        tm.assert_frame_equal(
            df.iloc[[], :],
            df.iloc[:0, :],
            check_index_type=True,
            check_column_type=True,
        )
        # horizontal empty
        tm.assert_frame_equal(
            df.iloc[[]], df.iloc[:0, :], check_index_type=True, check_column_type=True
        )

    def test_identity_slice_returns_new_object(self):
        # GH13873
        original_df = DataFrame({"a": [1, 2, 3]})
        sliced_df = original_df.iloc[:]
        assert sliced_df is not original_df

        # should be a shallow copy
        original_df["a"] = [4, 4, 4]
        assert (sliced_df["a"] == 4).all()

        original_series = Series([1, 2, 3, 4, 5, 6])
        sliced_series = original_series.iloc[:]
        assert sliced_series is not original_series

        # should also be a shallow copy
        original_series[:3] = [7, 8, 9]
        assert all(sliced_series[:3] == [7, 8, 9])

    def test_indexing_zerodim_np_array(self):
        # GH24919
        df = DataFrame([[1, 2], [3, 4]])
        result = df.iloc[np.array(0)]
        s = pd.Series([1, 2], name=0)
        tm.assert_series_equal(result, s)

    def test_series_indexing_zerodim_np_array(self):
        # GH24919
        s = Series([1, 2])
        result = s.iloc[np.array(0)]
        assert result == 1
