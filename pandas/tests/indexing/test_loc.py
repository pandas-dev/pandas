""" test label based indexing with loc """
from datetime import datetime, time, timedelta
from io import StringIO
import re

from dateutil.tz import gettz
import numpy as np
import pytest

import pandas.util._test_decorators as td

import pandas as pd
from pandas import (
    Categorical,
    CategoricalIndex,
    DataFrame,
    Index,
    MultiIndex,
    Series,
    SparseDtype,
    Timedelta,
    Timestamp,
    date_range,
    timedelta_range,
    to_datetime,
    to_timedelta,
)
import pandas._testing as tm
from pandas.api.types import is_scalar
from pandas.tests.indexing.common import Base


class TestLoc(Base):
    def test_loc_getitem_int(self):

        # int label
        self.check_result("loc", 2, typs=["labels"], fails=KeyError)

    def test_loc_getitem_label(self):

        # label
        self.check_result("loc", "c", typs=["empty"], fails=KeyError)

    def test_loc_getitem_label_out_of_range(self):

        # out of range label
        self.check_result(
            "loc", "f", typs=["ints", "uints", "labels", "mixed", "ts"], fails=KeyError
        )
        self.check_result("loc", "f", typs=["floats"], fails=KeyError)
        self.check_result("loc", "f", typs=["floats"], fails=KeyError)
        self.check_result("loc", 20, typs=["ints", "uints", "mixed"], fails=KeyError)
        self.check_result("loc", 20, typs=["labels"], fails=KeyError)
        self.check_result("loc", 20, typs=["ts"], axes=0, fails=KeyError)
        self.check_result("loc", 20, typs=["floats"], axes=0, fails=KeyError)

    def test_loc_getitem_label_list(self):
        # TODO: test something here?
        # list of labels
        pass

    def test_loc_getitem_label_list_with_missing(self):
        self.check_result("loc", [0, 1, 2], typs=["empty"], fails=KeyError)
        self.check_result(
            "loc", [0, 2, 10], typs=["ints", "uints", "floats"], axes=0, fails=KeyError
        )

        self.check_result(
            "loc", [3, 6, 7], typs=["ints", "uints", "floats"], axes=1, fails=KeyError
        )

        # GH 17758 - MultiIndex and missing keys
        self.check_result(
            "loc", [(1, 3), (1, 4), (2, 5)], typs=["multi"], axes=0, fails=KeyError
        )

    def test_loc_getitem_label_list_fails(self):
        # fails
        self.check_result(
            "loc", [20, 30, 40], typs=["ints", "uints"], axes=1, fails=KeyError
        )

    def test_loc_getitem_label_array_like(self):
        # TODO: test something?
        # array like
        pass

    def test_loc_getitem_bool(self):
        # boolean indexers
        b = [True, False, True, False]

        self.check_result("loc", b, typs=["empty"], fails=IndexError)

    def test_loc_getitem_label_slice(self):

        # label slices (with ints)

        # real label slices

        # GH 14316

        self.check_result(
            "loc",
            slice(1, 3),
            typs=["labels", "mixed", "empty", "ts", "floats"],
            fails=TypeError,
        )

        self.check_result(
            "loc", slice("20130102", "20130104"), typs=["ts"], axes=1, fails=TypeError
        )

        self.check_result("loc", slice(2, 8), typs=["mixed"], axes=0, fails=TypeError)
        self.check_result("loc", slice(2, 8), typs=["mixed"], axes=1, fails=KeyError)

        self.check_result(
            "loc", slice(2, 4, 2), typs=["mixed"], axes=0, fails=TypeError
        )

    def test_setitem_from_duplicate_axis(self):
        # GH#34034
        df = DataFrame(
            [[20, "a"], [200, "a"], [200, "a"]],
            columns=["col1", "col2"],
            index=[10, 1, 1],
        )
        df.loc[1, "col1"] = np.arange(2)
        expected = DataFrame(
            [[20, "a"], [0, "a"], [1, "a"]], columns=["col1", "col2"], index=[10, 1, 1]
        )
        tm.assert_frame_equal(df, expected)


class TestLoc2:
    # TODO: better name, just separating out things that rely on base class

    def test_loc_getitem_missing_unicode_key(self):
        df = DataFrame({"a": [1]})
        with pytest.raises(KeyError, match="\u05d0"):
            df.loc[:, "\u05d0"]  # should not raise UnicodeEncodeError

    def test_loc_getitem_dups(self):
        # GH 5678
        # repeated getitems on a dup index returning a ndarray
        df = DataFrame(
            np.random.random_sample((20, 5)), index=["ABCDE"[x % 5] for x in range(20)]
        )
        expected = df.loc["A", 0]
        result = df.loc[:, 0].loc["A"]
        tm.assert_series_equal(result, expected)

    def test_loc_getitem_dups2(self):

        # GH4726
        # dup indexing with iloc/loc
        df = DataFrame(
            [[1, 2, "foo", "bar", Timestamp("20130101")]],
            columns=["a", "a", "a", "a", "a"],
            index=[1],
        )
        expected = Series(
            [1, 2, "foo", "bar", Timestamp("20130101")],
            index=["a", "a", "a", "a", "a"],
            name=1,
        )

        result = df.iloc[0]
        tm.assert_series_equal(result, expected)

        result = df.loc[1]
        tm.assert_series_equal(result, expected)

    def test_loc_setitem_dups(self):

        # GH 6541
        df_orig = DataFrame(
            {
                "me": list("rttti"),
                "foo": list("aaade"),
                "bar": np.arange(5, dtype="float64") * 1.34 + 2,
                "bar2": np.arange(5, dtype="float64") * -0.34 + 2,
            }
        ).set_index("me")

        indexer = (
            "r",
            ["bar", "bar2"],
        )
        df = df_orig.copy()
        df.loc[indexer] *= 2.0
        tm.assert_series_equal(df.loc[indexer], 2.0 * df_orig.loc[indexer])

        indexer = (
            "r",
            "bar",
        )
        df = df_orig.copy()
        df.loc[indexer] *= 2.0
        assert df.loc[indexer] == 2.0 * df_orig.loc[indexer]

        indexer = (
            "t",
            ["bar", "bar2"],
        )
        df = df_orig.copy()
        df.loc[indexer] *= 2.0
        tm.assert_frame_equal(df.loc[indexer], 2.0 * df_orig.loc[indexer])

    def test_loc_setitem_slice(self):
        # GH10503

        # assigning the same type should not change the type
        df1 = DataFrame({"a": [0, 1, 1], "b": Series([100, 200, 300], dtype="uint32")})
        ix = df1["a"] == 1
        newb1 = df1.loc[ix, "b"] + 1
        df1.loc[ix, "b"] = newb1
        expected = DataFrame(
            {"a": [0, 1, 1], "b": Series([100, 201, 301], dtype="uint32")}
        )
        tm.assert_frame_equal(df1, expected)

        # assigning a new type should get the inferred type
        df2 = DataFrame({"a": [0, 1, 1], "b": [100, 200, 300]}, dtype="uint64")
        ix = df1["a"] == 1
        newb2 = df2.loc[ix, "b"]
        df1.loc[ix, "b"] = newb2
        expected = DataFrame({"a": [0, 1, 1], "b": [100, 200, 300]}, dtype="uint64")
        tm.assert_frame_equal(df2, expected)

    def test_loc_setitem_dtype(self):
        # GH31340
        df = DataFrame({"id": ["A"], "a": [1.2], "b": [0.0], "c": [-2.5]})
        cols = ["a", "b", "c"]
        df.loc[:, cols] = df.loc[:, cols].astype("float32")

        expected = DataFrame(
            {"id": ["A"], "a": [1.2], "b": [0.0], "c": [-2.5]}, dtype="float32"
        )  # id is inferred as object

        tm.assert_frame_equal(df, expected)

    def test_getitem_label_list_with_missing(self):
        s = Series(range(3), index=["a", "b", "c"])

        # consistency
        with pytest.raises(KeyError, match="with any missing labels"):
            s[["a", "d"]]

        s = Series(range(3))
        with pytest.raises(KeyError, match="with any missing labels"):
            s[[0, 3]]

    @pytest.mark.parametrize("index", [[True, False], [True, False, True, False]])
    def test_loc_getitem_bool_diff_len(self, index):
        # GH26658
        s = Series([1, 2, 3])
        msg = f"Boolean index has wrong length: {len(index)} instead of {len(s)}"
        with pytest.raises(IndexError, match=msg):
            _ = s.loc[index]

    def test_loc_getitem_int_slice(self):
        # TODO: test something here?
        pass

    def test_loc_to_fail(self):

        # GH3449
        df = DataFrame(
            np.random.random((3, 3)), index=["a", "b", "c"], columns=["e", "f", "g"]
        )

        # raise a KeyError?
        msg = (
            r"\"None of \[Int64Index\(\[1, 2\], dtype='int64'\)\] are "
            r"in the \[index\]\""
        )
        with pytest.raises(KeyError, match=msg):
            df.loc[[1, 2], [1, 2]]

        # GH  7496
        # loc should not fallback

        s = Series(dtype=object)
        s.loc[1] = 1
        s.loc["a"] = 2

        with pytest.raises(KeyError, match=r"^-1$"):
            s.loc[-1]

        msg = (
            r"\"None of \[Int64Index\(\[-1, -2\], dtype='int64'\)\] are "
            r"in the \[index\]\""
        )
        with pytest.raises(KeyError, match=msg):
            s.loc[[-1, -2]]

        msg = r"\"None of \[Index\(\['4'\], dtype='object'\)\] are in the \[index\]\""
        with pytest.raises(KeyError, match=msg):
            s.loc[["4"]]

        s.loc[-1] = 3
        with pytest.raises(KeyError, match="with any missing labels"):
            s.loc[[-1, -2]]

        s["a"] = 2
        msg = (
            r"\"None of \[Int64Index\(\[-2\], dtype='int64'\)\] are "
            r"in the \[index\]\""
        )
        with pytest.raises(KeyError, match=msg):
            s.loc[[-2]]

        del s["a"]

        with pytest.raises(KeyError, match=msg):
            s.loc[[-2]] = 0

        # inconsistency between .loc[values] and .loc[values,:]
        # GH 7999
        df = DataFrame([["a"], ["b"]], index=[1, 2], columns=["value"])

        msg = (
            r"\"None of \[Int64Index\(\[3\], dtype='int64'\)\] are "
            r"in the \[index\]\""
        )
        with pytest.raises(KeyError, match=msg):
            df.loc[[3], :]

        with pytest.raises(KeyError, match=msg):
            df.loc[[3]]

    def test_loc_getitem_list_with_fail(self):
        # 15747
        # should KeyError if *any* missing labels

        s = Series([1, 2, 3])

        s.loc[[2]]

        with pytest.raises(
            KeyError,
            match=re.escape(
                "\"None of [Int64Index([3], dtype='int64')] are in the [index]\""
            ),
        ):
            s.loc[[3]]

        # a non-match and a match
        with pytest.raises(KeyError, match="with any missing labels"):
            s.loc[[2, 3]]

    def test_loc_index(self):
        # gh-17131
        # a boolean index should index like a boolean numpy array

        df = DataFrame(
            np.random.random(size=(5, 10)),
            index=["alpha_0", "alpha_1", "alpha_2", "beta_0", "beta_1"],
        )

        mask = df.index.map(lambda x: "alpha" in x)
        expected = df.loc[np.array(mask)]

        result = df.loc[mask]
        tm.assert_frame_equal(result, expected)

        result = df.loc[mask.values]
        tm.assert_frame_equal(result, expected)

        result = df.loc[pd.array(mask, dtype="boolean")]
        tm.assert_frame_equal(result, expected)

    def test_loc_general(self):

        df = DataFrame(
            np.random.rand(4, 4),
            columns=["A", "B", "C", "D"],
            index=["A", "B", "C", "D"],
        )

        # want this to work
        result = df.loc[:, "A":"B"].iloc[0:2, :]
        assert (result.columns == ["A", "B"]).all()
        assert (result.index == ["A", "B"]).all()

        # mixed type
        result = DataFrame({"a": [Timestamp("20130101")], "b": [1]}).iloc[0]
        expected = Series([Timestamp("20130101"), 1], index=["a", "b"], name=0)
        tm.assert_series_equal(result, expected)
        assert result.dtype == object

    def test_loc_setitem_consistency(self):
        # GH 6149
        # coerce similarly for setitem and loc when rows have a null-slice
        expected = DataFrame(
            {
                "date": Series(0, index=range(5), dtype=np.int64),
                "val": Series(range(5), dtype=np.int64),
            }
        )

        df = DataFrame(
            {
                "date": date_range("2000-01-01", "2000-01-5"),
                "val": Series(range(5), dtype=np.int64),
            }
        )
        df.loc[:, "date"] = 0
        tm.assert_frame_equal(df, expected)

        df = DataFrame(
            {
                "date": date_range("2000-01-01", "2000-01-5"),
                "val": Series(range(5), dtype=np.int64),
            }
        )
        df.loc[:, "date"] = np.array(0, dtype=np.int64)
        tm.assert_frame_equal(df, expected)

        df = DataFrame(
            {
                "date": date_range("2000-01-01", "2000-01-5"),
                "val": Series(range(5), dtype=np.int64),
            }
        )
        df.loc[:, "date"] = np.array([0, 0, 0, 0, 0], dtype=np.int64)
        tm.assert_frame_equal(df, expected)

        expected = DataFrame(
            {
                "date": Series("foo", index=range(5)),
                "val": Series(range(5), dtype=np.int64),
            }
        )
        df = DataFrame(
            {
                "date": date_range("2000-01-01", "2000-01-5"),
                "val": Series(range(5), dtype=np.int64),
            }
        )
        df.loc[:, "date"] = "foo"
        tm.assert_frame_equal(df, expected)

        expected = DataFrame(
            {
                "date": Series(1.0, index=range(5)),
                "val": Series(range(5), dtype=np.int64),
            }
        )
        df = DataFrame(
            {
                "date": date_range("2000-01-01", "2000-01-5"),
                "val": Series(range(5), dtype=np.int64),
            }
        )
        df.loc[:, "date"] = 1.0
        tm.assert_frame_equal(df, expected)

        # GH 15494
        # setting on frame with single row
        df = DataFrame({"date": Series([Timestamp("20180101")])})
        df.loc[:, "date"] = "string"
        expected = DataFrame({"date": Series(["string"])})
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_consistency_empty(self):
        # empty (essentially noops)
        expected = DataFrame(columns=["x", "y"])
        expected["x"] = expected["x"].astype(np.int64)
        df = DataFrame(columns=["x", "y"])
        df.loc[:, "x"] = 1
        tm.assert_frame_equal(df, expected)

        df = DataFrame(columns=["x", "y"])
        df["x"] = 1
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_consistency_slice_column_len(self):
        # .loc[:,column] setting with slice == len of the column
        # GH10408
        data = """Level_0,,,Respondent,Respondent,Respondent,OtherCat,OtherCat
Level_1,,,Something,StartDate,EndDate,Yes/No,SomethingElse
Region,Site,RespondentID,,,,,
Region_1,Site_1,3987227376,A,5/25/2015 10:59,5/25/2015 11:22,Yes,
Region_1,Site_1,3980680971,A,5/21/2015 9:40,5/21/2015 9:52,Yes,Yes
Region_1,Site_2,3977723249,A,5/20/2015 8:27,5/20/2015 8:41,Yes,
Region_1,Site_2,3977723089,A,5/20/2015 8:33,5/20/2015 9:09,Yes,No"""

        df = pd.read_csv(StringIO(data), header=[0, 1], index_col=[0, 1, 2])
        df.loc[:, ("Respondent", "StartDate")] = pd.to_datetime(
            df.loc[:, ("Respondent", "StartDate")]
        )
        df.loc[:, ("Respondent", "EndDate")] = pd.to_datetime(
            df.loc[:, ("Respondent", "EndDate")]
        )
        df.loc[:, ("Respondent", "Duration")] = (
            df.loc[:, ("Respondent", "EndDate")]
            - df.loc[:, ("Respondent", "StartDate")]
        )

        df.loc[:, ("Respondent", "Duration")] = df.loc[
            :, ("Respondent", "Duration")
        ].astype("timedelta64[s]")
        expected = Series(
            [1380, 720, 840, 2160.0], index=df.index, name=("Respondent", "Duration")
        )
        tm.assert_series_equal(df[("Respondent", "Duration")], expected)

    @pytest.mark.parametrize("unit", ["Y", "M", "D", "h", "m", "s", "ms", "us"])
    def test_loc_assign_non_ns_datetime(self, unit):
        # GH 27395, non-ns dtype assignment via .loc should work
        # and return the same result when using simple assignment
        df = DataFrame(
            {
                "timestamp": [
                    np.datetime64("2017-02-11 12:41:29"),
                    np.datetime64("1991-11-07 04:22:37"),
                ]
            }
        )

        df.loc[:, unit] = df.loc[:, "timestamp"].values.astype(f"datetime64[{unit}]")
        df["expected"] = df.loc[:, "timestamp"].values.astype(f"datetime64[{unit}]")
        expected = Series(df.loc[:, "expected"], name=unit)
        tm.assert_series_equal(df.loc[:, unit], expected)

    def test_loc_modify_datetime(self):
        # see gh-28837
        df = DataFrame.from_dict(
            {"date": [1485264372711, 1485265925110, 1540215845888, 1540282121025]}
        )

        df["date_dt"] = pd.to_datetime(df["date"], unit="ms", cache=True)

        df.loc[:, "date_dt_cp"] = df.loc[:, "date_dt"]
        df.loc[[2, 3], "date_dt_cp"] = df.loc[[2, 3], "date_dt"]

        expected = DataFrame(
            [
                [1485264372711, "2017-01-24 13:26:12.711", "2017-01-24 13:26:12.711"],
                [1485265925110, "2017-01-24 13:52:05.110", "2017-01-24 13:52:05.110"],
                [1540215845888, "2018-10-22 13:44:05.888", "2018-10-22 13:44:05.888"],
                [1540282121025, "2018-10-23 08:08:41.025", "2018-10-23 08:08:41.025"],
            ],
            columns=["date", "date_dt", "date_dt_cp"],
        )

        columns = ["date_dt", "date_dt_cp"]
        expected[columns] = expected[columns].apply(pd.to_datetime)

        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_frame(self):
        df = DataFrame(np.random.randn(4, 4), index=list("abcd"), columns=list("ABCD"))

        result = df.iloc[0, 0]

        df.loc["a", "A"] = 1
        result = df.loc["a", "A"]
        assert result == 1

        result = df.iloc[0, 0]
        assert result == 1

        df.loc[:, "B":"D"] = 0
        expected = df.loc[:, "B":"D"]
        result = df.iloc[:, 1:]
        tm.assert_frame_equal(result, expected)

        # GH 6254
        # setting issue
        df = DataFrame(index=[3, 5, 4], columns=["A"])
        df.loc[[4, 3, 5], "A"] = np.array([1, 2, 3], dtype="int64")
        expected = DataFrame({"A": Series([1, 2, 3], index=[4, 3, 5])}).reindex(
            index=[3, 5, 4]
        )
        tm.assert_frame_equal(df, expected)

        # GH 6252
        # setting with an empty frame
        keys1 = ["@" + str(i) for i in range(5)]
        val1 = np.arange(5, dtype="int64")

        keys2 = ["@" + str(i) for i in range(4)]
        val2 = np.arange(4, dtype="int64")

        index = list(set(keys1).union(keys2))
        df = DataFrame(index=index)
        df["A"] = np.nan
        df.loc[keys1, "A"] = val1

        df["B"] = np.nan
        df.loc[keys2, "B"] = val2

        expected = DataFrame(
            {"A": Series(val1, index=keys1), "B": Series(val2, index=keys2)}
        ).reindex(index=index)
        tm.assert_frame_equal(df, expected)

        # GH 8669
        # invalid coercion of nan -> int
        df = DataFrame({"A": [1, 2, 3], "B": np.nan})
        df.loc[df.B > df.A, "B"] = df.A
        expected = DataFrame({"A": [1, 2, 3], "B": np.nan})
        tm.assert_frame_equal(df, expected)

        # GH 6546
        # setting with mixed labels
        df = DataFrame({1: [1, 2], 2: [3, 4], "a": ["a", "b"]})

        result = df.loc[0, [1, 2]]
        expected = Series([1, 3], index=[1, 2], dtype=object, name=0)
        tm.assert_series_equal(result, expected)

        expected = DataFrame({1: [5, 2], 2: [6, 4], "a": ["a", "b"]})
        df.loc[0, [1, 2]] = [5, 6]
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_frame_multiples(self):
        # multiple setting
        df = DataFrame(
            {"A": ["foo", "bar", "baz"], "B": Series(range(3), dtype=np.int64)}
        )
        rhs = df.loc[1:2]
        rhs.index = df.index[0:2]
        df.loc[0:1] = rhs
        expected = DataFrame(
            {"A": ["bar", "baz", "baz"], "B": Series([1, 2, 2], dtype=np.int64)}
        )
        tm.assert_frame_equal(df, expected)

        # multiple setting with frame on rhs (with M8)
        df = DataFrame(
            {
                "date": date_range("2000-01-01", "2000-01-5"),
                "val": Series(range(5), dtype=np.int64),
            }
        )
        expected = DataFrame(
            {
                "date": [
                    Timestamp("20000101"),
                    Timestamp("20000102"),
                    Timestamp("20000101"),
                    Timestamp("20000102"),
                    Timestamp("20000103"),
                ],
                "val": Series([0, 1, 0, 1, 2], dtype=np.int64),
            }
        )
        rhs = df.loc[0:2]
        rhs.index = df.index[2:5]
        df.loc[2:4] = rhs
        tm.assert_frame_equal(df, expected)

    @pytest.mark.parametrize(
        "indexer", [["A"], slice(None, "A", None), np.array(["A"])]
    )
    @pytest.mark.parametrize("value", [["Z"], np.array(["Z"])])
    def test_loc_setitem_with_scalar_index(self, indexer, value):
        # GH #19474
        # assigning like "df.loc[0, ['A']] = ['Z']" should be evaluated
        # elementwisely, not using "setter('A', ['Z'])".

        df = DataFrame([[1, 2], [3, 4]], columns=["A", "B"])
        df.loc[0, indexer] = value
        result = df.loc[0, "A"]

        assert is_scalar(result) and result == "Z"

    @pytest.mark.parametrize(
        "index,box,expected",
        [
            (
                ([0, 2], ["A", "B", "C", "D"]),
                7,
                DataFrame(
                    [[7, 7, 7, 7], [3, 4, np.nan, np.nan], [7, 7, 7, 7]],
                    columns=["A", "B", "C", "D"],
                ),
            ),
            (
                (1, ["C", "D"]),
                [7, 8],
                DataFrame(
                    [[1, 2, np.nan, np.nan], [3, 4, 7, 8], [5, 6, np.nan, np.nan]],
                    columns=["A", "B", "C", "D"],
                ),
            ),
            (
                (1, ["A", "B", "C"]),
                np.array([7, 8, 9], dtype=np.int64),
                DataFrame(
                    [[1, 2, np.nan], [7, 8, 9], [5, 6, np.nan]], columns=["A", "B", "C"]
                ),
            ),
            (
                (slice(1, 3, None), ["B", "C", "D"]),
                [[7, 8, 9], [10, 11, 12]],
                DataFrame(
                    [[1, 2, np.nan, np.nan], [3, 7, 8, 9], [5, 10, 11, 12]],
                    columns=["A", "B", "C", "D"],
                ),
            ),
            (
                (slice(1, 3, None), ["C", "A", "D"]),
                np.array([[7, 8, 9], [10, 11, 12]], dtype=np.int64),
                DataFrame(
                    [[1, 2, np.nan, np.nan], [8, 4, 7, 9], [11, 6, 10, 12]],
                    columns=["A", "B", "C", "D"],
                ),
            ),
            (
                (slice(None, None, None), ["A", "C"]),
                DataFrame([[7, 8], [9, 10], [11, 12]], columns=["A", "C"]),
                DataFrame(
                    [[7, 2, 8], [9, 4, 10], [11, 6, 12]], columns=["A", "B", "C"]
                ),
            ),
        ],
    )
    def test_loc_setitem_missing_columns(self, index, box, expected):
        # GH 29334
        df = DataFrame([[1, 2], [3, 4], [5, 6]], columns=["A", "B"])
        df.loc[index] = box
        tm.assert_frame_equal(df, expected)

    def test_loc_coercion(self):

        # 12411
        df = DataFrame({"date": [Timestamp("20130101").tz_localize("UTC"), pd.NaT]})
        expected = df.dtypes

        result = df.iloc[[0]]
        tm.assert_series_equal(result.dtypes, expected)

        result = df.iloc[[1]]
        tm.assert_series_equal(result.dtypes, expected)

        # 12045
        import datetime

        df = DataFrame(
            {"date": [datetime.datetime(2012, 1, 1), datetime.datetime(1012, 1, 2)]}
        )
        expected = df.dtypes

        result = df.iloc[[0]]
        tm.assert_series_equal(result.dtypes, expected)

        result = df.iloc[[1]]
        tm.assert_series_equal(result.dtypes, expected)

        # 11594
        df = DataFrame({"text": ["some words"] + [None] * 9})
        expected = df.dtypes

        result = df.iloc[0:2]
        tm.assert_series_equal(result.dtypes, expected)

        result = df.iloc[3:]
        tm.assert_series_equal(result.dtypes, expected)

    def test_setitem_new_key_tz(self):
        # GH#12862 should not raise on assigning the second value
        vals = [
            pd.to_datetime(42).tz_localize("UTC"),
            pd.to_datetime(666).tz_localize("UTC"),
        ]
        expected = Series(vals, index=["foo", "bar"])

        ser = Series(dtype=object)
        ser["foo"] = vals[0]
        ser["bar"] = vals[1]

        tm.assert_series_equal(ser, expected)

        ser = Series(dtype=object)
        ser.loc["foo"] = vals[0]
        ser.loc["bar"] = vals[1]

        tm.assert_series_equal(ser, expected)

    def test_loc_non_unique(self):
        # GH3659
        # non-unique indexer with loc slice
        # https://groups.google.com/forum/?fromgroups#!topic/pydata/zTm2No0crYs

        # these are going to raise because the we are non monotonic
        df = DataFrame(
            {"A": [1, 2, 3, 4, 5, 6], "B": [3, 4, 5, 6, 7, 8]}, index=[0, 1, 0, 1, 2, 3]
        )
        msg = "'Cannot get left slice bound for non-unique label: 1'"
        with pytest.raises(KeyError, match=msg):
            df.loc[1:]
        msg = "'Cannot get left slice bound for non-unique label: 0'"
        with pytest.raises(KeyError, match=msg):
            df.loc[0:]
        msg = "'Cannot get left slice bound for non-unique label: 1'"
        with pytest.raises(KeyError, match=msg):
            df.loc[1:2]

        # monotonic are ok
        df = DataFrame(
            {"A": [1, 2, 3, 4, 5, 6], "B": [3, 4, 5, 6, 7, 8]}, index=[0, 1, 0, 1, 2, 3]
        ).sort_index(axis=0)
        result = df.loc[1:]
        expected = DataFrame({"A": [2, 4, 5, 6], "B": [4, 6, 7, 8]}, index=[1, 1, 2, 3])
        tm.assert_frame_equal(result, expected)

        result = df.loc[0:]
        tm.assert_frame_equal(result, df)

        result = df.loc[1:2]
        expected = DataFrame({"A": [2, 4, 5], "B": [4, 6, 7]}, index=[1, 1, 2])
        tm.assert_frame_equal(result, expected)

    @pytest.mark.arm_slow
    def test_loc_non_unique_memory_error(self):

        # GH 4280
        # non_unique index with a large selection triggers a memory error

        columns = list("ABCDEFG")

        def gen_test(length, l2):
            return pd.concat(
                [
                    DataFrame(
                        np.random.randn(length, len(columns)),
                        index=np.arange(length),
                        columns=columns,
                    ),
                    DataFrame(
                        np.ones((l2, len(columns))), index=[0] * l2, columns=columns
                    ),
                ]
            )

        def gen_expected(df, mask):
            len_mask = len(mask)
            return pd.concat(
                [
                    df.take([0]),
                    DataFrame(
                        np.ones((len_mask, len(columns))),
                        index=[0] * len_mask,
                        columns=columns,
                    ),
                    df.take(mask[1:]),
                ]
            )

        df = gen_test(900, 100)
        assert df.index.is_unique is False

        mask = np.arange(100)
        result = df.loc[mask]
        expected = gen_expected(df, mask)
        tm.assert_frame_equal(result, expected)

        df = gen_test(900000, 100000)
        assert df.index.is_unique is False

        mask = np.arange(100000)
        result = df.loc[mask]
        expected = gen_expected(df, mask)
        tm.assert_frame_equal(result, expected)

    def test_loc_name(self):
        # GH 3880
        df = DataFrame([[1, 1], [1, 1]])
        df.index.name = "index_name"
        result = df.iloc[[0, 1]].index.name
        assert result == "index_name"

        result = df.loc[[0, 1]].index.name
        assert result == "index_name"

    def test_loc_empty_list_indexer_is_ok(self):

        df = tm.makeCustomDataframe(5, 2)
        # vertical empty
        tm.assert_frame_equal(
            df.loc[:, []], df.iloc[:, :0], check_index_type=True, check_column_type=True
        )
        # horizontal empty
        tm.assert_frame_equal(
            df.loc[[], :], df.iloc[:0, :], check_index_type=True, check_column_type=True
        )
        # horizontal empty
        tm.assert_frame_equal(
            df.loc[[]], df.iloc[:0, :], check_index_type=True, check_column_type=True
        )

    def test_identity_slice_returns_new_object(self):
        # GH13873
        original_df = DataFrame({"a": [1, 2, 3]})
        sliced_df = original_df.loc[:]
        assert sliced_df is not original_df
        assert original_df[:] is not original_df

        # should be a shallow copy
        original_df["a"] = [4, 4, 4]
        assert (sliced_df["a"] == 4).all()

        # These should not return copies
        assert original_df is original_df.loc[:, :]
        df = DataFrame(np.random.randn(10, 4))
        assert df[0] is df.loc[:, 0]

        # Same tests for Series
        original_series = Series([1, 2, 3, 4, 5, 6])
        sliced_series = original_series.loc[:]
        assert sliced_series is not original_series
        assert original_series[:] is not original_series

        original_series[:3] = [7, 8, 9]
        assert all(sliced_series[:3] == [7, 8, 9])

    @pytest.mark.xfail(reason="accidental fix reverted - GH37497")
    def test_loc_copy_vs_view(self):
        # GH 15631
        x = DataFrame(zip(range(3), range(3)), columns=["a", "b"])

        y = x.copy()
        q = y.loc[:, "a"]
        q += 2

        tm.assert_frame_equal(x, y)

        z = x.copy()
        q = z.loc[x.index, "a"]
        q += 2

        tm.assert_frame_equal(x, z)

    def test_loc_uint64(self):
        # GH20722
        # Test whether loc accept uint64 max value as index.
        s = Series([1, 2], index=[np.iinfo("uint64").max - 1, np.iinfo("uint64").max])

        result = s.loc[np.iinfo("uint64").max - 1]
        expected = s.iloc[0]
        assert result == expected

        result = s.loc[[np.iinfo("uint64").max - 1]]
        expected = s.iloc[[0]]
        tm.assert_series_equal(result, expected)

        result = s.loc[[np.iinfo("uint64").max - 1, np.iinfo("uint64").max]]
        tm.assert_series_equal(result, s)

    def test_loc_setitem_empty_append_expands_rows(self):
        # GH6173, various appends to an empty dataframe

        data = [1, 2, 3]
        expected = DataFrame({"x": data, "y": [None] * len(data)})

        # appends to fit length of data
        df = DataFrame(columns=["x", "y"])
        df.loc[:, "x"] = data
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_empty_append_expands_rows_mixed_dtype(self):
        # GH#37932 same as test_loc_setitem_empty_append_expands_rows
        #  but with mixed dtype so we go through take_split_path
        data = [1, 2, 3]
        expected = DataFrame({"x": data, "y": [None] * len(data)})

        df = DataFrame(columns=["x", "y"])
        df["x"] = df["x"].astype(np.int64)
        df.loc[:, "x"] = data
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_empty_append_single_value(self):
        # only appends one value
        expected = DataFrame({"x": [1.0], "y": [np.nan]})
        df = DataFrame(columns=["x", "y"], dtype=float)
        df.loc[0, "x"] = expected.loc[0, "x"]
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_empty_append_raises(self):
        # GH6173, various appends to an empty dataframe

        data = [1, 2]
        df = DataFrame(columns=["x", "y"])
        df.index = df.index.astype(np.int64)
        msg = (
            r"None of \[Int64Index\(\[0, 1\], dtype='int64'\)\] "
            r"are in the \[index\]"
        )
        with pytest.raises(KeyError, match=msg):
            df.loc[[0, 1], "x"] = data

        msg = "|".join(
            [
                "cannot copy sequence with size 2 to array axis with dimension 0",
                r"could not broadcast input array from shape \(2,\) into shape \(0,\)",
            ]
        )
        with pytest.raises(ValueError, match=msg):
            df.loc[0:2, "x"] = data

    def test_indexing_zerodim_np_array(self):
        # GH24924
        df = DataFrame([[1, 2], [3, 4]])
        result = df.loc[np.array(0)]
        s = Series([1, 2], name=0)
        tm.assert_series_equal(result, s)

    def test_series_indexing_zerodim_np_array(self):
        # GH24924
        s = Series([1, 2])
        result = s.loc[np.array(0)]
        assert result == 1

    def test_loc_reverse_assignment(self):
        # GH26939
        data = [1, 2, 3, 4, 5, 6] + [None] * 4
        expected = Series(data, index=range(2010, 2020))

        result = Series(index=range(2010, 2020), dtype=np.float64)
        result.loc[2015:2010:-1] = [6, 5, 4, 3, 2, 1]

        tm.assert_series_equal(result, expected)

    def test_loc_setitem_str_to_small_float_conversion_type(self):
        # GH#20388
        np.random.seed(13)
        col_data = [str(np.random.random() * 1e-12) for _ in range(5)]
        result = DataFrame(col_data, columns=["A"])
        expected = DataFrame(col_data, columns=["A"], dtype=object)
        tm.assert_frame_equal(result, expected)

        # change the dtype of the elements from object to float one by one
        result.loc[result.index, "A"] = [float(x) for x in col_data]
        expected = DataFrame(col_data, columns=["A"], dtype=float)
        tm.assert_frame_equal(result, expected)

    def test_loc_getitem_time_object(self, frame_or_series):
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        mask = (rng.hour == 9) & (rng.minute == 30)

        obj = DataFrame(np.random.randn(len(rng), 3), index=rng)
        if frame_or_series is Series:
            obj = obj[0]

        result = obj.loc[time(9, 30)]
        exp = obj.loc[mask]
        tm.assert_equal(result, exp)

        chunk = obj.loc["1/4/2000":]
        result = chunk.loc[time(9, 30)]
        expected = result[-1:]

        # Without resetting the freqs, these are 5 min and 1440 min, respectively
        result.index = result.index._with_freq(None)
        expected.index = expected.index._with_freq(None)
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize("spmatrix_t", ["coo_matrix", "csc_matrix", "csr_matrix"])
    @pytest.mark.parametrize("dtype", [np.int64, np.float64, complex])
    @td.skip_if_no_scipy
    def test_loc_getitem_range_from_spmatrix(self, spmatrix_t, dtype):
        import scipy.sparse

        spmatrix_t = getattr(scipy.sparse, spmatrix_t)

        # The bug is triggered by a sparse matrix with purely sparse columns.  So the
        # recipe below generates a rectangular matrix of dimension (5, 7) where all the
        # diagonal cells are ones, meaning the last two columns are purely sparse.
        rows, cols = 5, 7
        spmatrix = spmatrix_t(np.eye(rows, cols, dtype=dtype), dtype=dtype)
        df = DataFrame.sparse.from_spmatrix(spmatrix)

        # regression test for GH#34526
        itr_idx = range(2, rows)
        result = df.loc[itr_idx].values
        expected = spmatrix.toarray()[itr_idx]
        tm.assert_numpy_array_equal(result, expected)

        # regression test for GH#34540
        result = df.loc[itr_idx].dtypes.values
        expected = np.full(cols, SparseDtype(dtype, fill_value=0))
        tm.assert_numpy_array_equal(result, expected)

    def test_loc_getitem_listlike_all_retains_sparse(self):
        df = DataFrame({"A": pd.array([0, 0], dtype=SparseDtype("int64"))})
        result = df.loc[[0, 1]]
        tm.assert_frame_equal(result, df)

    @pytest.mark.parametrize("key_type", [iter, np.array, Series, Index])
    def test_loc_getitem_iterable(self, float_frame, key_type):
        idx = key_type(["A", "B", "C"])
        result = float_frame.loc[:, idx]
        expected = float_frame.loc[:, ["A", "B", "C"]]
        tm.assert_frame_equal(result, expected)

    def test_loc_getitem_timedelta_0seconds(self):
        # GH#10583
        df = DataFrame(np.random.normal(size=(10, 4)))
        df.index = timedelta_range(start="0s", periods=10, freq="s")
        expected = df.loc[Timedelta("0s") :, :]
        result = df.loc["0s":, :]
        tm.assert_frame_equal(expected, result)

    @pytest.mark.parametrize(
        "val,expected", [(2 ** 63 - 1, Series([1])), (2 ** 63, Series([2]))]
    )
    def test_loc_getitem_uint64_scalar(self, val, expected):
        # see GH#19399
        df = DataFrame([1, 2], index=[2 ** 63 - 1, 2 ** 63])
        result = df.loc[val]

        expected.name = val
        tm.assert_series_equal(result, expected)

    def test_loc_setitem_int_label_with_float64index(self):
        # note labels are floats
        ser = Series(["a", "b", "c"], index=[0, 0.5, 1])
        tmp = ser.copy()

        ser.loc[1] = "zoo"
        tmp.iloc[2] = "zoo"

        tm.assert_series_equal(ser, tmp)

    @pytest.mark.parametrize(
        "indexer, expected",
        [
            # The test name is a misnomer in the 0 case as df.index[indexer]
            #  is a scalar.
            (0, [20, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
            (slice(4, 8), [0, 1, 2, 3, 20, 20, 20, 20, 8, 9]),
            ([3, 5], [0, 1, 2, 20, 4, 20, 6, 7, 8, 9]),
        ],
    )
    def test_loc_setitem_listlike_with_timedelta64index(self, indexer, expected):
        # GH#16637
        tdi = to_timedelta(range(10), unit="s")
        df = DataFrame({"x": range(10)}, dtype="int64", index=tdi)

        df.loc[df.index[indexer], "x"] = 20

        expected = DataFrame(
            expected,
            index=tdi,
            columns=["x"],
            dtype="int64",
        )

        tm.assert_frame_equal(expected, df)


class TestLocWithMultiIndex:
    @pytest.mark.parametrize(
        "keys, expected",
        [
            (["b", "a"], [["b", "b", "a", "a"], [1, 2, 1, 2]]),
            (["a", "b"], [["a", "a", "b", "b"], [1, 2, 1, 2]]),
            ((["a", "b"], [1, 2]), [["a", "a", "b", "b"], [1, 2, 1, 2]]),
            ((["a", "b"], [2, 1]), [["a", "a", "b", "b"], [2, 1, 2, 1]]),
            ((["b", "a"], [2, 1]), [["b", "b", "a", "a"], [2, 1, 2, 1]]),
            ((["b", "a"], [1, 2]), [["b", "b", "a", "a"], [1, 2, 1, 2]]),
            ((["c", "a"], [2, 1]), [["c", "a", "a"], [1, 2, 1]]),
        ],
    )
    @pytest.mark.parametrize("dim", ["index", "columns"])
    def test_loc_getitem_multilevel_index_order(self, dim, keys, expected):
        # GH#22797
        # Try to respect order of keys given for MultiIndex.loc
        kwargs = {dim: [["c", "a", "a", "b", "b"], [1, 1, 2, 1, 2]]}
        df = DataFrame(np.arange(25).reshape(5, 5), **kwargs)
        exp_index = MultiIndex.from_arrays(expected)
        if dim == "index":
            res = df.loc[keys, :]
            tm.assert_index_equal(res.index, exp_index)
        elif dim == "columns":
            res = df.loc[:, keys]
            tm.assert_index_equal(res.columns, exp_index)

    def test_loc_preserve_names(self, multiindex_year_month_day_dataframe_random_data):
        ymd = multiindex_year_month_day_dataframe_random_data

        result = ymd.loc[2000]
        result2 = ymd["A"].loc[2000]
        assert result.index.names == ymd.index.names[1:]
        assert result2.index.names == ymd.index.names[1:]

        result = ymd.loc[2000, 2]
        result2 = ymd["A"].loc[2000, 2]
        assert result.index.name == ymd.index.names[2]
        assert result2.index.name == ymd.index.names[2]

    def test_loc_getitem_multiindex_nonunique_len_zero(self):
        # GH#13691
        mi = MultiIndex.from_product([[0], [1, 1]])
        ser = Series(0, index=mi)

        res = ser.loc[[]]

        expected = ser[:0]
        tm.assert_series_equal(res, expected)

        res2 = ser.loc[ser.iloc[0:0]]
        tm.assert_series_equal(res2, expected)

    def test_loc_getitem_access_none_value_in_multiindex(self):
        # GH#34318: test that you can access a None value using .loc
        #  through a Multiindex

        ser = Series([None], pd.MultiIndex.from_arrays([["Level1"], ["Level2"]]))
        result = ser.loc[("Level1", "Level2")]
        assert result is None

        midx = MultiIndex.from_product([["Level1"], ["Level2_a", "Level2_b"]])
        ser = Series([None] * len(midx), dtype=object, index=midx)
        result = ser.loc[("Level1", "Level2_a")]
        assert result is None

        ser = Series([1] * len(midx), dtype=object, index=midx)
        result = ser.loc[("Level1", "Level2_a")]
        assert result == 1

    def test_loc_setitem_multiindex_slice(self):
        # GH 34870

        index = pd.MultiIndex.from_tuples(
            zip(
                ["bar", "bar", "baz", "baz", "foo", "foo", "qux", "qux"],
                ["one", "two", "one", "two", "one", "two", "one", "two"],
            ),
            names=["first", "second"],
        )

        result = Series([1, 1, 1, 1, 1, 1, 1, 1], index=index)
        result.loc[("baz", "one"):("foo", "two")] = 100

        expected = Series([1, 1, 100, 100, 100, 100, 1, 1], index=index)

        tm.assert_series_equal(result, expected)

    def test_loc_getitem_slice_datetime_objs_with_datetimeindex(self):
        times = date_range("2000-01-01", freq="10min", periods=100000)
        ser = Series(range(100000), times)
        result = ser.loc[datetime(1900, 1, 1) : datetime(2100, 1, 1)]
        tm.assert_series_equal(result, ser)

    def test_loc_getitem_sorted_index_level_with_duplicates(self):
        # GH#4516 sorting a MultiIndex with duplicates and multiple dtypes
        mi = MultiIndex.from_tuples(
            [
                ("foo", "bar"),
                ("foo", "bar"),
                ("bah", "bam"),
                ("bah", "bam"),
                ("foo", "bar"),
                ("bah", "bam"),
            ],
            names=["A", "B"],
        )
        df = DataFrame(
            [
                [1.0, 1],
                [2.0, 2],
                [3.0, 3],
                [4.0, 4],
                [5.0, 5],
                [6.0, 6],
            ],
            index=mi,
            columns=["C", "D"],
        )
        df = df.sort_index(level=0)

        expected = DataFrame(
            [[1.0, 1], [2.0, 2], [5.0, 5]], columns=["C", "D"], index=mi.take([0, 1, 4])
        )

        result = df.loc[("foo", "bar")]
        tm.assert_frame_equal(result, expected)


class TestLocSetitemWithExpansion:
    @pytest.mark.slow
    def test_loc_setitem_with_expansion_large_dataframe(self):
        # GH#10692
        result = DataFrame({"x": range(10 ** 6)}, dtype="int64")
        result.loc[len(result)] = len(result) + 1
        expected = DataFrame({"x": range(10 ** 6 + 1)}, dtype="int64")
        tm.assert_frame_equal(result, expected)

    def test_loc_setitem_empty_series(self):
        # GH#5226

        # partially set with an empty object series
        ser = Series(dtype=object)
        ser.loc[1] = 1
        tm.assert_series_equal(ser, Series([1], index=[1]))
        ser.loc[3] = 3
        tm.assert_series_equal(ser, Series([1, 3], index=[1, 3]))

        ser = Series(dtype=object)
        ser.loc[1] = 1.0
        tm.assert_series_equal(ser, Series([1.0], index=[1]))
        ser.loc[3] = 3.0
        tm.assert_series_equal(ser, Series([1.0, 3.0], index=[1, 3]))

        ser = Series(dtype=object)
        ser.loc["foo"] = 1
        tm.assert_series_equal(ser, Series([1], index=["foo"]))
        ser.loc["bar"] = 3
        tm.assert_series_equal(ser, Series([1, 3], index=["foo", "bar"]))
        ser.loc[3] = 4
        tm.assert_series_equal(ser, Series([1, 3, 4], index=["foo", "bar", 3]))

    def test_loc_setitem_incremental_with_dst(self):
        # GH#20724
        base = datetime(2015, 11, 1, tzinfo=gettz("US/Pacific"))
        idxs = [base + timedelta(seconds=i * 900) for i in range(16)]
        result = Series([0], index=[idxs[0]])
        for ts in idxs:
            result.loc[ts] = 1
        expected = Series(1, index=idxs)
        tm.assert_series_equal(result, expected)

    def test_loc_setitem_datetime_keys_cast(self):
        # GH#9516
        dt1 = Timestamp("20130101 09:00:00")
        dt2 = Timestamp("20130101 10:00:00")

        for conv in [
            lambda x: x,
            lambda x: x.to_datetime64(),
            lambda x: x.to_pydatetime(),
            lambda x: np.datetime64(x),
        ]:

            df = DataFrame()
            df.loc[conv(dt1), "one"] = 100
            df.loc[conv(dt2), "one"] = 200

            expected = DataFrame({"one": [100.0, 200.0]}, index=[dt1, dt2])
            tm.assert_frame_equal(df, expected)

    def test_loc_setitem_categorical_column_retains_dtype(self, ordered):
        # GH16360
        result = DataFrame({"A": [1]})
        result.loc[:, "B"] = Categorical(["b"], ordered=ordered)
        expected = DataFrame({"A": [1], "B": Categorical(["b"], ordered=ordered)})
        tm.assert_frame_equal(result, expected)


class TestLocCallable:
    def test_frame_loc_getitem_callable(self):
        # GH#11485
        df = DataFrame({"A": [1, 2, 3, 4], "B": list("aabb"), "C": [1, 2, 3, 4]})
        # iloc cannot use boolean Series (see GH3635)

        # return bool indexer
        res = df.loc[lambda x: x.A > 2]
        tm.assert_frame_equal(res, df.loc[df.A > 2])

        res = df.loc[lambda x: x.A > 2]
        tm.assert_frame_equal(res, df.loc[df.A > 2])

        res = df.loc[lambda x: x.A > 2]
        tm.assert_frame_equal(res, df.loc[df.A > 2])

        res = df.loc[lambda x: x.A > 2]
        tm.assert_frame_equal(res, df.loc[df.A > 2])

        res = df.loc[lambda x: x.B == "b", :]
        tm.assert_frame_equal(res, df.loc[df.B == "b", :])

        res = df.loc[lambda x: x.B == "b", :]
        tm.assert_frame_equal(res, df.loc[df.B == "b", :])

        res = df.loc[lambda x: x.A > 2, lambda x: x.columns == "B"]
        tm.assert_frame_equal(res, df.loc[df.A > 2, [False, True, False]])

        res = df.loc[lambda x: x.A > 2, lambda x: x.columns == "B"]
        tm.assert_frame_equal(res, df.loc[df.A > 2, [False, True, False]])

        res = df.loc[lambda x: x.A > 2, lambda x: "B"]
        tm.assert_series_equal(res, df.loc[df.A > 2, "B"])

        res = df.loc[lambda x: x.A > 2, lambda x: "B"]
        tm.assert_series_equal(res, df.loc[df.A > 2, "B"])

        res = df.loc[lambda x: x.A > 2, lambda x: ["A", "B"]]
        tm.assert_frame_equal(res, df.loc[df.A > 2, ["A", "B"]])

        res = df.loc[lambda x: x.A > 2, lambda x: ["A", "B"]]
        tm.assert_frame_equal(res, df.loc[df.A > 2, ["A", "B"]])

        res = df.loc[lambda x: x.A == 2, lambda x: ["A", "B"]]
        tm.assert_frame_equal(res, df.loc[df.A == 2, ["A", "B"]])

        res = df.loc[lambda x: x.A == 2, lambda x: ["A", "B"]]
        tm.assert_frame_equal(res, df.loc[df.A == 2, ["A", "B"]])

        # scalar
        res = df.loc[lambda x: 1, lambda x: "A"]
        assert res == df.loc[1, "A"]

        res = df.loc[lambda x: 1, lambda x: "A"]
        assert res == df.loc[1, "A"]

    def test_frame_loc_getitem_callable_mixture(self):
        # GH#11485
        df = DataFrame({"A": [1, 2, 3, 4], "B": list("aabb"), "C": [1, 2, 3, 4]})

        res = df.loc[lambda x: x.A > 2, ["A", "B"]]
        tm.assert_frame_equal(res, df.loc[df.A > 2, ["A", "B"]])

        res = df.loc[lambda x: x.A > 2, ["A", "B"]]
        tm.assert_frame_equal(res, df.loc[df.A > 2, ["A", "B"]])

        res = df.loc[[2, 3], lambda x: ["A", "B"]]
        tm.assert_frame_equal(res, df.loc[[2, 3], ["A", "B"]])

        res = df.loc[[2, 3], lambda x: ["A", "B"]]
        tm.assert_frame_equal(res, df.loc[[2, 3], ["A", "B"]])

        res = df.loc[3, lambda x: ["A", "B"]]
        tm.assert_series_equal(res, df.loc[3, ["A", "B"]])

        res = df.loc[3, lambda x: ["A", "B"]]
        tm.assert_series_equal(res, df.loc[3, ["A", "B"]])

    def test_frame_loc_getitem_callable_labels(self):
        # GH#11485
        df = DataFrame({"X": [1, 2, 3, 4], "Y": list("aabb")}, index=list("ABCD"))

        # return label
        res = df.loc[lambda x: ["A", "C"]]
        tm.assert_frame_equal(res, df.loc[["A", "C"]])

        res = df.loc[lambda x: ["A", "C"]]
        tm.assert_frame_equal(res, df.loc[["A", "C"]])

        res = df.loc[lambda x: ["A", "C"], :]
        tm.assert_frame_equal(res, df.loc[["A", "C"], :])

        res = df.loc[lambda x: ["A", "C"], lambda x: "X"]
        tm.assert_series_equal(res, df.loc[["A", "C"], "X"])

        res = df.loc[lambda x: ["A", "C"], lambda x: ["X"]]
        tm.assert_frame_equal(res, df.loc[["A", "C"], ["X"]])

        # mixture
        res = df.loc[["A", "C"], lambda x: "X"]
        tm.assert_series_equal(res, df.loc[["A", "C"], "X"])

        res = df.loc[["A", "C"], lambda x: ["X"]]
        tm.assert_frame_equal(res, df.loc[["A", "C"], ["X"]])

        res = df.loc[lambda x: ["A", "C"], "X"]
        tm.assert_series_equal(res, df.loc[["A", "C"], "X"])

        res = df.loc[lambda x: ["A", "C"], ["X"]]
        tm.assert_frame_equal(res, df.loc[["A", "C"], ["X"]])

    def test_frame_loc_setitem_callable(self):
        # GH#11485
        df = DataFrame({"X": [1, 2, 3, 4], "Y": list("aabb")}, index=list("ABCD"))

        # return label
        res = df.copy()
        res.loc[lambda x: ["A", "C"]] = -20
        exp = df.copy()
        exp.loc[["A", "C"]] = -20
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.loc[lambda x: ["A", "C"], :] = 20
        exp = df.copy()
        exp.loc[["A", "C"], :] = 20
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.loc[lambda x: ["A", "C"], lambda x: "X"] = -1
        exp = df.copy()
        exp.loc[["A", "C"], "X"] = -1
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.loc[lambda x: ["A", "C"], lambda x: ["X"]] = [5, 10]
        exp = df.copy()
        exp.loc[["A", "C"], ["X"]] = [5, 10]
        tm.assert_frame_equal(res, exp)

        # mixture
        res = df.copy()
        res.loc[["A", "C"], lambda x: "X"] = np.array([-1, -2])
        exp = df.copy()
        exp.loc[["A", "C"], "X"] = np.array([-1, -2])
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.loc[["A", "C"], lambda x: ["X"]] = 10
        exp = df.copy()
        exp.loc[["A", "C"], ["X"]] = 10
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.loc[lambda x: ["A", "C"], "X"] = -2
        exp = df.copy()
        exp.loc[["A", "C"], "X"] = -2
        tm.assert_frame_equal(res, exp)

        res = df.copy()
        res.loc[lambda x: ["A", "C"], ["X"]] = -4
        exp = df.copy()
        exp.loc[["A", "C"], ["X"]] = -4
        tm.assert_frame_equal(res, exp)


class TestPartialStringSlicing:
    def test_loc_getitem_partial_string_slicing_datetimeindex(self):
        # GH#35509
        df = DataFrame(
            {"col1": ["a", "b", "c"], "col2": [1, 2, 3]},
            index=to_datetime(["2020-08-01", "2020-07-02", "2020-08-05"]),
        )
        expected = DataFrame(
            {"col1": ["a", "c"], "col2": [1, 3]},
            index=to_datetime(["2020-08-01", "2020-08-05"]),
        )
        result = df.loc["2020-08"]
        tm.assert_frame_equal(result, expected)

    def test_loc_getitem_partial_string_slicing_with_periodindex(self):
        pi = pd.period_range(start="2017-01-01", end="2018-01-01", freq="M")
        ser = pi.to_series()
        result = ser.loc[:"2017-12"]
        expected = ser.iloc[:-1]

        tm.assert_series_equal(result, expected)

    def test_loc_getitem_partial_string_slicing_with_timedeltaindex(self):
        ix = timedelta_range(start="1 day", end="2 days", freq="1H")
        ser = ix.to_series()
        result = ser.loc[:"1 days"]
        expected = ser.iloc[:-1]

        tm.assert_series_equal(result, expected)

    def test_loc_getitem_str_timedeltaindex(self):
        # GH#16896
        df = DataFrame({"x": range(3)}, index=to_timedelta(range(3), unit="days"))
        expected = df.iloc[0]
        sliced = df.loc["0 days"]
        tm.assert_series_equal(sliced, expected)


class TestLabelSlicing:
    def test_loc_getitem_label_slice_across_dst(self):
        # GH#21846
        idx = date_range(
            "2017-10-29 01:30:00", tz="Europe/Berlin", periods=5, freq="30 min"
        )
        series2 = Series([0, 1, 2, 3, 4], index=idx)

        t_1 = Timestamp("2017-10-29 02:30:00+02:00", tz="Europe/Berlin", freq="30min")
        t_2 = Timestamp("2017-10-29 02:00:00+01:00", tz="Europe/Berlin", freq="30min")
        result = series2.loc[t_1:t_2]
        expected = Series([2, 3], index=idx[2:4])
        tm.assert_series_equal(result, expected)

        result = series2[t_1]
        expected = 2
        assert result == expected

    def test_loc_getitem_label_slice_period(self):
        ix = pd.period_range(start="2017-01-01", end="2018-01-01", freq="M")
        ser = ix.to_series()
        result = ser.loc[: ix[-2]]
        expected = ser.iloc[:-1]

        tm.assert_series_equal(result, expected)

    def test_loc_getitem_label_slice_timedelta64(self):
        ix = timedelta_range(start="1 day", end="2 days", freq="1H")
        ser = ix.to_series()
        result = ser.loc[: ix[-2]]
        expected = ser.iloc[:-1]

        tm.assert_series_equal(result, expected)

    def test_loc_getitem_slice_floats_inexact(self):
        index = [52195.504153, 52196.303147, 52198.369883]
        df = DataFrame(np.random.rand(3, 2), index=index)

        s1 = df.loc[52195.1:52196.5]
        assert len(s1) == 2

        s1 = df.loc[52195.1:52196.6]
        assert len(s1) == 2

        s1 = df.loc[52195.1:52198.9]
        assert len(s1) == 3

    def test_loc_getitem_float_slice_float64index(self):
        ser = Series(np.random.rand(10), index=np.arange(10, 20, dtype=float))

        assert len(ser.loc[12.0:]) == 8
        assert len(ser.loc[12.5:]) == 7

        idx = np.arange(10, 20, dtype=float)
        idx[2] = 12.2
        ser.index = idx
        assert len(ser.loc[12.0:]) == 8
        assert len(ser.loc[12.5:]) == 7

    @pytest.mark.parametrize(
        "start,stop, expected_slice",
        [
            [np.timedelta64(0, "ns"), None, slice(0, 11)],
            [np.timedelta64(1, "D"), np.timedelta64(6, "D"), slice(1, 7)],
            [None, np.timedelta64(4, "D"), slice(0, 5)],
        ],
    )
    def test_loc_getitem_slice_label_td64obj(self, start, stop, expected_slice):
        # GH#20393
        ser = Series(range(11), timedelta_range("0 days", "10 days"))
        result = ser.loc[slice(start, stop)]
        expected = ser.iloc[expected_slice]
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("start", ["2018", "2020"])
    def test_loc_getitem_slice_unordered_dt_index(self, frame_or_series, start):
        obj = frame_or_series(
            [1, 2, 3],
            index=[Timestamp("2016"), Timestamp("2019"), Timestamp("2017")],
        )
        with tm.assert_produces_warning(FutureWarning):
            obj.loc[start:"2022"]

    @pytest.mark.parametrize("value", [1, 1.5])
    def test_loc_getitem_slice_labels_int_in_object_index(self, frame_or_series, value):
        # GH: 26491
        obj = frame_or_series(range(4), index=[value, "first", 2, "third"])
        result = obj.loc[value:"third"]
        expected = frame_or_series(range(4), index=[value, "first", 2, "third"])
        tm.assert_equal(result, expected)


class TestLocBooleanMask:
    def test_loc_setitem_bool_mask_timedeltaindex(self):
        # GH#14946
        df = DataFrame({"x": range(10)})
        df.index = to_timedelta(range(10), unit="s")
        conditions = [df["x"] > 3, df["x"] == 3, df["x"] < 3]
        expected_data = [
            [0, 1, 2, 3, 10, 10, 10, 10, 10, 10],
            [0, 1, 2, 10, 4, 5, 6, 7, 8, 9],
            [10, 10, 10, 3, 4, 5, 6, 7, 8, 9],
        ]
        for cond, data in zip(conditions, expected_data):
            result = df.copy()
            result.loc[cond, "x"] = 10

            expected = DataFrame(
                data,
                index=to_timedelta(range(10), unit="s"),
                columns=["x"],
                dtype="int64",
            )
            tm.assert_frame_equal(expected, result)

    def test_loc_setitem_mask_with_datetimeindex_tz(self):
        # GH#16889
        # support .loc with alignment and tz-aware DatetimeIndex
        mask = np.array([True, False, True, False])

        idx = date_range("20010101", periods=4, tz="UTC")
        df = DataFrame({"a": np.arange(4)}, index=idx).astype("float64")

        result = df.copy()
        result.loc[mask, :] = df.loc[mask, :]
        tm.assert_frame_equal(result, df)

        result = df.copy()
        result.loc[mask] = df.loc[mask]
        tm.assert_frame_equal(result, df)

        idx = date_range("20010101", periods=4)
        df = DataFrame({"a": np.arange(4)}, index=idx).astype("float64")

        result = df.copy()
        result.loc[mask, :] = df.loc[mask, :]
        tm.assert_frame_equal(result, df)

        result = df.copy()
        result.loc[mask] = df.loc[mask]
        tm.assert_frame_equal(result, df)

    def test_loc_setitem_mask_and_label_with_datetimeindex(self):
        # GH#9478
        # a datetimeindex alignment issue with partial setting
        df = DataFrame(
            np.arange(6.0).reshape(3, 2),
            columns=list("AB"),
            index=date_range("1/1/2000", periods=3, freq="1H"),
        )
        expected = df.copy()
        expected["C"] = [expected.index[0]] + [pd.NaT, pd.NaT]

        mask = df.A < 1
        df.loc[mask, "C"] = df.loc[mask].index
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_mask_td64_series_value(self):
        # GH#23462 key list of bools, value is a Series
        td1 = Timedelta(0)
        td2 = Timedelta(28767471428571405)
        df = DataFrame({"col": Series([td1, td2])})
        df_copy = df.copy()
        ser = Series([td1])

        expected = df["col"].iloc[1].value
        df.loc[[True, False]] = ser
        result = df["col"].iloc[1].value

        assert expected == result
        tm.assert_frame_equal(df, df_copy)


class TestLocListlike:
    @pytest.mark.parametrize("box", [lambda x: x, np.asarray, list])
    def test_loc_getitem_list_of_labels_categoricalindex_with_na(self, box):
        # passing a list can include valid categories _or_ NA values
        ci = CategoricalIndex(["A", "B", np.nan])
        ser = Series(range(3), index=ci)

        result = ser.loc[box(ci)]
        tm.assert_series_equal(result, ser)

        result = ser[box(ci)]
        tm.assert_series_equal(result, ser)

        result = ser.to_frame().loc[box(ci)]
        tm.assert_frame_equal(result, ser.to_frame())

        ser2 = ser[:-1]
        ci2 = ci[1:]
        # but if there are no NAs present, this should raise KeyError
        msg = (
            r"Passing list-likes to .loc or \[\] with any missing labels is no "
            "longer supported. The following labels were missing: "
            r"(Categorical)?Index\(\[nan\], .*\). "
            "See https"
        )
        with pytest.raises(KeyError, match=msg):
            ser2.loc[box(ci2)]

        with pytest.raises(KeyError, match=msg):
            ser2[box(ci2)]

        with pytest.raises(KeyError, match=msg):
            ser2.to_frame().loc[box(ci2)]


def test_series_loc_getitem_label_list_missing_values():
    # gh-11428
    key = np.array(
        ["2001-01-04", "2001-01-02", "2001-01-04", "2001-01-14"], dtype="datetime64"
    )
    s = Series([2, 5, 8, 11], date_range("2001-01-01", freq="D", periods=4))
    with pytest.raises(KeyError, match="with any missing labels"):
        s.loc[key]


def test_series_getitem_label_list_missing_integer_values():
    # GH: 25927
    s = Series(
        index=np.array([9730701000001104, 10049011000001109]),
        data=np.array([999000011000001104, 999000011000001104]),
    )
    with pytest.raises(KeyError, match="with any missing labels"):
        s.loc[np.array([9730701000001104, 10047311000001102])]


@pytest.mark.parametrize(
    "columns, column_key, expected_columns",
    [
        ([2011, 2012, 2013], [2011, 2012], [0, 1]),
        ([2011, 2012, "All"], [2011, 2012], [0, 1]),
        ([2011, 2012, "All"], [2011, "All"], [0, 2]),
    ],
)
def test_loc_getitem_label_list_integer_labels(columns, column_key, expected_columns):
    # gh-14836
    df = DataFrame(np.random.rand(3, 3), columns=columns, index=list("ABC"))
    expected = df.iloc[:, expected_columns]
    result = df.loc[["A", "B", "C"], column_key]

    if df.columns.is_object() and all(isinstance(x, int) for x in column_key):
        expected.columns = expected.columns.astype(int)

    tm.assert_frame_equal(result, expected, check_column_type=True)


def test_loc_setitem_float_intindex():
    # GH 8720
    rand_data = np.random.randn(8, 4)
    result = DataFrame(rand_data)
    result.loc[:, 0.5] = np.nan
    expected_data = np.hstack((rand_data, np.array([np.nan] * 8).reshape(8, 1)))
    expected = DataFrame(expected_data, columns=[0.0, 1.0, 2.0, 3.0, 0.5])
    tm.assert_frame_equal(result, expected)

    result = DataFrame(rand_data)
    result.loc[:, 0.5] = np.nan
    tm.assert_frame_equal(result, expected)


def test_loc_axis_1_slice():
    # GH 10586
    cols = [(yr, m) for yr in [2014, 2015] for m in [7, 8, 9, 10]]
    df = DataFrame(
        np.ones((10, 8)),
        index=tuple("ABCDEFGHIJ"),
        columns=pd.MultiIndex.from_tuples(cols),
    )
    result = df.loc(axis=1)[(2014, 9):(2015, 8)]
    expected = DataFrame(
        np.ones((10, 4)),
        index=tuple("ABCDEFGHIJ"),
        columns=pd.MultiIndex.from_tuples(
            [(2014, 9), (2014, 10), (2015, 7), (2015, 8)]
        ),
    )
    tm.assert_frame_equal(result, expected)


def test_loc_set_dataframe_multiindex():
    # GH 14592
    expected = DataFrame(
        "a", index=range(2), columns=pd.MultiIndex.from_product([range(2), range(2)])
    )
    result = expected.copy()
    result.loc[0, [(0, 1)]] = result.loc[0, [(0, 1)]]
    tm.assert_frame_equal(result, expected)


def test_loc_mixed_int_float():
    # GH#19456
    ser = Series(range(2), pd.Index([1, 2.0], dtype=object))

    result = ser.loc[1]
    assert result == 0


def test_loc_with_positional_slice_deprecation():
    # GH#31840
    ser = Series(range(4), index=["A", "B", "C", "D"])

    with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
        ser.loc[:3] = 2

    expected = Series([2, 2, 2, 3], index=["A", "B", "C", "D"])
    tm.assert_series_equal(ser, expected)


def test_loc_slice_disallows_positional():
    # GH#16121, GH#24612, GH#31810
    dti = pd.date_range("2016-01-01", periods=3)
    df = DataFrame(np.random.random((3, 2)), index=dti)

    ser = df[0]

    msg = (
        "cannot do slice indexing on DatetimeIndex with these "
        r"indexers \[1\] of type int"
    )

    for obj in [df, ser]:
        with pytest.raises(TypeError, match=msg):
            obj.loc[1:3]

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            # GH#31840 deprecated incorrect behavior
            obj.loc[1:3] = 1

    with pytest.raises(TypeError, match=msg):
        df.loc[1:3, 1]

    with tm.assert_produces_warning(FutureWarning):
        # GH#31840 deprecated incorrect behavior
        df.loc[1:3, 1] = 2


def test_loc_datetimelike_mismatched_dtypes():
    # GH#32650 dont mix and match datetime/timedelta/period dtypes

    df = DataFrame(
        np.random.randn(5, 3),
        columns=["a", "b", "c"],
        index=pd.date_range("2012", freq="H", periods=5),
    )
    # create dataframe with non-unique DatetimeIndex
    df = df.iloc[[0, 2, 2, 3]].copy()

    dti = df.index
    tdi = pd.TimedeltaIndex(dti.asi8)  # matching i8 values

    msg = r"None of \[TimedeltaIndex.* are in the \[index\]"
    with pytest.raises(KeyError, match=msg):
        df.loc[tdi]

    with pytest.raises(KeyError, match=msg):
        df["a"].loc[tdi]


def test_loc_with_period_index_indexer():
    # GH#4125
    idx = pd.period_range("2002-01", "2003-12", freq="M")
    df = DataFrame(np.random.randn(24, 10), index=idx)
    tm.assert_frame_equal(df, df.loc[idx])
    tm.assert_frame_equal(df, df.loc[list(idx)])
    tm.assert_frame_equal(df, df.loc[list(idx)])
    tm.assert_frame_equal(df.iloc[0:5], df.loc[idx[0:5]])
    tm.assert_frame_equal(df, df.loc[list(idx)])


class TestLocSeries:
    @pytest.mark.parametrize("val,expected", [(2 ** 63 - 1, 3), (2 ** 63, 4)])
    def test_loc_uint64(self, val, expected):
        # see GH#19399
        ser = Series({2 ** 63 - 1: 3, 2 ** 63: 4})
        assert ser.loc[val] == expected

    def test_loc_getitem(self, string_series, datetime_series):
        inds = string_series.index[[3, 4, 7]]
        tm.assert_series_equal(string_series.loc[inds], string_series.reindex(inds))
        tm.assert_series_equal(string_series.iloc[5::2], string_series[5::2])

        # slice with indices
        d1, d2 = datetime_series.index[[5, 15]]
        result = datetime_series.loc[d1:d2]
        expected = datetime_series.truncate(d1, d2)
        tm.assert_series_equal(result, expected)

        # boolean
        mask = string_series > string_series.median()
        tm.assert_series_equal(string_series.loc[mask], string_series[mask])

        # ask for index value
        assert datetime_series.loc[d1] == datetime_series[d1]
        assert datetime_series.loc[d2] == datetime_series[d2]

    def test_loc_getitem_not_monotonic(self, datetime_series):
        d1, d2 = datetime_series.index[[5, 15]]

        ts2 = datetime_series[::2][[1, 2, 0]]

        msg = r"Timestamp\('2000-01-10 00:00:00'\)"
        with pytest.raises(KeyError, match=msg):
            ts2.loc[d1:d2]
        with pytest.raises(KeyError, match=msg):
            ts2.loc[d1:d2] = 0

    def test_loc_getitem_setitem_integer_slice_keyerrors(self):
        ser = Series(np.random.randn(10), index=list(range(0, 20, 2)))

        # this is OK
        cp = ser.copy()
        cp.iloc[4:10] = 0
        assert (cp.iloc[4:10] == 0).all()

        # so is this
        cp = ser.copy()
        cp.iloc[3:11] = 0
        assert (cp.iloc[3:11] == 0).values.all()

        result = ser.iloc[2:6]
        result2 = ser.loc[3:11]
        expected = ser.reindex([4, 6, 8, 10])

        tm.assert_series_equal(result, expected)
        tm.assert_series_equal(result2, expected)

        # non-monotonic, raise KeyError
        s2 = ser.iloc[list(range(5)) + list(range(9, 4, -1))]
        with pytest.raises(KeyError, match=r"^3$"):
            s2.loc[3:11]
        with pytest.raises(KeyError, match=r"^3$"):
            s2.loc[3:11] = 0

    def test_loc_getitem_iterator(self, string_series):
        idx = iter(string_series.index[:10])
        result = string_series.loc[idx]
        tm.assert_series_equal(result, string_series[:10])

    def test_loc_setitem_boolean(self, string_series):
        mask = string_series > string_series.median()

        result = string_series.copy()
        result.loc[mask] = 0
        expected = string_series
        expected[mask] = 0
        tm.assert_series_equal(result, expected)

    def test_loc_setitem_corner(self, string_series):
        inds = list(string_series.index[[5, 8, 12]])
        string_series.loc[inds] = 5
        msg = r"\['foo'\] not in index"
        with pytest.raises(KeyError, match=msg):
            string_series.loc[inds + ["foo"]] = 5

    def test_basic_setitem_with_labels(self, datetime_series):
        indices = datetime_series.index[[5, 10, 15]]

        cp = datetime_series.copy()
        exp = datetime_series.copy()
        cp[indices] = 0
        exp.loc[indices] = 0
        tm.assert_series_equal(cp, exp)

        cp = datetime_series.copy()
        exp = datetime_series.copy()
        cp[indices[0] : indices[2]] = 0
        exp.loc[indices[0] : indices[2]] = 0
        tm.assert_series_equal(cp, exp)

    def test_loc_setitem_listlike_of_ints(self):

        # integer indexes, be careful
        ser = Series(np.random.randn(10), index=list(range(0, 20, 2)))
        inds = [0, 4, 6]
        arr_inds = np.array([0, 4, 6])

        cp = ser.copy()
        exp = ser.copy()
        ser[inds] = 0
        ser.loc[inds] = 0
        tm.assert_series_equal(cp, exp)

        cp = ser.copy()
        exp = ser.copy()
        ser[arr_inds] = 0
        ser.loc[arr_inds] = 0
        tm.assert_series_equal(cp, exp)

        inds_notfound = [0, 4, 5, 6]
        arr_inds_notfound = np.array([0, 4, 5, 6])
        msg = r"\[5\] not in index"
        with pytest.raises(KeyError, match=msg):
            ser[inds_notfound] = 0
        with pytest.raises(Exception, match=msg):
            ser[arr_inds_notfound] = 0

    def test_loc_setitem_dt64tz_values(self):
        # GH#12089
        ser = Series(
            date_range("2011-01-01", periods=3, tz="US/Eastern"),
            index=["a", "b", "c"],
        )
        s2 = ser.copy()
        expected = Timestamp("2011-01-03", tz="US/Eastern")
        s2.loc["a"] = expected
        result = s2.loc["a"]
        assert result == expected

        s2 = ser.copy()
        s2.iloc[0] = expected
        result = s2.iloc[0]
        assert result == expected

        s2 = ser.copy()
        s2["a"] = expected
        result = s2["a"]
        assert result == expected

    @pytest.mark.parametrize("array_fn", [np.array, pd.array, list, tuple])
    @pytest.mark.parametrize("size", [0, 4, 5, 6])
    def test_loc_iloc_setitem_with_listlike(self, size, array_fn):
        # GH37748
        # testing insertion, in a Series of size N (here 5), of a listlike object
        # of size  0, N-1, N, N+1

        arr = array_fn([0] * size)
        expected = Series([arr, 0, 0, 0, 0], index=list("abcde"), dtype=object)

        ser = Series(0, index=list("abcde"), dtype=object)
        ser.loc["a"] = arr
        tm.assert_series_equal(ser, expected)

        ser = Series(0, index=list("abcde"), dtype=object)
        ser.iloc[0] = arr
        tm.assert_series_equal(ser, expected)
