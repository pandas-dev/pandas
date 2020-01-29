""" test label based indexing with loc """
from io import StringIO
import re

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series, Timestamp, date_range
import pandas._testing as tm
from pandas.api.types import is_scalar
from pandas.tests.indexing.common import Base


class TestLoc(Base):
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

        indexer = tuple(["r", ["bar", "bar2"]])
        df = df_orig.copy()
        df.loc[indexer] *= 2.0
        tm.assert_series_equal(df.loc[indexer], 2.0 * df_orig.loc[indexer])

        indexer = tuple(["r", "bar"])
        df = df_orig.copy()
        df.loc[indexer] *= 2.0
        assert df.loc[indexer] == 2.0 * df_orig.loc[indexer]

        indexer = tuple(["t", ["bar", "bar2"]])
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

    def test_loc_getitem_int(self):

        # int label
        self.check_result("loc", 2, "loc", 2, typs=["label"], fails=KeyError)

    def test_loc_getitem_label(self):

        # label
        self.check_result("loc", "c", "loc", "c", typs=["empty"], fails=KeyError)

    def test_loc_getitem_label_out_of_range(self):

        # out of range label
        self.check_result(
            "loc",
            "f",
            "loc",
            "f",
            typs=["ints", "uints", "labels", "mixed", "ts"],
            fails=KeyError,
        )
        self.check_result("loc", "f", "ix", "f", typs=["floats"], fails=KeyError)
        self.check_result("loc", "f", "loc", "f", typs=["floats"], fails=KeyError)
        self.check_result(
            "loc", 20, "loc", 20, typs=["ints", "uints", "mixed"], fails=KeyError,
        )
        self.check_result("loc", 20, "loc", 20, typs=["labels"], fails=TypeError)
        self.check_result("loc", 20, "loc", 20, typs=["ts"], axes=0, fails=TypeError)
        self.check_result("loc", 20, "loc", 20, typs=["floats"], axes=0, fails=KeyError)

    def test_loc_getitem_label_list(self):
        # TODO: test something here?
        # list of labels
        pass

    def test_loc_getitem_label_list_with_missing(self):
        self.check_result(
            "loc", [0, 1, 2], "loc", [0, 1, 2], typs=["empty"], fails=KeyError,
        )
        self.check_result(
            "loc",
            [0, 2, 10],
            "ix",
            [0, 2, 10],
            typs=["ints", "uints", "floats"],
            axes=0,
            fails=KeyError,
        )

        self.check_result(
            "loc",
            [3, 6, 7],
            "ix",
            [3, 6, 7],
            typs=["ints", "uints", "floats"],
            axes=1,
            fails=KeyError,
        )

        # GH 17758 - MultiIndex and missing keys
        self.check_result(
            "loc",
            [(1, 3), (1, 4), (2, 5)],
            "ix",
            [(1, 3), (1, 4), (2, 5)],
            typs=["multi"],
            axes=0,
            fails=KeyError,
        )

    def test_getitem_label_list_with_missing(self):
        s = Series(range(3), index=["a", "b", "c"])

        # consistency
        with pytest.raises(KeyError, match="with any missing labels"):
            s[["a", "d"]]

        s = Series(range(3))
        with pytest.raises(KeyError, match="with any missing labels"):
            s[[0, 3]]

    def test_loc_getitem_label_list_fails(self):
        # fails
        self.check_result(
            "loc",
            [20, 30, 40],
            "loc",
            [20, 30, 40],
            typs=["ints", "uints"],
            axes=1,
            fails=KeyError,
        )

    def test_loc_getitem_label_array_like(self):
        # TODO: test something?
        # array like
        pass

    def test_loc_getitem_bool(self):
        # boolean indexers
        b = [True, False, True, False]

        self.check_result("loc", b, "loc", b, typs=["empty"], fails=IndexError)

    @pytest.mark.parametrize("index", [[True, False], [True, False, True, False]])
    def test_loc_getitem_bool_diff_len(self, index):
        # GH26658
        s = Series([1, 2, 3])
        msg = "Boolean index has wrong length: {} instead of {}".format(
            len(index), len(s)
        )
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
            r"\"None of \[Int64Index\(\[1, 2\], dtype='int64'\)\] are"
            r" in the \[index\]\""
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
            r"\"None of \[Int64Index\(\[-1, -2\], dtype='int64'\)\] are"
            r" in the \[index\]\""
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
            r"\"None of \[Int64Index\(\[-2\], dtype='int64'\)\] are"
            r" in the \[index\]\""
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
            r"\"None of \[Int64Index\(\[3\], dtype='int64'\)\] are"
            r" in the \[index\]\""
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

    def test_loc_getitem_label_slice(self):

        # label slices (with ints)

        # real label slices

        # GH 14316

        self.check_result(
            "loc",
            slice(1, 3),
            "loc",
            slice(1, 3),
            typs=["labels", "mixed", "empty", "ts", "floats"],
            fails=TypeError,
        )

        self.check_result(
            "loc",
            slice("20130102", "20130104"),
            "loc",
            slice("20130102", "20130104"),
            typs=["ts"],
            axes=1,
            fails=TypeError,
        )

        self.check_result(
            "loc",
            slice(2, 8),
            "loc",
            slice(2, 8),
            typs=["mixed"],
            axes=0,
            fails=TypeError,
        )
        self.check_result(
            "loc",
            slice(2, 8),
            "loc",
            slice(2, 8),
            typs=["mixed"],
            axes=1,
            fails=KeyError,
        )

        self.check_result(
            "loc",
            slice(2, 4, 2),
            "loc",
            slice(2, 4, 2),
            typs=["mixed"],
            axes=0,
            fails=TypeError,
        )

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

        df.loc[:, unit] = df.loc[:, "timestamp"].values.astype(
            "datetime64[{unit}]".format(unit=unit)
        )
        df["expected"] = df.loc[:, "timestamp"].values.astype(
            "datetime64[{unit}]".format(unit=unit)
        )
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
        df = self.frame_labels

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
        expected = DataFrame(dict(A=Series([1, 2, 3], index=[4, 3, 5]))).reindex(
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
            dict(A=Series(val1, index=keys1), B=Series(val2, index=keys2))
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

        df = pd.DataFrame([[1, 2], [3, 4]], columns=["A", "B"])
        df.loc[0, indexer] = value
        result = df.loc[0, "A"]

        assert is_scalar(result) and result == "Z"

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
        expected = pd.Series(vals, index=["foo", "bar"])

        ser = pd.Series(dtype=object)
        ser["foo"] = vals[0]
        ser["bar"] = vals[1]

        tm.assert_series_equal(ser, expected)

        ser = pd.Series(dtype=object)
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

    def test_loc_non_unique_memory_error(self):

        # GH 4280
        # non_unique index with a large selection triggers a memory error

        columns = list("ABCDEFG")

        def gen_test(l, l2):
            return pd.concat(
                [
                    DataFrame(
                        np.random.randn(l, len(columns)),
                        index=np.arange(l),
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

    def test_loc_uint64(self):
        # GH20722
        # Test whether loc accept uint64 max value as index.
        s = pd.Series(
            [1, 2], index=[np.iinfo("uint64").max - 1, np.iinfo("uint64").max]
        )

        result = s.loc[np.iinfo("uint64").max - 1]
        expected = s.iloc[0]
        assert result == expected

        result = s.loc[[np.iinfo("uint64").max - 1]]
        expected = s.iloc[[0]]
        tm.assert_series_equal(result, expected)

        result = s.loc[[np.iinfo("uint64").max - 1, np.iinfo("uint64").max]]
        tm.assert_series_equal(result, s)

    def test_loc_setitem_empty_append(self):
        # GH6173, various appends to an empty dataframe

        data = [1, 2, 3]
        expected = DataFrame({"x": data, "y": [None] * len(data)})

        # appends to fit length of data
        df = DataFrame(columns=["x", "y"])
        df.loc[:, "x"] = data
        tm.assert_frame_equal(df, expected)

        # only appends one value
        expected = DataFrame({"x": [1.0], "y": [np.nan]})
        df = DataFrame(columns=["x", "y"], dtype=np.float)
        df.loc[0, "x"] = expected.loc[0, "x"]
        tm.assert_frame_equal(df, expected)

    def test_loc_setitem_empty_append_raises(self):
        # GH6173, various appends to an empty dataframe

        data = [1, 2]
        df = DataFrame(columns=["x", "y"])
        msg = (
            r"None of \[Int64Index\(\[0, 1\], dtype='int64'\)\] "
            r"are in the \[index\]"
        )
        with pytest.raises(KeyError, match=msg):
            df.loc[[0, 1], "x"] = data

        msg = "cannot copy sequence with size 2 to array axis with dimension 0"
        with pytest.raises(ValueError, match=msg):
            df.loc[0:2, "x"] = data

    def test_indexing_zerodim_np_array(self):
        # GH24924
        df = DataFrame([[1, 2], [3, 4]])
        result = df.loc[np.array(0)]
        s = pd.Series([1, 2], name=0)
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

        result = pd.Series(index=range(2010, 2020), dtype=np.float64)
        result.loc[2015:2010:-1] = [6, 5, 4, 3, 2, 1]

        tm.assert_series_equal(result, expected)


def test_series_loc_getitem_label_list_missing_values():
    # gh-11428
    key = np.array(
        ["2001-01-04", "2001-01-02", "2001-01-04", "2001-01-14"], dtype="datetime64"
    )
    s = Series([2, 5, 8, 11], date_range("2001-01-01", freq="D", periods=4))
    with pytest.raises(KeyError, match="with any missing labels"):
        s.loc[key]


@pytest.mark.parametrize(
    "columns, column_key, expected_columns, check_column_type",
    [
        ([2011, 2012, 2013], [2011, 2012], [0, 1], True),
        ([2011, 2012, "All"], [2011, 2012], [0, 1], False),
        ([2011, 2012, "All"], [2011, "All"], [0, 2], True),
    ],
)
def test_loc_getitem_label_list_integer_labels(
    columns, column_key, expected_columns, check_column_type
):
    # gh-14836
    df = DataFrame(np.random.rand(3, 3), columns=columns, index=list("ABC"))
    expected = df.iloc[:, expected_columns]
    result = df.loc[["A", "B", "C"], column_key]
    tm.assert_frame_equal(result, expected, check_column_type=check_column_type)


def test_loc_setitem_float_intindex():
    # GH 8720
    rand_data = np.random.randn(8, 4)
    result = pd.DataFrame(rand_data)
    result.loc[:, 0.5] = np.nan
    expected_data = np.hstack((rand_data, np.array([np.nan] * 8).reshape(8, 1)))
    expected = pd.DataFrame(expected_data, columns=[0.0, 1.0, 2.0, 3.0, 0.5])
    tm.assert_frame_equal(result, expected)

    result = pd.DataFrame(rand_data)
    result.loc[:, 0.5] = np.nan
    tm.assert_frame_equal(result, expected)


def test_loc_axis_1_slice():
    # GH 10586
    cols = [(yr, m) for yr in [2014, 2015] for m in [7, 8, 9, 10]]
    df = pd.DataFrame(
        np.ones((10, 8)),
        index=tuple("ABCDEFGHIJ"),
        columns=pd.MultiIndex.from_tuples(cols),
    )
    result = df.loc(axis=1)[(2014, 9):(2015, 8)]
    expected = pd.DataFrame(
        np.ones((10, 4)),
        index=tuple("ABCDEFGHIJ"),
        columns=pd.MultiIndex.from_tuples(
            [(2014, 9), (2014, 10), (2015, 7), (2015, 8)]
        ),
    )
    tm.assert_frame_equal(result, expected)
