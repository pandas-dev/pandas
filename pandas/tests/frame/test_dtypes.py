from collections import OrderedDict
from datetime import timedelta

import numpy as np
import pytest

from pandas.core.dtypes.dtypes import CategoricalDtype, DatetimeTZDtype, IntervalDtype

import pandas as pd
from pandas import (
    Categorical,
    DataFrame,
    Series,
    Timedelta,
    Timestamp,
    _np_version_under1p14,
    concat,
    date_range,
    option_context,
)
import pandas._testing as tm
from pandas.core.arrays import integer_array


def _check_cast(df, v):
    """
    Check if all dtypes of df are equal to v
    """
    assert all(s.dtype.name == v for _, s in df.items())


class TestDataFrameDataTypes:
    def test_concat_empty_dataframe_dtypes(self):
        df = DataFrame(columns=list("abc"))
        df["a"] = df["a"].astype(np.bool_)
        df["b"] = df["b"].astype(np.int32)
        df["c"] = df["c"].astype(np.float64)

        result = pd.concat([df, df])
        assert result["a"].dtype == np.bool_
        assert result["b"].dtype == np.int32
        assert result["c"].dtype == np.float64

        result = pd.concat([df, df.astype(np.float64)])
        assert result["a"].dtype == np.object_
        assert result["b"].dtype == np.float64
        assert result["c"].dtype == np.float64

    def test_empty_frame_dtypes(self):
        empty_df = pd.DataFrame()
        tm.assert_series_equal(empty_df.dtypes, pd.Series(dtype=np.object))

        nocols_df = pd.DataFrame(index=[1, 2, 3])
        tm.assert_series_equal(nocols_df.dtypes, pd.Series(dtype=np.object))

        norows_df = pd.DataFrame(columns=list("abc"))
        tm.assert_series_equal(
            norows_df.dtypes, pd.Series(np.object, index=list("abc"))
        )

        norows_int_df = pd.DataFrame(columns=list("abc")).astype(np.int32)
        tm.assert_series_equal(
            norows_int_df.dtypes, pd.Series(np.dtype("int32"), index=list("abc"))
        )

        odict = OrderedDict
        df = pd.DataFrame(odict([("a", 1), ("b", True), ("c", 1.0)]), index=[1, 2, 3])
        ex_dtypes = pd.Series(
            odict([("a", np.int64), ("b", np.bool), ("c", np.float64)])
        )
        tm.assert_series_equal(df.dtypes, ex_dtypes)

        # same but for empty slice of df
        tm.assert_series_equal(df[:0].dtypes, ex_dtypes)

    def test_datetime_with_tz_dtypes(self):
        tzframe = DataFrame(
            {
                "A": date_range("20130101", periods=3),
                "B": date_range("20130101", periods=3, tz="US/Eastern"),
                "C": date_range("20130101", periods=3, tz="CET"),
            }
        )
        tzframe.iloc[1, 1] = pd.NaT
        tzframe.iloc[1, 2] = pd.NaT
        result = tzframe.dtypes.sort_index()
        expected = Series(
            [
                np.dtype("datetime64[ns]"),
                DatetimeTZDtype("ns", "US/Eastern"),
                DatetimeTZDtype("ns", "CET"),
            ],
            ["A", "B", "C"],
        )

        tm.assert_series_equal(result, expected)

    def test_dtypes_are_correct_after_column_slice(self):
        # GH6525
        df = pd.DataFrame(index=range(5), columns=list("abc"), dtype=np.float_)
        odict = OrderedDict
        tm.assert_series_equal(
            df.dtypes,
            pd.Series(odict([("a", np.float_), ("b", np.float_), ("c", np.float_)])),
        )
        tm.assert_series_equal(
            df.iloc[:, 2:].dtypes, pd.Series(odict([("c", np.float_)]))
        )
        tm.assert_series_equal(
            df.dtypes,
            pd.Series(odict([("a", np.float_), ("b", np.float_), ("c", np.float_)])),
        )

    def test_select_dtypes_include_using_list_like(self):
        df = DataFrame(
            {
                "a": list("abc"),
                "b": list(range(1, 4)),
                "c": np.arange(3, 6).astype("u1"),
                "d": np.arange(4.0, 7.0, dtype="float64"),
                "e": [True, False, True],
                "f": pd.Categorical(list("abc")),
                "g": pd.date_range("20130101", periods=3),
                "h": pd.date_range("20130101", periods=3, tz="US/Eastern"),
                "i": pd.date_range("20130101", periods=3, tz="CET"),
                "j": pd.period_range("2013-01", periods=3, freq="M"),
                "k": pd.timedelta_range("1 day", periods=3),
            }
        )

        ri = df.select_dtypes(include=[np.number])
        ei = df[["b", "c", "d", "k"]]
        tm.assert_frame_equal(ri, ei)

        ri = df.select_dtypes(include=[np.number], exclude=["timedelta"])
        ei = df[["b", "c", "d"]]
        tm.assert_frame_equal(ri, ei)

        ri = df.select_dtypes(include=[np.number, "category"], exclude=["timedelta"])
        ei = df[["b", "c", "d", "f"]]
        tm.assert_frame_equal(ri, ei)

        ri = df.select_dtypes(include=["datetime"])
        ei = df[["g"]]
        tm.assert_frame_equal(ri, ei)

        ri = df.select_dtypes(include=["datetime64"])
        ei = df[["g"]]
        tm.assert_frame_equal(ri, ei)

        ri = df.select_dtypes(include=["datetimetz"])
        ei = df[["h", "i"]]
        tm.assert_frame_equal(ri, ei)

        with pytest.raises(NotImplementedError, match=r"^$"):
            df.select_dtypes(include=["period"])

    def test_select_dtypes_exclude_using_list_like(self):
        df = DataFrame(
            {
                "a": list("abc"),
                "b": list(range(1, 4)),
                "c": np.arange(3, 6).astype("u1"),
                "d": np.arange(4.0, 7.0, dtype="float64"),
                "e": [True, False, True],
            }
        )
        re = df.select_dtypes(exclude=[np.number])
        ee = df[["a", "e"]]
        tm.assert_frame_equal(re, ee)

    def test_select_dtypes_exclude_include_using_list_like(self):
        df = DataFrame(
            {
                "a": list("abc"),
                "b": list(range(1, 4)),
                "c": np.arange(3, 6).astype("u1"),
                "d": np.arange(4.0, 7.0, dtype="float64"),
                "e": [True, False, True],
                "f": pd.date_range("now", periods=3).values,
            }
        )
        exclude = (np.datetime64,)
        include = np.bool_, "integer"
        r = df.select_dtypes(include=include, exclude=exclude)
        e = df[["b", "c", "e"]]
        tm.assert_frame_equal(r, e)

        exclude = ("datetime",)
        include = "bool", "int64", "int32"
        r = df.select_dtypes(include=include, exclude=exclude)
        e = df[["b", "e"]]
        tm.assert_frame_equal(r, e)

    def test_select_dtypes_include_using_scalars(self):
        df = DataFrame(
            {
                "a": list("abc"),
                "b": list(range(1, 4)),
                "c": np.arange(3, 6).astype("u1"),
                "d": np.arange(4.0, 7.0, dtype="float64"),
                "e": [True, False, True],
                "f": pd.Categorical(list("abc")),
                "g": pd.date_range("20130101", periods=3),
                "h": pd.date_range("20130101", periods=3, tz="US/Eastern"),
                "i": pd.date_range("20130101", periods=3, tz="CET"),
                "j": pd.period_range("2013-01", periods=3, freq="M"),
                "k": pd.timedelta_range("1 day", periods=3),
            }
        )

        ri = df.select_dtypes(include=np.number)
        ei = df[["b", "c", "d", "k"]]
        tm.assert_frame_equal(ri, ei)

        ri = df.select_dtypes(include="datetime")
        ei = df[["g"]]
        tm.assert_frame_equal(ri, ei)

        ri = df.select_dtypes(include="datetime64")
        ei = df[["g"]]
        tm.assert_frame_equal(ri, ei)

        ri = df.select_dtypes(include="category")
        ei = df[["f"]]
        tm.assert_frame_equal(ri, ei)

        with pytest.raises(NotImplementedError, match=r"^$"):
            df.select_dtypes(include="period")

    def test_select_dtypes_exclude_using_scalars(self):
        df = DataFrame(
            {
                "a": list("abc"),
                "b": list(range(1, 4)),
                "c": np.arange(3, 6).astype("u1"),
                "d": np.arange(4.0, 7.0, dtype="float64"),
                "e": [True, False, True],
                "f": pd.Categorical(list("abc")),
                "g": pd.date_range("20130101", periods=3),
                "h": pd.date_range("20130101", periods=3, tz="US/Eastern"),
                "i": pd.date_range("20130101", periods=3, tz="CET"),
                "j": pd.period_range("2013-01", periods=3, freq="M"),
                "k": pd.timedelta_range("1 day", periods=3),
            }
        )

        ri = df.select_dtypes(exclude=np.number)
        ei = df[["a", "e", "f", "g", "h", "i", "j"]]
        tm.assert_frame_equal(ri, ei)

        ri = df.select_dtypes(exclude="category")
        ei = df[["a", "b", "c", "d", "e", "g", "h", "i", "j", "k"]]
        tm.assert_frame_equal(ri, ei)

        with pytest.raises(NotImplementedError, match=r"^$"):
            df.select_dtypes(exclude="period")

    def test_select_dtypes_include_exclude_using_scalars(self):
        df = DataFrame(
            {
                "a": list("abc"),
                "b": list(range(1, 4)),
                "c": np.arange(3, 6).astype("u1"),
                "d": np.arange(4.0, 7.0, dtype="float64"),
                "e": [True, False, True],
                "f": pd.Categorical(list("abc")),
                "g": pd.date_range("20130101", periods=3),
                "h": pd.date_range("20130101", periods=3, tz="US/Eastern"),
                "i": pd.date_range("20130101", periods=3, tz="CET"),
                "j": pd.period_range("2013-01", periods=3, freq="M"),
                "k": pd.timedelta_range("1 day", periods=3),
            }
        )

        ri = df.select_dtypes(include=np.number, exclude="floating")
        ei = df[["b", "c", "k"]]
        tm.assert_frame_equal(ri, ei)

    def test_select_dtypes_include_exclude_mixed_scalars_lists(self):
        df = DataFrame(
            {
                "a": list("abc"),
                "b": list(range(1, 4)),
                "c": np.arange(3, 6).astype("u1"),
                "d": np.arange(4.0, 7.0, dtype="float64"),
                "e": [True, False, True],
                "f": pd.Categorical(list("abc")),
                "g": pd.date_range("20130101", periods=3),
                "h": pd.date_range("20130101", periods=3, tz="US/Eastern"),
                "i": pd.date_range("20130101", periods=3, tz="CET"),
                "j": pd.period_range("2013-01", periods=3, freq="M"),
                "k": pd.timedelta_range("1 day", periods=3),
            }
        )

        ri = df.select_dtypes(include=np.number, exclude=["floating", "timedelta"])
        ei = df[["b", "c"]]
        tm.assert_frame_equal(ri, ei)

        ri = df.select_dtypes(include=[np.number, "category"], exclude="floating")
        ei = df[["b", "c", "f", "k"]]
        tm.assert_frame_equal(ri, ei)

    def test_select_dtypes_duplicate_columns(self):
        # GH20839
        odict = OrderedDict
        df = DataFrame(
            odict(
                [
                    ("a", list("abc")),
                    ("b", list(range(1, 4))),
                    ("c", np.arange(3, 6).astype("u1")),
                    ("d", np.arange(4.0, 7.0, dtype="float64")),
                    ("e", [True, False, True]),
                    ("f", pd.date_range("now", periods=3).values),
                ]
            )
        )
        df.columns = ["a", "a", "b", "b", "b", "c"]

        expected = DataFrame(
            {"a": list(range(1, 4)), "b": np.arange(3, 6).astype("u1")}
        )

        result = df.select_dtypes(include=[np.number], exclude=["floating"])
        tm.assert_frame_equal(result, expected)

    def test_select_dtypes_not_an_attr_but_still_valid_dtype(self):
        df = DataFrame(
            {
                "a": list("abc"),
                "b": list(range(1, 4)),
                "c": np.arange(3, 6).astype("u1"),
                "d": np.arange(4.0, 7.0, dtype="float64"),
                "e": [True, False, True],
                "f": pd.date_range("now", periods=3).values,
            }
        )
        df["g"] = df.f.diff()
        assert not hasattr(np, "u8")
        r = df.select_dtypes(include=["i8", "O"], exclude=["timedelta"])
        e = df[["a", "b"]]
        tm.assert_frame_equal(r, e)

        r = df.select_dtypes(include=["i8", "O", "timedelta64[ns]"])
        e = df[["a", "b", "g"]]
        tm.assert_frame_equal(r, e)

    def test_select_dtypes_empty(self):
        df = DataFrame({"a": list("abc"), "b": list(range(1, 4))})
        msg = "at least one of include or exclude must be nonempty"
        with pytest.raises(ValueError, match=msg):
            df.select_dtypes()

    def test_select_dtypes_bad_datetime64(self):
        df = DataFrame(
            {
                "a": list("abc"),
                "b": list(range(1, 4)),
                "c": np.arange(3, 6).astype("u1"),
                "d": np.arange(4.0, 7.0, dtype="float64"),
                "e": [True, False, True],
                "f": pd.date_range("now", periods=3).values,
            }
        )
        with pytest.raises(ValueError, match=".+ is too specific"):
            df.select_dtypes(include=["datetime64[D]"])

        with pytest.raises(ValueError, match=".+ is too specific"):
            df.select_dtypes(exclude=["datetime64[as]"])

    def test_select_dtypes_datetime_with_tz(self):

        df2 = DataFrame(
            dict(
                A=Timestamp("20130102", tz="US/Eastern"),
                B=Timestamp("20130603", tz="CET"),
            ),
            index=range(5),
        )
        df3 = pd.concat([df2.A.to_frame(), df2.B.to_frame()], axis=1)
        result = df3.select_dtypes(include=["datetime64[ns]"])
        expected = df3.reindex(columns=[])
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype", [str, "str", np.string_, "S1", "unicode", np.unicode_, "U1"]
    )
    @pytest.mark.parametrize("arg", ["include", "exclude"])
    def test_select_dtypes_str_raises(self, dtype, arg):
        df = DataFrame(
            {
                "a": list("abc"),
                "g": list("abc"),
                "b": list(range(1, 4)),
                "c": np.arange(3, 6).astype("u1"),
                "d": np.arange(4.0, 7.0, dtype="float64"),
                "e": [True, False, True],
                "f": pd.date_range("now", periods=3).values,
            }
        )
        msg = "string dtypes are not allowed"
        kwargs = {arg: [dtype]}

        with pytest.raises(TypeError, match=msg):
            df.select_dtypes(**kwargs)

    def test_select_dtypes_bad_arg_raises(self):
        df = DataFrame(
            {
                "a": list("abc"),
                "g": list("abc"),
                "b": list(range(1, 4)),
                "c": np.arange(3, 6).astype("u1"),
                "d": np.arange(4.0, 7.0, dtype="float64"),
                "e": [True, False, True],
                "f": pd.date_range("now", periods=3).values,
            }
        )

        msg = "data type.*not understood"
        with pytest.raises(TypeError, match=msg):
            df.select_dtypes(["blargy, blarg, blarg"])

    def test_select_dtypes_typecodes(self):
        # GH 11990
        df = tm.makeCustomDataframe(30, 3, data_gen_f=lambda x, y: np.random.random())
        expected = df
        FLOAT_TYPES = list(np.typecodes["AllFloat"])
        tm.assert_frame_equal(df.select_dtypes(FLOAT_TYPES), expected)

    def test_dtypes_gh8722(self, float_string_frame):
        float_string_frame["bool"] = float_string_frame["A"] > 0
        result = float_string_frame.dtypes
        expected = Series(
            {k: v.dtype for k, v in float_string_frame.items()}, index=result.index
        )
        tm.assert_series_equal(result, expected)

        # compat, GH 8722
        with option_context("use_inf_as_na", True):
            df = DataFrame([[1]])
            result = df.dtypes
            tm.assert_series_equal(result, Series({0: np.dtype("int64")}))

    def test_astype_float(self, float_frame):
        casted = float_frame.astype(int)
        expected = DataFrame(
            float_frame.values.astype(int),
            index=float_frame.index,
            columns=float_frame.columns,
        )
        tm.assert_frame_equal(casted, expected)

        casted = float_frame.astype(np.int32)
        expected = DataFrame(
            float_frame.values.astype(np.int32),
            index=float_frame.index,
            columns=float_frame.columns,
        )
        tm.assert_frame_equal(casted, expected)

        float_frame["foo"] = "5"
        casted = float_frame.astype(int)
        expected = DataFrame(
            float_frame.values.astype(int),
            index=float_frame.index,
            columns=float_frame.columns,
        )
        tm.assert_frame_equal(casted, expected)

    def test_astype_mixed_float(self, mixed_float_frame):
        # mixed casting
        casted = mixed_float_frame.reindex(columns=["A", "B"]).astype("float32")
        _check_cast(casted, "float32")

        casted = mixed_float_frame.reindex(columns=["A", "B"]).astype("float16")
        _check_cast(casted, "float16")

    def test_astype_mixed_type(self, mixed_type_frame):
        # mixed casting
        mn = mixed_type_frame._get_numeric_data().copy()
        mn["little_float"] = np.array(12345.0, dtype="float16")
        mn["big_float"] = np.array(123456789101112.0, dtype="float64")

        casted = mn.astype("float64")
        _check_cast(casted, "float64")

        casted = mn.astype("int64")
        _check_cast(casted, "int64")

        casted = mn.reindex(columns=["little_float"]).astype("float16")
        _check_cast(casted, "float16")

        casted = mn.astype("float32")
        _check_cast(casted, "float32")

        casted = mn.astype("int32")
        _check_cast(casted, "int32")

        # to object
        casted = mn.astype("O")
        _check_cast(casted, "object")

    def test_astype_with_exclude_string(self, float_frame):
        df = float_frame.copy()
        expected = float_frame.astype(int)
        df["string"] = "foo"
        casted = df.astype(int, errors="ignore")

        expected["string"] = "foo"
        tm.assert_frame_equal(casted, expected)

        df = float_frame.copy()
        expected = float_frame.astype(np.int32)
        df["string"] = "foo"
        casted = df.astype(np.int32, errors="ignore")

        expected["string"] = "foo"
        tm.assert_frame_equal(casted, expected)

    def test_astype_with_view_float(self, float_frame):

        # this is the only real reason to do it this way
        tf = np.round(float_frame).astype(np.int32)
        casted = tf.astype(np.float32, copy=False)

        # TODO(wesm): verification?
        tf = float_frame.astype(np.float64)
        casted = tf.astype(np.int64, copy=False)  # noqa

    def test_astype_with_view_mixed_float(self, mixed_float_frame):

        tf = mixed_float_frame.reindex(columns=["A", "B", "C"])

        casted = tf.astype(np.int64)
        casted = tf.astype(np.float32)  # noqa

    @pytest.mark.parametrize("dtype", [np.int32, np.int64])
    @pytest.mark.parametrize("val", [np.nan, np.inf])
    def test_astype_cast_nan_inf_int(self, val, dtype):
        # see gh-14265
        #
        # Check NaN and inf --> raise error when converting to int.
        msg = "Cannot convert non-finite values \\(NA or inf\\) to integer"
        df = DataFrame([val])

        with pytest.raises(ValueError, match=msg):
            df.astype(dtype)

    def test_astype_str(self):
        # see gh-9757
        a = Series(date_range("2010-01-04", periods=5))
        b = Series(date_range("3/6/2012 00:00", periods=5, tz="US/Eastern"))
        c = Series([Timedelta(x, unit="d") for x in range(5)])
        d = Series(range(5))
        e = Series([0.0, 0.2, 0.4, 0.6, 0.8])

        df = DataFrame({"a": a, "b": b, "c": c, "d": d, "e": e})

        # Datetime-like
        result = df.astype(str)

        expected = DataFrame(
            {
                "a": list(map(str, map(lambda x: Timestamp(x)._date_repr, a._values))),
                "b": list(map(str, map(Timestamp, b._values))),
                "c": list(
                    map(
                        str,
                        map(lambda x: Timedelta(x)._repr_base(format="all"), c._values),
                    )
                ),
                "d": list(map(str, d._values)),
                "e": list(map(str, e._values)),
            }
        )

        tm.assert_frame_equal(result, expected)

    def test_astype_str_float(self):
        # see gh-11302
        result = DataFrame([np.NaN]).astype(str)
        expected = DataFrame(["nan"])

        tm.assert_frame_equal(result, expected)
        result = DataFrame([1.12345678901234567890]).astype(str)

        # < 1.14 truncates
        # >= 1.14 preserves the full repr
        val = "1.12345678901" if _np_version_under1p14 else "1.1234567890123457"
        expected = DataFrame([val])
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("dtype_class", [dict, Series])
    def test_astype_dict_like(self, dtype_class):
        # GH7271 & GH16717
        a = Series(date_range("2010-01-04", periods=5))
        b = Series(range(5))
        c = Series([0.0, 0.2, 0.4, 0.6, 0.8])
        d = Series(["1.0", "2", "3.14", "4", "5.4"])
        df = DataFrame({"a": a, "b": b, "c": c, "d": d})
        original = df.copy(deep=True)

        # change type of a subset of columns
        dt1 = dtype_class({"b": "str", "d": "float32"})
        result = df.astype(dt1)
        expected = DataFrame(
            {
                "a": a,
                "b": Series(["0", "1", "2", "3", "4"]),
                "c": c,
                "d": Series([1.0, 2.0, 3.14, 4.0, 5.4], dtype="float32"),
            }
        )
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(df, original)

        dt2 = dtype_class({"b": np.float32, "c": "float32", "d": np.float64})
        result = df.astype(dt2)
        expected = DataFrame(
            {
                "a": a,
                "b": Series([0.0, 1.0, 2.0, 3.0, 4.0], dtype="float32"),
                "c": Series([0.0, 0.2, 0.4, 0.6, 0.8], dtype="float32"),
                "d": Series([1.0, 2.0, 3.14, 4.0, 5.4], dtype="float64"),
            }
        )
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(df, original)

        # change all columns
        dt3 = dtype_class({"a": str, "b": str, "c": str, "d": str})
        tm.assert_frame_equal(df.astype(dt3), df.astype(str))
        tm.assert_frame_equal(df, original)

        # error should be raised when using something other than column labels
        # in the keys of the dtype dict
        dt4 = dtype_class({"b": str, 2: str})
        dt5 = dtype_class({"e": str})
        msg = "Only a column name can be used for the key in a dtype mappings argument"
        with pytest.raises(KeyError, match=msg):
            df.astype(dt4)
        with pytest.raises(KeyError, match=msg):
            df.astype(dt5)
        tm.assert_frame_equal(df, original)

        # if the dtypes provided are the same as the original dtypes, the
        # resulting DataFrame should be the same as the original DataFrame
        dt6 = dtype_class({col: df[col].dtype for col in df.columns})
        equiv = df.astype(dt6)
        tm.assert_frame_equal(df, equiv)
        tm.assert_frame_equal(df, original)

        # GH 16717
        # if dtypes provided is empty, the resulting DataFrame
        # should be the same as the original DataFrame
        dt7 = dtype_class({}) if dtype_class is dict else dtype_class({}, dtype=object)
        equiv = df.astype(dt7)
        tm.assert_frame_equal(df, equiv)
        tm.assert_frame_equal(df, original)

    def test_astype_duplicate_col(self):
        a1 = Series([1, 2, 3, 4, 5], name="a")
        b = Series([0.1, 0.2, 0.4, 0.6, 0.8], name="b")
        a2 = Series([0, 1, 2, 3, 4], name="a")
        df = concat([a1, b, a2], axis=1)

        result = df.astype(str)
        a1_str = Series(["1", "2", "3", "4", "5"], dtype="str", name="a")
        b_str = Series(["0.1", "0.2", "0.4", "0.6", "0.8"], dtype=str, name="b")
        a2_str = Series(["0", "1", "2", "3", "4"], dtype="str", name="a")
        expected = concat([a1_str, b_str, a2_str], axis=1)
        tm.assert_frame_equal(result, expected)

        result = df.astype({"a": "str"})
        expected = concat([a1_str, b, a2_str], axis=1)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype",
        [
            "category",
            CategoricalDtype(),
            CategoricalDtype(ordered=True),
            CategoricalDtype(ordered=False),
            CategoricalDtype(categories=list("abcdef")),
            CategoricalDtype(categories=list("edba"), ordered=False),
            CategoricalDtype(categories=list("edcb"), ordered=True),
        ],
        ids=repr,
    )
    def test_astype_categorical(self, dtype):
        # GH 18099
        d = {"A": list("abbc"), "B": list("bccd"), "C": list("cdde")}
        df = DataFrame(d)
        result = df.astype(dtype)
        expected = DataFrame({k: Categorical(d[k], dtype=dtype) for k in d})
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("cls", [CategoricalDtype, DatetimeTZDtype, IntervalDtype])
    def test_astype_categoricaldtype_class_raises(self, cls):
        df = DataFrame({"A": ["a", "a", "b", "c"]})
        xpr = "Expected an instance of {}".format(cls.__name__)
        with pytest.raises(TypeError, match=xpr):
            df.astype({"A": cls})

        with pytest.raises(TypeError, match=xpr):
            df["A"].astype(cls)

    def test_singlerow_slice_categoricaldtype_gives_series(self):
        # GH29521
        df = pd.DataFrame({"x": pd.Categorical("a b c d e".split())})
        result = df.iloc[0]
        raw_cat = pd.Categorical(["a"], categories=["a", "b", "c", "d", "e"])
        expected = pd.Series(raw_cat, index=["x"], name=0, dtype="category")

        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("dtype", ["Int64", "Int32", "Int16"])
    def test_astype_extension_dtypes(self, dtype):
        # GH 22578
        df = pd.DataFrame([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], columns=["a", "b"])

        expected1 = pd.DataFrame(
            {
                "a": integer_array([1, 3, 5], dtype=dtype),
                "b": integer_array([2, 4, 6], dtype=dtype),
            }
        )
        tm.assert_frame_equal(df.astype(dtype), expected1)
        tm.assert_frame_equal(df.astype("int64").astype(dtype), expected1)
        tm.assert_frame_equal(df.astype(dtype).astype("float64"), df)

        df = pd.DataFrame([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], columns=["a", "b"])
        df["b"] = df["b"].astype(dtype)
        expected2 = pd.DataFrame(
            {"a": [1.0, 3.0, 5.0], "b": integer_array([2, 4, 6], dtype=dtype)}
        )
        tm.assert_frame_equal(df, expected2)

        tm.assert_frame_equal(df.astype(dtype), expected1)
        tm.assert_frame_equal(df.astype("int64").astype(dtype), expected1)

    @pytest.mark.parametrize("dtype", ["Int64", "Int32", "Int16"])
    def test_astype_extension_dtypes_1d(self, dtype):
        # GH 22578
        df = pd.DataFrame({"a": [1.0, 2.0, 3.0]})

        expected1 = pd.DataFrame({"a": integer_array([1, 2, 3], dtype=dtype)})
        tm.assert_frame_equal(df.astype(dtype), expected1)
        tm.assert_frame_equal(df.astype("int64").astype(dtype), expected1)

        df = pd.DataFrame({"a": [1.0, 2.0, 3.0]})
        df["a"] = df["a"].astype(dtype)
        expected2 = pd.DataFrame({"a": integer_array([1, 2, 3], dtype=dtype)})
        tm.assert_frame_equal(df, expected2)

        tm.assert_frame_equal(df.astype(dtype), expected1)
        tm.assert_frame_equal(df.astype("int64").astype(dtype), expected1)

    @pytest.mark.parametrize("dtype", ["category", "Int64"])
    def test_astype_extension_dtypes_duplicate_col(self, dtype):
        # GH 24704
        a1 = Series([0, np.nan, 4], name="a")
        a2 = Series([np.nan, 3, 5], name="a")
        df = concat([a1, a2], axis=1)

        result = df.astype(dtype)
        expected = concat([a1.astype(dtype), a2.astype(dtype)], axis=1)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("kwargs", [dict(), dict(other=None)])
    def test_df_where_with_category(self, kwargs):
        # GH 16979
        df = DataFrame(np.arange(2 * 3).reshape(2, 3), columns=list("ABC"))
        mask = np.array([[True, False, True], [False, True, True]])

        # change type to category
        df.A = df.A.astype("category")
        df.B = df.B.astype("category")
        df.C = df.C.astype("category")

        result = df.A.where(mask[:, 0], **kwargs)
        expected = Series(pd.Categorical([0, np.nan], categories=[0, 3]), name="A")

        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype", [{100: "float64", 200: "uint64"}, "category", "float64"]
    )
    def test_astype_column_metadata(self, dtype):
        # GH 19920
        columns = pd.UInt64Index([100, 200, 300], name="foo")
        df = DataFrame(np.arange(15).reshape(5, 3), columns=columns)
        df = df.astype(dtype)
        tm.assert_index_equal(df.columns, columns)

    def test_df_where_change_dtype(self):
        # GH 16979
        df = DataFrame(np.arange(2 * 3).reshape(2, 3), columns=list("ABC"))
        mask = np.array([[True, False, False], [False, False, True]])

        result = df.where(mask)
        expected = DataFrame(
            [[0, np.nan, np.nan], [np.nan, np.nan, 5]], columns=list("ABC")
        )

        tm.assert_frame_equal(result, expected)

        # change type to category
        df.A = df.A.astype("category")
        df.B = df.B.astype("category")
        df.C = df.C.astype("category")

        result = df.where(mask)
        A = pd.Categorical([0, np.nan], categories=[0, 3])
        B = pd.Categorical([np.nan, np.nan], categories=[1, 4])
        C = pd.Categorical([np.nan, 5], categories=[2, 5])
        expected = DataFrame({"A": A, "B": B, "C": C})

        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("dtype", ["M8", "m8"])
    @pytest.mark.parametrize("unit", ["ns", "us", "ms", "s", "h", "m", "D"])
    def test_astype_from_datetimelike_to_objectt(self, dtype, unit):
        # tests astype to object dtype
        # gh-19223 / gh-12425
        dtype = "{}[{}]".format(dtype, unit)
        arr = np.array([[1, 2, 3]], dtype=dtype)
        df = DataFrame(arr)
        result = df.astype(object)
        assert (result.dtypes == object).all()

        if dtype.startswith("M8"):
            assert result.iloc[0, 0] == pd.to_datetime(1, unit=unit)
        else:
            assert result.iloc[0, 0] == pd.to_timedelta(1, unit=unit)

    @pytest.mark.parametrize("arr_dtype", [np.int64, np.float64])
    @pytest.mark.parametrize("dtype", ["M8", "m8"])
    @pytest.mark.parametrize("unit", ["ns", "us", "ms", "s", "h", "m", "D"])
    def test_astype_to_datetimelike_unit(self, arr_dtype, dtype, unit):
        # tests all units from numeric origination
        # gh-19223 / gh-12425
        dtype = "{}[{}]".format(dtype, unit)
        arr = np.array([[1, 2, 3]], dtype=arr_dtype)
        df = DataFrame(arr)
        result = df.astype(dtype)
        expected = DataFrame(arr.astype(dtype))

        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("unit", ["ns", "us", "ms", "s", "h", "m", "D"])
    def test_astype_to_datetime_unit(self, unit):
        # tests all units from datetime origination
        # gh-19223
        dtype = "M8[{}]".format(unit)
        arr = np.array([[1, 2, 3]], dtype=dtype)
        df = DataFrame(arr)
        result = df.astype(dtype)
        expected = DataFrame(arr.astype(dtype))

        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("unit", ["ns"])
    def test_astype_to_timedelta_unit_ns(self, unit):
        # preserver the timedelta conversion
        # gh-19223
        dtype = "m8[{}]".format(unit)
        arr = np.array([[1, 2, 3]], dtype=dtype)
        df = DataFrame(arr)
        result = df.astype(dtype)
        expected = DataFrame(arr.astype(dtype))

        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("unit", ["us", "ms", "s", "h", "m", "D"])
    def test_astype_to_timedelta_unit(self, unit):
        # coerce to float
        # gh-19223
        dtype = "m8[{}]".format(unit)
        arr = np.array([[1, 2, 3]], dtype=dtype)
        df = DataFrame(arr)
        result = df.astype(dtype)
        expected = DataFrame(df.values.astype(dtype).astype(float))

        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("unit", ["ns", "us", "ms", "s", "h", "m", "D"])
    def test_astype_to_incorrect_datetimelike(self, unit):
        # trying to astype a m to a M, or vice-versa
        # gh-19224
        dtype = "M8[{}]".format(unit)
        other = "m8[{}]".format(unit)

        df = DataFrame(np.array([[1, 2, 3]], dtype=dtype))
        msg = (
            r"cannot astype a datetimelike from \[datetime64\[ns\]\] to "
            r"\[timedelta64\[{}\]\]"
        ).format(unit)
        with pytest.raises(TypeError, match=msg):
            df.astype(other)

        msg = (
            r"cannot astype a timedelta from \[timedelta64\[ns\]\] to "
            r"\[datetime64\[{}\]\]"
        ).format(unit)
        df = DataFrame(np.array([[1, 2, 3]], dtype=other))
        with pytest.raises(TypeError, match=msg):
            df.astype(dtype)

    def test_timedeltas(self):
        df = DataFrame(
            dict(
                A=Series(date_range("2012-1-1", periods=3, freq="D")),
                B=Series([timedelta(days=i) for i in range(3)]),
            )
        )
        result = df.dtypes
        expected = Series(
            [np.dtype("datetime64[ns]"), np.dtype("timedelta64[ns]")], index=list("AB")
        )
        tm.assert_series_equal(result, expected)

        df["C"] = df["A"] + df["B"]
        result = df.dtypes
        expected = Series(
            [
                np.dtype("datetime64[ns]"),
                np.dtype("timedelta64[ns]"),
                np.dtype("datetime64[ns]"),
            ],
            index=list("ABC"),
        )
        tm.assert_series_equal(result, expected)

        # mixed int types
        df["D"] = 1
        result = df.dtypes
        expected = Series(
            [
                np.dtype("datetime64[ns]"),
                np.dtype("timedelta64[ns]"),
                np.dtype("datetime64[ns]"),
                np.dtype("int64"),
            ],
            index=list("ABCD"),
        )
        tm.assert_series_equal(result, expected)

    def test_arg_for_errors_in_astype(self):
        # issue #14878

        df = DataFrame([1, 2, 3])

        with pytest.raises(ValueError):
            df.astype(np.float64, errors=True)

        df.astype(np.int8, errors="ignore")

    def test_arg_for_errors_in_astype_dictlist(self):
        # GH-25905
        df = pd.DataFrame(
            [
                {"a": "1", "b": "16.5%", "c": "test"},
                {"a": "2.2", "b": "15.3", "c": "another_test"},
            ]
        )
        expected = pd.DataFrame(
            [
                {"a": 1.0, "b": "16.5%", "c": "test"},
                {"a": 2.2, "b": "15.3", "c": "another_test"},
            ]
        )
        type_dict = {"a": "float64", "b": "float64", "c": "object"}

        result = df.astype(dtype=type_dict, errors="ignore")

        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "input_vals",
        [
            ([1, 2]),
            (["1", "2"]),
            (list(pd.date_range("1/1/2011", periods=2, freq="H"))),
            (list(pd.date_range("1/1/2011", periods=2, freq="H", tz="US/Eastern"))),
            ([pd.Interval(left=0, right=5)]),
        ],
    )
    def test_constructor_list_str(self, input_vals, string_dtype):
        # GH 16605
        # Ensure that data elements are converted to strings when
        # dtype is str, 'str', or 'U'

        result = DataFrame({"A": input_vals}, dtype=string_dtype)
        expected = DataFrame({"A": input_vals}).astype({"A": string_dtype})
        tm.assert_frame_equal(result, expected)

    def test_constructor_list_str_na(self, string_dtype):

        result = DataFrame({"A": [1.0, 2.0, None]}, dtype=string_dtype)
        expected = DataFrame({"A": ["1.0", "2.0", None]}, dtype=object)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "data, expected",
        [
            # empty
            (DataFrame(), True),
            # multi-same
            (DataFrame({"A": [1, 2], "B": [1, 2]}), True),
            # multi-object
            (
                DataFrame(
                    {
                        "A": np.array([1, 2], dtype=object),
                        "B": np.array(["a", "b"], dtype=object),
                    }
                ),
                True,
            ),
            # multi-extension
            (
                DataFrame(
                    {"A": pd.Categorical(["a", "b"]), "B": pd.Categorical(["a", "b"])}
                ),
                True,
            ),
            # differ types
            (DataFrame({"A": [1, 2], "B": [1.0, 2.0]}), False),
            # differ sizes
            (
                DataFrame(
                    {
                        "A": np.array([1, 2], dtype=np.int32),
                        "B": np.array([1, 2], dtype=np.int64),
                    }
                ),
                False,
            ),
            # multi-extension differ
            (
                DataFrame(
                    {"A": pd.Categorical(["a", "b"]), "B": pd.Categorical(["b", "c"])}
                ),
                False,
            ),
        ],
    )
    def test_is_homogeneous_type(self, data, expected):
        assert data._is_homogeneous_type is expected

    def test_asarray_homogenous(self):
        df = pd.DataFrame({"A": pd.Categorical([1, 2]), "B": pd.Categorical([1, 2])})
        result = np.asarray(df)
        # may change from object in the future
        expected = np.array([[1, 1], [2, 2]], dtype="object")
        tm.assert_numpy_array_equal(result, expected)

    def test_str_to_small_float_conversion_type(self):
        # GH 20388
        np.random.seed(13)
        col_data = [str(np.random.random() * 1e-12) for _ in range(5)]
        result = pd.DataFrame(col_data, columns=["A"])
        expected = pd.DataFrame(col_data, columns=["A"], dtype=object)
        tm.assert_frame_equal(result, expected)
        # change the dtype of the elements from object to float one by one
        result.loc[result.index, "A"] = [float(x) for x in col_data]
        expected = pd.DataFrame(col_data, columns=["A"], dtype=float)
        tm.assert_frame_equal(result, expected)


class TestDataFrameDatetimeWithTZ:
    def test_interleave(self, timezone_frame):

        # interleave with object
        result = timezone_frame.assign(D="foo").values
        expected = np.array(
            [
                [
                    Timestamp("2013-01-01 00:00:00"),
                    Timestamp("2013-01-02 00:00:00"),
                    Timestamp("2013-01-03 00:00:00"),
                ],
                [
                    Timestamp("2013-01-01 00:00:00-0500", tz="US/Eastern"),
                    pd.NaT,
                    Timestamp("2013-01-03 00:00:00-0500", tz="US/Eastern"),
                ],
                [
                    Timestamp("2013-01-01 00:00:00+0100", tz="CET"),
                    pd.NaT,
                    Timestamp("2013-01-03 00:00:00+0100", tz="CET"),
                ],
                ["foo", "foo", "foo"],
            ],
            dtype=object,
        ).T
        tm.assert_numpy_array_equal(result, expected)

        # interleave with only datetime64[ns]
        result = timezone_frame.values
        expected = np.array(
            [
                [
                    Timestamp("2013-01-01 00:00:00"),
                    Timestamp("2013-01-02 00:00:00"),
                    Timestamp("2013-01-03 00:00:00"),
                ],
                [
                    Timestamp("2013-01-01 00:00:00-0500", tz="US/Eastern"),
                    pd.NaT,
                    Timestamp("2013-01-03 00:00:00-0500", tz="US/Eastern"),
                ],
                [
                    Timestamp("2013-01-01 00:00:00+0100", tz="CET"),
                    pd.NaT,
                    Timestamp("2013-01-03 00:00:00+0100", tz="CET"),
                ],
            ],
            dtype=object,
        ).T
        tm.assert_numpy_array_equal(result, expected)

    def test_astype(self, timezone_frame):
        # astype
        expected = np.array(
            [
                [
                    Timestamp("2013-01-01 00:00:00"),
                    Timestamp("2013-01-02 00:00:00"),
                    Timestamp("2013-01-03 00:00:00"),
                ],
                [
                    Timestamp("2013-01-01 00:00:00-0500", tz="US/Eastern"),
                    pd.NaT,
                    Timestamp("2013-01-03 00:00:00-0500", tz="US/Eastern"),
                ],
                [
                    Timestamp("2013-01-01 00:00:00+0100", tz="CET"),
                    pd.NaT,
                    Timestamp("2013-01-03 00:00:00+0100", tz="CET"),
                ],
            ],
            dtype=object,
        ).T
        expected = DataFrame(
            expected,
            index=timezone_frame.index,
            columns=timezone_frame.columns,
            dtype=object,
        )
        result = timezone_frame.astype(object)
        tm.assert_frame_equal(result, expected)

        result = timezone_frame.astype("datetime64[ns]")
        expected = DataFrame(
            {
                "A": date_range("20130101", periods=3),
                "B": (
                    date_range("20130101", periods=3, tz="US/Eastern")
                    .tz_convert("UTC")
                    .tz_localize(None)
                ),
                "C": (
                    date_range("20130101", periods=3, tz="CET")
                    .tz_convert("UTC")
                    .tz_localize(None)
                ),
            }
        )
        expected.iloc[1, 1] = pd.NaT
        expected.iloc[1, 2] = pd.NaT
        tm.assert_frame_equal(result, expected)

    def test_astype_str(self, timezone_frame):
        # str formatting
        result = timezone_frame.astype(str)
        expected = DataFrame(
            [
                [
                    "2013-01-01",
                    "2013-01-01 00:00:00-05:00",
                    "2013-01-01 00:00:00+01:00",
                ],
                ["2013-01-02", "NaT", "NaT"],
                [
                    "2013-01-03",
                    "2013-01-03 00:00:00-05:00",
                    "2013-01-03 00:00:00+01:00",
                ],
            ],
            columns=timezone_frame.columns,
        )
        tm.assert_frame_equal(result, expected)

        with option_context("display.max_columns", 20):
            result = str(timezone_frame)
            assert (
                "0 2013-01-01 2013-01-01 00:00:00-05:00 2013-01-01 00:00:00+01:00"
            ) in result
            assert (
                "1 2013-01-02                       NaT                       NaT"
            ) in result
            assert (
                "2 2013-01-03 2013-01-03 00:00:00-05:00 2013-01-03 00:00:00+01:00"
            ) in result
