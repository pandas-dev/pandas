from collections import OrderedDict

import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Timestamp
import pandas._testing as tm


class TestSelectDtypes:
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
