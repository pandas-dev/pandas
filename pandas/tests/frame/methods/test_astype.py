import re

import numpy as np
import pytest

from pandas.errors import Pandas4Warning
import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm


def _check_cast(df, v):
    """
    Check if all dtypes of df are equal to v
    """
    assert all(s.dtype.name == v for _, s in df.items())


class TestAstype:
    def test_astype_float(self, float_frame):
        casted = float_frame.astype(int)
        expected = pd.DataFrame(
            float_frame.values.astype(int),
            index=float_frame.index,
            columns=float_frame.columns,
        )
        tm.assert_frame_equal(casted, expected)

        casted = float_frame.astype(np.int32)
        expected = pd.DataFrame(
            float_frame.values.astype(np.int32),
            index=float_frame.index,
            columns=float_frame.columns,
        )
        tm.assert_frame_equal(casted, expected)

        float_frame["foo"] = "5"
        casted = float_frame.astype(int)
        expected = pd.DataFrame(
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

    def test_astype_mixed_type(self):
        # mixed casting
        df = pd.DataFrame(
            {
                "a": 1.0,
                "b": 2,
                "c": "foo",
                "float32": np.array([1.0] * 10, dtype="float32"),
                "int32": np.array([1] * 10, dtype="int32"),
            },
            index=np.arange(10),
        )
        mn = df._get_numeric_data().copy()
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
        tf.astype(np.float32)

        # TODO(wesm): verification?
        tf = float_frame.astype(np.float64)
        tf.astype(np.int64)

    def test_astype_with_view_mixed_float(self, mixed_float_frame):
        tf = mixed_float_frame.reindex(columns=["A", "B", "C"])

        tf.astype(np.int64)
        tf.astype(np.float32)

    @pytest.mark.parametrize("val", [np.nan, np.inf])
    def test_astype_cast_nan_inf_int(self, val, any_int_numpy_dtype):
        # see GH#14265
        #
        # Check NaN and inf --> raise error when converting to int.
        msg = "Cannot convert non-finite values \\(NA or inf\\) to integer"
        df = pd.DataFrame([val])

        with pytest.raises(ValueError, match=msg):
            df.astype(any_int_numpy_dtype)

    def test_astype_str(self):
        # see GH#9757
        a = pd.Series(pd.date_range("2010-01-04", periods=5))
        b = pd.Series(pd.date_range("3/6/2012 00:00", periods=5, tz="US/Eastern"))
        c = pd.Series([pd.Timedelta(x, unit="D") for x in range(5)])
        d = pd.Series(range(5))
        e = pd.Series([0.0, 0.2, 0.4, 0.6, 0.8])

        df = pd.DataFrame({"a": a, "b": b, "c": c, "d": d, "e": e})

        # Datetime-like
        result = df.astype(str)

        expected = pd.DataFrame(
            {
                "a": list(map(str, (pd.Timestamp(x)._date_repr for x in a._values))),
                "b": list(map(str, map(pd.Timestamp, b._values))),
                "c": [pd.Timedelta(x)._repr_base() for x in c._values],
                "d": list(map(str, d._values)),
                "e": list(map(str, e._values)),
            },
            dtype="str",
        )

        tm.assert_frame_equal(result, expected)

    def test_astype_str_float(self, using_infer_string):
        # see GH#11302
        result = pd.DataFrame([np.nan]).astype(str)
        expected = pd.DataFrame([np.nan if using_infer_string else "nan"], dtype="str")

        tm.assert_frame_equal(result, expected)
        result = pd.DataFrame([1.12345678901234567890]).astype(str)

        val = "1.1234567890123457"
        expected = pd.DataFrame([val], dtype="str")
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("dtype_class", [dict, pd.Series])
    def test_astype_dict_like(self, dtype_class):
        # GH7271 & GH16717
        a = pd.Series(pd.date_range("2010-01-04", periods=5))
        b = pd.Series(range(5))
        c = pd.Series([0.0, 0.2, 0.4, 0.6, 0.8])
        d = pd.Series(["1.0", "2", "3.14", "4", "5.4"])
        df = pd.DataFrame({"a": a, "b": b, "c": c, "d": d})
        original = df.copy(deep=True)

        # change type of a subset of columns
        dt1 = dtype_class({"b": "str", "d": "float32"})
        result = df.astype(dt1)
        expected = pd.DataFrame(
            {
                "a": a,
                "b": pd.Series(["0", "1", "2", "3", "4"], dtype="str"),
                "c": c,
                "d": pd.Series([1.0, 2.0, 3.14, 4.0, 5.4], dtype="float32"),
            }
        )
        tm.assert_frame_equal(result, expected)
        tm.assert_frame_equal(df, original)

        dt2 = dtype_class({"b": np.float32, "c": "float32", "d": np.float64})
        result = df.astype(dt2)
        expected = pd.DataFrame(
            {
                "a": a,
                "b": pd.Series([0.0, 1.0, 2.0, 3.0, 4.0], dtype="float32"),
                "c": pd.Series([0.0, 0.2, 0.4, 0.6, 0.8], dtype="float32"),
                "d": pd.Series([1.0, 2.0, 3.14, 4.0, 5.4], dtype="float64"),
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
        msg_frame = (
            "Only a column name can be used for the key in a dtype mappings argument. "
            "'{}' not found in columns."
        )
        with pytest.raises(KeyError, match=msg_frame.format(2)):
            df.astype(dt4)
        with pytest.raises(KeyError, match=msg_frame.format("e")):
            df.astype(dt5)
        tm.assert_frame_equal(df, original)

        # if the dtypes provided are the same as the original dtypes, the
        # resulting DataFrame should be the same as the original DataFrame
        dt6 = dtype_class({col: df[col].dtype for col in df.columns})
        equiv = df.astype(dt6)
        tm.assert_frame_equal(df, equiv)
        tm.assert_frame_equal(df, original)

        # GH#16717
        # if dtypes provided is empty, the resulting DataFrame
        # should be the same as the original DataFrame
        dt7 = dtype_class({}) if dtype_class is dict else dtype_class({}, dtype=object)
        equiv = df.astype(dt7)
        tm.assert_frame_equal(df, equiv)
        tm.assert_frame_equal(df, original)

    def test_astype_duplicate_col(self):
        a1 = pd.Series([1, 2, 3, 4, 5], name="a")
        b = pd.Series([0.1, 0.2, 0.4, 0.6, 0.8], name="b")
        a2 = pd.Series([0, 1, 2, 3, 4], name="a")
        df = pd.concat([a1, b, a2], axis=1)

        result = df.astype("str")
        a1_str = pd.Series(["1", "2", "3", "4", "5"], dtype="str", name="a")
        b_str = pd.Series(["0.1", "0.2", "0.4", "0.6", "0.8"], dtype="str", name="b")
        a2_str = pd.Series(["0", "1", "2", "3", "4"], dtype="str", name="a")
        expected = pd.concat([a1_str, b_str, a2_str], axis=1)
        tm.assert_frame_equal(result, expected)

        result = df.astype({"a": "str"})
        expected = pd.concat([a1_str, b, a2_str], axis=1)
        tm.assert_frame_equal(result, expected)

    def test_astype_duplicate_col_series_arg(self):
        # GH#44417
        vals = np.random.default_rng(2).standard_normal((3, 4))
        df = pd.DataFrame(vals, columns=["A", "B", "C", "A"])
        dtypes = df.dtypes
        dtypes.iloc[0] = str
        dtypes.iloc[2] = "Float64"

        result = df.astype(dtypes)
        expected = pd.DataFrame(
            {
                0: pd.Series(vals[:, 0].astype(str), dtype="str"),
                1: vals[:, 1],
                2: pd.array(vals[:, 2], dtype="Float64"),
                3: vals[:, 3],
            }
        )
        expected.columns = df.columns
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype",
        [
            "category",
            pd.CategoricalDtype(),
            pd.CategoricalDtype(ordered=True),
            pd.CategoricalDtype(ordered=False),
            pd.CategoricalDtype(categories=list("abcdef")),
            pd.CategoricalDtype(categories=list("edba"), ordered=False),
            pd.CategoricalDtype(categories=list("edcb"), ordered=True),
        ],
        ids=repr,
    )
    @pytest.mark.filterwarnings(
        "ignore:Constructing a Categorical with a dtype and values"
    )
    def test_astype_categorical(self, dtype):
        # GH#18099
        d = {"A": list("abbc"), "B": list("bccd"), "C": list("cdde")}
        df = pd.DataFrame(d)
        result = df.astype(dtype)
        expected = pd.DataFrame(
            {k: pd.Categorical(v, dtype=dtype) for k, v in d.items()}
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "cls", [pd.CategoricalDtype, pd.DatetimeTZDtype, pd.IntervalDtype]
    )
    def test_astype_categoricaldtype_class_raises(self, cls):
        df = pd.DataFrame({"A": ["a", "a", "b", "c"]})
        xpr = f"Expected an instance of {cls.__name__}"
        with pytest.raises(TypeError, match=xpr):
            df.astype({"A": cls})

        with pytest.raises(TypeError, match=xpr):
            df["A"].astype(cls)

    def test_astype_extension_dtypes(self, any_int_ea_dtype):
        # GH#22578
        dtype = any_int_ea_dtype
        df = pd.DataFrame([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], columns=["a", "b"])

        expected1 = pd.DataFrame(
            {
                "a": pd.array([1, 3, 5], dtype=dtype),
                "b": pd.array([2, 4, 6], dtype=dtype),
            }
        )
        tm.assert_frame_equal(df.astype(dtype), expected1)
        tm.assert_frame_equal(df.astype("int64").astype(dtype), expected1)
        tm.assert_frame_equal(df.astype(dtype).astype("float64"), df)

        df = pd.DataFrame([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], columns=["a", "b"])
        df["b"] = df["b"].astype(dtype)
        expected2 = pd.DataFrame(
            {"a": [1.0, 3.0, 5.0], "b": pd.array([2, 4, 6], dtype=dtype)}
        )
        tm.assert_frame_equal(df, expected2)

        tm.assert_frame_equal(df.astype(dtype), expected1)
        tm.assert_frame_equal(df.astype("int64").astype(dtype), expected1)

    def test_astype_extension_dtypes_1d(self, any_int_ea_dtype):
        # GH#22578
        dtype = any_int_ea_dtype
        df = pd.DataFrame({"a": [1.0, 2.0, 3.0]})

        expected1 = pd.DataFrame({"a": pd.array([1, 2, 3], dtype=dtype)})
        tm.assert_frame_equal(df.astype(dtype), expected1)
        tm.assert_frame_equal(df.astype("int64").astype(dtype), expected1)

        df = pd.DataFrame({"a": [1.0, 2.0, 3.0]})
        df["a"] = df["a"].astype(dtype)
        expected2 = pd.DataFrame({"a": pd.array([1, 2, 3], dtype=dtype)})
        tm.assert_frame_equal(df, expected2)

        tm.assert_frame_equal(df.astype(dtype), expected1)
        tm.assert_frame_equal(df.astype("int64").astype(dtype), expected1)

    @pytest.mark.parametrize("dtype", ["category", "Int64"])
    def test_astype_extension_dtypes_duplicate_col(self, dtype, using_nan_is_na):
        # GH#24704
        a1 = pd.Series([0, np.nan, 4], name="a")
        a2 = pd.Series([np.nan, 3, 5], name="a")
        df = pd.concat([a1, a2], axis=1)

        if dtype == "Int64" and not using_nan_is_na:
            msg = "Cannot cast NaN value to Integer dtype"
            with pytest.raises(ValueError, match=msg):
                df.astype(dtype)
            with pytest.raises(ValueError, match=msg):
                a1.astype(dtype)
            with pytest.raises(ValueError, match=msg):
                a2.astype(dtype)
            return

        result = df.astype(dtype)
        expected = pd.concat([a1.astype(dtype), a2.astype(dtype)], axis=1)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype", [{100: "float64", 200: "uint64"}, "category", "float64"]
    )
    def test_astype_column_metadata(self, dtype):
        # GH#19920
        columns = pd.Index([100, 200, 300], dtype=np.uint64, name="foo")
        df = pd.DataFrame(np.arange(15).reshape(5, 3), columns=columns)
        df = df.astype(dtype)
        tm.assert_index_equal(df.columns, columns)

    @pytest.mark.parametrize("unit", ["Y", "M", "W", "D", "h", "m"])
    def test_astype_from_object_to_datetime_unit(self, unit):
        vals = [
            ["2015-01-01", "2015-01-02", "2015-01-03"],
            ["2017-01-01", "2017-01-02", "2017-02-03"],
        ]
        df = pd.DataFrame(vals, dtype=object)
        msg = (
            rf"Unexpected value for 'dtype': 'datetime64\[{unit}\]'. "
            r"Must be 'datetime64\[s\]', 'datetime64\[ms\]', 'datetime64\[us\]', "
            r"'datetime64\[ns\]' or DatetimeTZDtype"
        )
        with pytest.raises(ValueError, match=msg):
            df.astype(f"M8[{unit}]")

    @pytest.mark.parametrize("unit", ["Y", "M", "W", "D", "h", "m"])
    def test_astype_from_object_to_timedelta_unit(self, unit):
        vals = [
            ["1 Day", "2 Days", "3 Days"],
            ["4 Days", "5 Days", "6 Days"],
        ]
        df = pd.DataFrame(vals, dtype=object)
        msg = (
            r"Cannot convert from timedelta64\[us\] to timedelta64\[.*\]. "
            "Supported resolutions are 's', 'ms', 'us', 'ns'"
        )
        with pytest.raises(ValueError, match=msg):
            # TODO: this is ValueError while for DatetimeArray it is TypeError;
            #  get these consistent
            df.astype(f"m8[{unit}]")

    @pytest.mark.parametrize("dtype", ["M8", "m8"])
    @pytest.mark.parametrize("unit", ["ns", "us", "ms", "s", "h", "m", "D"])
    def test_astype_from_datetimelike_to_object(self, dtype, unit):
        # tests astype to object dtype
        # GH#19223 / GH#12425
        dtype = f"{dtype}[{unit}]"
        arr = np.array([[1, 2, 3]], dtype=dtype)
        df = pd.DataFrame(arr)
        result = df.astype(object)
        assert (result.dtypes == object).all()

        if dtype.startswith("M8"):
            assert result.iloc[0, 0] == pd.Timestamp(1, unit=unit)
        else:
            assert result.iloc[0, 0] == pd.Timedelta(1, unit=unit)

    @pytest.mark.parametrize("dtype", ["M8", "m8"])
    @pytest.mark.parametrize("unit", ["ns", "us", "ms", "s", "h", "m", "D"])
    def test_astype_to_datetimelike_unit(self, any_real_numpy_dtype, dtype, unit):
        # tests all units from numeric origination
        # GH#19223 / GH#12425
        dtype = f"{dtype}[{unit}]"
        arr = np.array([[1, 2, 3]], dtype=any_real_numpy_dtype)
        df = pd.DataFrame(arr)
        result = df.astype(dtype)
        expected = pd.DataFrame(arr.astype(dtype))

        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("unit", ["ns", "us", "ms", "s", "h", "m", "D"])
    def test_astype_to_datetime_unit(self, unit):
        # tests all units from datetime origination
        # GH#19223
        dtype = f"M8[{unit}]"
        arr = np.array([[1, 2, 3]], dtype=dtype)
        df = pd.DataFrame(arr)
        ser = df.iloc[:, 0]
        idx = pd.Index(ser)
        dta = ser._values

        if unit in ["ns", "us", "ms", "s"]:
            # GH#48928
            result = df.astype(dtype)
        else:
            # we use the nearest supported dtype (i.e. M8[s])
            msg = rf"Cannot cast DatetimeArray to dtype datetime64\[{unit}\]"
            with pytest.raises(TypeError, match=msg):
                df.astype(dtype)

            with pytest.raises(TypeError, match=msg):
                ser.astype(dtype)

            with pytest.raises(TypeError, match=msg.replace("Array", "Index")):
                idx.astype(dtype)

            with pytest.raises(TypeError, match=msg):
                dta.astype(dtype)

            return

        exp_df = pd.DataFrame(arr.astype(dtype))
        assert (exp_df.dtypes == dtype).all()
        tm.assert_frame_equal(result, exp_df)

        res_ser = ser.astype(dtype)
        exp_ser = exp_df.iloc[:, 0]
        assert exp_ser.dtype == dtype
        tm.assert_series_equal(res_ser, exp_ser)

        exp_dta = exp_ser._values

        res_index = idx.astype(dtype)
        exp_index = pd.Index(exp_ser)
        assert exp_index.dtype == dtype
        tm.assert_index_equal(res_index, exp_index)

        res_dta = dta.astype(dtype)
        assert exp_dta.dtype == dtype
        tm.assert_extension_array_equal(res_dta, exp_dta)

    def test_astype_to_timedelta_unit_ns(self):
        # preserver the timedelta conversion
        # GH#19223
        dtype = "m8[ns]"
        arr = np.array([[1, 2, 3]], dtype=dtype)
        df = pd.DataFrame(arr)
        result = df.astype(dtype)
        expected = pd.DataFrame(arr.astype(dtype))

        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("unit", ["us", "ms", "s", "h", "m", "D"])
    def test_astype_to_timedelta_unit(self, unit):
        # coerce to float
        # GH#19223 until 2.0 used to coerce to float
        dtype = f"m8[{unit}]"
        arr = np.array([[1, 2, 3]], dtype=dtype)
        df = pd.DataFrame(arr)
        ser = df.iloc[:, 0]
        tdi = pd.Index(ser)
        tda = tdi._values

        if unit in ["us", "ms", "s"]:
            assert (df.dtypes == dtype).all()
            result = df.astype(dtype)
        else:
            # We get the nearest supported unit, i.e. "s"
            assert (df.dtypes == "m8[s]").all()

            msg = (
                rf"Cannot convert from timedelta64\[s\] to timedelta64\[{unit}\]. "
                "Supported resolutions are 's', 'ms', 'us', 'ns'"
            )
            with pytest.raises(ValueError, match=msg):
                df.astype(dtype)
            with pytest.raises(ValueError, match=msg):
                ser.astype(dtype)
            with pytest.raises(ValueError, match=msg):
                tdi.astype(dtype)
            with pytest.raises(ValueError, match=msg):
                tda.astype(dtype)

            return

        result = df.astype(dtype)
        # The conversion is a no-op, so we just get a copy
        expected = df
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("unit", ["ns", "us", "ms", "s", "h", "m", "D"])
    def test_astype_to_incorrect_datetimelike(self, unit):
        # trying to astype an m to an M, or vice-versa
        # GH#19224
        dtype = f"M8[{unit}]"
        other = f"m8[{unit}]"

        df = pd.DataFrame(np.array([[1, 2, 3]], dtype=dtype))
        msg = rf"Cannot cast DatetimeArray to dtype timedelta64\[{unit}\]"
        with pytest.raises(TypeError, match=msg):
            df.astype(other)

        msg = rf"Cannot cast TimedeltaArray to dtype datetime64\[{unit}\]"
        df = pd.DataFrame(np.array([[1, 2, 3]], dtype=other))
        with pytest.raises(TypeError, match=msg):
            df.astype(dtype)

    def test_astype_arg_for_errors(self):
        # GH#14878

        df = pd.DataFrame([1, 2, 3])

        msg = (
            "Expected value of kwarg 'errors' to be one of "
            "['raise', 'ignore']. Supplied value is 'True'"
        )
        with pytest.raises(ValueError, match=re.escape(msg)):
            df.astype(np.float64, errors=True)

        df.astype(np.int8, errors="ignore")

    def test_astype_invalid_conversion(self):
        # GH#47571
        df = pd.DataFrame({"a": [1, 2, "text"], "b": [1, 2, 3]})

        msg = (
            "invalid literal for int() with base 10: 'text': "
            "Error while type casting for column 'a'"
        )

        with pytest.raises(ValueError, match=re.escape(msg)):
            df.astype({"a": int})

    def test_astype_arg_for_errors_dictlist(self):
        # GH#25905
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
        expected["c"] = expected["c"].astype("object")
        type_dict = {"a": "float64", "b": "float64", "c": "object"}

        result = df.astype(dtype=type_dict, errors="ignore")

        tm.assert_frame_equal(result, expected)

    def test_astype_dt64tz(self, timezone_frame):
        # astype
        expected = np.array(
            [
                [
                    pd.Timestamp("2013-01-01 00:00:00"),
                    pd.Timestamp("2013-01-02 00:00:00"),
                    pd.Timestamp("2013-01-03 00:00:00"),
                ],
                [
                    pd.Timestamp("2013-01-01 00:00:00-0500", tz="US/Eastern"),
                    pd.NaT,
                    pd.Timestamp("2013-01-03 00:00:00-0500", tz="US/Eastern"),
                ],
                [
                    pd.Timestamp("2013-01-01 00:00:00+0100", tz="CET"),
                    pd.NaT,
                    pd.Timestamp("2013-01-03 00:00:00+0100", tz="CET"),
                ],
            ],
            dtype=object,
        ).T
        expected = pd.DataFrame(
            expected,
            index=timezone_frame.index,
            columns=timezone_frame.columns,
            dtype=object,
        )
        result = timezone_frame.astype(object)
        tm.assert_frame_equal(result, expected)

        msg = "Cannot use .astype to convert from timezone-aware dtype to timezone-"
        with pytest.raises(TypeError, match=msg):
            # dt64tz->dt64 deprecated
            timezone_frame.astype("datetime64[ns]")

    def test_astype_dt64tz_to_str(self, timezone_frame, using_infer_string):
        # str formatting
        result = timezone_frame.astype(str)
        na_value = np.nan if using_infer_string else "NaT"
        expected = pd.DataFrame(
            [
                [
                    "2013-01-01",
                    "2013-01-01 00:00:00-05:00",
                    "2013-01-01 00:00:00+01:00",
                ],
                ["2013-01-02", na_value, na_value],
                [
                    "2013-01-03",
                    "2013-01-03 00:00:00-05:00",
                    "2013-01-03 00:00:00+01:00",
                ],
            ],
            columns=timezone_frame.columns,
            dtype="str",
        )
        tm.assert_frame_equal(result, expected)

        with pd.option_context("display.max_columns", 20):
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

    def test_astype_empty_dtype_dict(self):
        # issue mentioned further down in the following issue's thread
        # https://github.com/pandas-dev/pandas/issues/33113
        df = pd.DataFrame()
        result = df.astype({})
        tm.assert_frame_equal(result, df)
        assert result is not df

    @pytest.mark.parametrize(
        "data, dtype",
        [
            (["x", "y", "z"], "string[python]"),
            pytest.param(
                ["x", "y", "z"],
                "string[pyarrow]",
                marks=td.skip_if_no("pyarrow"),
            ),
            (["x", "y", "z"], "category"),
            (3 * [pd.Timestamp("2020-01-01", tz="UTC")], None),
            (3 * [pd.Interval(0, 1)], None),
        ],
    )
    @pytest.mark.parametrize("errors", ["raise", "ignore"])
    def test_astype_ignores_errors_for_extension_dtypes(self, data, dtype, errors):
        # https://github.com/pandas-dev/pandas/issues/35471
        df = pd.DataFrame(pd.Series(data, dtype=dtype))
        if errors == "ignore":
            expected = df
            result = df.astype(float, errors=errors)
            tm.assert_frame_equal(result, expected)
        else:
            msg = "|".join(["Cannot cast", "could not convert"])
            with pytest.raises((ValueError, TypeError), match=msg):
                df.astype(float, errors=errors)

    def test_astype_tz_conversion(self):
        # GH 35973, GH#58998
        msg = "'d' is deprecated and will be removed in a future version."
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            val = {
                "tz": pd.date_range(
                    "2020-08-30", freq="d", periods=2, tz="Europe/London", unit="ns"
                )
            }
        df = pd.DataFrame(val)
        result = df.astype({"tz": "datetime64[ns, Europe/Berlin]"})

        expected = df
        expected["tz"] = expected["tz"].dt.tz_convert("Europe/Berlin")
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("tz", ["UTC", "Europe/Berlin"])
    def test_astype_tz_object_conversion(self, tz):
        # GH 35973
        val = {
            "tz": pd.date_range(
                "2020-08-30", freq="D", periods=2, tz="Europe/London", unit="ns"
            )
        }
        expected = pd.DataFrame(val)

        # convert expected to object dtype from other tz str (independently tested)
        result = expected.astype({"tz": f"datetime64[ns, {tz}]"})
        result = result.astype({"tz": "object"})

        # do real test: object dtype to a specified tz, different from construction tz.
        result = result.astype({"tz": "datetime64[ns, Europe/London]"})
        tm.assert_frame_equal(result, expected)

    def test_astype_dt64_to_string(self, frame_or_series, tz_naive_fixture):
        # GH#41409
        tz = tz_naive_fixture

        dti = pd.date_range("2016-01-01", periods=3, tz=tz)
        dta = dti._data
        dta[0] = pd.NaT

        obj = frame_or_series(dta)
        result = obj.astype("string")

        # Check that Series/DataFrame.astype matches DatetimeArray.astype
        expected = frame_or_series(dta.astype("string"))
        tm.assert_equal(result, expected)

        item = result.iloc[0]
        if frame_or_series is pd.DataFrame:
            item = item.iloc[0]
        assert item is pd.NA

        # For non-NA values, we should match what we get for non-EA str
        alt = obj.astype(str)
        assert np.all(alt.iloc[1:] == result.iloc[1:])

    def test_astype_td64_to_string(self, frame_or_series):
        # GH#41409
        tdi = pd.timedelta_range("1 Day", periods=3)
        obj = frame_or_series(tdi)

        expected = frame_or_series(["1 days", "2 days", "3 days"], dtype="string")
        result = obj.astype("string")
        tm.assert_equal(result, expected)

    def test_astype_bytes(self):
        # GH#39474
        result = pd.DataFrame(["foo", "bar", "baz"]).astype(bytes)
        assert result.dtypes[0] == np.dtype("S3")

    @pytest.mark.parametrize(
        "index_slice",
        [
            np.s_[:2, :2],
            np.s_[:1, :2],
            np.s_[:2, :1],
            np.s_[::2, ::2],
            np.s_[::1, ::2],
            np.s_[::2, ::1],
        ],
    )
    def test_astype_noncontiguous(self, index_slice):
        # GH#42396
        data = np.arange(16).reshape(4, 4)
        df = pd.DataFrame(data)

        result = df.iloc[index_slice].astype("int16")
        expected = df.iloc[index_slice]
        tm.assert_frame_equal(result, expected, check_dtype=False)

    def test_astype_retain_attrs(self, any_numpy_dtype):
        # GH#44414
        df = pd.DataFrame({"a": [0, 1, 2], "b": [3, 4, 5]})
        df.attrs["Location"] = "Michigan"

        result = df.astype({"a": any_numpy_dtype}).attrs
        expected = df.attrs

        tm.assert_dict_equal(expected, result)


class TestAstypeCategorical:
    def test_astype_from_categorical3(self):
        df = pd.DataFrame({"cats": [1, 2, 3, 4, 5, 6], "vals": [1, 2, 3, 4, 5, 6]})
        cats = pd.Categorical([1, 2, 3, 4, 5, 6])
        exp_df = pd.DataFrame({"cats": cats, "vals": [1, 2, 3, 4, 5, 6]})
        df["cats"] = df["cats"].astype("category")
        tm.assert_frame_equal(exp_df, df)

    def test_astype_from_categorical4(self):
        df = pd.DataFrame(
            {"cats": ["a", "b", "b", "a", "a", "d"], "vals": [1, 2, 3, 4, 5, 6]}
        )
        cats = pd.Categorical(["a", "b", "b", "a", "a", "d"])
        exp_df = pd.DataFrame({"cats": cats, "vals": [1, 2, 3, 4, 5, 6]})
        df["cats"] = df["cats"].astype("category")
        tm.assert_frame_equal(exp_df, df)

    def test_categorical_astype_to_int(self, any_int_dtype):
        # GH#39402

        df = pd.DataFrame(data={"col1": pd.array([2.0, 1.0, 3.0])})
        df.col1 = df.col1.astype("category")
        df.col1 = df.col1.astype(any_int_dtype)
        expected = pd.DataFrame({"col1": pd.array([2, 1, 3], dtype=any_int_dtype)})
        tm.assert_frame_equal(df, expected)

    def test_astype_categorical_to_string_missing(self):
        # https://github.com/pandas-dev/pandas/issues/41797
        df = pd.DataFrame(["a", "b", np.nan])
        expected = df.astype(str)
        cat = df.astype("category")
        result = cat.astype(str)
        tm.assert_frame_equal(result, expected)


class IntegerArrayNoCopy(pd.core.arrays.IntegerArray):
    # GH 42501

    def copy(self):
        raise NotImplementedError


class Int16DtypeNoCopy(pd.Int16Dtype):
    # GH 42501

    def construct_array_type(self):
        return IntegerArrayNoCopy


def test_frame_astype_no_copy():
    # GH 42501
    df = pd.DataFrame({"a": [1, 4, None, 5], "b": [6, 7, 8, 9]}, dtype=object)
    result = df.astype({"a": Int16DtypeNoCopy()})

    assert result.a.dtype == pd.Int16Dtype()
    assert np.shares_memory(df.b.values, result.b.values)


@pytest.mark.parametrize("dtype", ["int64", "Int64"])
def test_astype_copies(dtype):
    # GH#50984
    pytest.importorskip("pyarrow")
    df = pd.DataFrame({"a": [1, 2, 3]}, dtype=dtype)
    result = df.astype("int64[pyarrow]")
    df.iloc[0, 0] = 100
    expected = pd.DataFrame({"a": [1, 2, 3]}, dtype="int64[pyarrow]")
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize("val", [None, 1, 1.5, np.nan, pd.NaT])
def test_astype_to_string_not_modifying_input(string_storage, val):
    # GH#51073
    df = pd.DataFrame({"a": ["a", "b", val]})
    expected = df.copy()
    with pd.option_context("mode.string_storage", string_storage):
        df.astype("string")
    tm.assert_frame_equal(df, expected)


@pytest.mark.parametrize("val", [None, 1, 1.5, np.nan, pd.NaT])
def test_astype_to_string_dtype_not_modifying_input(any_string_dtype, val):
    # GH#51073 - variant of the above test with explicit dtype instances
    df = pd.DataFrame({"a": ["a", "b", val]})
    expected = df.copy()
    df.astype(any_string_dtype)
    tm.assert_frame_equal(df, expected)
