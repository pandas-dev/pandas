from datetime import timedelta

import numpy as np
import pytest

from pandas import Categorical, DataFrame, NaT, Period, Series, Timedelta, Timestamp
import pandas._testing as tm


class TestSeriesFillNA:
    def test_fillna_pytimedelta(self):
        # GH#8209
        ser = Series([np.nan, Timedelta("1 days")], index=["A", "B"])

        result = ser.fillna(timedelta(1))
        expected = Series(Timedelta("1 days"), index=["A", "B"])
        tm.assert_series_equal(result, expected)

    def test_fillna_period(self):
        # GH#13737
        ser = Series([Period("2011-01", freq="M"), Period("NaT", freq="M")])

        res = ser.fillna(Period("2012-01", freq="M"))
        exp = Series([Period("2011-01", freq="M"), Period("2012-01", freq="M")])
        tm.assert_series_equal(res, exp)
        assert res.dtype == "Period[M]"

    def test_fillna_dt64_timestamp(self):
        ser = Series(
            [
                Timestamp("20130101"),
                Timestamp("20130101"),
                Timestamp("20130102"),
                Timestamp("20130103 9:01:01"),
            ]
        )
        ser[2] = np.nan

        # reg fillna
        result = ser.fillna(Timestamp("20130104"))
        expected = Series(
            [
                Timestamp("20130101"),
                Timestamp("20130101"),
                Timestamp("20130104"),
                Timestamp("20130103 9:01:01"),
            ]
        )
        tm.assert_series_equal(result, expected)

        result = ser.fillna(NaT)
        expected = ser
        tm.assert_series_equal(result, expected)

    def test_fillna_dt64_non_nao(self):
        # GH#27419
        ser = Series([Timestamp("2010-01-01"), NaT, Timestamp("2000-01-01")])
        val = np.datetime64("1975-04-05", "ms")

        result = ser.fillna(val)
        expected = Series(
            [Timestamp("2010-01-01"), Timestamp("1975-04-05"), Timestamp("2000-01-01")]
        )
        tm.assert_series_equal(result, expected)

    def test_fillna_numeric_inplace(self):
        x = Series([np.nan, 1.0, np.nan, 3.0, np.nan], ["z", "a", "b", "c", "d"])
        y = x.copy()

        return_value = y.fillna(value=0, inplace=True)
        assert return_value is None

        expected = x.fillna(value=0)
        tm.assert_series_equal(y, expected)

    # ---------------------------------------------------------------
    # CategoricalDtype

    @pytest.mark.parametrize(
        "fill_value, expected_output",
        [
            ("a", ["a", "a", "b", "a", "a"]),
            ({1: "a", 3: "b", 4: "b"}, ["a", "a", "b", "b", "b"]),
            ({1: "a"}, ["a", "a", "b", np.nan, np.nan]),
            ({1: "a", 3: "b"}, ["a", "a", "b", "b", np.nan]),
            (Series("a"), ["a", np.nan, "b", np.nan, np.nan]),
            (Series("a", index=[1]), ["a", "a", "b", np.nan, np.nan]),
            (Series({1: "a", 3: "b"}), ["a", "a", "b", "b", np.nan]),
            (Series(["a", "b"], index=[3, 4]), ["a", np.nan, "b", "a", "b"]),
        ],
    )
    def test_fillna_categorical(self, fill_value, expected_output):
        # GH#17033
        # Test fillna for a Categorical series
        data = ["a", np.nan, "b", np.nan, np.nan]
        ser = Series(Categorical(data, categories=["a", "b"]))
        exp = Series(Categorical(expected_output, categories=["a", "b"]))
        result = ser.fillna(fill_value)
        tm.assert_series_equal(result, exp)

    @pytest.mark.parametrize(
        "fill_value, expected_output",
        [
            (Series(["a", "b", "c", "d", "e"]), ["a", "b", "b", "d", "e"]),
            (Series(["b", "d", "a", "d", "a"]), ["a", "d", "b", "d", "a"]),
            (
                Series(
                    Categorical(
                        ["b", "d", "a", "d", "a"], categories=["b", "c", "d", "e", "a"]
                    )
                ),
                ["a", "d", "b", "d", "a"],
            ),
        ],
    )
    def test_fillna_categorical_with_new_categories(self, fill_value, expected_output):
        # GH#26215
        data = ["a", np.nan, "b", np.nan, np.nan]
        ser = Series(Categorical(data, categories=["a", "b", "c", "d", "e"]))
        exp = Series(Categorical(expected_output, categories=["a", "b", "c", "d", "e"]))
        result = ser.fillna(fill_value)
        tm.assert_series_equal(result, exp)

    def test_fillna_categorical_raises(self):
        data = ["a", np.nan, "b", np.nan, np.nan]
        ser = Series(Categorical(data, categories=["a", "b"]))

        with pytest.raises(ValueError, match="fill value must be in categories"):
            ser.fillna("d")

        with pytest.raises(ValueError, match="fill value must be in categories"):
            ser.fillna(Series("d"))

        with pytest.raises(ValueError, match="fill value must be in categories"):
            ser.fillna({1: "d", 3: "a"})

        msg = '"value" parameter must be a scalar or dict, but you passed a "list"'
        with pytest.raises(TypeError, match=msg):
            ser.fillna(["a", "b"])

        msg = '"value" parameter must be a scalar or dict, but you passed a "tuple"'
        with pytest.raises(TypeError, match=msg):
            ser.fillna(("a", "b"))

        msg = (
            '"value" parameter must be a scalar, dict '
            'or Series, but you passed a "DataFrame"'
        )
        with pytest.raises(TypeError, match=msg):
            ser.fillna(DataFrame({1: ["a"], 3: ["b"]}))

    # ---------------------------------------------------------------
    # Invalid Usages

    def test_fillna_listlike_invalid(self):
        ser = Series(np.random.randint(-100, 100, 50))
        msg = '"value" parameter must be a scalar or dict, but you passed a "list"'
        with pytest.raises(TypeError, match=msg):
            ser.fillna([1, 2])

        msg = '"value" parameter must be a scalar or dict, but you passed a "tuple"'
        with pytest.raises(TypeError, match=msg):
            ser.fillna((1, 2))

    def test_fillna_method_and_limit_invalid(self):

        # related GH#9217, make sure limit is an int and greater than 0
        ser = Series([1, 2, 3, None])
        msg = (
            r"Cannot specify both 'value' and 'method'\.|"
            r"Limit must be greater than 0|"
            "Limit must be an integer"
        )
        for limit in [-1, 0, 1.0, 2.0]:
            for method in ["backfill", "bfill", "pad", "ffill", None]:
                with pytest.raises(ValueError, match=msg):
                    ser.fillna(1, limit=limit, method=method)
