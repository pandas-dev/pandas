import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestSeriesValueCounts:
    def test_value_counts_datetime(self, unit):
        # most dtypes are tested in tests/base
        values = [
            pd.Timestamp("2011-01-01 09:00"),
            pd.Timestamp("2011-01-01 10:00"),
            pd.Timestamp("2011-01-01 11:00"),
            pd.Timestamp("2011-01-01 09:00"),
            pd.Timestamp("2011-01-01 09:00"),
            pd.Timestamp("2011-01-01 11:00"),
        ]

        exp_idx = pd.DatetimeIndex(
            ["2011-01-01 09:00", "2011-01-01 11:00", "2011-01-01 10:00"],
            name="xxx",
        ).as_unit(unit)
        exp = pd.Series([3, 2, 1], index=exp_idx, name="count")

        ser = pd.Series(values, name="xxx").dt.as_unit(unit)
        tm.assert_series_equal(ser.value_counts(), exp)
        # check DatetimeIndex outputs the same result
        idx = pd.DatetimeIndex(values, name="xxx").as_unit(unit)
        tm.assert_series_equal(idx.value_counts(), exp)

        # normalize
        exp = pd.Series(np.array([3.0, 2.0, 1]) / 6.0, index=exp_idx, name="proportion")
        tm.assert_series_equal(ser.value_counts(normalize=True), exp)
        tm.assert_series_equal(idx.value_counts(normalize=True), exp)

    def test_value_counts_datetime_tz(self, unit):
        values = [
            pd.Timestamp("2011-01-01 09:00", tz="US/Eastern"),
            pd.Timestamp("2011-01-01 10:00", tz="US/Eastern"),
            pd.Timestamp("2011-01-01 11:00", tz="US/Eastern"),
            pd.Timestamp("2011-01-01 09:00", tz="US/Eastern"),
            pd.Timestamp("2011-01-01 09:00", tz="US/Eastern"),
            pd.Timestamp("2011-01-01 11:00", tz="US/Eastern"),
        ]

        exp_idx = pd.DatetimeIndex(
            ["2011-01-01 09:00", "2011-01-01 11:00", "2011-01-01 10:00"],
            tz="US/Eastern",
            name="xxx",
        ).as_unit(unit)
        exp = pd.Series([3, 2, 1], index=exp_idx, name="count")

        ser = pd.Series(values, name="xxx").dt.as_unit(unit)
        tm.assert_series_equal(ser.value_counts(), exp)
        idx = pd.DatetimeIndex(values, name="xxx").as_unit(unit)
        tm.assert_series_equal(idx.value_counts(), exp)

        exp = pd.Series(np.array([3.0, 2.0, 1]) / 6.0, index=exp_idx, name="proportion")
        tm.assert_series_equal(ser.value_counts(normalize=True), exp)
        tm.assert_series_equal(idx.value_counts(normalize=True), exp)

    def test_value_counts_period(self):
        values = [
            pd.Period("2011-01", freq="M"),
            pd.Period("2011-02", freq="M"),
            pd.Period("2011-03", freq="M"),
            pd.Period("2011-01", freq="M"),
            pd.Period("2011-01", freq="M"),
            pd.Period("2011-03", freq="M"),
        ]

        exp_idx = pd.PeriodIndex(
            ["2011-01", "2011-03", "2011-02"], freq="M", name="xxx"
        )
        exp = pd.Series([3, 2, 1], index=exp_idx, name="count")

        ser = pd.Series(values, name="xxx")
        tm.assert_series_equal(ser.value_counts(), exp)
        # check DatetimeIndex outputs the same result
        idx = pd.PeriodIndex(values, name="xxx")
        tm.assert_series_equal(idx.value_counts(), exp)

        # normalize
        exp = pd.Series(np.array([3.0, 2.0, 1]) / 6.0, index=exp_idx, name="proportion")
        tm.assert_series_equal(ser.value_counts(normalize=True), exp)
        tm.assert_series_equal(idx.value_counts(normalize=True), exp)

    def test_value_counts_categorical_ordered(self):
        # most dtypes are tested in tests/base
        values = pd.Categorical([1, 2, 3, 1, 1, 3], ordered=True)

        exp_idx = pd.CategoricalIndex(
            [1, 3, 2], categories=[1, 2, 3], ordered=True, name="xxx"
        )
        exp = pd.Series([3, 2, 1], index=exp_idx, name="count")

        ser = pd.Series(values, name="xxx")
        tm.assert_series_equal(ser.value_counts(), exp)
        # check CategoricalIndex outputs the same result
        idx = pd.CategoricalIndex(values, name="xxx")
        tm.assert_series_equal(idx.value_counts(), exp)

        # normalize
        exp = pd.Series(np.array([3.0, 2.0, 1]) / 6.0, index=exp_idx, name="proportion")
        tm.assert_series_equal(ser.value_counts(normalize=True), exp)
        tm.assert_series_equal(idx.value_counts(normalize=True), exp)

    def test_value_counts_categorical_not_ordered(self):
        values = pd.Categorical([1, 2, 3, 1, 1, 3], ordered=False)

        exp_idx = pd.CategoricalIndex(
            [1, 3, 2], categories=[1, 2, 3], ordered=False, name="xxx"
        )
        exp = pd.Series([3, 2, 1], index=exp_idx, name="count")

        ser = pd.Series(values, name="xxx")
        tm.assert_series_equal(ser.value_counts(), exp)
        # check CategoricalIndex outputs the same result
        idx = pd.CategoricalIndex(values, name="xxx")
        tm.assert_series_equal(idx.value_counts(), exp)

        # normalize
        exp = pd.Series(np.array([3.0, 2.0, 1]) / 6.0, index=exp_idx, name="proportion")
        tm.assert_series_equal(ser.value_counts(normalize=True), exp)
        tm.assert_series_equal(idx.value_counts(normalize=True), exp)

    def test_value_counts_categorical(self):
        # GH#12835
        cats = pd.Categorical(list("abcccb"), categories=list("cabd"))
        ser = pd.Series(cats, name="xxx")
        res = ser.value_counts(sort=False)

        exp_index = pd.CategoricalIndex(
            list("cabd"), categories=cats.categories, name="xxx"
        )
        exp = pd.Series([3, 1, 2, 0], name="count", index=exp_index)
        tm.assert_series_equal(res, exp)

        res = ser.value_counts(sort=True)

        exp_index = pd.CategoricalIndex(
            list("cbad"), categories=cats.categories, name="xxx"
        )
        exp = pd.Series([3, 2, 1, 0], name="count", index=exp_index)
        tm.assert_series_equal(res, exp)

        # check object dtype handles the Series.name as the same
        # (tested in tests/base)
        ser = pd.Series(["a", "b", "c", "c", "c", "b"], name="xxx")
        res = ser.value_counts()
        exp = pd.Series(
            [3, 2, 1], name="count", index=pd.Index(["c", "b", "a"], name="xxx")
        )
        tm.assert_series_equal(res, exp)

    def test_value_counts_categorical_with_nan(self):
        # see GH#9443

        # sanity check
        ser = pd.Series(["a", "b", "a"], dtype="category")
        exp = pd.Series([2, 1], index=pd.CategoricalIndex(["a", "b"]), name="count")

        res = ser.value_counts(dropna=True)
        tm.assert_series_equal(res, exp)

        res = ser.value_counts(dropna=True)
        tm.assert_series_equal(res, exp)

        # same Series via two different constructions --> same behaviour
        series = [
            pd.Series(["a", "b", None, "a", None, None], dtype="category"),
            pd.Series(
                pd.Categorical(["a", "b", None, "a", None, None], categories=["a", "b"])
            ),
        ]

        for ser in series:
            # None is a NaN value, so we exclude its count here
            exp = pd.Series([2, 1], index=pd.CategoricalIndex(["a", "b"]), name="count")
            res = ser.value_counts(dropna=True)
            tm.assert_series_equal(res, exp)

            # we don't exclude the count of None and sort by counts
            exp = pd.Series(
                [3, 2, 1], index=pd.CategoricalIndex([np.nan, "a", "b"]), name="count"
            )
            res = ser.value_counts(dropna=False)
            tm.assert_series_equal(res, exp)

            # When we aren't sorting by counts, and np.nan isn't a
            # category, it should be last.
            exp = pd.Series(
                [2, 1, 3], index=pd.CategoricalIndex(["a", "b", np.nan]), name="count"
            )
            res = ser.value_counts(dropna=False, sort=False)
            tm.assert_series_equal(res, exp)

    @pytest.mark.parametrize(
        "ser, dropna, exp",
        [
            (
                pd.Series([False, True, True, pd.NA]),
                False,
                pd.Series([2, 1, 1], index=[True, False, pd.NA], name="count"),
            ),
            (
                pd.Series([False, True, True, pd.NA]),
                True,
                pd.Series(
                    [2, 1], index=pd.Index([True, False], dtype=object), name="count"
                ),
            ),
            (
                pd.Series(range(3), index=[True, False, np.nan]).index,
                False,
                pd.Series([1, 1, 1], index=[True, False, np.nan], name="count"),
            ),
        ],
    )
    def test_value_counts_bool_with_nan(self, ser, dropna, exp):
        # GH32146
        out = ser.value_counts(dropna=dropna)
        tm.assert_series_equal(out, exp)

    @pytest.mark.parametrize("dtype", [np.complex128, np.complex64])
    def test_value_counts_complex_numbers(self, dtype):
        # GH 17927
        input_array = np.array([1 + 1j, 1 + 1j, 1, 3j, 3j, 3j], dtype=dtype)
        result = pd.Series(input_array).value_counts()
        expected = pd.Series(
            [3, 2, 1], index=pd.Index([3j, 1 + 1j, 1], dtype=dtype), name="count"
        )
        tm.assert_series_equal(result, expected)

    def test_value_counts_masked(self):
        # GH#54984
        dtype = "Int64"
        ser = pd.Series([1, 2, None, 2, None, 3], dtype=dtype)
        result = ser.value_counts(dropna=False)
        expected = pd.Series(
            [2, 2, 1, 1],
            index=pd.Index([2, None, 1, 3], dtype=dtype),
            dtype=dtype,
            name="count",
        )
        tm.assert_series_equal(result, expected)

        result = ser.value_counts(dropna=True)
        expected = pd.Series(
            [2, 1, 1], index=pd.Index([2, 1, 3], dtype=dtype), dtype=dtype, name="count"
        )
        tm.assert_series_equal(result, expected)
