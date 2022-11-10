import numpy as np
import pytest

import pandas as pd
from pandas import (
    Categorical,
    CategoricalIndex,
    Series,
)
import pandas._testing as tm

VALUE_COUNTS_NAME_MSG = (
    r"In pandas 2.0.0, the name of the resulting Series will be 'count' "
    r"\(or 'proportion' if `normalize=True`\)"
)


class TestSeriesValueCounts:
    def test_value_counts_datetime(self):
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
            ["2011-01-01 09:00", "2011-01-01 11:00", "2011-01-01 10:00"]
        )
        exp = Series([3, 2, 1], index=exp_idx, name="xxx")

        ser = Series(values, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            tm.assert_series_equal(ser.value_counts(), exp)
        # check DatetimeIndex outputs the same result
        idx = pd.DatetimeIndex(values, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            tm.assert_series_equal(idx.value_counts(), exp)

        # normalize
        exp = Series(np.array([3.0, 2.0, 1]) / 6.0, index=exp_idx, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            tm.assert_series_equal(ser.value_counts(normalize=True), exp)
            tm.assert_series_equal(idx.value_counts(normalize=True), exp)

    def test_value_counts_datetime_tz(self):
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
        )
        exp = Series([3, 2, 1], index=exp_idx, name="xxx")

        ser = Series(values, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            tm.assert_series_equal(ser.value_counts(), exp)
        idx = pd.DatetimeIndex(values, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            tm.assert_series_equal(idx.value_counts(), exp)

        exp = Series(np.array([3.0, 2.0, 1]) / 6.0, index=exp_idx, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
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

        exp_idx = pd.PeriodIndex(["2011-01", "2011-03", "2011-02"], freq="M")
        exp = Series([3, 2, 1], index=exp_idx, name="xxx")

        ser = Series(values, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            tm.assert_series_equal(ser.value_counts(), exp)
        # check DatetimeIndex outputs the same result
        idx = pd.PeriodIndex(values, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            tm.assert_series_equal(idx.value_counts(), exp)

        # normalize
        exp = Series(np.array([3.0, 2.0, 1]) / 6.0, index=exp_idx, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            tm.assert_series_equal(ser.value_counts(normalize=True), exp)
            tm.assert_series_equal(idx.value_counts(normalize=True), exp)

    def test_value_counts_categorical_ordered(self):
        # most dtypes are tested in tests/base
        values = Categorical([1, 2, 3, 1, 1, 3], ordered=True)

        exp_idx = CategoricalIndex([1, 3, 2], categories=[1, 2, 3], ordered=True)
        exp = Series([3, 2, 1], index=exp_idx, name="xxx")

        ser = Series(values, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            tm.assert_series_equal(ser.value_counts(), exp)
        # check CategoricalIndex outputs the same result
        idx = CategoricalIndex(values, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            tm.assert_series_equal(idx.value_counts(), exp)

        # normalize
        exp = Series(np.array([3.0, 2.0, 1]) / 6.0, index=exp_idx, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            tm.assert_series_equal(ser.value_counts(normalize=True), exp)
            tm.assert_series_equal(idx.value_counts(normalize=True), exp)

    def test_value_counts_categorical_not_ordered(self):
        values = Categorical([1, 2, 3, 1, 1, 3], ordered=False)

        exp_idx = CategoricalIndex([1, 3, 2], categories=[1, 2, 3], ordered=False)
        exp = Series([3, 2, 1], index=exp_idx, name="xxx")

        ser = Series(values, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            tm.assert_series_equal(ser.value_counts(), exp)
        # check CategoricalIndex outputs the same result
        idx = CategoricalIndex(values, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            tm.assert_series_equal(idx.value_counts(), exp)

        # normalize
        exp = Series(np.array([3.0, 2.0, 1]) / 6.0, index=exp_idx, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            tm.assert_series_equal(ser.value_counts(normalize=True), exp)
            tm.assert_series_equal(idx.value_counts(normalize=True), exp)

    def test_value_counts_categorical(self):
        # GH#12835
        cats = Categorical(list("abcccb"), categories=list("cabd"))
        ser = Series(cats, name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            res = ser.value_counts(sort=False)

        exp_index = CategoricalIndex(list("cabd"), categories=cats.categories)
        exp = Series([3, 1, 2, 0], name="xxx", index=exp_index)
        tm.assert_series_equal(res, exp)

        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            res = ser.value_counts(sort=True)

        exp_index = CategoricalIndex(list("cbad"), categories=cats.categories)
        exp = Series([3, 2, 1, 0], name="xxx", index=exp_index)
        tm.assert_series_equal(res, exp)

        # check object dtype handles the Series.name as the same
        # (tested in tests/base)
        ser = Series(["a", "b", "c", "c", "c", "b"], name="xxx")
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            res = ser.value_counts()
        exp = Series([3, 2, 1], name="xxx", index=["c", "b", "a"])
        tm.assert_series_equal(res, exp)

    def test_value_counts_categorical_with_nan(self):
        # see GH#9443

        # sanity check
        ser = Series(["a", "b", "a"], dtype="category")
        exp = Series([2, 1], index=CategoricalIndex(["a", "b"]))

        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            res = ser.value_counts(dropna=True)
        tm.assert_series_equal(res, exp)

        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            res = ser.value_counts(dropna=True)
        tm.assert_series_equal(res, exp)

        # same Series via two different constructions --> same behaviour
        series = [
            Series(["a", "b", None, "a", None, None], dtype="category"),
            Series(
                Categorical(["a", "b", None, "a", None, None], categories=["a", "b"])
            ),
        ]

        for ser in series:
            # None is a NaN value, so we exclude its count here
            exp = Series([2, 1], index=CategoricalIndex(["a", "b"]))
            with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
                res = ser.value_counts(dropna=True)
            tm.assert_series_equal(res, exp)

            # we don't exclude the count of None and sort by counts
            exp = Series([3, 2, 1], index=CategoricalIndex([np.nan, "a", "b"]))
            with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
                res = ser.value_counts(dropna=False)
            tm.assert_series_equal(res, exp)

            # When we aren't sorting by counts, and np.nan isn't a
            # category, it should be last.
            exp = Series([2, 1, 3], index=CategoricalIndex(["a", "b", np.nan]))
            with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
                res = ser.value_counts(dropna=False, sort=False)
            tm.assert_series_equal(res, exp)

    @pytest.mark.parametrize(
        "ser, dropna, exp",
        [
            (
                Series([False, True, True, pd.NA]),
                False,
                Series([2, 1, 1], index=[True, False, pd.NA]),
            ),
            (
                Series([False, True, True, pd.NA]),
                True,
                Series([2, 1], index=pd.Index([True, False], dtype=object)),
            ),
            (
                Series(range(3), index=[True, False, np.nan]).index,
                False,
                Series([1, 1, 1], index=[True, False, np.nan]),
            ),
        ],
    )
    def test_value_counts_bool_with_nan(self, ser, dropna, exp):
        # GH32146
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            out = ser.value_counts(dropna=dropna)
        tm.assert_series_equal(out, exp)

    @pytest.mark.parametrize(
        "input_array,expected",
        [
            (
                [1 + 1j, 1 + 1j, 1, 3j, 3j, 3j],
                Series([3, 2, 1], index=pd.Index([3j, 1 + 1j, 1], dtype=np.complex128)),
            ),
            (
                np.array([1 + 1j, 1 + 1j, 1, 3j, 3j, 3j], dtype=np.complex64),
                Series([3, 2, 1], index=pd.Index([3j, 1 + 1j, 1], dtype=np.complex64)),
            ),
        ],
    )
    def test_value_counts_complex_numbers(self, input_array, expected):
        # GH 17927
        with tm.assert_produces_warning(FutureWarning, match=VALUE_COUNTS_NAME_MSG):
            result = Series(input_array).value_counts()
        tm.assert_series_equal(result, expected)

    def test_value_counts_name(self):
        # https://github.com/pandas-dev/pandas/issues/49497
        ser = Series([1, 2, 3], name="foo")
        result = ser.value_counts(name="count")
        expected = Series([1, 1, 1], index=[1, 2, 3], name="count")
        tm.assert_series_equal(result, expected)
