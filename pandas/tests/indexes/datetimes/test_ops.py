from datetime import datetime

import numpy as np
import pytest

import pandas as pd
from pandas import (
    DateOffset,
    DatetimeIndex,
    Index,
    Series,
    Timestamp,
    bdate_range,
    date_range,
)
import pandas._testing as tm

from pandas.tseries.offsets import BDay, Day, Hour

START, END = datetime(2009, 1, 1), datetime(2010, 1, 1)


class TestDatetimeIndexOps:
    def test_ops_properties_basic(self, datetime_series):

        # sanity check that the behavior didn't change
        # GH#7206
        for op in ["year", "day", "second", "weekday"]:
            msg = f"'Series' object has no attribute '{op}'"
            with pytest.raises(AttributeError, match=msg):
                getattr(datetime_series, op)

        # attribute access should still work!
        s = Series(dict(year=2000, month=1, day=10))
        assert s.year == 2000
        assert s.month == 1
        assert s.day == 10
        msg = "'Series' object has no attribute 'weekday'"
        with pytest.raises(AttributeError, match=msg):
            s.weekday

    def test_repeat_range(self, tz_naive_fixture):
        tz = tz_naive_fixture
        rng = date_range("1/1/2000", "1/1/2001")

        result = rng.repeat(5)
        assert result.freq is None
        assert len(result) == 5 * len(rng)

        index = pd.date_range("2001-01-01", periods=2, freq="D", tz=tz)
        exp = pd.DatetimeIndex(
            ["2001-01-01", "2001-01-01", "2001-01-02", "2001-01-02"], tz=tz
        )
        for res in [index.repeat(2), np.repeat(index, 2)]:
            tm.assert_index_equal(res, exp)
            assert res.freq is None

        index = pd.date_range("2001-01-01", periods=2, freq="2D", tz=tz)
        exp = pd.DatetimeIndex(
            ["2001-01-01", "2001-01-01", "2001-01-03", "2001-01-03"], tz=tz
        )
        for res in [index.repeat(2), np.repeat(index, 2)]:
            tm.assert_index_equal(res, exp)
            assert res.freq is None

        index = pd.DatetimeIndex(["2001-01-01", "NaT", "2003-01-01"], tz=tz)
        exp = pd.DatetimeIndex(
            [
                "2001-01-01",
                "2001-01-01",
                "2001-01-01",
                "NaT",
                "NaT",
                "NaT",
                "2003-01-01",
                "2003-01-01",
                "2003-01-01",
            ],
            tz=tz,
        )
        for res in [index.repeat(3), np.repeat(index, 3)]:
            tm.assert_index_equal(res, exp)
            assert res.freq is None

    def test_repeat(self, tz_naive_fixture):
        tz = tz_naive_fixture
        reps = 2
        msg = "the 'axis' parameter is not supported"

        rng = pd.date_range(start="2016-01-01", periods=2, freq="30Min", tz=tz)

        expected_rng = DatetimeIndex(
            [
                Timestamp("2016-01-01 00:00:00", tz=tz, freq="30T"),
                Timestamp("2016-01-01 00:00:00", tz=tz, freq="30T"),
                Timestamp("2016-01-01 00:30:00", tz=tz, freq="30T"),
                Timestamp("2016-01-01 00:30:00", tz=tz, freq="30T"),
            ]
        )

        res = rng.repeat(reps)
        tm.assert_index_equal(res, expected_rng)
        assert res.freq is None

        tm.assert_index_equal(np.repeat(rng, reps), expected_rng)
        with pytest.raises(ValueError, match=msg):
            np.repeat(rng, reps, axis=1)

    def test_resolution(self, tz_naive_fixture):
        tz = tz_naive_fixture
        for freq, expected in zip(
            ["A", "Q", "M", "D", "H", "T", "S", "L", "U"],
            [
                "day",
                "day",
                "day",
                "day",
                "hour",
                "minute",
                "second",
                "millisecond",
                "microsecond",
            ],
        ):
            idx = pd.date_range(start="2013-04-01", periods=30, freq=freq, tz=tz)
            assert idx.resolution == expected

    def test_value_counts_unique(self, tz_naive_fixture):
        tz = tz_naive_fixture
        # GH 7735
        idx = pd.date_range("2011-01-01 09:00", freq="H", periods=10)
        # create repeated values, 'n'th element is repeated by n+1 times
        idx = DatetimeIndex(np.repeat(idx.values, range(1, len(idx) + 1)), tz=tz)

        exp_idx = pd.date_range("2011-01-01 18:00", freq="-1H", periods=10, tz=tz)
        expected = Series(range(10, 0, -1), index=exp_idx, dtype="int64")
        expected.index = expected.index._with_freq(None)

        for obj in [idx, Series(idx)]:

            tm.assert_series_equal(obj.value_counts(), expected)

        expected = pd.date_range("2011-01-01 09:00", freq="H", periods=10, tz=tz)
        expected = expected._with_freq(None)
        tm.assert_index_equal(idx.unique(), expected)

        idx = DatetimeIndex(
            [
                "2013-01-01 09:00",
                "2013-01-01 09:00",
                "2013-01-01 09:00",
                "2013-01-01 08:00",
                "2013-01-01 08:00",
                pd.NaT,
            ],
            tz=tz,
        )

        exp_idx = DatetimeIndex(["2013-01-01 09:00", "2013-01-01 08:00"], tz=tz)
        expected = Series([3, 2], index=exp_idx)

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        exp_idx = DatetimeIndex(["2013-01-01 09:00", "2013-01-01 08:00", pd.NaT], tz=tz)
        expected = Series([3, 2, 1], index=exp_idx)

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(dropna=False), expected)

        tm.assert_index_equal(idx.unique(), exp_idx)

    @pytest.mark.parametrize(
        "idx",
        [
            DatetimeIndex(
                ["2011-01-01", "2011-01-02", "2011-01-03"], freq="D", name="idx"
            ),
            DatetimeIndex(
                ["2011-01-01 09:00", "2011-01-01 10:00", "2011-01-01 11:00"],
                freq="H",
                name="tzidx",
                tz="Asia/Tokyo",
            ),
        ],
    )
    def test_order_with_freq(self, idx):
        ordered = idx.sort_values()
        tm.assert_index_equal(ordered, idx)
        assert ordered.freq == idx.freq

        ordered = idx.sort_values(ascending=False)
        expected = idx[::-1]
        tm.assert_index_equal(ordered, expected)
        assert ordered.freq == expected.freq
        assert ordered.freq.n == -1

        ordered, indexer = idx.sort_values(return_indexer=True)
        tm.assert_index_equal(ordered, idx)
        tm.assert_numpy_array_equal(indexer, np.array([0, 1, 2]), check_dtype=False)
        assert ordered.freq == idx.freq

        ordered, indexer = idx.sort_values(return_indexer=True, ascending=False)
        expected = idx[::-1]
        tm.assert_index_equal(ordered, expected)
        tm.assert_numpy_array_equal(indexer, np.array([2, 1, 0]), check_dtype=False)
        assert ordered.freq == expected.freq
        assert ordered.freq.n == -1

    @pytest.mark.parametrize(
        "index_dates,expected_dates",
        [
            (
                ["2011-01-01", "2011-01-03", "2011-01-05", "2011-01-02", "2011-01-01"],
                ["2011-01-01", "2011-01-01", "2011-01-02", "2011-01-03", "2011-01-05"],
            ),
            (
                ["2011-01-01", "2011-01-03", "2011-01-05", "2011-01-02", "2011-01-01"],
                ["2011-01-01", "2011-01-01", "2011-01-02", "2011-01-03", "2011-01-05"],
            ),
            (
                [pd.NaT, "2011-01-03", "2011-01-05", "2011-01-02", pd.NaT],
                [pd.NaT, pd.NaT, "2011-01-02", "2011-01-03", "2011-01-05"],
            ),
        ],
    )
    def test_order_without_freq(self, index_dates, expected_dates, tz_naive_fixture):
        tz = tz_naive_fixture

        # without freq
        index = DatetimeIndex(index_dates, tz=tz, name="idx")
        expected = DatetimeIndex(expected_dates, tz=tz, name="idx")

        ordered = index.sort_values()
        tm.assert_index_equal(ordered, expected)
        assert ordered.freq is None

        ordered = index.sort_values(ascending=False)
        tm.assert_index_equal(ordered, expected[::-1])
        assert ordered.freq is None

        ordered, indexer = index.sort_values(return_indexer=True)
        tm.assert_index_equal(ordered, expected)

        exp = np.array([0, 4, 3, 1, 2])
        tm.assert_numpy_array_equal(indexer, exp, check_dtype=False)
        assert ordered.freq is None

        ordered, indexer = index.sort_values(return_indexer=True, ascending=False)
        tm.assert_index_equal(ordered, expected[::-1])

        exp = np.array([2, 1, 3, 4, 0])
        tm.assert_numpy_array_equal(indexer, exp, check_dtype=False)
        assert ordered.freq is None

    def test_drop_duplicates_metadata(self, freq_sample):
        # GH 10115
        idx = pd.date_range("2011-01-01", freq=freq_sample, periods=10, name="idx")
        result = idx.drop_duplicates()
        tm.assert_index_equal(idx, result)
        assert idx.freq == result.freq

        idx_dup = idx.append(idx)
        assert idx_dup.freq is None  # freq is reset
        result = idx_dup.drop_duplicates()
        expected = idx._with_freq(None)
        tm.assert_index_equal(result, expected)
        assert result.freq is None

    @pytest.mark.parametrize(
        "keep, expected, index",
        [
            ("first", np.concatenate(([False] * 10, [True] * 5)), np.arange(0, 10)),
            ("last", np.concatenate(([True] * 5, [False] * 10)), np.arange(5, 15)),
            (
                False,
                np.concatenate(([True] * 5, [False] * 5, [True] * 5)),
                np.arange(5, 10),
            ),
        ],
    )
    def test_drop_duplicates(self, freq_sample, keep, expected, index):
        # to check Index/Series compat
        idx = pd.date_range("2011-01-01", freq=freq_sample, periods=10, name="idx")
        idx = idx.append(idx[:5])

        tm.assert_numpy_array_equal(idx.duplicated(keep=keep), expected)
        expected = idx[~expected]

        result = idx.drop_duplicates(keep=keep)
        tm.assert_index_equal(result, expected)

        result = Series(idx).drop_duplicates(keep=keep)
        tm.assert_series_equal(result, Series(expected, index=index))

    def test_infer_freq(self, freq_sample):
        # GH 11018
        idx = pd.date_range("2011-01-01 09:00:00", freq=freq_sample, periods=10)
        result = pd.DatetimeIndex(idx.asi8, freq="infer")
        tm.assert_index_equal(idx, result)
        assert result.freq == freq_sample

    def test_nat(self, tz_naive_fixture):
        tz = tz_naive_fixture
        assert pd.DatetimeIndex._na_value is pd.NaT
        assert pd.DatetimeIndex([])._na_value is pd.NaT

        idx = pd.DatetimeIndex(["2011-01-01", "2011-01-02"], tz=tz)
        assert idx._can_hold_na

        tm.assert_numpy_array_equal(idx._isnan, np.array([False, False]))
        assert idx.hasnans is False
        tm.assert_numpy_array_equal(idx._nan_idxs, np.array([], dtype=np.intp))

        idx = pd.DatetimeIndex(["2011-01-01", "NaT"], tz=tz)
        assert idx._can_hold_na

        tm.assert_numpy_array_equal(idx._isnan, np.array([False, True]))
        assert idx.hasnans is True
        tm.assert_numpy_array_equal(idx._nan_idxs, np.array([1], dtype=np.intp))

    def test_equals(self):
        # GH 13107
        idx = pd.DatetimeIndex(["2011-01-01", "2011-01-02", "NaT"])
        assert idx.equals(idx)
        assert idx.equals(idx.copy())
        assert idx.equals(idx.astype(object))
        assert idx.astype(object).equals(idx)
        assert idx.astype(object).equals(idx.astype(object))
        assert not idx.equals(list(idx))
        assert not idx.equals(pd.Series(idx))

        idx2 = pd.DatetimeIndex(["2011-01-01", "2011-01-02", "NaT"], tz="US/Pacific")
        assert not idx.equals(idx2)
        assert not idx.equals(idx2.copy())
        assert not idx.equals(idx2.astype(object))
        assert not idx.astype(object).equals(idx2)
        assert not idx.equals(list(idx2))
        assert not idx.equals(pd.Series(idx2))

        # same internal, different tz
        idx3 = pd.DatetimeIndex(idx.asi8, tz="US/Pacific")
        tm.assert_numpy_array_equal(idx.asi8, idx3.asi8)
        assert not idx.equals(idx3)
        assert not idx.equals(idx3.copy())
        assert not idx.equals(idx3.astype(object))
        assert not idx.astype(object).equals(idx3)
        assert not idx.equals(list(idx3))
        assert not idx.equals(pd.Series(idx3))

        # check that we do not raise when comparing with OutOfBounds objects
        oob = pd.Index([datetime(2500, 1, 1)] * 3, dtype=object)
        assert not idx.equals(oob)
        assert not idx2.equals(oob)
        assert not idx3.equals(oob)

        # check that we do not raise when comparing with OutOfBounds dt64
        oob2 = oob.map(np.datetime64)
        assert not idx.equals(oob2)
        assert not idx2.equals(oob2)
        assert not idx3.equals(oob2)

    @pytest.mark.parametrize("values", [["20180101", "20180103", "20180105"], []])
    @pytest.mark.parametrize("freq", ["2D", Day(2), "2B", BDay(2), "48H", Hour(48)])
    @pytest.mark.parametrize("tz", [None, "US/Eastern"])
    def test_freq_setter(self, values, freq, tz):
        # GH 20678
        idx = DatetimeIndex(values, tz=tz)

        # can set to an offset, converting from string if necessary
        idx._data.freq = freq
        assert idx.freq == freq
        assert isinstance(idx.freq, DateOffset)

        # can reset to None
        idx._data.freq = None
        assert idx.freq is None

    def test_freq_setter_errors(self):
        # GH 20678
        idx = DatetimeIndex(["20180101", "20180103", "20180105"])

        # setting with an incompatible freq
        msg = (
            "Inferred frequency 2D from passed values does not conform to "
            "passed frequency 5D"
        )
        with pytest.raises(ValueError, match=msg):
            idx._data.freq = "5D"

        # setting with non-freq string
        with pytest.raises(ValueError, match="Invalid frequency"):
            idx._data.freq = "foo"

    def test_freq_view_safe(self):
        # Setting the freq for one DatetimeIndex shouldn't alter the freq
        #  for another that views the same data

        dti = pd.date_range("2016-01-01", periods=5)
        dta = dti._data

        dti2 = DatetimeIndex(dta)._with_freq(None)
        assert dti2.freq is None

        # Original was not altered
        assert dti.freq == "D"
        assert dta.freq == "D"


class TestBusinessDatetimeIndex:
    def setup_method(self, method):
        self.rng = bdate_range(START, END)

    def test_comparison(self):
        d = self.rng[10]

        comp = self.rng > d
        assert comp[11]
        assert not comp[9]

    def test_copy(self):
        cp = self.rng.copy()
        repr(cp)
        tm.assert_index_equal(cp, self.rng)

    def test_equals(self):
        assert not self.rng.equals(list(self.rng))

    def test_identical(self):
        t1 = self.rng.copy()
        t2 = self.rng.copy()
        assert t1.identical(t2)

        # name
        t1 = t1.rename("foo")
        assert t1.equals(t2)
        assert not t1.identical(t2)
        t2 = t2.rename("foo")
        assert t1.identical(t2)

        # freq
        t2v = Index(t2.values)
        assert t1.equals(t2v)
        assert not t1.identical(t2v)


class TestCustomDatetimeIndex:
    def setup_method(self, method):
        self.rng = bdate_range(START, END, freq="C")

    def test_comparison(self):
        d = self.rng[10]

        comp = self.rng > d
        assert comp[11]
        assert not comp[9]

    def test_copy(self):
        cp = self.rng.copy()
        repr(cp)
        tm.assert_index_equal(cp, self.rng)

    def test_equals(self):
        assert not self.rng.equals(list(self.rng))
