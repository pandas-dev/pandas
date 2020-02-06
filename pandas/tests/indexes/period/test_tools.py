from datetime import datetime

import numpy as np
import pytest

from pandas._libs.tslibs import IncompatibleFrequency

from pandas import (
    DatetimeIndex,
    NaT,
    Period,
    PeriodIndex,
    Timedelta,
    Timestamp,
    date_range,
    period_range,
)
import pandas._testing as tm


class TestPeriodRepresentation:
    """
    Wish to match NumPy units
    """

    def _check_freq(self, freq, base_date):
        rng = period_range(start=base_date, periods=10, freq=freq)
        exp = np.arange(10, dtype=np.int64)

        tm.assert_numpy_array_equal(rng.asi8, exp)

    def test_annual(self):
        self._check_freq("A", 1970)

    def test_monthly(self):
        self._check_freq("M", "1970-01")

    @pytest.mark.parametrize("freq", ["W-THU", "D", "B", "H", "T", "S", "L", "U", "N"])
    def test_freq(self, freq):
        self._check_freq(freq, "1970-01-01")


class TestSearchsorted:
    @pytest.mark.parametrize("freq", ["D", "2D"])
    def test_searchsorted(self, freq):
        pidx = PeriodIndex(
            ["2014-01-01", "2014-01-02", "2014-01-03", "2014-01-04", "2014-01-05"],
            freq=freq,
        )

        p1 = Period("2014-01-01", freq=freq)
        assert pidx.searchsorted(p1) == 0

        p2 = Period("2014-01-04", freq=freq)
        assert pidx.searchsorted(p2) == 3

        assert pidx.searchsorted(NaT) == 0

        msg = "Input has different freq=H from PeriodArray"
        with pytest.raises(IncompatibleFrequency, match=msg):
            pidx.searchsorted(Period("2014-01-01", freq="H"))

        msg = "Input has different freq=5D from PeriodArray"
        with pytest.raises(IncompatibleFrequency, match=msg):
            pidx.searchsorted(Period("2014-01-01", freq="5D"))

    def test_searchsorted_invalid(self):
        pidx = PeriodIndex(
            ["2014-01-01", "2014-01-02", "2014-01-03", "2014-01-04", "2014-01-05"],
            freq="D",
        )

        other = np.array([0, 1], dtype=np.int64)

        msg = "|".join(
            [
                "searchsorted requires compatible dtype or scalar",
                "Unexpected type for 'value'",
            ]
        )
        with pytest.raises(TypeError, match=msg):
            pidx.searchsorted(other)

        with pytest.raises(TypeError, match=msg):
            pidx.searchsorted(other.astype("timedelta64[ns]"))

        with pytest.raises(TypeError, match=msg):
            pidx.searchsorted(np.timedelta64(4))

        with pytest.raises(TypeError, match=msg):
            pidx.searchsorted(np.timedelta64("NaT", "ms"))

        with pytest.raises(TypeError, match=msg):
            pidx.searchsorted(np.datetime64(4, "ns"))

        with pytest.raises(TypeError, match=msg):
            pidx.searchsorted(np.datetime64("NaT", "ns"))


class TestPeriodIndexConversion:
    def test_tolist(self):
        index = period_range(freq="A", start="1/1/2001", end="12/1/2009")
        rs = index.tolist()
        for x in rs:
            assert isinstance(x, Period)

        recon = PeriodIndex(rs)
        tm.assert_index_equal(index, recon)


class TestToTimestamp:
    def test_to_timestamp_freq(self):
        idx = period_range("2017", periods=12, freq="A-DEC")
        result = idx.to_timestamp()
        expected = date_range("2017", periods=12, freq="AS-JAN")
        tm.assert_index_equal(result, expected)

    def test_to_timestamp_pi_nat(self):
        # GH#7228
        index = PeriodIndex(["NaT", "2011-01", "2011-02"], freq="M", name="idx")

        result = index.to_timestamp("D")
        expected = DatetimeIndex(
            [NaT, datetime(2011, 1, 1), datetime(2011, 2, 1)], name="idx"
        )
        tm.assert_index_equal(result, expected)
        assert result.name == "idx"

        result2 = result.to_period(freq="M")
        tm.assert_index_equal(result2, index)
        assert result2.name == "idx"

        result3 = result.to_period(freq="3M")
        exp = PeriodIndex(["NaT", "2011-01", "2011-02"], freq="3M", name="idx")
        tm.assert_index_equal(result3, exp)
        assert result3.freqstr == "3M"

        msg = "Frequency must be positive, because it represents span: -2A"
        with pytest.raises(ValueError, match=msg):
            result.to_period(freq="-2A")

    def test_to_timestamp_preserve_name(self):
        index = period_range(freq="A", start="1/1/2001", end="12/1/2009", name="foo")
        assert index.name == "foo"

        conv = index.to_timestamp("D")
        assert conv.name == "foo"

    def test_to_timestamp_quarterly_bug(self):
        years = np.arange(1960, 2000).repeat(4)
        quarters = np.tile(list(range(1, 5)), 40)

        pindex = PeriodIndex(year=years, quarter=quarters)

        stamps = pindex.to_timestamp("D", "end")
        expected = DatetimeIndex([x.to_timestamp("D", "end") for x in pindex])
        tm.assert_index_equal(stamps, expected)

    def test_to_timestamp_pi_mult(self):
        idx = PeriodIndex(["2011-01", "NaT", "2011-02"], freq="2M", name="idx")

        result = idx.to_timestamp()
        expected = DatetimeIndex(["2011-01-01", "NaT", "2011-02-01"], name="idx")
        tm.assert_index_equal(result, expected)

        result = idx.to_timestamp(how="E")
        expected = DatetimeIndex(["2011-02-28", "NaT", "2011-03-31"], name="idx")
        expected = expected + Timedelta(1, "D") - Timedelta(1, "ns")
        tm.assert_index_equal(result, expected)

    def test_to_timestamp_pi_combined(self):
        idx = period_range(start="2011", periods=2, freq="1D1H", name="idx")

        result = idx.to_timestamp()
        expected = DatetimeIndex(["2011-01-01 00:00", "2011-01-02 01:00"], name="idx")
        tm.assert_index_equal(result, expected)

        result = idx.to_timestamp(how="E")
        expected = DatetimeIndex(
            ["2011-01-02 00:59:59", "2011-01-03 01:59:59"], name="idx"
        )
        expected = expected + Timedelta(1, "s") - Timedelta(1, "ns")
        tm.assert_index_equal(result, expected)

        result = idx.to_timestamp(how="E", freq="H")
        expected = DatetimeIndex(["2011-01-02 00:00", "2011-01-03 01:00"], name="idx")
        expected = expected + Timedelta(1, "h") - Timedelta(1, "ns")
        tm.assert_index_equal(result, expected)

    def test_to_timestamp_1703(self):
        index = period_range("1/1/2012", periods=4, freq="D")

        result = index.to_timestamp()
        assert result[0] == Timestamp("1/1/2012")
