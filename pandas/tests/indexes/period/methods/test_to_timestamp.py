from datetime import datetime

import numpy as np
import pytest

from pandas import (
    DatetimeIndex,
    NaT,
    PeriodIndex,
    Timedelta,
    Timestamp,
    date_range,
    period_range,
)
import pandas._testing as tm


class TestToTimestamp:
    def test_to_timestamp_non_contiguous(self):
        # GH#44100
        dti = date_range("2021-10-18", periods=9, freq="D", unit="ns")
        pi = dti.to_period()

        result = pi[::2].to_timestamp()
        expected = dti[::2].as_unit("us")
        tm.assert_index_equal(result, expected)

        result = pi._data[::2].to_timestamp()
        expected = dti._data[::2].as_unit("us")
        tm.assert_datetime_array_equal(result, expected)

        result = pi[::-1].to_timestamp()
        expected = dti[::-1].as_unit("us")
        tm.assert_index_equal(result, expected)

        result = pi._data[::-1].to_timestamp()
        expected = dti._data[::-1].as_unit("us")
        tm.assert_datetime_array_equal(result, expected)

        result = pi[::2][::-1].to_timestamp()
        expected = dti[::2][::-1].as_unit("us")
        tm.assert_index_equal(result, expected)

        result = pi._data[::2][::-1].to_timestamp()
        expected = dti._data[::2][::-1].as_unit("us")
        tm.assert_datetime_array_equal(result, expected)

    def test_to_timestamp_freq(self):
        idx = period_range("2017", periods=12, freq="Y-DEC")
        result = idx.to_timestamp()
        expected = date_range("2017", periods=12, freq="YS-JAN")
        tm.assert_index_equal(result, expected)

    def test_to_timestamp_pi_nat(self):
        # GH#7228
        index = PeriodIndex(["NaT", "2011-01", "2011-02"], freq="M", name="idx")

        result = index.to_timestamp("D")
        expected = DatetimeIndex(
            [NaT, datetime(2011, 1, 1), datetime(2011, 2, 1)],
            dtype="M8[us]",
            name="idx",
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

        msg = "Frequency must be positive, because it represents span: -2Y"
        with pytest.raises(ValueError, match=msg):
            result.to_period(freq="-2Y")

    def test_to_timestamp_preserve_name(self):
        index = period_range(freq="Y", start="1/1/2001", end="12/1/2009", name="foo")
        assert index.name == "foo"

        conv = index.to_timestamp("D")
        assert conv.name == "foo"

    def test_to_timestamp_quarterly_bug(self):
        years = np.arange(1960, 2000).repeat(4)
        quarters = np.tile(list(range(1, 5)), 40)

        pindex = PeriodIndex.from_fields(year=years, quarter=quarters)

        stamps = pindex.to_timestamp("D", "end")
        expected = DatetimeIndex([x.to_timestamp("D", "end") for x in pindex])
        tm.assert_index_equal(stamps, expected)
        assert stamps.freq == expected.freq

    def test_to_timestamp_pi_mult(self):
        idx = PeriodIndex(["2011-01", "NaT", "2011-02"], freq="2M", name="idx")

        result = idx.to_timestamp()
        expected = DatetimeIndex(
            ["2011-01-01", "NaT", "2011-02-01"], dtype="M8[us]", name="idx"
        )
        tm.assert_index_equal(result, expected)

        result = idx.to_timestamp(how="E")
        expected = DatetimeIndex(
            ["2011-02-28", "NaT", "2011-03-31"], dtype="M8[us]", name="idx"
        )
        expected = expected + Timedelta(1, "D") - Timedelta(1, "us")
        tm.assert_index_equal(result, expected)

    def test_to_timestamp_pi_combined(self):
        idx = period_range(start="2011", periods=2, freq="1D1h", name="idx")

        result = idx.to_timestamp()
        expected = DatetimeIndex(
            ["2011-01-01 00:00", "2011-01-02 01:00"], dtype="M8[us]", name="idx"
        )
        tm.assert_index_equal(result, expected)

        result = idx.to_timestamp(how="E")
        expected = DatetimeIndex(
            ["2011-01-02 00:59:59", "2011-01-03 01:59:59"], name="idx", dtype="M8[us]"
        )
        expected = expected + Timedelta(1, "s") - Timedelta(1, "us")
        tm.assert_index_equal(result, expected)

        result = idx.to_timestamp(how="E", freq="h")
        expected = DatetimeIndex(
            ["2011-01-02 00:00", "2011-01-03 01:00"], dtype="M8[us]", name="idx"
        )
        expected = expected + Timedelta(1, "h") - Timedelta(1, "us")
        tm.assert_index_equal(result, expected)

    def test_to_timestamp_1703(self):
        index = period_range("1/1/2012", periods=4, freq="D")

        result = index.to_timestamp()
        assert result[0] == Timestamp("1/1/2012")


def test_ms_to_timestamp_error_message():
    # https://github.com/pandas-dev/pandas/issues/58974#issuecomment-2164265446
    ix = period_range("2000", periods=3, freq="M")
    with pytest.raises(ValueError, match="for Period, please use 'M' instead of 'MS'"):
        ix.to_timestamp("MS")


@pytest.mark.parametrize("freq", ["ns", "1ns"])
def test_to_timestamp_nanosecond_target(freq):
    # GH#63760 a target freq that normalizes to nanoseconds must give a
    #  nanosecond-unit result rather than misinterpreting the ordinal
    pi = period_range("2020-01-01", periods=3, freq="D")
    result = pi.to_timestamp(freq)
    expected = date_range("2020-01-01", periods=3, freq="D", unit="ns")
    tm.assert_index_equal(result, expected)


@pytest.mark.parametrize("freq", ["D", "s", "us"])
def test_to_timestamp_from_nanosecond_period(freq):
    # GH#63760 converting a nanosecond PeriodIndex with a coarser target must
    #  yield a microsecond DatetimeIndex, not reinterpret the value as ns
    pi = period_range("2020-01-01", periods=3, freq="ns")
    result = pi.to_timestamp(freq)
    expected = DatetimeIndex(["2020-01-01"] * 3, dtype="M8[us]")
    tm.assert_index_equal(result, expected)


@pytest.mark.parametrize("freq", ["ns", "1ns"])
def test_to_timestamp_end_nanosecond_target(freq):
    # GH#63760 how="end" keeps nanosecond precision for a target freq that
    #  normalizes to nanoseconds
    pi = period_range("2020-01-01", periods=2, freq="D")
    result = pi.to_timestamp(freq, how="E")
    expected = DatetimeIndex(
        ["2020-01-01 23:59:59.999999999", "2020-01-02 23:59:59.999999999"],
        dtype="M8[ns]",
    )
    tm.assert_index_equal(result, expected)


def test_to_timestamp_end_from_nanosecond_period():
    # GH#63760 nanosecond periods keep their nanosecond end bounds with a
    #  coarser target freq
    pi = period_range("2020-01-01", periods=2, freq="ns")
    result = pi.to_timestamp("D", how="E")
    expected = DatetimeIndex(
        ["2020-01-01 00:00:00", "2020-01-01 00:00:00.000000001"], dtype="M8[ns]"
    )
    tm.assert_index_equal(result, expected)


def test_to_timestamp_end_bday_nanosecond_target():
    # GH#63760 the B roll-forward path keeps nanosecond precision too
    msg = "Period with BDay freq is deprecated"
    with tm.assert_produces_warning(FutureWarning, match=msg):
        pi = period_range("2020-01-02", periods=2, freq="B")
        result = pi.to_timestamp("ns", how="E")
    expected = DatetimeIndex(
        ["2020-01-02 23:59:59.999999999", "2020-01-03 23:59:59.999999999"],
        dtype="M8[ns]",
    )
    tm.assert_index_equal(result, expected)
