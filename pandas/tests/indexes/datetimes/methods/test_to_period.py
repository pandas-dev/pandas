from datetime import timezone

import dateutil.tz
from dateutil.tz import tzlocal
import pytest

from pandas._libs.tslibs.ccalendar import MONTHS
from pandas._libs.tslibs.offsets import MonthEnd
from pandas._libs.tslibs.period import INVALID_FREQ_ERR_MSG

from pandas import (
    DatetimeIndex,
    Period,
    PeriodIndex,
    Timestamp,
    date_range,
    period_range,
)
import pandas._testing as tm


class TestToPeriod:
    def test_dti_to_period(self):
        dti = date_range(start="1/1/2005", end="12/1/2005", freq="ME")
        pi1 = dti.to_period()
        pi2 = dti.to_period(freq="D")
        pi3 = dti.to_period(freq="3D")

        assert pi1[0] == Period("Jan 2005", freq="M")
        assert pi2[0] == Period("1/31/2005", freq="D")
        assert pi3[0] == Period("1/31/2005", freq="3D")

        assert pi1[-1] == Period("Nov 2005", freq="M")
        assert pi2[-1] == Period("11/30/2005", freq="D")
        assert pi3[-1], Period("11/30/2005", freq="3D")

        tm.assert_index_equal(pi1, period_range("1/1/2005", "11/1/2005", freq="M"))
        tm.assert_index_equal(
            pi2, period_range("1/1/2005", "11/1/2005", freq="M").asfreq("D")
        )
        tm.assert_index_equal(
            pi3, period_range("1/1/2005", "11/1/2005", freq="M").asfreq("3D")
        )

    @pytest.mark.parametrize("month", MONTHS)
    def test_to_period_quarterly(self, month):
        # make sure we can make the round trip
        freq = f"Q-{month}"
        rng = period_range("1989Q3", "1991Q3", freq=freq)
        stamps = rng.to_timestamp()
        result = stamps.to_period(freq)
        tm.assert_index_equal(rng, result)

    @pytest.mark.parametrize("off", ["BQE", "QS", "BQS"])
    def test_to_period_quarterlyish(self, off):
        rng = date_range("01-Jan-2012", periods=8, freq=off)
        prng = rng.to_period()
        assert prng.freq == "QE-DEC"

    @pytest.mark.parametrize("off", ["BYE", "YS", "BYS"])
    def test_to_period_annualish(self, off):
        rng = date_range("01-Jan-2012", periods=8, freq=off)
        prng = rng.to_period()
        assert prng.freq == "YE-DEC"

    @pytest.mark.parametrize(
        "off, expected_freq",
        [
            ("QS-APR", "Q-MAR"),
            ("BQS-APR", "Q-MAR"),
            ("BQE-APR", "Q-APR"),
            # start-offset month is shifted back one; JAN wraps to DEC
            ("QS-JAN", "Q-DEC"),
            ("QS-DEC", "Q-NOV"),
        ],
    )
    def test_to_period_quarterly_anchored(self, off, expected_freq):
        # GH#36939 - anchored quarter offsets should preserve the anchor
        rng = date_range("01-Apr-2012", periods=4, freq=off)
        result = rng.to_period()
        assert result.freqstr == expected_freq

    @pytest.mark.parametrize(
        "dates, expected_freq, expected_periods",
        [
            (
                ["2012-01-01", "2012-04-01", "2012-07-01", "2012-10-01"],
                "Q-DEC",
                ["2012Q1", "2012Q2", "2012Q3", "2012Q4"],
            ),
            (
                ["2012-02-01", "2012-05-01", "2012-08-01", "2012-11-01"],
                "Q-JAN",
                ["2013Q1", "2013Q2", "2013Q3", "2013Q4"],
            ),
            (
                ["2012-03-01", "2012-06-01", "2012-09-01", "2012-12-01"],
                "Q-FEB",
                ["2013Q1", "2013Q2", "2013Q3", "2013Q4"],
            ),
        ],
    )
    def test_to_period_quarter_start_inferred(
        self, dates, expected_freq, expected_periods
    ):
        # GH#36939 - without a set freq the anchor is only inferred up to
        #  mod 3; calendar quarter starts must keep the calendar Q-DEC labels
        dti = DatetimeIndex(dates)
        result = dti.to_period()
        expected = PeriodIndex(expected_periods, freq=expected_freq)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize(
        "off, expected_freq",
        [
            ("YS-APR", "Y-MAR"),
            ("BYS-APR", "Y-MAR"),
            ("BYE-APR", "Y-APR"),
            # start-offset month is shifted back one; JAN wraps to DEC
            ("YS-JAN", "Y-DEC"),
            ("YS-DEC", "Y-NOV"),
        ],
    )
    def test_to_period_annual_anchored(self, off, expected_freq):
        # GH#36939 - anchored annual offsets should preserve the anchor
        rng = date_range("01-Apr-2012", periods=3, freq=off)
        result = rng.to_period()
        assert result.freqstr == expected_freq

    def test_to_period_monthish(self):
        offsets = ["MS", "BME"]
        for off in offsets:
            rng = date_range("01-Jan-2012", periods=8, freq=off)
            prng = rng.to_period()
            assert prng.freqstr == "M"

        rng = date_range("01-Jan-2012", periods=8, freq="ME")
        prng = rng.to_period()
        assert prng.freqstr == "M"

        with pytest.raises(ValueError, match=INVALID_FREQ_ERR_MSG):
            date_range("01-Jan-2012", periods=8, freq="EOM")

    @pytest.mark.parametrize(
        "freq_offset, freq_period",
        [
            ("2ME", "2M"),
            (MonthEnd(2), MonthEnd(2)),
        ],
    )
    def test_dti_to_period_2monthish(self, freq_offset, freq_period):
        dti = date_range("2020-01-01", periods=3, freq=freq_offset)
        pi = dti.to_period()

        tm.assert_index_equal(pi, period_range("2020-01", "2020-05", freq=freq_period))

    @pytest.mark.parametrize(
        "freq", ["2ME", "1me", "2QE", "2QE-SEP", "1YE", "ye", "2YE-MAR"]
    )
    def test_to_period_frequency_M_Q_Y_raises(self, freq):
        msg = f"Invalid frequency: {freq}"

        rng = date_range("01-Jan-2012", periods=8, freq="ME")
        with pytest.raises(ValueError, match=msg):
            rng.to_period(freq)

    def test_to_period_infer(self):
        # https://github.com/pandas-dev/pandas/issues/33358
        rng = date_range(
            start="2019-12-22 06:40:00+00:00",
            end="2019-12-22 08:45:00+00:00",
            freq="5min",
        )

        with tm.assert_produces_warning(UserWarning, match="drop timezone info"):
            pi1 = rng.to_period("5min")

        with tm.assert_produces_warning(UserWarning, match="drop timezone info"):
            pi2 = rng.to_period()

        tm.assert_index_equal(pi1, pi2)

    @pytest.mark.filterwarnings(r"ignore:PeriodDtype\[B\] is deprecated:FutureWarning")
    def test_period_dt64_round_trip(self):
        dti = date_range("1/1/2000", "1/7/2002", freq="B", unit="ns")
        pi = dti.to_period()
        tm.assert_index_equal(pi.to_timestamp(), dti.as_unit("us"))

        dti = date_range("1/1/2000", "1/7/2002", freq="B", unit="ns")
        pi = dti.to_period(freq="h")
        tm.assert_index_equal(pi.to_timestamp(), dti.as_unit("us"))

    def test_to_period_millisecond(self):
        index = DatetimeIndex(
            [
                Timestamp("2007-01-01 10:11:12.123456Z"),
                Timestamp("2007-01-01 10:11:13.789123Z"),
            ]
        )

        with tm.assert_produces_warning(UserWarning, match="drop timezone info"):
            period = index.to_period(freq="ms")
        assert 2 == len(period)
        assert period[0] == Period("2007-01-01 10:11:12.123Z", "ms")
        assert period[1] == Period("2007-01-01 10:11:13.789Z", "ms")

    def test_to_period_microsecond(self):
        index = DatetimeIndex(
            [
                Timestamp("2007-01-01 10:11:12.123456Z"),
                Timestamp("2007-01-01 10:11:13.789123Z"),
            ]
        )

        with tm.assert_produces_warning(UserWarning, match="drop timezone info"):
            period = index.to_period(freq="us")
        assert 2 == len(period)
        assert period[0] == Period("2007-01-01 10:11:12.123456Z", "us")
        assert period[1] == Period("2007-01-01 10:11:13.789123Z", "us")

    @pytest.mark.parametrize(
        "tz",
        [
            "US/Eastern",
            timezone.utc,
            tzlocal(),
            "dateutil/US/Eastern",
            dateutil.tz.tzutc(),
        ],
    )
    def test_to_period_tz(self, tz):
        ts = date_range("1/1/2000", "2/1/2000", tz=tz)

        with tm.assert_produces_warning(UserWarning, match="drop timezone info"):
            result = ts.to_period()[0]
            expected = ts[0].to_period(ts.freq)

        assert result == expected

        expected = date_range("1/1/2000", "2/1/2000").to_period()

        with tm.assert_produces_warning(UserWarning, match="drop timezone info"):
            result = ts.to_period(ts.freq)

        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("tz", ["Etc/GMT-1", "Etc/GMT+1"])
    def test_to_period_tz_utc_offset_consistency(self, tz):
        # GH#22905
        ts = date_range("1/1/2000", "2/1/2000", tz="Etc/GMT-1")
        with tm.assert_produces_warning(UserWarning, match="drop timezone info"):
            result = ts.to_period()[0]
            expected = ts[0].to_period(ts.freq)
            assert result == expected

    def test_to_period_nofreq(self):
        idx = DatetimeIndex(["2000-01-01", "2000-01-02", "2000-01-04"])
        msg = "You must pass a freq argument as current index has none."
        with pytest.raises(ValueError, match=msg):
            idx.to_period()

        idx = DatetimeIndex(["2000-01-01", "2000-01-02", "2000-01-03"], freq="infer")
        assert idx.freqstr == "D"
        expected = PeriodIndex(["2000-01-01", "2000-01-02", "2000-01-03"], freq="D")
        tm.assert_index_equal(idx.to_period(), expected)

        # GH#7606
        idx = DatetimeIndex(["2000-01-01", "2000-01-02", "2000-01-03"])
        assert idx.freqstr is None
        tm.assert_index_equal(idx.to_period(), expected)

    @pytest.mark.parametrize("freq", ["2BME", "SME-15", "2BMS"])
    def test_to_period_offsets_not_supported(self, freq):
        # GH#56243
        msg = "|".join(
            [
                f"Invalid frequency: {freq}",
                f"{freq} is not supported as period frequency",
            ]
        )

        ts = date_range("1/1/2012", periods=4, freq=freq)
        with pytest.raises(ValueError, match=msg):
            ts.to_period()
