import pytest

from pandas._libs.tslibs.ccalendar import MONTHS
from pandas._libs.tslibs.frequencies import INVALID_FREQ_ERR_MSG

from pandas import Period, date_range, period_range
import pandas._testing as tm


class TestToPeriod:
    def test_dti_to_period(self):
        dti = date_range(start="1/1/2005", end="12/1/2005", freq="M")
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
        freq = "Q-{month}".format(month=month)
        rng = period_range("1989Q3", "1991Q3", freq=freq)
        stamps = rng.to_timestamp()
        result = stamps.to_period(freq)
        tm.assert_index_equal(rng, result)

    @pytest.mark.parametrize("off", ["BQ", "QS", "BQS"])
    def test_to_period_quarterlyish(self, off):
        rng = date_range("01-Jan-2012", periods=8, freq=off)
        prng = rng.to_period()
        assert prng.freq == "Q-DEC"

    @pytest.mark.parametrize("off", ["BA", "AS", "BAS"])
    def test_to_period_annualish(self, off):
        rng = date_range("01-Jan-2012", periods=8, freq=off)
        prng = rng.to_period()
        assert prng.freq == "A-DEC"

    def test_to_period_monthish(self):
        offsets = ["MS", "BM"]
        for off in offsets:
            rng = date_range("01-Jan-2012", periods=8, freq=off)
            prng = rng.to_period()
            assert prng.freq == "M"

        rng = date_range("01-Jan-2012", periods=8, freq="M")
        prng = rng.to_period()
        assert prng.freq == "M"

        with pytest.raises(ValueError, match=INVALID_FREQ_ERR_MSG):
            date_range("01-Jan-2012", periods=8, freq="EOM")

    def test_period_dt64_round_trip(self):
        dti = date_range("1/1/2000", "1/7/2002", freq="B")
        pi = dti.to_period()
        tm.assert_index_equal(pi.to_timestamp(), dti)

        dti = date_range("1/1/2000", "1/7/2002", freq="B")
        pi = dti.to_period(freq="H")
        tm.assert_index_equal(pi.to_timestamp(), dti)
