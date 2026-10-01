import pytest

import pandas as pd
import pandas._testing as tm

from pandas.tseries.offsets import (
    BaseOffset,
    BDay,
    Day,
    Hour,
)


class TestFreq:
    def test_freq_setter_errors(self):
        # GH#20678
        idx = pd.DatetimeIndex(["20180101", "20180103", "20180105"])

        # setting with an incompatible freq
        msg = (
            "Inferred frequency 2D from passed values does not conform to "
            "passed frequency 5D"
        )
        with pytest.raises(ValueError, match=msg):
            idx.freq = "5D"

        # setting with non-freq string
        with pytest.raises(ValueError, match="Invalid frequency"):
            idx.freq = "foo"

    @pytest.mark.parametrize("values", [["20180101", "20180103", "20180105"], []])
    @pytest.mark.parametrize("freq", ["2D", Day(2), "2B", BDay(2), "48h", Hour(48)])
    @pytest.mark.parametrize("tz", [None, "US/Eastern"])
    def test_freq_setter(self, values, freq, tz):
        # GH#20678
        idx = pd.DatetimeIndex(values, tz=tz)

        # can set to an offset, converting from string if necessary
        idx.freq = freq
        assert idx.freq == freq
        assert isinstance(idx.freq, BaseOffset)

        # can reset to None
        idx.freq = None
        assert idx.freq is None

    def test_freq_view_safe(self):
        # Setting the freq for one DatetimeIndex shouldn't alter the freq
        #  for another that views the same data

        dti = pd.date_range("2016-01-01", periods=5)
        dta = dti._data

        dti2 = pd.DatetimeIndex(dta)._with_freq(None)
        assert dti2.freq is None

        # Original was not altered. freq is now Index-level state.
        assert dti.freq == "D"

    @pytest.mark.parametrize("freq", ["2D", Day(2), "48h", Hour(48)])
    def test_set_freq(self, freq):
        # GH#61094
        idx = pd.DatetimeIndex(["20180101", "20180103", "20180105"])
        result = idx.set_freq(freq)
        assert result.freq == freq
        assert isinstance(result.freq, BaseOffset)
        assert idx.freq is None
        tm.assert_index_equal(result, idx, check_freq=False)

        assert result.set_freq(None).freq is None
        assert result.freq == freq

    def test_set_freq_errors(self):
        # GH#61094
        idx = pd.DatetimeIndex(["20180101", "20180103", "20180105"])
        msg = (
            "Inferred frequency 2D from passed values does not conform to "
            "passed frequency 5D"
        )
        with pytest.raises(ValueError, match=msg):
            idx.set_freq("5D")

        with pytest.raises(ValueError, match="Invalid frequency"):
            idx.set_freq("foo")

    def test_set_freq_not_cached_from_original(self):
        # GH#61094 the repr depends on freq, so it must not reuse the
        # original's cached formatting
        idx = pd.DatetimeIndex(["2020-01-01"])
        repr(idx)
        result = idx.set_freq("h")
        expected = pd.DatetimeIndex(["2020-01-01"], freq="h")
        assert repr(result) == repr(expected)
