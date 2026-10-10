import pytest

import pandas as pd
import pandas._testing as tm

from pandas.tseries.offsets import (
    BaseOffset,
    Day,
    Hour,
    MonthEnd,
)


class TestFreq:
    @pytest.mark.parametrize("values", [["0 days", "2 days", "4 days"], []])
    @pytest.mark.parametrize("freq", ["2D", Day(2), "48h", Hour(48)])
    def test_freq_setter(self, values, freq):
        # GH#20678
        idx = pd.TimedeltaIndex(values)

        # can set to an offset, converting from string if necessary
        idx.freq = freq
        assert idx.freq == freq
        assert isinstance(idx.freq, BaseOffset)

        # can reset to None
        idx.freq = None
        assert idx.freq is None

    def test_with_freq_empty_requires_tick(self):
        idx = pd.TimedeltaIndex([])

        off = MonthEnd(1)
        msg = "TimedeltaArray/Index freq must be a Tick"
        with pytest.raises(TypeError, match=msg):
            idx._with_freq(off)

    def test_freq_setter_errors(self):
        # GH#20678
        idx = pd.TimedeltaIndex(["0 days", "2 days", "4 days"])

        # setting with an incompatible freq
        msg = (
            "Inferred frequency 2D from passed values does not conform to "
            "passed frequency 5D"
        )
        with pytest.raises(ValueError, match=msg):
            idx.freq = "5D"

        # setting with a non-fixed frequency
        msg = r"<2 \* BusinessDays> is a non-fixed frequency"
        with pytest.raises(ValueError, match=msg):
            idx.freq = "2B"

        # setting with non-freq string
        with pytest.raises(ValueError, match="Invalid frequency"):
            idx.freq = "foo"

    def test_freq_view_safe(self):
        # Setting the freq for one TimedeltaIndex shouldn't alter the freq
        #  for another that views the same data

        tdi = pd.TimedeltaIndex(["0 days", "2 days", "4 days"], freq="2D")
        tda = tdi._data

        tdi2 = pd.TimedeltaIndex(tda)._with_freq(None)
        assert tdi2.freq is None

        # Original was not altered. freq is now Index-level state.
        assert tdi.freq == "2D"

    @pytest.mark.parametrize("freq", ["2D", Day(2), "48h", Hour(48)])
    def test_set_freq(self, freq):
        # GH#61094
        idx = pd.TimedeltaIndex(["0 days", "2 days", "4 days"])
        result = idx.set_freq(freq)
        assert result.freq == freq
        assert isinstance(result.freq, BaseOffset)
        assert idx.freq is None
        tm.assert_index_equal(result, idx, check_freq=False)

        assert result.set_freq(None).freq is None
        assert result.freq == freq

    def test_set_freq_errors(self):
        # GH#61094
        idx = pd.TimedeltaIndex(["0 days", "2 days", "4 days"])
        msg = r"<2 \* BusinessDays> is a non-fixed frequency"
        with pytest.raises(ValueError, match=msg):
            idx.set_freq("2B")

        msg = "TimedeltaArray/Index freq must be a Tick"
        with pytest.raises(TypeError, match=msg):
            pd.TimedeltaIndex([]).set_freq(MonthEnd(1))
