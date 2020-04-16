import numpy as np
import pytest

from pandas import Timedelta, Timestamp, date_range, timedelta_range, to_timedelta
import pandas._testing as tm

from pandas.tseries.offsets import Day, Second


class TestTimedeltas:
    def test_timedelta_range(self):

        expected = to_timedelta(np.arange(5), unit="D")
        result = timedelta_range("0 days", periods=5, freq="D")
        tm.assert_index_equal(result, expected)

        expected = to_timedelta(np.arange(11), unit="D")
        result = timedelta_range("0 days", "10 days", freq="D")
        tm.assert_index_equal(result, expected)

        expected = to_timedelta(np.arange(5), unit="D") + Second(2) + Day()
        result = timedelta_range("1 days, 00:00:02", "5 days, 00:00:02", freq="D")
        tm.assert_index_equal(result, expected)

        expected = to_timedelta([1, 3, 5, 7, 9], unit="D") + Second(2)
        result = timedelta_range("1 days, 00:00:02", periods=5, freq="2D")
        tm.assert_index_equal(result, expected)

        expected = to_timedelta(np.arange(50), unit="T") * 30
        result = timedelta_range("0 days", freq="30T", periods=50)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize(
        "periods, freq", [(3, "2D"), (5, "D"), (6, "19H12T"), (7, "16H"), (9, "12H")]
    )
    def test_linspace_behavior(self, periods, freq):
        # GH 20976
        result = timedelta_range(start="0 days", end="4 days", periods=periods)
        expected = timedelta_range(start="0 days", end="4 days", freq=freq)
        tm.assert_index_equal(result, expected)
        assert result.freq == freq

    def test_errors(self):
        # not enough params
        msg = (
            "Of the four parameters: start, end, periods, and freq, "
            "exactly three must be specified"
        )
        with pytest.raises(ValueError, match=msg):
            timedelta_range(start="0 days")

        with pytest.raises(ValueError, match=msg):
            timedelta_range(end="5 days")

        with pytest.raises(ValueError, match=msg):
            timedelta_range(periods=2)

        with pytest.raises(ValueError, match=msg):
            timedelta_range()

        # too many params
        with pytest.raises(ValueError, match=msg):
            timedelta_range(start="0 days", end="5 days", periods=10, freq="H")

    @pytest.mark.parametrize(
        "start, end, freq",
        [
            ("1D", "10D", "2D"),
            ("2D", "30D", "3D"),
            ("2s", "50s", "5s"),
            # tests that worked before GH 33498:
            ("4D", "16D", "3D"),
            ("8D", "16D", "40s"),
        ],
    )
    def test_timedelta_range_freq_divide_end(self, start, end, freq):
        # GH 33498 only the cases where `(end % freq) == 0` used to fail

        def mock_timedelta_range(start=None, end=None, **kwargs):
            epoch = Timestamp(0)
            if start is not None:
                start = epoch + Timedelta(start)
            if end is not None:
                end = epoch + Timedelta(end)
            result = date_range(start=start, end=end, **kwargs) - epoch
            result.freq = freq
            return result

        res = timedelta_range(start=start, end=end, freq=freq)
        exp = mock_timedelta_range(start=start, end=end, freq=freq)

        tm.assert_index_equal(res, exp)
        assert res.freq == exp.freq
