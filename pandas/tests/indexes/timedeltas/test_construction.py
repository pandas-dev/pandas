from datetime import timedelta

import numpy as np
import pytest

import pandas as pd
from pandas import Timedelta, TimedeltaIndex, timedelta_range, to_timedelta
from pandas.core.arrays import TimedeltaArray
import pandas.util.testing as tm


class TestTimedeltaIndex:
    def test_int64_nocopy(self):
        # GH#23539 check that a copy isn't made when we pass int64 data
        #  and copy=False
        arr = np.arange(10, dtype=np.int64)
        tdi = TimedeltaIndex(arr, copy=False)
        assert tdi._data._data.base is arr

    def test_infer_from_tdi(self):
        # GH#23539
        # fast-path for inferring a frequency if the passed data already
        #  has one
        tdi = pd.timedelta_range("1 second", periods=10 ** 7, freq="1s")

        result = pd.TimedeltaIndex(tdi, freq="infer")
        assert result.freq == tdi.freq

        # check that inferred_freq was not called by checking that the
        #  value has not been cached
        assert "inferred_freq" not in getattr(result, "_cache", {})

    def test_infer_from_tdi_mismatch(self):
        # GH#23539
        # fast-path for invalidating a frequency if the passed data already
        #  has one and it does not match the `freq` input
        tdi = pd.timedelta_range("1 second", periods=100, freq="1s")

        msg = (
            "Inferred frequency .* from passed values does "
            "not conform to passed frequency"
        )
        with pytest.raises(ValueError, match=msg):
            TimedeltaIndex(tdi, freq="D")

        with pytest.raises(ValueError, match=msg):
            # GH#23789
            TimedeltaArray(tdi, freq="D")

    def test_dt64_data_invalid(self):
        # GH#23539
        # passing tz-aware DatetimeIndex raises, naive or ndarray[datetime64]
        #  raise as of GH#29794
        dti = pd.date_range("2016-01-01", periods=3)

        msg = "cannot be converted to timedelta64"
        with pytest.raises(TypeError, match=msg):
            TimedeltaIndex(dti.tz_localize("Europe/Brussels"))

        with pytest.raises(TypeError, match=msg):
            TimedeltaIndex(dti)

        with pytest.raises(TypeError, match=msg):
            TimedeltaIndex(np.asarray(dti))

    def test_float64_ns_rounded(self):
        # GH#23539 without specifying a unit, floats are regarded as nanos,
        #  and fractional portions are truncated
        tdi = TimedeltaIndex([2.3, 9.7])
        expected = TimedeltaIndex([2, 9])
        tm.assert_index_equal(tdi, expected)

        # integral floats are non-lossy
        tdi = TimedeltaIndex([2.0, 9.0])
        expected = TimedeltaIndex([2, 9])
        tm.assert_index_equal(tdi, expected)

        # NaNs get converted to NaT
        tdi = TimedeltaIndex([2.0, np.nan])
        expected = TimedeltaIndex([pd.Timedelta(nanoseconds=2), pd.NaT])
        tm.assert_index_equal(tdi, expected)

    def test_float64_unit_conversion(self):
        # GH#23539
        tdi = TimedeltaIndex([1.5, 2.25], unit="D")
        expected = TimedeltaIndex([Timedelta(days=1.5), Timedelta(days=2.25)])
        tm.assert_index_equal(tdi, expected)

    def test_construction_base_constructor(self):
        arr = [pd.Timedelta("1 days"), pd.NaT, pd.Timedelta("3 days")]
        tm.assert_index_equal(pd.Index(arr), pd.TimedeltaIndex(arr))
        tm.assert_index_equal(pd.Index(np.array(arr)), pd.TimedeltaIndex(np.array(arr)))

        arr = [np.nan, pd.NaT, pd.Timedelta("1 days")]
        tm.assert_index_equal(pd.Index(arr), pd.TimedeltaIndex(arr))
        tm.assert_index_equal(pd.Index(np.array(arr)), pd.TimedeltaIndex(np.array(arr)))

    def test_constructor(self):
        expected = TimedeltaIndex(
            [
                "1 days",
                "1 days 00:00:05",
                "2 days",
                "2 days 00:00:02",
                "0 days 00:00:03",
            ]
        )
        result = TimedeltaIndex(
            [
                "1 days",
                "1 days, 00:00:05",
                np.timedelta64(2, "D"),
                timedelta(days=2, seconds=2),
                pd.offsets.Second(3),
            ]
        )
        tm.assert_index_equal(result, expected)

        # unicode
        result = TimedeltaIndex(
            [
                "1 days",
                "1 days, 00:00:05",
                np.timedelta64(2, "D"),
                timedelta(days=2, seconds=2),
                pd.offsets.Second(3),
            ]
        )

        expected = TimedeltaIndex(
            ["0 days 00:00:00", "0 days 00:00:01", "0 days 00:00:02"]
        )
        tm.assert_index_equal(TimedeltaIndex(range(3), unit="s"), expected)
        expected = TimedeltaIndex(
            ["0 days 00:00:00", "0 days 00:00:05", "0 days 00:00:09"]
        )
        tm.assert_index_equal(TimedeltaIndex([0, 5, 9], unit="s"), expected)
        expected = TimedeltaIndex(
            ["0 days 00:00:00.400", "0 days 00:00:00.450", "0 days 00:00:01.200"]
        )
        tm.assert_index_equal(TimedeltaIndex([400, 450, 1200], unit="ms"), expected)

    def test_constructor_iso(self):
        # GH #21877
        expected = timedelta_range("1s", periods=9, freq="s")
        durations = ["P0DT0H0M{}S".format(i) for i in range(1, 10)]
        result = to_timedelta(durations)
        tm.assert_index_equal(result, expected)

    def test_constructor_coverage(self):
        rng = timedelta_range("1 days", periods=10.5)
        exp = timedelta_range("1 days", periods=10)
        tm.assert_index_equal(rng, exp)

        msg = "periods must be a number, got foo"
        with pytest.raises(TypeError, match=msg):
            timedelta_range(start="1 days", periods="foo", freq="D")

        with pytest.raises(TypeError):
            TimedeltaIndex("1 days")

        # generator expression
        gen = (timedelta(i) for i in range(10))
        result = TimedeltaIndex(gen)
        expected = TimedeltaIndex([timedelta(i) for i in range(10)])
        tm.assert_index_equal(result, expected)

        # NumPy string array
        strings = np.array(["1 days", "2 days", "3 days"])
        result = TimedeltaIndex(strings)
        expected = to_timedelta([1, 2, 3], unit="d")
        tm.assert_index_equal(result, expected)

        from_ints = TimedeltaIndex(expected.asi8)
        tm.assert_index_equal(from_ints, expected)

        # non-conforming freq
        msg = (
            "Inferred frequency None from passed values does not conform to"
            " passed frequency D"
        )
        with pytest.raises(ValueError, match=msg):
            TimedeltaIndex(["1 days", "2 days", "4 days"], freq="D")

        msg = (
            "Of the four parameters: start, end, periods, and freq, exactly"
            " three must be specified"
        )
        with pytest.raises(ValueError, match=msg):
            timedelta_range(periods=10, freq="D")

    def test_constructor_name(self):
        idx = timedelta_range(start="1 days", periods=1, freq="D", name="TEST")
        assert idx.name == "TEST"

        # GH10025
        idx2 = TimedeltaIndex(idx, name="something else")
        assert idx2.name == "something else"

    def test_constructor_no_precision_raises(self):
        # GH-24753, GH-24739

        msg = "with no precision is not allowed"
        with pytest.raises(ValueError, match=msg):
            pd.TimedeltaIndex(["2000"], dtype="timedelta64")

        with pytest.raises(ValueError, match=msg):
            pd.Index(["2000"], dtype="timedelta64")

    def test_constructor_wrong_precision_raises(self):
        with pytest.raises(ValueError):
            pd.TimedeltaIndex(["2000"], dtype="timedelta64[us]")
