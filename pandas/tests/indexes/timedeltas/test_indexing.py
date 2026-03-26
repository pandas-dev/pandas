from datetime import datetime
import re

import numpy as np
import pytest

from pandas.errors import Pandas4Warning

from pandas import (
    DataFrame,
    Index,
    NaT,
    Series,
    Timedelta,
    TimedeltaIndex,
    Timestamp,
    notna,
    offsets,
    timedelta_range,
    to_timedelta,
)
import pandas._testing as tm


class TestGetItem:
    def test_getitem_slice_keeps_name(self):
        # GH#4226, GH#59051
        msg = "'d' is deprecated and will be removed in a future version."
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            tdi = timedelta_range("1d", "5d", freq="h", name="timebucket")
        assert tdi[1:].name == tdi.name

    def test_getitem(self):
        idx1 = timedelta_range("1 day", "31 day", freq="D", name="idx")

        for idx in [idx1]:
            result = idx[0]
            assert result == Timedelta("1 day")

            result = idx[0:5]
            expected = timedelta_range("1 day", "5 day", freq="D", name="idx")
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

            result = idx[0:10:2]
            expected = timedelta_range("1 day", "9 day", freq="2D", name="idx")
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

            result = idx[-20:-5:3]
            expected = timedelta_range("12 day", "24 day", freq="3D", name="idx")
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

            result = idx[4::-1]
            expected = TimedeltaIndex(
                ["5 day", "4 day", "3 day", "2 day", "1 day"], freq="-1D", name="idx"
            )
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

    @pytest.mark.parametrize(
        "key",
        [
            Timestamp("1970-01-01"),
            Timestamp("1970-01-02"),
            datetime(1970, 1, 1),
            Timestamp("1970-01-03").to_datetime64(),
            # non-matching NA values
            np.datetime64("NaT", "ns"),
        ],
    )
    def test_timestamp_invalid_key(self, key):
        # GH#20464
        tdi = timedelta_range(0, periods=10)
        with pytest.raises(KeyError, match=re.escape(repr(key))):
            tdi.get_loc(key)


class TestGetLoc:
    def test_get_loc_key_unit_mismatch(self):
        idx = to_timedelta(["0 days", "1 days", "2 days"])
        key = idx[1].as_unit("ms")
        loc = idx.get_loc(key)
        assert loc == 1

    def test_get_loc_key_unit_mismatch_not_castable(self):
        tdi = to_timedelta(["0 days", "1 days", "2 days"]).astype("m8[s]")
        assert tdi.dtype == "m8[s]"
        key = tdi[0].as_unit("ns") + Timedelta(1)

        with pytest.raises(KeyError, match=r"Timedelta\('0 days 00:00:00.000000001'\)"):
            tdi.get_loc(key)

        assert key not in tdi

    def test_get_loc(self):
        idx = to_timedelta(["0 days", "1 days", "2 days"])

        # GH 16909
        assert idx.get_loc(idx[1].to_timedelta64()) == 1

        # GH 16896
        assert idx.get_loc("0 days") == 0

    def test_get_loc_nat(self):
        tidx = TimedeltaIndex(["1 days 01:00:00", "NaT", "2 days 01:00:00"])

        assert tidx.get_loc(NaT) == 1
        assert tidx.get_loc(None) == 1
        assert tidx.get_loc(float("nan")) == 1
        assert tidx.get_loc(np.nan) == 1


class TestGetIndexer:
    def test_get_indexer(self):
        idx = to_timedelta(["0 days", "1 days", "2 days"])
        tm.assert_numpy_array_equal(
            idx.get_indexer(idx), np.array([0, 1, 2], dtype=np.intp)
        )

        target = to_timedelta(["-1 hour", "12 hours", "1 day 1 hour"])
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, "pad"), np.array([-1, 0, 1], dtype=np.intp)
        )
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, "backfill"), np.array([0, 1, 2], dtype=np.intp)
        )
        tm.assert_numpy_array_equal(
            idx.get_indexer(target, "nearest"), np.array([0, 1, 1], dtype=np.intp)
        )

        res = idx.get_indexer(target, "nearest", tolerance=Timedelta("1 hour"))
        tm.assert_numpy_array_equal(res, np.array([0, -1, 1], dtype=np.intp))

    @pytest.mark.parametrize("method", ["pad", "backfill", "nearest"])
    def test_get_indexer_nat_target(self, method):
        # GH#32572 NaT in the target should not be matched
        tdi = to_timedelta(["0 days", "1 days", "2 days"])
        target = TimedeltaIndex([NaT])
        result = tdi.get_indexer(target, method=method)
        expected = np.array([-1], dtype=np.intp)
        tm.assert_numpy_array_equal(result, expected)


class TestWhere:
    def test_where_doesnt_retain_freq(self):
        tdi = timedelta_range("1 day", periods=3, freq="D", name="idx")
        cond = [True, True, False]
        expected = TimedeltaIndex([tdi[0], tdi[1], tdi[0]], freq=None, name="idx")

        result = tdi.where(cond, tdi[::-1])
        tm.assert_index_equal(result, expected)

    def test_where_invalid_dtypes(self, fixed_now_ts):
        tdi = timedelta_range("1 day", periods=3, freq="D", name="idx")

        tail = tdi[2:].tolist()
        i2 = Index([NaT, NaT, *tail])
        mask = notna(i2)

        expected = Index([NaT._value, NaT._value, *tail], dtype=object, name="idx")
        assert isinstance(expected[0], int)
        result = tdi.where(mask, i2.asi8)
        tm.assert_index_equal(result, expected)

        ts = i2 + fixed_now_ts
        expected = Index([ts[0], ts[1], *tail], dtype=object, name="idx")
        result = tdi.where(mask, ts)
        tm.assert_index_equal(result, expected)

        per = (i2 + fixed_now_ts).to_period("D")
        expected = Index([per[0], per[1], *tail], dtype=object, name="idx")
        result = tdi.where(mask, per)
        tm.assert_index_equal(result, expected)

        ts = fixed_now_ts
        expected = Index([ts, ts, *tail], dtype=object, name="idx")
        result = tdi.where(mask, ts)
        tm.assert_index_equal(result, expected)

    def test_where_mismatched_nat(self):
        tdi = timedelta_range("1 day", periods=3, freq="D", name="idx")
        cond = np.array([True, False, False])

        dtnat = np.datetime64("NaT", "ns")
        expected = Index([tdi[0], dtnat, dtnat], dtype=object, name="idx")
        assert expected[2] is dtnat
        result = tdi.where(cond, dtnat)
        tm.assert_index_equal(result, expected)


class TestTake:
    def test_take(self):
        # GH 10295
        idx1 = timedelta_range("1 day", "31 day", freq="D", name="idx")

        for idx in [idx1]:
            result = idx.take([0])
            assert result == Timedelta("1 day")

            result = idx.take([-1])
            assert result == Timedelta("31 day")

            result = idx.take([0, 1, 2])
            expected = timedelta_range("1 day", "3 day", freq="D", name="idx")
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

            result = idx.take([0, 2, 4])
            expected = timedelta_range("1 day", "5 day", freq="2D", name="idx")
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

            result = idx.take([7, 4, 1])
            expected = timedelta_range("8 day", "2 day", freq="-3D", name="idx")
            tm.assert_index_equal(result, expected)
            assert result.freq == expected.freq

            result = idx.take([3, 2, 5])
            expected = TimedeltaIndex(["4 day", "3 day", "6 day"], name="idx")
            tm.assert_index_equal(result, expected)
            assert result.freq is None

            result = idx.take([-3, 2, 5])
            expected = TimedeltaIndex(["29 day", "3 day", "6 day"], name="idx")
            tm.assert_index_equal(result, expected)
            assert result.freq is None

    def test_take_invalid_kwargs(self):
        idx = timedelta_range("1 day", "31 day", freq="D", name="idx")
        indices = [1, 6, 5, 9, 10, 13, 15, 3]

        msg = r"take\(\) got an unexpected keyword argument 'foo'"
        with pytest.raises(TypeError, match=msg):
            idx.take(indices, foo=2)

        msg = "the 'out' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            idx.take(indices, out=indices)

        msg = "the 'mode' parameter is not supported"
        with pytest.raises(ValueError, match=msg):
            idx.take(indices, mode="clip")

    def test_take_equiv_getitem(self):
        tds = ["1day 02:00:00", "1 day 04:00:00", "1 day 10:00:00"]
        idx = timedelta_range(start="1D", end="2D", freq="h", name="idx")
        expected = TimedeltaIndex(tds, freq=None, name="idx")

        taken1 = idx.take([2, 4, 10])
        taken2 = idx[[2, 4, 10]]

        for taken in [taken1, taken2]:
            tm.assert_index_equal(taken, expected)
            assert isinstance(taken, TimedeltaIndex)
            assert taken.freq is None
            assert taken.name == expected.name

    def test_take_fill_value(self):
        # GH 12631
        idx = TimedeltaIndex(["1 days", "2 days", "3 days"], name="xxx")
        result = idx.take(np.array([1, 0, -1]))
        expected = TimedeltaIndex(["2 days", "1 days", "3 days"], name="xxx")
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=NaT)
        expected = TimedeltaIndex(["2 days", "1 days", "NaT"], name="xxx")
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False)
        expected = TimedeltaIndex(["2 days", "1 days", "3 days"], name="xxx")
        tm.assert_index_equal(result, expected)

        msg = "When allow_fill=True, all indices must be >= -1"
        with pytest.raises(ValueError, match=msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with pytest.raises(ValueError, match=msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        msg = "index -5 is out of bounds for (axis 0 with )?size 3"
        with pytest.raises(IndexError, match=msg):
            idx.take(np.array([1, -5]))


class TestMaybeCastSliceBound:
    @pytest.fixture(params=["increasing", "decreasing", None])
    def monotonic(self, request):
        return request.param

    @pytest.fixture
    def tdi(self, monotonic):
        tdi = timedelta_range("1 Day", periods=10)
        if monotonic == "decreasing":
            tdi = tdi[::-1]
        elif monotonic is None:
            taker = np.arange(10, dtype=np.intp)
            np.random.default_rng(2).shuffle(taker)
            tdi = tdi.take(taker)
        return tdi

    def test_maybe_cast_slice_bound_invalid_str(self, tdi):
        # test the low-level _maybe_cast_slice_bound and that we get the
        #  expected exception+message all the way up the stack
        msg = (
            "cannot do slice indexing on TimedeltaIndex with these "
            r"indexers \[foo\] of type str"
        )
        with pytest.raises(TypeError, match=msg):
            tdi._maybe_cast_slice_bound("foo", side="left")
        with pytest.raises(TypeError, match=msg):
            tdi.get_slice_bound("foo", side="left")
        with pytest.raises(TypeError, match=msg):
            tdi.slice_locs("foo", None, None)

    def test_slice_invalid_str_with_timedeltaindex(
        self, tdi, frame_or_series, indexer_sl
    ):
        obj = frame_or_series(range(10), index=tdi)

        msg = (
            "cannot do slice indexing on TimedeltaIndex with these "
            r"indexers \[foo\] of type str"
        )
        with pytest.raises(TypeError, match=msg):
            indexer_sl(obj)["foo":]
        with pytest.raises(TypeError, match=msg):
            indexer_sl(obj)["foo":-1]
        with pytest.raises(TypeError, match=msg):
            indexer_sl(obj)[:"foo"]
        with pytest.raises(TypeError, match=msg):
            indexer_sl(obj)[tdi[0] : "foo"]


class TestStringSliceResolution:
    """Tests for GH#33603 - string resolution for TimedeltaIndex slicing."""

    def test_string_slice_reso_seconds_not_minutes(self):
        # GH#33603 - "720s" should have second resolution, not minute
        # even though 720s = 12 minutes exactly
        fs = 50000
        idx = to_timedelta(np.arange(900 * fs, dtype=int) / fs * 1e9, unit="ns")
        df = DataFrame({"dummy": np.arange(len(idx))}, index=idx)

        result = df["710s":"720s"]
        assert np.isclose(len(result) / fs, 11.0, atol=1e-4)

    def test_string_slice_reso_value_not_elevated(self):
        # GH#33603 - "60s" should have second resolution, not minute
        # "3600s" should have second resolution, not hour
        tdi = timedelta_range("0s", "7200s", freq="100ms")
        ser = Series(np.arange(len(tdi)), index=tdi)

        result_60 = ser["59s":"60s"]
        expected_60 = ser[Timedelta("59s") : Timedelta("60s") + Timedelta("999999us")]
        tm.assert_series_equal(result_60, expected_60)

        result_3600 = ser["3599s":"3600s"]
        expected_3600 = ser[
            Timedelta("3599s") : Timedelta("3600s") + Timedelta("999999us")
        ]
        tm.assert_series_equal(result_3600, expected_3600)

    def test_string_slice_reso_day_with_zero_hours(self):
        # GH#33603 - "2D 0 hours" should have hour resolution, not day
        ser = Series(
            np.arange(100),
            index=timedelta_range("1 days", periods=100, freq="h"),
        )
        str_count = len(ser[:"2D 0 hours"])
        td_count = (ser.index <= Timedelta("2D 0 hours")).sum()
        # With hour resolution, the upper bound is 2D 0h 59min 59.999999s,
        # so str_count should equal td_count (both include only up to 2D 0h)
        assert str_count == td_count

    def test_string_slice_consistency_with_timedelta(self):
        # GH#33603 - string slicing and Timedelta slicing should be
        # consistent (string slicing includes the resolution window)
        tdi = timedelta_range(0, "10s", freq="100ms")
        ser = Series(np.arange(len(tdi)), index=tdi)

        # String "3s" has second resolution -> includes 3.0s through 3.999...s
        str_result = ser.loc[:"3s"]
        td_result = ser.loc[: Timedelta("3s") + Timedelta("999999us")]
        tm.assert_series_equal(str_result, td_result)

    def test_string_slice_reso_fractional_seconds(self):
        # "3.5s" has value-based resolution "ms", which is finer than
        # string resolution "s", so effective resolution should be "ms"
        tdi = timedelta_range("0s", "10s", freq="100us")
        ser = Series(np.arange(len(tdi)), index=tdi)

        result = ser.loc[:"3.5s"]
        # With ms resolution, upper bound = 3.500s + 1ms - 1us = 3.500999s
        assert result.index[-1] <= Timedelta("3.500999s")
        assert result.index[-1] >= Timedelta("3.5s")

    def test_string_slice_reso_hh_mm_ss(self):
        # hh:mm:ss format should have second resolution
        tdi = timedelta_range("1:00:00", periods=200, freq="100ms")
        ser = Series(np.arange(len(tdi)), index=tdi)

        result = ser.loc[:"1:00:05"]
        assert result.index[-1] <= Timedelta("1:00:05.999999")
        assert result.index[-1] >= Timedelta("1:00:05")

    def test_string_slice_reso_hh_mm_ss_fractional(self):
        # hh:mm:ss.fff should have ms resolution
        tdi = timedelta_range("1:00:00", periods=2000, freq="100us")
        ser = Series(np.arange(len(tdi)), index=tdi)

        result = ser.loc[:"1:00:00.100"]
        assert result.index[-1] <= Timedelta("1:00:00.100999")
        assert result.index[-1] >= Timedelta("1:00:00.100")

    def test_parsed_string_to_bounds_uses_reso(self):
        # Verify _parsed_string_to_bounds uses the reso parameter
        tdi = timedelta_range("0s", "1000s", freq="100ms")

        parsed, reso = tdi._parse_with_reso("720s")
        lower, upper = tdi._parsed_string_to_bounds(reso, parsed)
        assert lower == Timedelta("720s")
        assert upper == Timedelta("720s") + Timedelta("999999us")

        parsed_60, reso_60 = tdi._parse_with_reso("60s")
        lower_60, upper_60 = tdi._parsed_string_to_bounds(reso_60, parsed_60)
        assert lower_60 == Timedelta("60s")
        assert upper_60 == Timedelta("60s") + Timedelta("999999us")


class TestContains:
    def test_contains_nonunique(self):
        # GH#9512
        for vals in (
            [0, 1, 0],
            [0, 0, -1],
            [0, -1, -1],
            ["00:01:00", "00:01:00", "00:02:00"],
            ["00:01:00", "00:01:00", "00:00:01"],
        ):
            idx = TimedeltaIndex(vals)
            assert idx[0] in idx

    def test_contains(self):
        # Checking for any NaT-like objects
        # GH#13603, GH#59051
        msg = "'d' is deprecated and will be removed in a future version."
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            td = to_timedelta(range(5), unit="d") + offsets.Hour(1)
        for v in [NaT, None, float("nan"), np.nan]:
            assert v not in td

        td = to_timedelta([NaT])
        for v in [NaT, None, float("nan"), np.nan]:
            assert v in td
