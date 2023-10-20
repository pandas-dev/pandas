"""
Tests for DatetimeIndex timezone-related methods
"""
from datetime import (
    date,
    datetime,
    time,
    timedelta,
    timezone,
    tzinfo,
)

from dateutil.tz import gettz
import numpy as np
import pytest
import pytz

from pandas._libs.tslibs import (
    conversion,
    timezones,
)

import pandas as pd
from pandas import (
    DatetimeIndex,
    Index,
    Timestamp,
    bdate_range,
    date_range,
    isna,
    to_datetime,
)
import pandas._testing as tm


class FixedOffset(tzinfo):
    """Fixed offset in minutes east from UTC."""

    def __init__(self, offset, name) -> None:
        self.__offset = timedelta(minutes=offset)
        self.__name = name

    def utcoffset(self, dt):
        return self.__offset

    def tzname(self, dt):
        return self.__name

    def dst(self, dt):
        return timedelta(0)


fixed_off_no_name = FixedOffset(-330, None)


class TestDatetimeIndexTimezones:
    # ------------------------------------------------------------
    # DatetimeIndex.__new__

    @pytest.mark.parametrize("prefix", ["", "dateutil/"])
    def test_dti_constructor_static_tzinfo(self, prefix):
        # it works!
        index = DatetimeIndex([datetime(2012, 1, 1)], tz=prefix + "EST")
        index.hour
        index[0]

    def test_dti_constructor_with_fixed_tz(self):
        off = FixedOffset(420, "+07:00")
        start = datetime(2012, 3, 11, 5, 0, 0, tzinfo=off)
        end = datetime(2012, 6, 11, 5, 0, 0, tzinfo=off)
        rng = date_range(start=start, end=end)
        assert off == rng.tz

        rng2 = date_range(start, periods=len(rng), tz=off)
        tm.assert_index_equal(rng, rng2)

        rng3 = date_range("3/11/2012 05:00:00+07:00", "6/11/2012 05:00:00+07:00")
        assert (rng.values == rng3.values).all()

    @pytest.mark.parametrize("tzstr", ["US/Eastern", "dateutil/US/Eastern"])
    def test_dti_convert_datetime_list(self, tzstr):
        dr = date_range("2012-06-02", periods=10, tz=tzstr, name="foo")
        dr2 = DatetimeIndex(list(dr), name="foo", freq="D")
        tm.assert_index_equal(dr, dr2)

    def test_dti_construction_univalent(self):
        rng = date_range("03/12/2012 00:00", periods=10, freq="W-FRI", tz="US/Eastern")
        rng2 = DatetimeIndex(data=rng, tz="US/Eastern")
        tm.assert_index_equal(rng, rng2)

    @pytest.mark.parametrize("tz", [pytz.timezone("US/Eastern"), gettz("US/Eastern")])
    def test_dti_from_tzaware_datetime(self, tz):
        d = [datetime(2012, 8, 19, tzinfo=tz)]

        index = DatetimeIndex(d)
        assert timezones.tz_compare(index.tz, tz)

    @pytest.mark.parametrize("tzstr", ["US/Eastern", "dateutil/US/Eastern"])
    def test_dti_tz_constructors(self, tzstr):
        """Test different DatetimeIndex constructions with timezone
        Follow-up of GH#4229
        """
        arr = ["11/10/2005 08:00:00", "11/10/2005 09:00:00"]

        idx1 = to_datetime(arr).tz_localize(tzstr)
        idx2 = date_range(start="2005-11-10 08:00:00", freq="h", periods=2, tz=tzstr)
        idx2 = idx2._with_freq(None)  # the others all have freq=None
        idx3 = DatetimeIndex(arr, tz=tzstr)
        idx4 = DatetimeIndex(np.array(arr), tz=tzstr)

        for other in [idx2, idx3, idx4]:
            tm.assert_index_equal(idx1, other)

    # -------------------------------------------------------------
    # Unsorted

    @pytest.mark.parametrize(
        "dtype",
        [None, "datetime64[ns, CET]", "datetime64[ns, EST]", "datetime64[ns, UTC]"],
    )
    def test_date_accessor(self, dtype):
        # Regression test for GH#21230
        expected = np.array([date(2018, 6, 4), pd.NaT])

        index = DatetimeIndex(["2018-06-04 10:00:00", pd.NaT], dtype=dtype)
        result = index.date

        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype",
        [None, "datetime64[ns, CET]", "datetime64[ns, EST]", "datetime64[ns, UTC]"],
    )
    def test_time_accessor(self, dtype):
        # Regression test for GH#21267
        expected = np.array([time(10, 20, 30), pd.NaT])

        index = DatetimeIndex(["2018-06-04 10:20:30", pd.NaT], dtype=dtype)
        result = index.time

        tm.assert_numpy_array_equal(result, expected)

    def test_timetz_accessor(self, tz_naive_fixture):
        # GH21358
        tz = timezones.maybe_get_tz(tz_naive_fixture)

        expected = np.array([time(10, 20, 30, tzinfo=tz), pd.NaT])

        index = DatetimeIndex(["2018-06-04 10:20:30", pd.NaT], tz=tz)
        result = index.timetz

        tm.assert_numpy_array_equal(result, expected)

    def test_dti_drop_dont_lose_tz(self):
        # GH#2621
        ind = date_range("2012-12-01", periods=10, tz="utc")
        ind = ind.drop(ind[-1])

        assert ind.tz is not None

    def test_dti_tz_conversion_freq(self, tz_naive_fixture):
        # GH25241
        t3 = DatetimeIndex(["2019-01-01 10:00"], freq="h")
        assert t3.tz_localize(tz=tz_naive_fixture).freq == t3.freq
        t4 = DatetimeIndex(["2019-01-02 12:00"], tz="UTC", freq="min")
        assert t4.tz_convert(tz="UTC").freq == t4.freq

    def test_drop_dst_boundary(self):
        # see gh-18031
        tz = "Europe/Brussels"
        freq = "15min"

        start = Timestamp("201710290100", tz=tz)
        end = Timestamp("201710290300", tz=tz)
        index = date_range(start=start, end=end, freq=freq)

        expected = DatetimeIndex(
            [
                "201710290115",
                "201710290130",
                "201710290145",
                "201710290200",
                "201710290215",
                "201710290230",
                "201710290245",
                "201710290200",
                "201710290215",
                "201710290230",
                "201710290245",
                "201710290300",
            ],
            tz=tz,
            freq=freq,
            ambiguous=[
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                False,
                False,
                False,
                False,
                False,
            ],
        )
        result = index.drop(index[0])
        tm.assert_index_equal(result, expected)

    def test_date_range_localize(self):
        rng = date_range("3/11/2012 03:00", periods=15, freq="h", tz="US/Eastern")
        rng2 = DatetimeIndex(["3/11/2012 03:00", "3/11/2012 04:00"], tz="US/Eastern")
        rng3 = date_range("3/11/2012 03:00", periods=15, freq="h")
        rng3 = rng3.tz_localize("US/Eastern")

        tm.assert_index_equal(rng._with_freq(None), rng3)

        # DST transition time
        val = rng[0]
        exp = Timestamp("3/11/2012 03:00", tz="US/Eastern")

        assert val.hour == 3
        assert exp.hour == 3
        assert val == exp  # same UTC value
        tm.assert_index_equal(rng[:2], rng2)

        # Right before the DST transition
        rng = date_range("3/11/2012 00:00", periods=2, freq="h", tz="US/Eastern")
        rng2 = DatetimeIndex(
            ["3/11/2012 00:00", "3/11/2012 01:00"], tz="US/Eastern", freq="h"
        )
        tm.assert_index_equal(rng, rng2)
        exp = Timestamp("3/11/2012 00:00", tz="US/Eastern")
        assert exp.hour == 0
        assert rng[0] == exp
        exp = Timestamp("3/11/2012 01:00", tz="US/Eastern")
        assert exp.hour == 1
        assert rng[1] == exp

        rng = date_range("3/11/2012 00:00", periods=10, freq="h", tz="US/Eastern")
        assert rng[2].hour == 3

    def test_timestamp_equality_different_timezones(self):
        utc_range = date_range("1/1/2000", periods=20, tz="UTC")
        eastern_range = utc_range.tz_convert("US/Eastern")
        berlin_range = utc_range.tz_convert("Europe/Berlin")

        for a, b, c in zip(utc_range, eastern_range, berlin_range):
            assert a == b
            assert b == c
            assert a == c

        assert (utc_range == eastern_range).all()
        assert (utc_range == berlin_range).all()
        assert (berlin_range == eastern_range).all()

    def test_dti_equals_with_tz(self):
        left = date_range("1/1/2011", periods=100, freq="h", tz="utc")
        right = date_range("1/1/2011", periods=100, freq="h", tz="US/Eastern")

        assert not left.equals(right)

    @pytest.mark.parametrize("tzstr", ["US/Eastern", "dateutil/US/Eastern"])
    def test_dti_tz_nat(self, tzstr):
        idx = DatetimeIndex([Timestamp("2013-1-1", tz=tzstr), pd.NaT])

        assert isna(idx[1])
        assert idx[0].tzinfo is not None

    @pytest.mark.parametrize("tzstr", ["US/Eastern", "dateutil/US/Eastern"])
    def test_dti_with_timezone_repr(self, tzstr):
        rng = date_range("4/13/2010", "5/6/2010")

        rng_eastern = rng.tz_localize(tzstr)

        rng_repr = repr(rng_eastern)
        assert "2010-04-13 00:00:00" in rng_repr

    @pytest.mark.parametrize("tzstr", ["US/Eastern", "dateutil/US/Eastern"])
    def test_dti_take_dont_lose_meta(self, tzstr):
        rng = date_range("1/1/2000", periods=20, tz=tzstr)

        result = rng.take(range(5))
        assert result.tz == rng.tz
        assert result.freq == rng.freq

    @pytest.mark.parametrize("tzstr", ["US/Eastern", "dateutil/US/Eastern"])
    def test_utc_box_timestamp_and_localize(self, tzstr):
        tz = timezones.maybe_get_tz(tzstr)

        rng = date_range("3/11/2012", "3/12/2012", freq="h", tz="utc")
        rng_eastern = rng.tz_convert(tzstr)

        expected = rng[-1].astimezone(tz)

        stamp = rng_eastern[-1]
        assert stamp == expected
        assert stamp.tzinfo == expected.tzinfo

        # right tzinfo
        rng = date_range("3/13/2012", "3/14/2012", freq="h", tz="utc")
        rng_eastern = rng.tz_convert(tzstr)
        # test not valid for dateutil timezones.
        # assert 'EDT' in repr(rng_eastern[0].tzinfo)
        assert "EDT" in repr(rng_eastern[0].tzinfo) or "tzfile" in repr(
            rng_eastern[0].tzinfo
        )

    @pytest.mark.parametrize("tz", [pytz.timezone("US/Central"), gettz("US/Central")])
    def test_with_tz(self, tz):
        # just want it to work
        start = datetime(2011, 3, 12, tzinfo=pytz.utc)
        dr = bdate_range(start, periods=50, freq=pd.offsets.Hour())
        assert dr.tz is pytz.utc

        # DateRange with naive datetimes
        dr = bdate_range("1/1/2005", "1/1/2009", tz=pytz.utc)
        dr = bdate_range("1/1/2005", "1/1/2009", tz=tz)

        # normalized
        central = dr.tz_convert(tz)
        assert central.tz is tz
        naive = central[0].to_pydatetime().replace(tzinfo=None)
        comp = conversion.localize_pydatetime(naive, tz).tzinfo
        assert central[0].tz is comp

        # compare vs a localized tz
        naive = dr[0].to_pydatetime().replace(tzinfo=None)
        comp = conversion.localize_pydatetime(naive, tz).tzinfo
        assert central[0].tz is comp

        # datetimes with tzinfo set
        dr = bdate_range(
            datetime(2005, 1, 1, tzinfo=pytz.utc), datetime(2009, 1, 1, tzinfo=pytz.utc)
        )
        msg = "Start and end cannot both be tz-aware with different timezones"
        with pytest.raises(Exception, match=msg):
            bdate_range(datetime(2005, 1, 1, tzinfo=pytz.utc), "1/1/2009", tz=tz)

    @pytest.mark.parametrize("prefix", ["", "dateutil/"])
    def test_field_access_localize(self, prefix):
        strdates = ["1/1/2012", "3/1/2012", "4/1/2012"]
        rng = DatetimeIndex(strdates, tz=prefix + "US/Eastern")
        assert (rng.hour == 0).all()

        # a more unusual time zone, #1946
        dr = date_range(
            "2011-10-02 00:00", freq="h", periods=10, tz=prefix + "America/Atikokan"
        )

        expected = Index(np.arange(10, dtype=np.int32))
        tm.assert_index_equal(dr.hour, expected)

    @pytest.mark.parametrize("tz", [pytz.timezone("US/Eastern"), gettz("US/Eastern")])
    def test_dti_convert_tz_aware_datetime_datetime(self, tz):
        # GH#1581
        dates = [datetime(2000, 1, 1), datetime(2000, 1, 2), datetime(2000, 1, 3)]

        dates_aware = [conversion.localize_pydatetime(x, tz) for x in dates]
        result = DatetimeIndex(dates_aware)
        assert timezones.tz_compare(result.tz, tz)

        converted = to_datetime(dates_aware, utc=True)
        ex_vals = np.array([Timestamp(x).as_unit("ns")._value for x in dates_aware])
        tm.assert_numpy_array_equal(converted.asi8, ex_vals)
        assert converted.tz is timezone.utc
