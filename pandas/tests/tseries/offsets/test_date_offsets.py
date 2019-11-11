from datetime import date, datetime, time as dt_time, timedelta
from typing import Dict, List, Tuple, Type

import numpy as np
import pytest

from pandas._libs.tslibs import (
    NaT,
    OutOfBoundsDatetime,
    Timestamp,
    conversion,
    timezones,
)
from pandas._libs.tslibs.frequencies import (
    INVALID_FREQ_ERR_MSG,
    get_freq_code,
    get_freq_str,
)
import pandas._libs.tslibs.offsets as liboffsets
from pandas._libs.tslibs.offsets import ApplyTypeError
import pandas.compat as compat
from pandas.compat.numpy import np_datetime64_compat

from pandas.core.indexes.datetimes import DatetimeIndex, _to_M8, date_range
from pandas.core.series import Series
import pandas.util.testing as tm

from pandas.io.pickle import read_pickle
from pandas.tseries.frequencies import _offset_map, get_offset
from pandas.tseries.holiday import USFederalHolidayCalendar
import pandas.tseries.offsets as offsets
from pandas.tseries.offsets import (
    BaseOffset,
    DateOffset,
    Day,
    Easter,
    LastWeekOfMonth,
    MonthBegin,
    MonthEnd,
    Nano,
    QuarterBegin,
    QuarterEnd,
    SemiMonthBegin,
    SemiMonthEnd,
    Tick,
    Week,
    WeekOfMonth,
    YearBegin,
    YearEnd,
)

from .common import assert_offset_equal, assert_onOffset


class WeekDay:
    # TODO: Remove: This is not used outside of tests
    MON = 0
    TUE = 1
    WED = 2
    THU = 3
    FRI = 4
    SAT = 5
    SUN = 6


####
# Misc function tests
####


def test_to_M8():
    valb = datetime(2007, 10, 1)
    valu = _to_M8(valb)
    assert isinstance(valu, np.datetime64)


#####
# DateOffset Tests
#####
_ApplyCases = List[Tuple[BaseOffset, Dict[datetime, datetime]]]


class Base:
    _offset = None  # type: Type[DateOffset]
    d = Timestamp(datetime(2008, 1, 2))

    timezones = [
        None,
        "UTC",
        "Asia/Tokyo",
        "US/Eastern",
        "dateutil/Asia/Tokyo",
        "dateutil/US/Pacific",
    ]

    def _get_offset(self, klass, value=1, normalize=False):
        print(klass)
        # create instance from offset class
        if klass is LastWeekOfMonth:
            klass = klass(n=value, weekday=5, normalize=normalize)
        elif klass is WeekOfMonth:
            klass = klass(n=value, week=1, weekday=5, normalize=normalize)
        elif klass is Week:
            klass = klass(n=value, weekday=5, normalize=normalize)
        elif klass is DateOffset:
            klass = klass(days=value, normalize=normalize)
        else:
            klass = klass(value, normalize=normalize)
        return klass

    def test_apply_out_of_range(self, tz_naive_fixture):
        tz = tz_naive_fixture
        if self._offset is None:
            return

        # try to create an out-of-bounds result timestamp; if we can't create
        # the offset skip
        try:
            offset = self._get_offset(self._offset, value=10000)

            result = Timestamp("20080101") + offset
            assert isinstance(result, datetime)
            assert result.tzinfo is None

            # Check tz is preserved
            t = Timestamp("20080101", tz=tz)
            result = t + offset
            assert isinstance(result, datetime)
            assert t.tzinfo == result.tzinfo

        except OutOfBoundsDatetime:
            pass
        except (ValueError, KeyError):
            # we are creating an invalid offset
            # so ignore
            pass

    def test_offsets_compare_equal(self):
        # root cause of GH#456: __ne__ was not implemented
        if self._offset is None:
            return
        offset1 = self._offset()
        offset2 = self._offset()
        assert not offset1 != offset2
        assert offset1 == offset2

    def test_rsub(self):
        if self._offset is None or not hasattr(self, "offset2"):
            # i.e. skip for TestCommon and YQM subclasses that do not have
            # offset2 attr
            return
        assert self.d - self.offset2 == (-self.offset2).apply(self.d)

    def test_radd(self):
        if self._offset is None or not hasattr(self, "offset2"):
            # i.e. skip for TestCommon and YQM subclasses that do not have
            # offset2 attr
            return
        assert self.d + self.offset2 == self.offset2 + self.d

    def test_sub(self):
        if self._offset is None or not hasattr(self, "offset2"):
            # i.e. skip for TestCommon and YQM subclasses that do not have
            # offset2 attr
            return
        off = self.offset2
        msg = "Cannot subtract datetime from offset"
        with pytest.raises(TypeError, match=msg):
            off - self.d

        assert 2 * off - off == off
        assert self.d - self.offset2 == self.d + self._offset(-2)
        assert self.d - self.offset2 == self.d - (2 * off - off)

    def testMult1(self):
        if self._offset is None or not hasattr(self, "offset1"):
            # i.e. skip for TestCommon and YQM subclasses that do not have
            # offset1 attr
            return
        assert self.d + 10 * self.offset1 == self.d + self._offset(10)
        assert self.d + 5 * self.offset1 == self.d + self._offset(5)

    def testMult2(self):
        if self._offset is None:
            return
        assert self.d + (-5 * self._offset(-10)) == self.d + self._offset(50)
        assert self.d + (-3 * self._offset(-2)) == self.d + self._offset(6)

    def test_compare_str(self):
        # GH#23524
        # comparing to strings that cannot be cast to DateOffsets should
        #  not raise for __eq__ or __ne__
        if self._offset is None:
            return
        off = self._get_offset(self._offset)

        assert not off == "infer"
        assert off != "foo"
        # Note: inequalities are only implemented for Tick subclasses;
        #  tests for this are in test_ticks


class TestCommon(Base):
    # exected value created by Base._get_offset
    # are applied to 2011/01/01 09:00 (Saturday)
    # used for .apply and .rollforward
    expecteds = {
        "Day": Timestamp("2011-01-02 09:00:00"),
        "DateOffset": Timestamp("2011-01-02 09:00:00"),
        "MonthBegin": Timestamp("2011-02-01 09:00:00"),
        "MonthEnd": Timestamp("2011-01-31 09:00:00"),
        "SemiMonthEnd": Timestamp("2011-01-15 09:00:00"),
        "SemiMonthBegin": Timestamp("2011-01-15 09:00:00"),
        "YearBegin": Timestamp("2012-01-01 09:00:00"),
        "YearEnd": Timestamp("2011-12-31 09:00:00"),
        "QuarterBegin": Timestamp("2011-03-01 09:00:00"),
        "QuarterEnd": Timestamp("2011-03-31 09:00:00"),
        "WeekOfMonth": Timestamp("2011-01-08 09:00:00"),
        "LastWeekOfMonth": Timestamp("2011-01-29 09:00:00"),
        "Week": Timestamp("2011-01-08 09:00:00"),
        "Easter": Timestamp("2011-04-24 09:00:00"),
        "Hour": Timestamp("2011-01-01 10:00:00"),
        "Minute": Timestamp("2011-01-01 09:01:00"),
        "Second": Timestamp("2011-01-01 09:00:01"),
        "Milli": Timestamp("2011-01-01 09:00:00.001000"),
        "Micro": Timestamp("2011-01-01 09:00:00.000001"),
        "Nano": Timestamp(np_datetime64_compat("2011-01-01T09:00:00.000000001Z")),
    }

    def test_immutable(self, offset_types):
        # GH#21341 check that __setattr__ raises
        offset = self._get_offset(offset_types)
        with pytest.raises(AttributeError):
            offset.normalize = True
        with pytest.raises(AttributeError):
            offset.n = 91

    def test_return_type(self, offset_types):
        offset = self._get_offset(offset_types)

        # make sure that we are returning a Timestamp
        result = Timestamp("20080101") + offset
        assert isinstance(result, Timestamp)

        # make sure that we are returning NaT
        assert NaT + offset is NaT
        assert offset + NaT is NaT

        assert NaT - offset is NaT
        assert (-offset).apply(NaT) is NaT

    def test_offset_n(self, offset_types):
        offset = self._get_offset(offset_types)
        assert offset.n == 1

        neg_offset = offset * -1
        assert neg_offset.n == -1

        mul_offset = offset * 3
        assert mul_offset.n == 3

    def test_offset_timedelta64_arg(self, offset_types):
        # check that offset._validate_n raises TypeError on a timedelt64
        #  object
        off = self._get_offset(offset_types)

        td64 = np.timedelta64(4567, "s")
        with pytest.raises(TypeError, match="argument must be an integer"):
            type(off)(n=td64, **off.kwds)

    def test_offset_mul_ndarray(self, offset_types):
        off = self._get_offset(offset_types)

        expected = np.array([[off, off * 2], [off * 3, off * 4]])

        result = np.array([[1, 2], [3, 4]]) * off
        tm.assert_numpy_array_equal(result, expected)

        result = off * np.array([[1, 2], [3, 4]])
        tm.assert_numpy_array_equal(result, expected)

    def test_offset_freqstr(self, offset_types):
        offset = self._get_offset(offset_types)

        freqstr = offset.freqstr
        if freqstr not in ("<Easter>", "<DateOffset: days=1>", "LWOM-SAT"):
            code = get_offset(freqstr)
            assert offset.rule_code == code

    def _check_offsetfunc_works(self, offset, funcname, dt, expected, normalize=False):

        if normalize and issubclass(offset, Tick):
            # normalize=True disallowed for Tick subclasses GH#21427
            return

        offset_s = self._get_offset(offset, normalize=normalize)
        func = getattr(offset_s, funcname)

        result = func(dt)
        assert isinstance(result, Timestamp)
        assert result == expected

        result = func(Timestamp(dt))
        assert isinstance(result, Timestamp)
        assert result == expected

        # see gh-14101
        exp_warning = None
        ts = Timestamp(dt) + Nano(5)

        if (
            offset_s.__class__.__name__ == "DateOffset"
            and (funcname == "apply" or normalize)
            and ts.nanosecond > 0
        ):
            exp_warning = UserWarning

        # test nanosecond is preserved
        with tm.assert_produces_warning(exp_warning, check_stacklevel=False):
            result = func(ts)
        assert isinstance(result, Timestamp)
        if normalize is False:
            assert result == expected + Nano(5)
        else:
            assert result == expected

        if isinstance(dt, np.datetime64):
            # test tz when input is datetime or Timestamp
            return

        for tz in self.timezones:
            expected_localize = expected.tz_localize(tz)
            tz_obj = timezones.maybe_get_tz(tz)
            dt_tz = conversion.localize_pydatetime(dt, tz_obj)

            result = func(dt_tz)
            assert isinstance(result, Timestamp)
            assert result == expected_localize

            result = func(Timestamp(dt, tz=tz))
            assert isinstance(result, Timestamp)
            assert result == expected_localize

            # see gh-14101
            exp_warning = None
            ts = Timestamp(dt, tz=tz) + Nano(5)

            if (
                offset_s.__class__.__name__ == "DateOffset"
                and (funcname == "apply" or normalize)
                and ts.nanosecond > 0
            ):
                exp_warning = UserWarning

            # test nanosecond is preserved
            with tm.assert_produces_warning(exp_warning, check_stacklevel=False):
                result = func(ts)
            assert isinstance(result, Timestamp)
            if normalize is False:
                assert result == expected_localize + Nano(5)
            else:
                assert result == expected_localize

    def test_apply(self, offset_types):
        sdt = datetime(2011, 1, 1, 9, 0)
        ndt = np_datetime64_compat("2011-01-01 09:00Z")

        for dt in [sdt, ndt]:
            expected = self.expecteds[offset_types.__name__]
            self._check_offsetfunc_works(offset_types, "apply", dt, expected)

            expected = Timestamp(expected.date())
            self._check_offsetfunc_works(
                offset_types, "apply", dt, expected, normalize=True
            )

    def test_rollforward(self, offset_types):
        expecteds = self.expecteds.copy()

        # result will not be changed if the target is on the offset
        no_changes = [
            "Day",
            "MonthBegin",
            "SemiMonthBegin",
            "YearBegin",
            "Week",
            "Hour",
            "Minute",
            "Second",
            "Milli",
            "Micro",
            "Nano",
            "DateOffset",
        ]
        for n in no_changes:
            expecteds[n] = Timestamp("2011/01/01 09:00")

        # but be changed when normalize=True
        norm_expected = expecteds.copy()
        for k in norm_expected:
            norm_expected[k] = Timestamp(norm_expected[k].date())

        normalized = {
            "Day": Timestamp("2011-01-02 00:00:00"),
            "DateOffset": Timestamp("2011-01-02 00:00:00"),
            "MonthBegin": Timestamp("2011-02-01 00:00:00"),
            "SemiMonthBegin": Timestamp("2011-01-15 00:00:00"),
            "YearBegin": Timestamp("2012-01-01 00:00:00"),
            "Week": Timestamp("2011-01-08 00:00:00"),
            "Hour": Timestamp("2011-01-01 00:00:00"),
            "Minute": Timestamp("2011-01-01 00:00:00"),
            "Second": Timestamp("2011-01-01 00:00:00"),
            "Milli": Timestamp("2011-01-01 00:00:00"),
            "Micro": Timestamp("2011-01-01 00:00:00"),
        }
        norm_expected.update(normalized)

        sdt = datetime(2011, 1, 1, 9, 0)
        ndt = np_datetime64_compat("2011-01-01 09:00Z")

        for dt in [sdt, ndt]:
            expected = expecteds[offset_types.__name__]
            self._check_offsetfunc_works(offset_types, "rollforward", dt, expected)
            expected = norm_expected[offset_types.__name__]
            self._check_offsetfunc_works(
                offset_types, "rollforward", dt, expected, normalize=True
            )

    def test_rollback(self, offset_types):
        expecteds = {
            "MonthEnd": Timestamp("2010-12-31 09:00:00"),
            "SemiMonthEnd": Timestamp("2010-12-31 09:00:00"),
            "YearEnd": Timestamp("2010-12-31 09:00:00"),
            "QuarterBegin": Timestamp("2010-12-01 09:00:00"),
            "QuarterEnd": Timestamp("2010-12-31 09:00:00"),
            "WeekOfMonth": Timestamp("2010-12-11 09:00:00"),
            "LastWeekOfMonth": Timestamp("2010-12-25 09:00:00"),
            "Easter": Timestamp("2010-04-04 09:00:00"),
        }

        # result will not be changed if the target is on the offset
        for n in [
            "Day",
            "MonthBegin",
            "SemiMonthBegin",
            "YearBegin",
            "Week",
            "Hour",
            "Minute",
            "Second",
            "Milli",
            "Micro",
            "Nano",
            "DateOffset",
        ]:
            expecteds[n] = Timestamp("2011/01/01 09:00")

        # but be changed when normalize=True
        norm_expected = expecteds.copy()
        for k in norm_expected:
            norm_expected[k] = Timestamp(norm_expected[k].date())

        normalized = {
            "Day": Timestamp("2010-12-31 00:00:00"),
            "DateOffset": Timestamp("2010-12-31 00:00:00"),
            "MonthBegin": Timestamp("2010-12-01 00:00:00"),
            "SemiMonthBegin": Timestamp("2010-12-15 00:00:00"),
            "YearBegin": Timestamp("2010-01-01 00:00:00"),
            "Week": Timestamp("2010-12-25 00:00:00"),
            "Hour": Timestamp("2011-01-01 00:00:00"),
            "Minute": Timestamp("2011-01-01 00:00:00"),
            "Second": Timestamp("2011-01-01 00:00:00"),
            "Milli": Timestamp("2011-01-01 00:00:00"),
            "Micro": Timestamp("2011-01-01 00:00:00"),
        }
        norm_expected.update(normalized)

        sdt = datetime(2011, 1, 1, 9, 0)
        ndt = np_datetime64_compat("2011-01-01 09:00Z")

        for dt in [sdt, ndt]:
            expected = expecteds[offset_types.__name__]
            self._check_offsetfunc_works(offset_types, "rollback", dt, expected)

            expected = norm_expected[offset_types.__name__]
            self._check_offsetfunc_works(
                offset_types, "rollback", dt, expected, normalize=True
            )

    def test_onOffset(self, offset_types):
        dt = self.expecteds[offset_types.__name__]
        offset_s = self._get_offset(offset_types)
        assert offset_s.onOffset(dt)

        # when normalize=True, onOffset checks time is 00:00:00
        if issubclass(offset_types, Tick):
            # normalize=True disallowed for Tick subclasses GH#21427
            return
        offset_n = self._get_offset(offset_types, normalize=True)
        assert not offset_n.onOffset(dt)

        date = datetime(dt.year, dt.month, dt.day)
        assert offset_n.onOffset(date)

    def test_add(self, offset_types, tz_naive_fixture):
        tz = tz_naive_fixture
        dt = datetime(2011, 1, 1, 9, 0)

        offset_s = self._get_offset(offset_types)
        expected = self.expecteds[offset_types.__name__]

        result_dt = dt + offset_s
        result_ts = Timestamp(dt) + offset_s
        for result in [result_dt, result_ts]:
            assert isinstance(result, Timestamp)
            assert result == expected

        expected_localize = expected.tz_localize(tz)
        result = Timestamp(dt, tz=tz) + offset_s
        assert isinstance(result, Timestamp)
        assert result == expected_localize

        # normalize=True, disallowed for Tick subclasses GH#21427
        if issubclass(offset_types, Tick):
            return
        offset_s = self._get_offset(offset_types, normalize=True)
        expected = Timestamp(expected.date())

        result_dt = dt + offset_s
        result_ts = Timestamp(dt) + offset_s
        for result in [result_dt, result_ts]:
            assert isinstance(result, Timestamp)
            assert result == expected

        expected_localize = expected.tz_localize(tz)
        result = Timestamp(dt, tz=tz) + offset_s
        assert isinstance(result, Timestamp)
        assert result == expected_localize

    def test_pickle_v0_15_2(self, datapath):
        offsets = {
            "DateOffset": DateOffset(years=1),
            "MonthBegin": MonthBegin(1),
            "Day": Day(1),
            "YearBegin": YearBegin(1),
            "Week": Week(1),
        }

        pickle_path = datapath("tseries", "offsets", "data", "dateoffset_0_15_2.pickle")
        # This code was executed once on v0.15.2 to generate the pickle:
        # with open(pickle_path, 'wb') as f: pickle.dump(offsets, f)
        #
        tm.assert_dict_equal(offsets, read_pickle(pickle_path))


class TestDateOffset(Base):
    def setup_method(self, method):
        self.d = Timestamp(datetime(2008, 1, 2))
        _offset_map.clear()

    def test_repr(self):
        repr(DateOffset())
        repr(DateOffset(2))
        repr(2 * DateOffset())
        repr(2 * DateOffset(months=2))

    def test_mul(self):
        assert DateOffset(2) == 2 * DateOffset(1)
        assert DateOffset(2) == DateOffset(1) * 2

    def test_constructor(self):

        assert (self.d + DateOffset(months=2)) == datetime(2008, 3, 2)
        assert (self.d - DateOffset(months=2)) == datetime(2007, 11, 2)

        assert (self.d + DateOffset(2)) == datetime(2008, 1, 4)

        assert not DateOffset(2).isAnchored()
        assert DateOffset(1).isAnchored()

        d = datetime(2008, 1, 31)
        assert (d + DateOffset(months=1)) == datetime(2008, 2, 29)

    def test_copy(self):
        assert DateOffset(months=2).copy() == DateOffset(months=2)

    def test_eq(self):
        offset1 = DateOffset(days=1)
        offset2 = DateOffset(days=365)

        assert offset1 != offset2


class TestWeek(Base):
    _offset = Week
    d = Timestamp(datetime(2008, 1, 2))
    offset1 = _offset()
    offset2 = _offset(2)

    def test_repr(self):
        assert repr(Week(weekday=0)) == "<Week: weekday=0>"
        assert repr(Week(n=-1, weekday=0)) == "<-1 * Week: weekday=0>"
        assert repr(Week(n=-2, weekday=0)) == "<-2 * Weeks: weekday=0>"

    def test_corner(self):
        with pytest.raises(ValueError):
            Week(weekday=7)

        with pytest.raises(ValueError, match="Day must be"):
            Week(weekday=-1)

    def test_isAnchored(self):
        assert Week(weekday=0).isAnchored()
        assert not Week().isAnchored()
        assert not Week(2, weekday=2).isAnchored()
        assert not Week(2).isAnchored()

    offset_cases = []
    # not business week
    offset_cases.append(
        (
            Week(),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 8),
                datetime(2008, 1, 4): datetime(2008, 1, 11),
                datetime(2008, 1, 5): datetime(2008, 1, 12),
                datetime(2008, 1, 6): datetime(2008, 1, 13),
                datetime(2008, 1, 7): datetime(2008, 1, 14),
            },
        )
    )

    # Mon
    offset_cases.append(
        (
            Week(weekday=0),
            {
                datetime(2007, 12, 31): datetime(2008, 1, 7),
                datetime(2008, 1, 4): datetime(2008, 1, 7),
                datetime(2008, 1, 5): datetime(2008, 1, 7),
                datetime(2008, 1, 6): datetime(2008, 1, 7),
                datetime(2008, 1, 7): datetime(2008, 1, 14),
            },
        )
    )

    # n=0 -> roll forward. Mon
    offset_cases.append(
        (
            Week(0, weekday=0),
            {
                datetime(2007, 12, 31): datetime(2007, 12, 31),
                datetime(2008, 1, 4): datetime(2008, 1, 7),
                datetime(2008, 1, 5): datetime(2008, 1, 7),
                datetime(2008, 1, 6): datetime(2008, 1, 7),
                datetime(2008, 1, 7): datetime(2008, 1, 7),
            },
        )
    )

    # n=0 -> roll forward. Mon
    offset_cases.append(
        (
            Week(-2, weekday=1),
            {
                datetime(2010, 4, 6): datetime(2010, 3, 23),
                datetime(2010, 4, 8): datetime(2010, 3, 30),
                datetime(2010, 4, 5): datetime(2010, 3, 23),
            },
        )
    )

    @pytest.mark.parametrize("case", offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in cases.items():
            assert_offset_equal(offset, base, expected)

    @pytest.mark.parametrize("weekday", range(7))
    def test_onOffset(self, weekday):
        offset = Week(weekday=weekday)

        for day in range(1, 8):
            date = datetime(2008, 1, day)

            if day % 7 == weekday:
                expected = True
            else:
                expected = False
        assert_onOffset(offset, date, expected)


class TestWeekOfMonth(Base):
    _offset = WeekOfMonth
    offset1 = _offset()
    offset2 = _offset(2)

    def test_constructor(self):
        with pytest.raises(ValueError, match="^Week"):
            WeekOfMonth(n=1, week=4, weekday=0)

        with pytest.raises(ValueError, match="^Week"):
            WeekOfMonth(n=1, week=-1, weekday=0)

        with pytest.raises(ValueError, match="^Day"):
            WeekOfMonth(n=1, week=0, weekday=-1)

        with pytest.raises(ValueError, match="^Day"):
            WeekOfMonth(n=1, week=0, weekday=-7)

    def test_repr(self):
        assert (
            repr(WeekOfMonth(weekday=1, week=2)) == "<WeekOfMonth: week=2, weekday=1>"
        )

    def test_offset(self):
        date1 = datetime(2011, 1, 4)  # 1st Tuesday of Month
        date2 = datetime(2011, 1, 11)  # 2nd Tuesday of Month
        date3 = datetime(2011, 1, 18)  # 3rd Tuesday of Month
        date4 = datetime(2011, 1, 25)  # 4th Tuesday of Month

        # see for loop for structure
        test_cases = [
            (-2, 2, 1, date1, datetime(2010, 11, 16)),
            (-2, 2, 1, date2, datetime(2010, 11, 16)),
            (-2, 2, 1, date3, datetime(2010, 11, 16)),
            (-2, 2, 1, date4, datetime(2010, 12, 21)),
            (-1, 2, 1, date1, datetime(2010, 12, 21)),
            (-1, 2, 1, date2, datetime(2010, 12, 21)),
            (-1, 2, 1, date3, datetime(2010, 12, 21)),
            (-1, 2, 1, date4, datetime(2011, 1, 18)),
            (0, 0, 1, date1, datetime(2011, 1, 4)),
            (0, 0, 1, date2, datetime(2011, 2, 1)),
            (0, 0, 1, date3, datetime(2011, 2, 1)),
            (0, 0, 1, date4, datetime(2011, 2, 1)),
            (0, 1, 1, date1, datetime(2011, 1, 11)),
            (0, 1, 1, date2, datetime(2011, 1, 11)),
            (0, 1, 1, date3, datetime(2011, 2, 8)),
            (0, 1, 1, date4, datetime(2011, 2, 8)),
            (0, 0, 1, date1, datetime(2011, 1, 4)),
            (0, 1, 1, date2, datetime(2011, 1, 11)),
            (0, 2, 1, date3, datetime(2011, 1, 18)),
            (0, 3, 1, date4, datetime(2011, 1, 25)),
            (1, 0, 0, date1, datetime(2011, 2, 7)),
            (1, 0, 0, date2, datetime(2011, 2, 7)),
            (1, 0, 0, date3, datetime(2011, 2, 7)),
            (1, 0, 0, date4, datetime(2011, 2, 7)),
            (1, 0, 1, date1, datetime(2011, 2, 1)),
            (1, 0, 1, date2, datetime(2011, 2, 1)),
            (1, 0, 1, date3, datetime(2011, 2, 1)),
            (1, 0, 1, date4, datetime(2011, 2, 1)),
            (1, 0, 2, date1, datetime(2011, 1, 5)),
            (1, 0, 2, date2, datetime(2011, 2, 2)),
            (1, 0, 2, date3, datetime(2011, 2, 2)),
            (1, 0, 2, date4, datetime(2011, 2, 2)),
            (1, 2, 1, date1, datetime(2011, 1, 18)),
            (1, 2, 1, date2, datetime(2011, 1, 18)),
            (1, 2, 1, date3, datetime(2011, 2, 15)),
            (1, 2, 1, date4, datetime(2011, 2, 15)),
            (2, 2, 1, date1, datetime(2011, 2, 15)),
            (2, 2, 1, date2, datetime(2011, 2, 15)),
            (2, 2, 1, date3, datetime(2011, 3, 15)),
            (2, 2, 1, date4, datetime(2011, 3, 15)),
        ]

        for n, week, weekday, dt, expected in test_cases:
            offset = WeekOfMonth(n, week=week, weekday=weekday)
            assert_offset_equal(offset, dt, expected)

        # try subtracting
        result = datetime(2011, 2, 1) - WeekOfMonth(week=1, weekday=2)
        assert result == datetime(2011, 1, 12)

        result = datetime(2011, 2, 3) - WeekOfMonth(week=0, weekday=2)
        assert result == datetime(2011, 2, 2)

    on_offset_cases = [
        (0, 0, datetime(2011, 2, 7), True),
        (0, 0, datetime(2011, 2, 6), False),
        (0, 0, datetime(2011, 2, 14), False),
        (1, 0, datetime(2011, 2, 14), True),
        (0, 1, datetime(2011, 2, 1), True),
        (0, 1, datetime(2011, 2, 8), False),
    ]

    @pytest.mark.parametrize("case", on_offset_cases)
    def test_onOffset(self, case):
        week, weekday, dt, expected = case
        offset = WeekOfMonth(week=week, weekday=weekday)
        assert offset.onOffset(dt) == expected


class TestLastWeekOfMonth(Base):
    _offset = LastWeekOfMonth
    offset1 = _offset()
    offset2 = _offset(2)

    def test_constructor(self):
        with pytest.raises(ValueError, match="^N cannot be 0"):
            LastWeekOfMonth(n=0, weekday=1)

        with pytest.raises(ValueError, match="^Day"):
            LastWeekOfMonth(n=1, weekday=-1)

        with pytest.raises(ValueError, match="^Day"):
            LastWeekOfMonth(n=1, weekday=7)

    def test_offset(self):
        # Saturday
        last_sat = datetime(2013, 8, 31)
        next_sat = datetime(2013, 9, 28)
        offset_sat = LastWeekOfMonth(n=1, weekday=5)

        one_day_before = last_sat + timedelta(days=-1)
        assert one_day_before + offset_sat == last_sat

        one_day_after = last_sat + timedelta(days=+1)
        assert one_day_after + offset_sat == next_sat

        # Test On that day
        assert last_sat + offset_sat == next_sat

        # Thursday

        offset_thur = LastWeekOfMonth(n=1, weekday=3)
        last_thurs = datetime(2013, 1, 31)
        next_thurs = datetime(2013, 2, 28)

        one_day_before = last_thurs + timedelta(days=-1)
        assert one_day_before + offset_thur == last_thurs

        one_day_after = last_thurs + timedelta(days=+1)
        assert one_day_after + offset_thur == next_thurs

        # Test on that day
        assert last_thurs + offset_thur == next_thurs

        three_before = last_thurs + timedelta(days=-3)
        assert three_before + offset_thur == last_thurs

        two_after = last_thurs + timedelta(days=+2)
        assert two_after + offset_thur == next_thurs

        offset_sunday = LastWeekOfMonth(n=1, weekday=WeekDay.SUN)
        assert datetime(2013, 7, 31) + offset_sunday == datetime(2013, 8, 25)

    on_offset_cases = [
        (WeekDay.SUN, datetime(2013, 1, 27), True),
        (WeekDay.SAT, datetime(2013, 3, 30), True),
        (WeekDay.MON, datetime(2013, 2, 18), False),  # Not the last Mon
        (WeekDay.SUN, datetime(2013, 2, 25), False),  # Not a SUN
        (WeekDay.MON, datetime(2013, 2, 25), True),
        (WeekDay.SAT, datetime(2013, 11, 30), True),
        (WeekDay.SAT, datetime(2006, 8, 26), True),
        (WeekDay.SAT, datetime(2007, 8, 25), True),
        (WeekDay.SAT, datetime(2008, 8, 30), True),
        (WeekDay.SAT, datetime(2009, 8, 29), True),
        (WeekDay.SAT, datetime(2010, 8, 28), True),
        (WeekDay.SAT, datetime(2011, 8, 27), True),
        (WeekDay.SAT, datetime(2019, 8, 31), True),
    ]

    @pytest.mark.parametrize("case", on_offset_cases)
    def test_onOffset(self, case):
        weekday, dt, expected = case
        offset = LastWeekOfMonth(weekday=weekday)
        assert offset.onOffset(dt) == expected


class TestSemiMonthEnd(Base):
    _offset = SemiMonthEnd
    offset1 = _offset()
    offset2 = _offset(2)

    def test_offset_whole_year(self):
        dates = (
            datetime(2007, 12, 31),
            datetime(2008, 1, 15),
            datetime(2008, 1, 31),
            datetime(2008, 2, 15),
            datetime(2008, 2, 29),
            datetime(2008, 3, 15),
            datetime(2008, 3, 31),
            datetime(2008, 4, 15),
            datetime(2008, 4, 30),
            datetime(2008, 5, 15),
            datetime(2008, 5, 31),
            datetime(2008, 6, 15),
            datetime(2008, 6, 30),
            datetime(2008, 7, 15),
            datetime(2008, 7, 31),
            datetime(2008, 8, 15),
            datetime(2008, 8, 31),
            datetime(2008, 9, 15),
            datetime(2008, 9, 30),
            datetime(2008, 10, 15),
            datetime(2008, 10, 31),
            datetime(2008, 11, 15),
            datetime(2008, 11, 30),
            datetime(2008, 12, 15),
            datetime(2008, 12, 31),
        )

        for base, exp_date in zip(dates[:-1], dates[1:]):
            assert_offset_equal(SemiMonthEnd(), base, exp_date)

        # ensure .apply_index works as expected
        s = DatetimeIndex(dates[:-1])
        with tm.assert_produces_warning(None):
            # GH#22535 check that we don't get a FutureWarning from adding
            # an integer array to PeriodIndex
            result = SemiMonthEnd().apply_index(s)

        exp = DatetimeIndex(dates[1:])
        tm.assert_index_equal(result, exp)

        # ensure generating a range with DatetimeIndex gives same result
        result = date_range(start=dates[0], end=dates[-1], freq="SM")
        exp = DatetimeIndex(dates)
        tm.assert_index_equal(result, exp)

    offset_cases = []
    offset_cases.append(
        (
            SemiMonthEnd(),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 15),
                datetime(2008, 1, 15): datetime(2008, 1, 31),
                datetime(2008, 1, 31): datetime(2008, 2, 15),
                datetime(2006, 12, 14): datetime(2006, 12, 15),
                datetime(2006, 12, 29): datetime(2006, 12, 31),
                datetime(2006, 12, 31): datetime(2007, 1, 15),
                datetime(2007, 1, 1): datetime(2007, 1, 15),
                datetime(2006, 12, 1): datetime(2006, 12, 15),
                datetime(2006, 12, 15): datetime(2006, 12, 31),
            },
        )
    )

    offset_cases.append(
        (
            SemiMonthEnd(day_of_month=20),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 20),
                datetime(2008, 1, 15): datetime(2008, 1, 20),
                datetime(2008, 1, 21): datetime(2008, 1, 31),
                datetime(2008, 1, 31): datetime(2008, 2, 20),
                datetime(2006, 12, 14): datetime(2006, 12, 20),
                datetime(2006, 12, 29): datetime(2006, 12, 31),
                datetime(2006, 12, 31): datetime(2007, 1, 20),
                datetime(2007, 1, 1): datetime(2007, 1, 20),
                datetime(2006, 12, 1): datetime(2006, 12, 20),
                datetime(2006, 12, 15): datetime(2006, 12, 20),
            },
        )
    )

    offset_cases.append(
        (
            SemiMonthEnd(0),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 15),
                datetime(2008, 1, 16): datetime(2008, 1, 31),
                datetime(2008, 1, 15): datetime(2008, 1, 15),
                datetime(2008, 1, 31): datetime(2008, 1, 31),
                datetime(2006, 12, 29): datetime(2006, 12, 31),
                datetime(2006, 12, 31): datetime(2006, 12, 31),
                datetime(2007, 1, 1): datetime(2007, 1, 15),
            },
        )
    )

    offset_cases.append(
        (
            SemiMonthEnd(0, day_of_month=16),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 16),
                datetime(2008, 1, 16): datetime(2008, 1, 16),
                datetime(2008, 1, 15): datetime(2008, 1, 16),
                datetime(2008, 1, 31): datetime(2008, 1, 31),
                datetime(2006, 12, 29): datetime(2006, 12, 31),
                datetime(2006, 12, 31): datetime(2006, 12, 31),
                datetime(2007, 1, 1): datetime(2007, 1, 16),
            },
        )
    )

    offset_cases.append(
        (
            SemiMonthEnd(2),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 31),
                datetime(2008, 1, 31): datetime(2008, 2, 29),
                datetime(2006, 12, 29): datetime(2007, 1, 15),
                datetime(2006, 12, 31): datetime(2007, 1, 31),
                datetime(2007, 1, 1): datetime(2007, 1, 31),
                datetime(2007, 1, 16): datetime(2007, 2, 15),
                datetime(2006, 11, 1): datetime(2006, 11, 30),
            },
        )
    )

    offset_cases.append(
        (
            SemiMonthEnd(-1),
            {
                datetime(2007, 1, 1): datetime(2006, 12, 31),
                datetime(2008, 6, 30): datetime(2008, 6, 15),
                datetime(2008, 12, 31): datetime(2008, 12, 15),
                datetime(2006, 12, 29): datetime(2006, 12, 15),
                datetime(2006, 12, 30): datetime(2006, 12, 15),
                datetime(2007, 1, 1): datetime(2006, 12, 31),
            },
        )
    )

    offset_cases.append(
        (
            SemiMonthEnd(-1, day_of_month=4),
            {
                datetime(2007, 1, 1): datetime(2006, 12, 31),
                datetime(2007, 1, 4): datetime(2006, 12, 31),
                datetime(2008, 6, 30): datetime(2008, 6, 4),
                datetime(2008, 12, 31): datetime(2008, 12, 4),
                datetime(2006, 12, 5): datetime(2006, 12, 4),
                datetime(2006, 12, 30): datetime(2006, 12, 4),
                datetime(2007, 1, 1): datetime(2006, 12, 31),
            },
        )
    )

    offset_cases.append(
        (
            SemiMonthEnd(-2),
            {
                datetime(2007, 1, 1): datetime(2006, 12, 15),
                datetime(2008, 6, 30): datetime(2008, 5, 31),
                datetime(2008, 3, 15): datetime(2008, 2, 15),
                datetime(2008, 12, 31): datetime(2008, 11, 30),
                datetime(2006, 12, 29): datetime(2006, 11, 30),
                datetime(2006, 12, 14): datetime(2006, 11, 15),
                datetime(2007, 1, 1): datetime(2006, 12, 15),
            },
        )
    )

    @pytest.mark.parametrize("case", offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in cases.items():
            assert_offset_equal(offset, base, expected)

    @pytest.mark.parametrize("case", offset_cases)
    def test_apply_index(self, case):
        offset, cases = case
        s = DatetimeIndex(cases.keys())
        with tm.assert_produces_warning(None):
            # GH#22535 check that we don't get a FutureWarning from adding
            # an integer array to PeriodIndex
            result = offset.apply_index(s)

        exp = DatetimeIndex(cases.values())
        tm.assert_index_equal(result, exp)

    on_offset_cases = [
        (datetime(2007, 12, 31), True),
        (datetime(2007, 12, 15), True),
        (datetime(2007, 12, 14), False),
        (datetime(2007, 12, 1), False),
        (datetime(2008, 2, 29), True),
    ]

    @pytest.mark.parametrize("case", on_offset_cases)
    def test_onOffset(self, case):
        dt, expected = case
        assert_onOffset(SemiMonthEnd(), dt, expected)

    @pytest.mark.parametrize("klass", [Series, DatetimeIndex])
    def test_vectorized_offset_addition(self, klass):
        s = klass(
            [
                Timestamp("2000-01-15 00:15:00", tz="US/Central"),
                Timestamp("2000-02-15", tz="US/Central"),
            ],
            name="a",
        )

        with tm.assert_produces_warning(None):
            # GH#22535 check that we don't get a FutureWarning from adding
            # an integer array to PeriodIndex
            result = s + SemiMonthEnd()
            result2 = SemiMonthEnd() + s

        exp = klass(
            [
                Timestamp("2000-01-31 00:15:00", tz="US/Central"),
                Timestamp("2000-02-29", tz="US/Central"),
            ],
            name="a",
        )
        tm.assert_equal(result, exp)
        tm.assert_equal(result2, exp)

        s = klass(
            [
                Timestamp("2000-01-01 00:15:00", tz="US/Central"),
                Timestamp("2000-02-01", tz="US/Central"),
            ],
            name="a",
        )

        with tm.assert_produces_warning(None):
            # GH#22535 check that we don't get a FutureWarning from adding
            # an integer array to PeriodIndex
            result = s + SemiMonthEnd()
            result2 = SemiMonthEnd() + s

        exp = klass(
            [
                Timestamp("2000-01-15 00:15:00", tz="US/Central"),
                Timestamp("2000-02-15", tz="US/Central"),
            ],
            name="a",
        )
        tm.assert_equal(result, exp)
        tm.assert_equal(result2, exp)


class TestSemiMonthBegin(Base):
    _offset = SemiMonthBegin
    offset1 = _offset()
    offset2 = _offset(2)

    def test_offset_whole_year(self):
        dates = (
            datetime(2007, 12, 15),
            datetime(2008, 1, 1),
            datetime(2008, 1, 15),
            datetime(2008, 2, 1),
            datetime(2008, 2, 15),
            datetime(2008, 3, 1),
            datetime(2008, 3, 15),
            datetime(2008, 4, 1),
            datetime(2008, 4, 15),
            datetime(2008, 5, 1),
            datetime(2008, 5, 15),
            datetime(2008, 6, 1),
            datetime(2008, 6, 15),
            datetime(2008, 7, 1),
            datetime(2008, 7, 15),
            datetime(2008, 8, 1),
            datetime(2008, 8, 15),
            datetime(2008, 9, 1),
            datetime(2008, 9, 15),
            datetime(2008, 10, 1),
            datetime(2008, 10, 15),
            datetime(2008, 11, 1),
            datetime(2008, 11, 15),
            datetime(2008, 12, 1),
            datetime(2008, 12, 15),
        )

        for base, exp_date in zip(dates[:-1], dates[1:]):
            assert_offset_equal(SemiMonthBegin(), base, exp_date)

        # ensure .apply_index works as expected
        s = DatetimeIndex(dates[:-1])
        with tm.assert_produces_warning(None):
            # GH#22535 check that we don't get a FutureWarning from adding
            # an integer array to PeriodIndex
            result = SemiMonthBegin().apply_index(s)

        exp = DatetimeIndex(dates[1:])
        tm.assert_index_equal(result, exp)

        # ensure generating a range with DatetimeIndex gives same result
        result = date_range(start=dates[0], end=dates[-1], freq="SMS")
        exp = DatetimeIndex(dates)
        tm.assert_index_equal(result, exp)

    offset_cases = []
    offset_cases.append(
        (
            SemiMonthBegin(),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 15),
                datetime(2008, 1, 15): datetime(2008, 2, 1),
                datetime(2008, 1, 31): datetime(2008, 2, 1),
                datetime(2006, 12, 14): datetime(2006, 12, 15),
                datetime(2006, 12, 29): datetime(2007, 1, 1),
                datetime(2006, 12, 31): datetime(2007, 1, 1),
                datetime(2007, 1, 1): datetime(2007, 1, 15),
                datetime(2006, 12, 1): datetime(2006, 12, 15),
                datetime(2006, 12, 15): datetime(2007, 1, 1),
            },
        )
    )

    offset_cases.append(
        (
            SemiMonthBegin(day_of_month=20),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 20),
                datetime(2008, 1, 15): datetime(2008, 1, 20),
                datetime(2008, 1, 21): datetime(2008, 2, 1),
                datetime(2008, 1, 31): datetime(2008, 2, 1),
                datetime(2006, 12, 14): datetime(2006, 12, 20),
                datetime(2006, 12, 29): datetime(2007, 1, 1),
                datetime(2006, 12, 31): datetime(2007, 1, 1),
                datetime(2007, 1, 1): datetime(2007, 1, 20),
                datetime(2006, 12, 1): datetime(2006, 12, 20),
                datetime(2006, 12, 15): datetime(2006, 12, 20),
            },
        )
    )

    offset_cases.append(
        (
            SemiMonthBegin(0),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 1),
                datetime(2008, 1, 16): datetime(2008, 2, 1),
                datetime(2008, 1, 15): datetime(2008, 1, 15),
                datetime(2008, 1, 31): datetime(2008, 2, 1),
                datetime(2006, 12, 29): datetime(2007, 1, 1),
                datetime(2006, 12, 2): datetime(2006, 12, 15),
                datetime(2007, 1, 1): datetime(2007, 1, 1),
            },
        )
    )

    offset_cases.append(
        (
            SemiMonthBegin(0, day_of_month=16),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 1),
                datetime(2008, 1, 16): datetime(2008, 1, 16),
                datetime(2008, 1, 15): datetime(2008, 1, 16),
                datetime(2008, 1, 31): datetime(2008, 2, 1),
                datetime(2006, 12, 29): datetime(2007, 1, 1),
                datetime(2006, 12, 31): datetime(2007, 1, 1),
                datetime(2007, 1, 5): datetime(2007, 1, 16),
                datetime(2007, 1, 1): datetime(2007, 1, 1),
            },
        )
    )

    offset_cases.append(
        (
            SemiMonthBegin(2),
            {
                datetime(2008, 1, 1): datetime(2008, 2, 1),
                datetime(2008, 1, 31): datetime(2008, 2, 15),
                datetime(2006, 12, 1): datetime(2007, 1, 1),
                datetime(2006, 12, 29): datetime(2007, 1, 15),
                datetime(2006, 12, 15): datetime(2007, 1, 15),
                datetime(2007, 1, 1): datetime(2007, 2, 1),
                datetime(2007, 1, 16): datetime(2007, 2, 15),
                datetime(2006, 11, 1): datetime(2006, 12, 1),
            },
        )
    )

    offset_cases.append(
        (
            SemiMonthBegin(-1),
            {
                datetime(2007, 1, 1): datetime(2006, 12, 15),
                datetime(2008, 6, 30): datetime(2008, 6, 15),
                datetime(2008, 6, 14): datetime(2008, 6, 1),
                datetime(2008, 12, 31): datetime(2008, 12, 15),
                datetime(2006, 12, 29): datetime(2006, 12, 15),
                datetime(2006, 12, 15): datetime(2006, 12, 1),
                datetime(2007, 1, 1): datetime(2006, 12, 15),
            },
        )
    )

    offset_cases.append(
        (
            SemiMonthBegin(-1, day_of_month=4),
            {
                datetime(2007, 1, 1): datetime(2006, 12, 4),
                datetime(2007, 1, 4): datetime(2007, 1, 1),
                datetime(2008, 6, 30): datetime(2008, 6, 4),
                datetime(2008, 12, 31): datetime(2008, 12, 4),
                datetime(2006, 12, 5): datetime(2006, 12, 4),
                datetime(2006, 12, 30): datetime(2006, 12, 4),
                datetime(2006, 12, 2): datetime(2006, 12, 1),
                datetime(2007, 1, 1): datetime(2006, 12, 4),
            },
        )
    )

    offset_cases.append(
        (
            SemiMonthBegin(-2),
            {
                datetime(2007, 1, 1): datetime(2006, 12, 1),
                datetime(2008, 6, 30): datetime(2008, 6, 1),
                datetime(2008, 6, 14): datetime(2008, 5, 15),
                datetime(2008, 12, 31): datetime(2008, 12, 1),
                datetime(2006, 12, 29): datetime(2006, 12, 1),
                datetime(2006, 12, 15): datetime(2006, 11, 15),
                datetime(2007, 1, 1): datetime(2006, 12, 1),
            },
        )
    )

    @pytest.mark.parametrize("case", offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in cases.items():
            assert_offset_equal(offset, base, expected)

    @pytest.mark.parametrize("case", offset_cases)
    def test_apply_index(self, case):
        offset, cases = case
        s = DatetimeIndex(cases.keys())

        with tm.assert_produces_warning(None):
            # GH#22535 check that we don't get a FutureWarning from adding
            # an integer array to PeriodIndex
            result = offset.apply_index(s)

        exp = DatetimeIndex(cases.values())
        tm.assert_index_equal(result, exp)

    on_offset_cases = [
        (datetime(2007, 12, 1), True),
        (datetime(2007, 12, 15), True),
        (datetime(2007, 12, 14), False),
        (datetime(2007, 12, 31), False),
        (datetime(2008, 2, 15), True),
    ]

    @pytest.mark.parametrize("case", on_offset_cases)
    def test_onOffset(self, case):
        dt, expected = case
        assert_onOffset(SemiMonthBegin(), dt, expected)

    @pytest.mark.parametrize("klass", [Series, DatetimeIndex])
    def test_vectorized_offset_addition(self, klass):
        s = klass(
            [
                Timestamp("2000-01-15 00:15:00", tz="US/Central"),
                Timestamp("2000-02-15", tz="US/Central"),
            ],
            name="a",
        )
        with tm.assert_produces_warning(None):
            # GH#22535 check that we don't get a FutureWarning from adding
            # an integer array to PeriodIndex
            result = s + SemiMonthBegin()
            result2 = SemiMonthBegin() + s

        exp = klass(
            [
                Timestamp("2000-02-01 00:15:00", tz="US/Central"),
                Timestamp("2000-03-01", tz="US/Central"),
            ],
            name="a",
        )
        tm.assert_equal(result, exp)
        tm.assert_equal(result2, exp)

        s = klass(
            [
                Timestamp("2000-01-01 00:15:00", tz="US/Central"),
                Timestamp("2000-02-01", tz="US/Central"),
            ],
            name="a",
        )
        with tm.assert_produces_warning(None):
            # GH#22535 check that we don't get a FutureWarning from adding
            # an integer array to PeriodIndex
            result = s + SemiMonthBegin()
            result2 = SemiMonthBegin() + s

        exp = klass(
            [
                Timestamp("2000-01-15 00:15:00", tz="US/Central"),
                Timestamp("2000-02-15", tz="US/Central"),
            ],
            name="a",
        )
        tm.assert_equal(result, exp)
        tm.assert_equal(result2, exp)


def test_Easter():
    assert_offset_equal(Easter(), datetime(2010, 1, 1), datetime(2010, 4, 4))
    assert_offset_equal(Easter(), datetime(2010, 4, 5), datetime(2011, 4, 24))
    assert_offset_equal(Easter(2), datetime(2010, 1, 1), datetime(2011, 4, 24))

    assert_offset_equal(Easter(), datetime(2010, 4, 4), datetime(2011, 4, 24))
    assert_offset_equal(Easter(2), datetime(2010, 4, 4), datetime(2012, 4, 8))

    assert_offset_equal(-Easter(), datetime(2011, 1, 1), datetime(2010, 4, 4))
    assert_offset_equal(-Easter(), datetime(2010, 4, 5), datetime(2010, 4, 4))
    assert_offset_equal(-Easter(2), datetime(2011, 1, 1), datetime(2009, 4, 12))

    assert_offset_equal(-Easter(), datetime(2010, 4, 4), datetime(2009, 4, 12))
    assert_offset_equal(-Easter(2), datetime(2010, 4, 4), datetime(2008, 3, 23))


class TestOffsetNames:
    def test_get_offset_name(self):
        assert Week(weekday=0).freqstr == "W-MON"
        assert Week(weekday=1).freqstr == "W-TUE"
        assert Week(weekday=2).freqstr == "W-WED"
        assert Week(weekday=3).freqstr == "W-THU"
        assert Week(weekday=4).freqstr == "W-FRI"

        assert LastWeekOfMonth(weekday=WeekDay.SUN).freqstr == "LWOM-SUN"


def test_get_offset():
    with pytest.raises(ValueError, match=INVALID_FREQ_ERR_MSG):
        get_offset("gibberish")
    with pytest.raises(ValueError, match=INVALID_FREQ_ERR_MSG):
        get_offset("QS-JAN-B")

    pairs = [
        ("W-MON", Week(weekday=0)),
        ("W-TUE", Week(weekday=1)),
        ("W-WED", Week(weekday=2)),
        ("W-THU", Week(weekday=3)),
        ("W-FRI", Week(weekday=4)),
    ]

    for name, expected in pairs:
        offset = get_offset(name)
        assert (
            offset == expected
        ), "Expected {name!r} to yield {expected!r} (actual: {offset!r})".format(
            name=name, expected=expected, offset=offset
        )


def test_get_offset_legacy():
    pairs = [("w@Sat", Week(weekday=5))]
    for name, expected in pairs:
        with pytest.raises(ValueError, match=INVALID_FREQ_ERR_MSG):
            get_offset(name)


class TestOffsetAliases:
    def setup_method(self, method):
        _offset_map.clear()

    def test_alias_equality(self):
        for k, v in _offset_map.items():
            if v is None:
                continue
            assert k == v.copy()

    def test_rule_code(self):
        lst = ["M", "MS", "BM", "BMS", "D", "B", "H", "T", "S", "L", "U"]
        for k in lst:
            assert k == get_offset(k).rule_code
            # should be cached - this is kind of an internals test...
            assert k in _offset_map
            assert k == (get_offset(k) * 3).rule_code

        suffix_lst = ["MON", "TUE", "WED", "THU", "FRI", "SAT", "SUN"]
        base = "W"
        for v in suffix_lst:
            alias = "-".join([base, v])
            assert alias == get_offset(alias).rule_code
            assert alias == (get_offset(alias) * 5).rule_code

        suffix_lst = [
            "JAN",
            "FEB",
            "MAR",
            "APR",
            "MAY",
            "JUN",
            "JUL",
            "AUG",
            "SEP",
            "OCT",
            "NOV",
            "DEC",
        ]
        base_lst = ["A", "AS", "BA", "BAS", "Q", "QS", "BQ", "BQS"]
        for base in base_lst:
            for v in suffix_lst:
                alias = "-".join([base, v])
                assert alias == get_offset(alias).rule_code
                assert alias == (get_offset(alias) * 5).rule_code

        lst = ["M", "D", "B", "H", "T", "S", "L", "U"]
        for k in lst:
            code, stride = get_freq_code("3" + k)
            assert isinstance(code, int)
            assert stride == 3
            assert k == get_freq_str(code)


def test_dateoffset_misc():
    oset = offsets.DateOffset(months=2, days=4)
    # it works
    oset.freqstr

    assert not offsets.DateOffset(months=2) == 2


class TestReprNames:
    def test_str_for_named_is_name(self):
        # look at all the amazing combinations!
        month_prefixes = ["A", "AS", "BA", "BAS", "Q", "BQ", "BQS", "QS"]
        names = [
            prefix + "-" + month
            for prefix in month_prefixes
            for month in [
                "JAN",
                "FEB",
                "MAR",
                "APR",
                "MAY",
                "JUN",
                "JUL",
                "AUG",
                "SEP",
                "OCT",
                "NOV",
                "DEC",
            ]
        ]
        days = ["MON", "TUE", "WED", "THU", "FRI", "SAT", "SUN"]
        names += ["W-" + day for day in days]
        names += ["WOM-" + week + day for week in ("1", "2", "3", "4") for day in days]
        _offset_map.clear()
        for name in names:
            offset = get_offset(name)
            assert offset.freqstr == name


def get_utc_offset_hours(ts):
    # take a Timestamp and compute total hours of utc offset
    o = ts.utcoffset()
    return (o.days * 24 * 3600 + o.seconds) / 3600.0


class TestDST:
    """
    test DateOffset additions over Daylight Savings Time
    """

    # one microsecond before the DST transition
    ts_pre_fallback = "2013-11-03 01:59:59.999999"
    ts_pre_springfwd = "2013-03-10 01:59:59.999999"

    # test both basic names and dateutil timezones
    timezone_utc_offsets = {
        "US/Eastern": dict(utc_offset_daylight=-4, utc_offset_standard=-5),
        "dateutil/US/Pacific": dict(utc_offset_daylight=-7, utc_offset_standard=-8),
    }
    valid_date_offsets_singular = [
        "weekday",
        "day",
        "hour",
        "minute",
        "second",
        "microsecond",
    ]
    valid_date_offsets_plural = [
        "weeks",
        "days",
        "hours",
        "minutes",
        "seconds",
        "milliseconds",
        "microseconds",
    ]

    def _test_all_offsets(self, n, **kwds):
        valid_offsets = (
            self.valid_date_offsets_plural
            if n > 1
            else self.valid_date_offsets_singular
        )

        for name in valid_offsets:
            self._test_offset(offset_name=name, offset_n=n, **kwds)

    def _test_offset(self, offset_name, offset_n, tstart, expected_utc_offset):
        offset = DateOffset(**{offset_name: offset_n})

        t = tstart + offset
        if expected_utc_offset is not None:
            assert get_utc_offset_hours(t) == expected_utc_offset

        if offset_name == "weeks":
            # dates should match
            assert t.date() == timedelta(days=7 * offset.kwds["weeks"]) + tstart.date()
            # expect the same day of week, hour of day, minute, second, ...
            assert (
                t.dayofweek == tstart.dayofweek
                and t.hour == tstart.hour
                and t.minute == tstart.minute
                and t.second == tstart.second
            )
        elif offset_name == "days":
            # dates should match
            assert timedelta(offset.kwds["days"]) + tstart.date() == t.date()
            # expect the same hour of day, minute, second, ...
            assert (
                t.hour == tstart.hour
                and t.minute == tstart.minute
                and t.second == tstart.second
            )
        elif offset_name in self.valid_date_offsets_singular:
            # expect the singular offset value to match between tstart and t
            datepart_offset = getattr(
                t, offset_name if offset_name != "weekday" else "dayofweek"
            )
            assert datepart_offset == offset.kwds[offset_name]
        else:
            # the offset should be the same as if it was done in UTC
            assert t == (tstart.tz_convert("UTC") + offset).tz_convert("US/Pacific")

    def _make_timestamp(self, string, hrs_offset, tz):
        if hrs_offset >= 0:
            offset_string = "{hrs:02d}00".format(hrs=hrs_offset)
        else:
            offset_string = "-{hrs:02d}00".format(hrs=-1 * hrs_offset)
        return Timestamp(string + offset_string).tz_convert(tz)

    def test_springforward_plural(self):
        # test moving from standard to daylight savings
        for tz, utc_offsets in self.timezone_utc_offsets.items():
            hrs_pre = utc_offsets["utc_offset_standard"]
            hrs_post = utc_offsets["utc_offset_daylight"]
            self._test_all_offsets(
                n=3,
                tstart=self._make_timestamp(self.ts_pre_springfwd, hrs_pre, tz),
                expected_utc_offset=hrs_post,
            )

    def test_fallback_singular(self):
        # in the case of singular offsets, we don't necessarily know which utc
        # offset the new Timestamp will wind up in (the tz for 1 month may be
        # different from 1 second) so we don't specify an expected_utc_offset
        for tz, utc_offsets in self.timezone_utc_offsets.items():
            hrs_pre = utc_offsets["utc_offset_standard"]
            self._test_all_offsets(
                n=1,
                tstart=self._make_timestamp(self.ts_pre_fallback, hrs_pre, tz),
                expected_utc_offset=None,
            )

    def test_springforward_singular(self):
        for tz, utc_offsets in self.timezone_utc_offsets.items():
            hrs_pre = utc_offsets["utc_offset_standard"]
            self._test_all_offsets(
                n=1,
                tstart=self._make_timestamp(self.ts_pre_springfwd, hrs_pre, tz),
                expected_utc_offset=None,
            )

    offset_classes = {
        MonthBegin: ["11/2/2012", "12/1/2012"],
        MonthEnd: ["11/2/2012", "11/30/2012"],
        SemiMonthBegin: ["11/2/2012", "11/15/2012"],
        SemiMonthEnd: ["11/2/2012", "11/15/2012"],
        Week: ["11/2/2012", "11/9/2012"],
        YearBegin: ["11/2/2012", "1/1/2013"],
        YearEnd: ["11/2/2012", "12/31/2012"],
        QuarterBegin: ["11/2/2012", "12/1/2012"],
        QuarterEnd: ["11/2/2012", "12/31/2012"],
        Day: ["11/4/2012", "11/4/2012 23:00"],
    }.items()

    @pytest.mark.parametrize("tup", offset_classes)
    def test_all_offset_classes(self, tup):
        offset, test_values = tup

        first = Timestamp(test_values[0], tz="US/Eastern") + offset()
        second = Timestamp(test_values[1], tz="US/Eastern")
        assert first == second


# ---------------------------------------------------------------------
def test_get_offset_day_error():
    # subclass of _BaseOffset must override _day_opt attribute, or we should
    # get a NotImplementedError

    with pytest.raises(NotImplementedError):
        DateOffset()._get_offset_day(datetime.now())


def test_valid_default_arguments(offset_types):
    # GH#19142 check that the calling the constructors without passing
    # any keyword arguments produce valid offsets
    cls = offset_types
    cls()


@pytest.mark.parametrize("kwd", sorted(list(liboffsets.relativedelta_kwds)))
def test_valid_month_attributes(kwd, month_classes):
    # GH#18226
    cls = month_classes
    # check that we cannot create e.g. MonthEnd(weeks=3)
    with pytest.raises(TypeError):
        cls(**{kwd: 3})


@pytest.mark.parametrize("kwd", sorted(list(liboffsets.relativedelta_kwds)))
def test_valid_relativedelta_kwargs(kwd):
    # Check that all the arguments specified in liboffsets.relativedelta_kwds
    # are in fact valid relativedelta keyword args
    DateOffset(**{kwd: 1})


@pytest.mark.parametrize("kwd", sorted(list(liboffsets.relativedelta_kwds)))
def test_valid_tick_attributes(kwd, tick_classes):
    # GH#18226
    cls = tick_classes
    # check that we cannot create e.g. Hour(weeks=3)
    with pytest.raises(TypeError):
        cls(**{kwd: 3})


def test_validate_n_error():
    with pytest.raises(TypeError):
        DateOffset(n="Doh!")

    with pytest.raises(TypeError):
        MonthBegin(n=timedelta(1))


def test_require_integers(offset_types):
    cls = offset_types
    with pytest.raises(ValueError):
        cls(n=1.5)


def test_tick_normalize_raises(tick_classes):
    # check that trying to create a Tick object with normalize=True raises
    # GH#21427
    cls = tick_classes
    with pytest.raises(ValueError):
        cls(n=3, normalize=True)


def test_weeks_onoffset():
    # GH#18510 Week with weekday = None, normalize = False should always
    # be onOffset
    offset = Week(n=2, weekday=None)
    ts = Timestamp("1862-01-13 09:03:34.873477378+0210", tz="Africa/Lusaka")
    fast = offset.onOffset(ts)
    slow = (ts + offset) - offset == ts
    assert fast == slow

    # negative n
    offset = Week(n=2, weekday=None)
    ts = Timestamp("1856-10-24 16:18:36.556360110-0717", tz="Pacific/Easter")
    fast = offset.onOffset(ts)
    slow = (ts + offset) - offset == ts
    assert fast == slow


def test_weekofmonth_onoffset():
    # GH#18864
    # Make sure that nanoseconds don't trip up onOffset (and with it apply)
    offset = WeekOfMonth(n=2, week=2, weekday=0)
    ts = Timestamp("1916-05-15 01:14:49.583410462+0422", tz="Asia/Qyzylorda")
    fast = offset.onOffset(ts)
    slow = (ts + offset) - offset == ts
    assert fast == slow

    # negative n
    offset = WeekOfMonth(n=-3, week=1, weekday=0)
    ts = Timestamp("1980-12-08 03:38:52.878321185+0500", tz="Asia/Oral")
    fast = offset.onOffset(ts)
    slow = (ts + offset) - offset == ts
    assert fast == slow


def test_last_week_of_month_on_offset():
    # GH#19036, GH#18977 _adjust_dst was incorrect for LastWeekOfMonth
    offset = LastWeekOfMonth(n=4, weekday=6)
    ts = Timestamp("1917-05-27 20:55:27.084284178+0200", tz="Europe/Warsaw")
    slow = (ts + offset) - offset == ts
    fast = offset.onOffset(ts)
    assert fast == slow

    # negative n
    offset = LastWeekOfMonth(n=-4, weekday=5)
    ts = Timestamp("2005-08-27 05:01:42.799392561-0500", tz="America/Rainy_River")
    slow = (ts + offset) - offset == ts
    fast = offset.onOffset(ts)
    assert fast == slow


def test_week_add_invalid():
    # Week with weekday should raise TypeError and _not_ AttributeError
    #  when adding invalid offset
    offset = Week(weekday=1)
    other = Day()
    with pytest.raises(TypeError, match="Cannot add"):
        offset + other
