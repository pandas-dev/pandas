"""
Tests for DatetimeIndex methods behaving like their Timestamp counterparts
"""

import calendar
from datetime import (
    date,
    datetime,
    time,
)
import locale
import unicodedata

import numpy as np
import pytest

from pandas._libs.tslibs import timezones
from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import DatetimeArray


class TestDatetimeIndexOps:
    def test_dti_no_millisecond_field(self):
        msg = "type object 'DatetimeIndex' has no attribute 'millisecond'"
        with pytest.raises(AttributeError, match=msg):
            pd.DatetimeIndex.millisecond

        msg = "'DatetimeIndex' object has no attribute 'millisecond'"
        with pytest.raises(AttributeError, match=msg):
            pd.DatetimeIndex([]).millisecond

    def test_dti_time(self):
        rng = pd.date_range("1/1/2000", freq="12min", periods=10)
        result = pd.Index(rng).time
        expected = [t.time() for t in rng]
        assert (result == expected).all()

    def test_dti_date(self):
        rng = pd.date_range("1/1/2000", freq="12h", periods=10)
        result = pd.Index(rng).date
        expected = [t.date() for t in rng]
        assert (result == expected).all()

    @pytest.mark.parametrize(
        "dtype",
        [None, "datetime64[ns, CET]", "datetime64[ns, EST]", "datetime64[ns, UTC]"],
    )
    def test_dti_date2(self, dtype):
        # Regression test for GH#21230
        expected = np.array([date(2018, 6, 4), pd.NaT])

        index = pd.DatetimeIndex(["2018-06-04 10:00:00", pd.NaT], dtype=dtype)
        result = index.date

        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype",
        [None, "datetime64[ns, CET]", "datetime64[ns, EST]", "datetime64[ns, UTC]"],
    )
    def test_dti_time2(self, dtype):
        # Regression test for GH#21267
        expected = np.array([time(10, 20, 30), pd.NaT])

        index = pd.DatetimeIndex(["2018-06-04 10:20:30", pd.NaT], dtype=dtype)
        result = index.time

        tm.assert_numpy_array_equal(result, expected)

    def test_dti_timetz(self, tz_naive_fixture):
        # GH#21358
        tz = timezones.maybe_get_tz(tz_naive_fixture)

        expected = np.array([time(10, 20, 30, tzinfo=tz), pd.NaT])

        index = pd.DatetimeIndex(["2018-06-04 10:20:30", pd.NaT], tz=tz)
        result = index.timetz

        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize(
        "field",
        [
            "day_of_week",
            "day_of_year",
            "quarter",
            "days_in_month",
            "is_month_start",
            "is_month_end",
            "is_quarter_start",
            "is_quarter_end",
            "is_year_start",
            "is_year_end",
        ],
    )
    def test_dti_timestamp_fields(self, field):
        # extra fields from DatetimeIndex like quarter and week
        idx = pd.date_range("2020-01-01", periods=10)
        expected = getattr(idx, field)[-1]

        result = getattr(pd.Timestamp(idx[-1]), field)
        assert result == expected

    def test_dti_nanosecond(self):
        dti = pd.DatetimeIndex(np.arange(10))
        expected = pd.Index(np.arange(10, dtype=np.int32))

        tm.assert_index_equal(dti.nanosecond, expected)

    @pytest.mark.parametrize("prefix", ["", "dateutil/"])
    def test_dti_hour_tzaware(self, prefix):
        strdates = ["1/1/2012", "3/1/2012", "4/1/2012"]
        rng = pd.DatetimeIndex(strdates, tz=prefix + "US/Eastern")
        assert (rng.hour == 0).all()

        # a more unusual time zone, GH#1946
        dr = pd.date_range(
            "2011-10-02 00:00", freq="h", periods=10, tz=prefix + "America/Atikokan"
        )

        expected = pd.Index(np.arange(10, dtype=np.int32))
        tm.assert_index_equal(dr.hour, expected)

    # GH#12806
    @pytest.mark.parametrize(
        "time_locale",
        [None, *tm.get_locales()],
    )
    def test_day_name_month_name(self, time_locale):
        # Test Monday -> Sunday and January -> December, in that sequence
        if time_locale is None:
            # If the time_locale is None, day-name and month_name should
            # return the english attributes
            expected_days = [
                "Monday",
                "Tuesday",
                "Wednesday",
                "Thursday",
                "Friday",
                "Saturday",
                "Sunday",
            ]
            expected_months = [
                "January",
                "February",
                "March",
                "April",
                "May",
                "June",
                "July",
                "August",
                "September",
                "October",
                "November",
                "December",
            ]
        else:
            with tm.set_locale(time_locale, locale.LC_TIME):
                expected_days = calendar.day_name[:]
                expected_months = calendar.month_name[1:]

        # GH#11128
        dti = pd.date_range(freq="D", start=datetime(1998, 1, 1), periods=365)
        english_days = [
            "Monday",
            "Tuesday",
            "Wednesday",
            "Thursday",
            "Friday",
            "Saturday",
            "Sunday",
        ]
        for day, name, eng_name in zip(
            range(4, 11), expected_days, english_days, strict=True
        ):
            name = name.capitalize()
            assert dti.day_name(locale=time_locale)[day] == name
            assert dti.day_name(locale=None)[day] == eng_name
            ts = pd.Timestamp(datetime(2016, 4, day))
            assert ts.day_name(locale=time_locale) == name
        dti = dti.append(pd.DatetimeIndex([pd.NaT]))
        assert np.isnan(dti.day_name(locale=time_locale)[-1])
        ts = pd.Timestamp(pd.NaT)
        assert np.isnan(ts.day_name(locale=time_locale))

        # GH#12805
        dti = pd.date_range(freq="ME", start="2012", end="2013")
        result = dti.month_name(locale=time_locale)
        expected = pd.Index([month.capitalize() for month in expected_months])

        # work around different normalization schemes GH#22342
        result = result.str.normalize("NFD")
        expected = expected.str.normalize("NFD")

        tm.assert_index_equal(result, expected)

        for item, expected in zip(dti, expected_months, strict=True):
            result = item.month_name(locale=time_locale)
            expected = expected.capitalize()

            result = unicodedata.normalize("NFD", result)
            expected = unicodedata.normalize("NFD", result)

            assert result == expected
        dti = dti.append(pd.DatetimeIndex([pd.NaT]))
        assert np.isnan(dti.month_name(locale=time_locale)[-1])

    def test_dti_week(self):
        # GH#6538: Check that DatetimeIndex and its TimeStamp elements
        # return the same weekofyear accessor close to new year w/ tz
        dates = ["2013/12/29", "2013/12/30", "2013/12/31"]
        dates = pd.DatetimeIndex(dates, tz="Europe/Brussels")
        expected = [52, 1, 1]
        assert dates.isocalendar().week.tolist() == expected
        assert [d.weekofyear for d in dates] == expected

    @pytest.mark.parametrize("tz", [None, "US/Eastern"])
    def test_dti_fields(self, tz):
        # GH#13303
        dti = pd.date_range(
            freq="D", start=datetime(1998, 1, 1), periods=365, tz=tz, unit="ns"
        )
        assert dti.year[0] == 1998
        assert dti.month[0] == 1
        assert dti.day[0] == 1
        assert dti.hour[0] == 0
        assert dti.minute[0] == 0
        assert dti.second[0] == 0
        assert dti.microsecond[0] == 0
        assert dti.day_of_week[0] == 3

        assert dti.day_of_year[0] == 1
        assert dti.day_of_year[120] == 121

        assert dti.isocalendar().week.iloc[0] == 1
        assert dti.isocalendar().week.iloc[120] == 18

        assert dti.quarter[0] == 1
        assert dti.quarter[120] == 2

        assert dti.days_in_month[0] == 31
        assert dti.days_in_month[90] == 30

        assert dti.is_month_start[0]
        assert not dti.is_month_start[1]
        assert dti.is_month_start[31]
        assert dti.is_quarter_start[0]
        assert dti.is_quarter_start[90]
        assert dti.is_year_start[0]
        assert not dti.is_year_start[364]
        assert not dti.is_month_end[0]
        assert dti.is_month_end[30]
        assert not dti.is_month_end[31]
        assert dti.is_month_end[364]
        assert not dti.is_quarter_end[0]
        assert not dti.is_quarter_end[30]
        assert dti.is_quarter_end[89]
        assert dti.is_quarter_end[364]
        assert not dti.is_year_end[0]
        assert dti.is_year_end[364]

        assert len(dti.year) == 365
        assert len(dti.month) == 365
        assert len(dti.day) == 365
        assert len(dti.hour) == 365
        assert len(dti.minute) == 365
        assert len(dti.second) == 365
        assert len(dti.microsecond) == 365
        assert len(dti.day_of_week) == 365
        assert len(dti.day_of_year) == 365
        assert len(dti.isocalendar()) == 365
        assert len(dti.quarter) == 365
        assert len(dti.is_month_start) == 365
        assert len(dti.is_month_end) == 365
        assert len(dti.is_quarter_start) == 365
        assert len(dti.is_quarter_end) == 365
        assert len(dti.is_year_start) == 365
        assert len(dti.is_year_end) == 365

    @pytest.mark.parametrize(
        "old_attr, new_attr",
        [
            ("dayofweek", "day_of_week"),
            ("dayofyear", "day_of_year"),
            ("daysinmonth", "days_in_month"),
        ],
    )
    def test_deprecated_day_attrs(self, old_attr, new_attr):
        # GH#46768
        dti = pd.date_range(start="1/1/2005", end="12/1/2005", freq="ME")
        msg = f"DatetimeIndex.{old_attr} is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            old_val = getattr(dti, old_attr)
        tm.assert_index_equal(old_val, getattr(dti, new_attr))

    def test_field_ops(self):
        dti = pd.date_range(freq="D", start="1/1/1998", periods=365)
        dti.name = "name"

        # non boolean accessors -> return Index
        for accessor in DatetimeArray._field_ops:
            if accessor == "weekday":
                # GH#12816 weekday is deprecated
                with tm.assert_produces_warning(Pandas4Warning, match="weekday"):
                    res = getattr(dti, accessor)
            else:
                res = getattr(dti, accessor)
            assert len(res) == 365
            assert isinstance(res, pd.Index)
            assert res.name == "name"

        # boolean accessors -> return array
        for accessor in DatetimeArray._bool_ops:
            res = getattr(dti, accessor)
            assert len(res) == 365
            assert isinstance(res, np.ndarray)

        # test boolean indexing
        res = dti[dti.is_quarter_start]
        exp = dti[[0, 90, 181, 273]]
        tm.assert_index_equal(res, exp)
        res = dti[dti.is_leap_year]
        exp = pd.DatetimeIndex([], freq="D", tz=dti.tz, name="name").as_unit(dti.unit)
        tm.assert_index_equal(res, exp)

    def test_dti_is_year_quarter_start(self):
        dti = pd.date_range(freq="BQE-FEB", start=datetime(1998, 1, 1), periods=4)

        assert sum(dti.is_quarter_start) == 0
        assert sum(dti.is_quarter_end) == 4
        assert sum(dti.is_year_start) == 0
        assert sum(dti.is_year_end) == 1

    def test_dti_is_month_start(self):
        dti = pd.DatetimeIndex(["2000-01-01", "2000-01-02", "2000-01-03"])

        assert dti.is_month_start[0] == 1

    def test_dti_is_month_start_custom(self):
        # Ensure is_start/end accessors throw ValueError for CustomBusinessDay,
        bday_egypt = pd.offsets.CustomBusinessDay(weekmask="Sun Mon Tue Wed Thu")
        dti = pd.date_range(datetime(2013, 4, 30), periods=5, freq=bday_egypt)
        msg = "Custom business days is not supported by is_month_start"
        with pytest.raises(ValueError, match=msg):
            dti.is_month_start

    @pytest.mark.parametrize(
        "timestamp, freq, periods, expected_values",
        [
            ("2017-12-01", "MS", 3, np.array([False, True, False])),
            ("2017-12-01", "QS", 3, np.array([True, False, False])),
            ("2017-12-01", "YS", 3, np.array([True, True, True])),
        ],
    )
    def test_dti_dr_is_year_start(self, timestamp, freq, periods, expected_values):
        # GH57377
        result = pd.date_range(timestamp, freq=freq, periods=periods).is_year_start
        tm.assert_numpy_array_equal(result, expected_values)

    @pytest.mark.parametrize(
        "timestamp, freq, periods, expected_values",
        [
            ("2017-12-01", "ME", 3, np.array([True, False, False])),
            ("2017-12-01", "QE", 3, np.array([True, False, False])),
            ("2017-12-01", "YE", 3, np.array([True, True, True])),
        ],
    )
    def test_dti_dr_is_year_end(self, timestamp, freq, periods, expected_values):
        # GH57377
        result = pd.date_range(timestamp, freq=freq, periods=periods).is_year_end
        tm.assert_numpy_array_equal(result, expected_values)

    @pytest.mark.parametrize(
        "timestamp, freq, periods, expected_values",
        [
            ("2017-12-01", "MS", 3, np.array([False, True, False])),
            ("2017-12-01", "QS", 3, np.array([True, True, True])),
            ("2017-12-01", "YS", 3, np.array([True, True, True])),
        ],
    )
    def test_dti_dr_is_quarter_start(self, timestamp, freq, periods, expected_values):
        # GH57377
        result = pd.date_range(timestamp, freq=freq, periods=periods).is_quarter_start
        tm.assert_numpy_array_equal(result, expected_values)

    @pytest.mark.parametrize(
        "timestamp, freq, periods, expected_values",
        [
            ("2017-12-01", "ME", 3, np.array([True, False, False])),
            ("2017-12-01", "QE", 3, np.array([True, True, True])),
            ("2017-12-01", "YE", 3, np.array([True, True, True])),
        ],
    )
    def test_dti_dr_is_quarter_end(self, timestamp, freq, periods, expected_values):
        # GH57377
        result = pd.date_range(timestamp, freq=freq, periods=periods).is_quarter_end
        tm.assert_numpy_array_equal(result, expected_values)

    @pytest.mark.parametrize(
        "timestamp, freq, periods, expected_values",
        [
            ("2017-12-01", "MS", 3, np.array([True, True, True])),
            ("2017-12-01", "QS", 3, np.array([True, True, True])),
            ("2017-12-01", "YS", 3, np.array([True, True, True])),
        ],
    )
    def test_dti_dr_is_month_start(self, timestamp, freq, periods, expected_values):
        # GH57377
        result = pd.date_range(timestamp, freq=freq, periods=periods).is_month_start
        tm.assert_numpy_array_equal(result, expected_values)

    @pytest.mark.parametrize(
        "timestamp, freq, periods, expected_values",
        [
            ("2017-12-01", "ME", 3, np.array([True, True, True])),
            ("2017-12-01", "QE", 3, np.array([True, True, True])),
            ("2017-12-01", "YE", 3, np.array([True, True, True])),
        ],
    )
    def test_dti_dr_is_month_end(self, timestamp, freq, periods, expected_values):
        # GH57377
        result = pd.date_range(timestamp, freq=freq, periods=periods).is_month_end
        tm.assert_numpy_array_equal(result, expected_values)

    def test_dti_is_year_quarter_start_doubledigit_freq(self):
        # GH#58523
        dr = pd.date_range("2017-01-01", periods=2, freq="10YS")
        assert all(dr.is_year_start)

        dr = pd.date_range("2017-01-01", periods=2, freq="10QS")
        assert all(dr.is_quarter_start)

    def test_dti_is_year_start_freq_custom_business_day_with_digit(self):
        # GH#58664
        dr = pd.date_range("2020-01-01", periods=2, freq="2C")
        msg = "Custom business days is not supported by is_year_start"
        with pytest.raises(ValueError, match=msg):
            dr.is_year_start

    @pytest.mark.parametrize("freq", ["3BMS", pd.offsets.BusinessMonthBegin(3)])
    def test_dti_is_year_quarter_start_freq_business_month_begin(self, freq):
        # GH#58729
        dr = pd.date_range("2020-01-01", periods=5, freq=freq)
        result = [x.is_year_start for x in dr]
        assert result == [True, False, False, False, True]

        dr = pd.date_range("2020-01-01", periods=4, freq=freq)
        result = [x.is_quarter_start for x in dr]
        assert all(dr.is_quarter_start)

    def test_dti_is_year_start_freq_two_business_days(self):
        # GH#58524
        dr = pd.date_range("2017-01-01", periods=2, freq="2B")
        result = dr.is_year_start
        expected = np.array([True, False])
        tm.assert_numpy_array_equal(result, expected)


@pytest.mark.parametrize("freq", ["MS", "QS", "YS"])
@pytest.mark.parametrize("n", [1, 2, 10])
@pytest.mark.parametrize(
    "dt",
    [
        datetime(1960, 1, 1),
        datetime(1970, 1, 1),
        datetime(1979, 12, 31),
    ],
)
def test_against_scalar_parametric(freq, dt, n):
    # https://github.com/pandas-dev/pandas/issues/49606
    freq = f"{n}{freq}"
    d = pd.date_range(dt, periods=3, freq=freq)
    result = list(d.is_year_start)
    expected = [x.is_year_start for x in d]
    assert result == expected


@pytest.mark.parametrize(
    "value",
    [
        np.iinfo(np.int64).max,  # year 292277026596
        2**32 * 86400 * 365,
        -(2**31) * 86400 * 400,
        np.iinfo(np.int64).min + 1,
    ],
)
def test_day_of_week_year_beyond_int32(value):
    # GH#66549 the year was passed to ccalendar.dayofweek as a C int, so
    #  second-resolution dates beyond year 2**31 got a truncated year
    dti = pd.DatetimeIndex(np.array([value], dtype="M8[s]"))
    # 1970-01-01 was a Thursday, i.e. day_of_week 3
    expected = (value // 86400 + 3) % 7

    assert dti.day_of_week[0] == expected
    assert dti[0].day_of_week == expected
    assert dti[0].weekday() == expected
    assert pd.Series(dti).dt.day_of_week[0] == expected
    # scalar day_name goes through the same ccalendar helper
    assert dti[0].day_name() == calendar.day_name[expected]


@pytest.mark.parametrize(
    "value, expected",
    [
        (np.iinfo(np.int64).max, 292277026596),
        (2**32 * 86400 * 365, 4292117654),
        (np.iinfo(np.int64).min + 1, -292277022657),
    ],
)
def test_year_beyond_int32_raises(value, expected):
    # GH#66549 the vectorized year is int32, and used to come back truncated,
    #  e.g. 219250468 for year 292277026596
    dti = pd.DatetimeIndex(np.array([value], dtype="M8[s]"))
    assert dti[0].year == expected

    msg = f"year {expected} is out of range for the int32 result"
    with pytest.raises(ValueError, match=msg):
        dti.year
    with pytest.raises(ValueError, match=msg):
        pd.Series(dti).dt.year
    with pytest.raises(ValueError, match=f"ISO year {expected}"):
        dti.isocalendar()


def test_day_of_year_and_is_leap_year_beyond_int32():
    # GH#66549 ccalendar.get_day_of_year took the year as a C int, so both the
    #  scalar and the array path read the leap status off a truncated year.
    #  10143198492 is a leap year, its int32 truncation 1553263900 is not.
    dti = pd.DatetimeIndex(np.array(["10143198492-12-29"], dtype="M8[s]"))

    assert dti[0].day_of_year == 364
    assert dti.day_of_year[0] == 364
    assert pd.Series(dti).dt.day_of_year[0] == 364

    assert dti[0].is_leap_year
    assert dti.is_leap_year[0]
    assert pd.Series(dti).dt.is_leap_year[0]
