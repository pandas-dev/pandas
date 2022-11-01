from datetime import datetime

import pandas._testing as tm

from pandas.tseries.holiday import (
    AbstractHolidayCalendar,
    USFederalHolidayCalendar,
    USMartinLutherKingJr,
    USMemorialDay,
)


def test_no_mlk_before_1986():
    # see gh-10278
    class MLKCalendar(AbstractHolidayCalendar):
        rules = [USMartinLutherKingJr]

    holidays = MLKCalendar().holidays(start="1984", end="1988").to_pydatetime().tolist()

    # Testing to make sure holiday is not incorrectly observed before 1986.
    assert holidays == [datetime(1986, 1, 20, 0, 0), datetime(1987, 1, 19, 0, 0)]


def test_memorial_day():
    class MemorialDay(AbstractHolidayCalendar):
        rules = [USMemorialDay]

    holidays = MemorialDay().holidays(start="1971", end="1980").to_pydatetime().tolist()

    # Fixes 5/31 error and checked manually against Wikipedia.
    assert holidays == [
        datetime(1971, 5, 31, 0, 0),
        datetime(1972, 5, 29, 0, 0),
        datetime(1973, 5, 28, 0, 0),
        datetime(1974, 5, 27, 0, 0),
        datetime(1975, 5, 26, 0, 0),
        datetime(1976, 5, 31, 0, 0),
        datetime(1977, 5, 30, 0, 0),
        datetime(1978, 5, 29, 0, 0),
        datetime(1979, 5, 28, 0, 0),
    ]


def test_federal_holiday_inconsistent_returntype():
    # GH 49075 test case
    # Instantiate two calendars to rule out _cache
    cal1 = USFederalHolidayCalendar()
    cal2 = USFederalHolidayCalendar()

    results_2018 = cal1.holidays(start=datetime(2018, 8, 1), end=datetime(2018, 8, 31))
    results_2019 = cal2.holidays(start=datetime(2019, 8, 1), end=datetime(2019, 8, 31))
    expected_res = DatetimeIndex([], dtype='datetime64[ns]', freq=None)
    # Check to make sure the returned types are identical
    # per GH49075 example code 
    tm.assert_index_equal(results_2018, results_2019)
    # Check each against expected results to rule out if 
    # the half-open interval logic has failed for one case
    tm.assert_index_equal(results_2018, expected_res)
    tm.assert_index_equal(results_2019, expected_res)