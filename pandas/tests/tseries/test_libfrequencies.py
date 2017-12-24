# -*- coding: utf-8 -*-

import pandas.util.testing as tm

from pandas.tseries import offsets
from pandas._libs.tslibs.frequencies import (get_rule_month,
                                             _period_str_to_code,
                                             _INVALID_FREQ_ERROR,
                                             is_superperiod, is_subperiod)


def test_get_rule_month():
    result = get_rule_month('W')
    assert (result == 'DEC')
    result = get_rule_month(offsets.Week())
    assert (result == 'DEC')

    result = get_rule_month('D')
    assert (result == 'DEC')
    result = get_rule_month(offsets.Day())
    assert (result == 'DEC')

    result = get_rule_month('Q')
    assert (result == 'DEC')
    result = get_rule_month(offsets.QuarterEnd(startingMonth=12))
    print(result == 'DEC')

    result = get_rule_month('Q-JAN')
    assert (result == 'JAN')
    result = get_rule_month(offsets.QuarterEnd(startingMonth=1))
    assert (result == 'JAN')

    result = get_rule_month('A-DEC')
    assert (result == 'DEC')
    result = get_rule_month('Y-DEC')
    assert (result == 'DEC')
    result = get_rule_month(offsets.YearEnd())
    assert (result == 'DEC')

    result = get_rule_month('A-MAY')
    assert (result == 'MAY')
    result = get_rule_month('Y-MAY')
    assert (result == 'MAY')
    result = get_rule_month(offsets.YearEnd(month=5))
    assert (result == 'MAY')


def test_period_str_to_code():
    assert (_period_str_to_code('A') == 1000)
    assert (_period_str_to_code('A-DEC') == 1000)
    assert (_period_str_to_code('A-JAN') == 1001)
    assert (_period_str_to_code('Y') == 1000)
    assert (_period_str_to_code('Y-DEC') == 1000)
    assert (_period_str_to_code('Y-JAN') == 1001)

    assert (_period_str_to_code('Q') == 2000)
    assert (_period_str_to_code('Q-DEC') == 2000)
    assert (_period_str_to_code('Q-FEB') == 2002)

    def _assert_depr(freq, expected, aliases):
        assert isinstance(aliases, list)
        assert (_period_str_to_code(freq) == expected)

        msg = _INVALID_FREQ_ERROR
        for alias in aliases:
            with tm.assert_raises_regex(ValueError, msg):
                _period_str_to_code(alias)

    _assert_depr("M", 3000, ["MTH", "MONTH", "MONTHLY"])

    assert (_period_str_to_code('W') == 4000)
    assert (_period_str_to_code('W-SUN') == 4000)
    assert (_period_str_to_code('W-FRI') == 4005)

    _assert_depr("B", 5000, ["BUS", "BUSINESS", "BUSINESSLY", "WEEKDAY"])
    _assert_depr("D", 6000, ["DAY", "DLY", "DAILY"])
    _assert_depr("H", 7000, ["HR", "HOUR", "HRLY", "HOURLY"])

    _assert_depr("T", 8000, ["minute", "MINUTE", "MINUTELY"])
    assert (_period_str_to_code('Min') == 8000)

    _assert_depr("S", 9000, ["sec", "SEC", "SECOND", "SECONDLY"])
    _assert_depr("L", 10000, ["MILLISECOND", "MILLISECONDLY"])
    assert (_period_str_to_code('ms') == 10000)

    _assert_depr("U", 11000, ["MICROSECOND", "MICROSECONDLY"])
    assert (_period_str_to_code('US') == 11000)

    _assert_depr("N", 12000, ["NANOSECOND", "NANOSECONDLY"])
    assert (_period_str_to_code('NS') == 12000)


def test_is_superperiod_subperiod():

    # input validation
    assert not (is_superperiod(offsets.YearEnd(), None))
    assert not (is_subperiod(offsets.MonthEnd(), None))
    assert not (is_superperiod(None, offsets.YearEnd()))
    assert not (is_subperiod(None, offsets.MonthEnd()))
    assert not (is_superperiod(None, None))
    assert not (is_subperiod(None, None))

    assert (is_superperiod(offsets.YearEnd(), offsets.MonthEnd()))
    assert (is_subperiod(offsets.MonthEnd(), offsets.YearEnd()))

    assert (is_superperiod(offsets.Hour(), offsets.Minute()))
    assert (is_subperiod(offsets.Minute(), offsets.Hour()))

    assert (is_superperiod(offsets.Second(), offsets.Milli()))
    assert (is_subperiod(offsets.Milli(), offsets.Second()))

    assert (is_superperiod(offsets.Milli(), offsets.Micro()))
    assert (is_subperiod(offsets.Micro(), offsets.Milli()))

    assert (is_superperiod(offsets.Micro(), offsets.Nano()))
    assert (is_subperiod(offsets.Nano(), offsets.Micro()))
