"""
Assertion helpers and base class for offsets tests
"""
from __future__ import annotations

from datetime import datetime

from pandas._libs.tslibs import Timestamp
from pandas._libs.tslibs.offsets import (
    FY5253,
    BaseOffset,
    DateOffset,
    FY5253Quarter,
    LastWeekOfMonth,
    Week,
    WeekOfMonth,
)


def assert_offset_equal(offset, base, expected):
    actual = offset + base
    actual_swapped = base + offset
    actual_apply = offset._apply(base)
    try:
        assert actual == expected
        assert actual_swapped == expected
        assert actual_apply == expected
    except AssertionError as err:
        raise AssertionError(
            f"\nExpected: {expected}\nActual: {actual}\nFor Offset: {offset})"
            f"\nAt Date: {base}"
        ) from err


def assert_is_on_offset(offset, date, expected):
    actual = offset.is_on_offset(date)
    assert actual == expected, (
        f"\nExpected: {expected}\nActual: {actual}\nFor Offset: {offset})"
        f"\nAt Date: {date}"
    )


class WeekDay:
    MON = 0
    TUE = 1
    WED = 2
    THU = 3
    FRI = 4
    SAT = 5
    SUN = 6


class Base:
    _offset: type[BaseOffset] | None = None
    d = Timestamp(datetime(2008, 1, 2))

    def _get_offset(self, klass, value=1, normalize=False):
        # create instance from offset class
        if klass is FY5253:
            klass = klass(
                n=value,
                startingMonth=1,
                weekday=1,
                variation="last",
                normalize=normalize,
            )
        elif klass is FY5253Quarter:
            klass = klass(
                n=value,
                startingMonth=1,
                weekday=1,
                qtr_with_extra_week=1,
                variation="last",
                normalize=normalize,
            )
        elif klass is LastWeekOfMonth:
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
