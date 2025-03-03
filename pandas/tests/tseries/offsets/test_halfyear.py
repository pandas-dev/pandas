"""
Tests for the following offsets:
- HalfYearBegin
- HalfYearEnd
"""

from __future__ import annotations

from datetime import datetime

import pytest

from pandas.tests.tseries.offsets.common import (
    assert_is_on_offset,
    assert_offset_equal,
)

from pandas.tseries.offsets import (
    HalfYearBegin,
    HalfYearEnd,
)


@pytest.mark.parametrize("klass", (HalfYearBegin, HalfYearEnd))
def test_halfyearly_dont_normalize(klass):
    date = datetime(2012, 3, 31, 5, 30)
    result = date + klass()
    assert result.time() == date.time()


@pytest.mark.parametrize("offset", [HalfYearBegin(), HalfYearEnd()])
@pytest.mark.parametrize(
    "date",
    [
        datetime(2016, m, d)
        for m in [7, 8, 9, 10, 11, 12]
        for d in [1, 2, 3, 28, 29, 30, 31]
        if not (m in {9, 11} and d == 31)
    ],
)
def test_on_offset(offset, date):
    res = offset.is_on_offset(date)
    slow_version = date == (date + offset) - offset
    assert res == slow_version


class TestHalfYearBegin:
    def test_repr(self):
        expected = "<HalfYearBegin: startingMonth=1>"
        assert repr(HalfYearBegin()) == expected
        expected = "<HalfYearBegin: startingMonth=3>"
        assert repr(HalfYearBegin(startingMonth=3)) == expected
        expected = "<HalfYearBegin: startingMonth=1>"
        assert repr(HalfYearBegin(startingMonth=1)) == expected

    def test_offset_corner_case(self):
        # corner
        offset = HalfYearBegin(n=-1, startingMonth=1)
        assert datetime(2010, 2, 1) + offset == datetime(2010, 1, 1)

    offset_cases = []
    offset_cases.append(
        (
            HalfYearBegin(startingMonth=1),
            {
                datetime(2007, 12, 1): datetime(2008, 1, 1),
                datetime(2008, 1, 1): datetime(2008, 7, 1),
                datetime(2008, 2, 15): datetime(2008, 7, 1),
                datetime(2008, 2, 29): datetime(2008, 7, 1),
                datetime(2008, 3, 15): datetime(2008, 7, 1),
                datetime(2008, 3, 31): datetime(2008, 7, 1),
                datetime(2008, 4, 15): datetime(2008, 7, 1),
                datetime(2008, 4, 1): datetime(2008, 7, 1),
                datetime(2008, 7, 1): datetime(2009, 1, 1),
                datetime(2008, 7, 15): datetime(2009, 1, 1),
            },
        )
    )

    offset_cases.append(
        (
            HalfYearBegin(startingMonth=2),
            {
                datetime(2008, 1, 1): datetime(2008, 2, 1),
                datetime(2008, 1, 31): datetime(2008, 2, 1),
                datetime(2008, 1, 15): datetime(2008, 2, 1),
                datetime(2008, 2, 29): datetime(2008, 8, 1),
                datetime(2008, 3, 15): datetime(2008, 8, 1),
                datetime(2008, 3, 31): datetime(2008, 8, 1),
                datetime(2008, 4, 15): datetime(2008, 8, 1),
                datetime(2008, 4, 30): datetime(2008, 8, 1),
            },
        )
    )

    offset_cases.append(
        (
            HalfYearBegin(startingMonth=1, n=0),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 1),
                datetime(2008, 12, 1): datetime(2009, 1, 1),
                datetime(2008, 2, 15): datetime(2008, 7, 1),
                datetime(2008, 2, 29): datetime(2008, 7, 1),
                datetime(2008, 3, 15): datetime(2008, 7, 1),
                datetime(2008, 3, 31): datetime(2008, 7, 1),
                datetime(2008, 4, 15): datetime(2008, 7, 1),
                datetime(2008, 4, 30): datetime(2008, 7, 1),
                datetime(2008, 7, 1): datetime(2008, 7, 1),
                datetime(2008, 7, 15): datetime(2009, 1, 1),
            },
        )
    )

    offset_cases.append(
        (
            HalfYearBegin(startingMonth=1, n=-1),
            {
                datetime(2008, 1, 1): datetime(2007, 7, 1),
                datetime(2008, 1, 31): datetime(2008, 1, 1),
                datetime(2008, 2, 15): datetime(2008, 1, 1),
                datetime(2008, 2, 29): datetime(2008, 1, 1),
                datetime(2008, 3, 15): datetime(2008, 1, 1),
                datetime(2008, 3, 31): datetime(2008, 1, 1),
                datetime(2008, 4, 15): datetime(2008, 1, 1),
                datetime(2008, 4, 30): datetime(2008, 1, 1),
                datetime(2008, 7, 1): datetime(2008, 1, 1),
                datetime(2008, 7, 15): datetime(2008, 7, 1),
            },
        )
    )

    offset_cases.append(
        (
            HalfYearBegin(startingMonth=1, n=2),
            {
                datetime(2008, 1, 1): datetime(2009, 1, 1),
                datetime(2008, 2, 15): datetime(2009, 1, 1),
                datetime(2008, 2, 29): datetime(2009, 1, 1),
                datetime(2008, 3, 15): datetime(2009, 1, 1),
                datetime(2008, 3, 31): datetime(2009, 1, 1),
                datetime(2008, 4, 15): datetime(2009, 1, 1),
                datetime(2008, 4, 1): datetime(2009, 1, 1),
                datetime(2008, 7, 15): datetime(2009, 7, 1),
                datetime(2008, 7, 1): datetime(2009, 7, 1),
            },
        )
    )

    @pytest.mark.parametrize("case", offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in cases.items():
            assert_offset_equal(offset, base, expected)

    on_offset_cases = [
        (HalfYearBegin(1, startingMonth=1), datetime(2008, 1, 1), True),
        (HalfYearBegin(1, startingMonth=1), datetime(2007, 12, 1), False),
        (HalfYearBegin(1, startingMonth=1), datetime(2008, 2, 1), False),
        (HalfYearBegin(1, startingMonth=1), datetime(2007, 3, 1), False),
        (HalfYearBegin(1, startingMonth=1), datetime(2008, 4, 1), False),
        (HalfYearBegin(1, startingMonth=1), datetime(2008, 5, 1), False),
        (HalfYearBegin(1, startingMonth=1), datetime(2007, 6, 1), False),
        (HalfYearBegin(1, startingMonth=3), datetime(2008, 1, 1), False),
        (HalfYearBegin(1, startingMonth=3), datetime(2007, 12, 1), False),
        (HalfYearBegin(1, startingMonth=3), datetime(2008, 2, 1), False),
        (HalfYearBegin(1, startingMonth=3), datetime(2007, 3, 1), True),
        (HalfYearBegin(1, startingMonth=3), datetime(2008, 4, 1), False),
        (HalfYearBegin(1, startingMonth=3), datetime(2008, 5, 1), False),
        (HalfYearBegin(1, startingMonth=3), datetime(2008, 5, 2), False),
        (HalfYearBegin(1, startingMonth=3), datetime(2007, 6, 1), False),
        (HalfYearBegin(1, startingMonth=3), datetime(2007, 6, 2), False),
        (HalfYearBegin(1, startingMonth=6), datetime(2008, 1, 1), False),
        (HalfYearBegin(1, startingMonth=6), datetime(2007, 12, 1), True),
        (HalfYearBegin(1, startingMonth=6), datetime(2008, 2, 1), False),
        (HalfYearBegin(1, startingMonth=6), datetime(2007, 3, 1), False),
        (HalfYearBegin(1, startingMonth=6), datetime(2007, 3, 2), False),
        (HalfYearBegin(1, startingMonth=6), datetime(2008, 4, 1), False),
        (HalfYearBegin(1, startingMonth=6), datetime(2008, 5, 1), False),
        (HalfYearBegin(1, startingMonth=6), datetime(2008, 5, 2), False),
        (HalfYearBegin(1, startingMonth=6), datetime(2007, 6, 1), True),
    ]

    @pytest.mark.parametrize("case", on_offset_cases)
    def test_is_on_offset(self, case):
        offset, dt, expected = case
        assert_is_on_offset(offset, dt, expected)


class TestHalfYearEnd:
    def test_repr(self):
        expected = "<HalfYearEnd: startingMonth=6>"
        assert repr(HalfYearEnd()) == expected
        expected = "<HalfYearEnd: startingMonth=3>"
        assert repr(HalfYearEnd(startingMonth=3)) == expected
        expected = "<HalfYearEnd: startingMonth=1>"
        assert repr(HalfYearEnd(startingMonth=1)) == expected

    def test_offset_corner_case(self):
        # corner
        offset = HalfYearEnd(n=-1, startingMonth=1)
        assert datetime(2010, 2, 1) + offset == datetime(2010, 1, 31)

    offset_cases = []
    offset_cases.append(
        (
            HalfYearEnd(startingMonth=1),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 31),
                datetime(2008, 1, 31): datetime(2008, 7, 31),
                datetime(2008, 2, 15): datetime(2008, 7, 31),
                datetime(2008, 2, 29): datetime(2008, 7, 31),
                datetime(2008, 3, 15): datetime(2008, 7, 31),
                datetime(2008, 3, 31): datetime(2008, 7, 31),
                datetime(2008, 4, 15): datetime(2008, 7, 31),
                datetime(2008, 7, 31): datetime(2009, 1, 31),
            },
        )
    )

    offset_cases.append(
        (
            HalfYearEnd(startingMonth=2),
            {
                datetime(2008, 1, 1): datetime(2008, 2, 29),
                datetime(2008, 1, 31): datetime(2008, 2, 29),
                datetime(2008, 2, 15): datetime(2008, 2, 29),
                datetime(2008, 2, 29): datetime(2008, 8, 31),
                datetime(2008, 3, 15): datetime(2008, 8, 31),
                datetime(2008, 3, 31): datetime(2008, 8, 31),
                datetime(2008, 4, 15): datetime(2008, 8, 31),
                datetime(2008, 8, 30): datetime(2008, 8, 31),
                datetime(2008, 8, 31): datetime(2009, 2, 28),
            },
        )
    )

    offset_cases.append(
        (
            HalfYearEnd(startingMonth=1, n=0),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 31),
                datetime(2008, 1, 31): datetime(2008, 1, 31),
                datetime(2008, 2, 15): datetime(2008, 7, 31),
                datetime(2008, 2, 29): datetime(2008, 7, 31),
                datetime(2008, 3, 15): datetime(2008, 7, 31),
                datetime(2008, 3, 31): datetime(2008, 7, 31),
                datetime(2008, 4, 15): datetime(2008, 7, 31),
                datetime(2008, 7, 31): datetime(2008, 7, 31),
            },
        )
    )

    offset_cases.append(
        (
            HalfYearEnd(startingMonth=1, n=-1),
            {
                datetime(2008, 1, 1): datetime(2007, 7, 31),
                datetime(2008, 1, 31): datetime(2007, 7, 31),
                datetime(2008, 2, 15): datetime(2008, 1, 31),
                datetime(2008, 2, 29): datetime(2008, 1, 31),
                datetime(2008, 3, 15): datetime(2008, 1, 31),
                datetime(2008, 3, 31): datetime(2008, 1, 31),
                datetime(2008, 7, 15): datetime(2008, 1, 31),
                datetime(2008, 7, 30): datetime(2008, 1, 31),
                datetime(2008, 7, 31): datetime(2008, 1, 31),
                datetime(2008, 8, 1): datetime(2008, 7, 31),
            },
        )
    )

    offset_cases.append(
        (
            HalfYearEnd(startingMonth=6, n=2),
            {
                datetime(2008, 1, 31): datetime(2008, 12, 31),
                datetime(2008, 2, 15): datetime(2008, 12, 31),
                datetime(2008, 2, 29): datetime(2008, 12, 31),
                datetime(2008, 3, 15): datetime(2008, 12, 31),
                datetime(2008, 3, 31): datetime(2008, 12, 31),
                datetime(2008, 4, 15): datetime(2008, 12, 31),
                datetime(2008, 4, 30): datetime(2008, 12, 31),
                datetime(2008, 6, 30): datetime(2009, 6, 30),
            },
        )
    )

    @pytest.mark.parametrize("case", offset_cases)
    def test_offset(self, case):
        offset, cases = case
        for base, expected in cases.items():
            assert_offset_equal(offset, base, expected)

    on_offset_cases = [
        (HalfYearEnd(1, startingMonth=1), datetime(2008, 1, 31), True),
        (HalfYearEnd(1, startingMonth=1), datetime(2007, 12, 31), False),
        (HalfYearEnd(1, startingMonth=1), datetime(2008, 2, 29), False),
        (HalfYearEnd(1, startingMonth=1), datetime(2007, 3, 30), False),
        (HalfYearEnd(1, startingMonth=1), datetime(2007, 3, 31), False),
        (HalfYearEnd(1, startingMonth=1), datetime(2008, 4, 30), False),
        (HalfYearEnd(1, startingMonth=1), datetime(2008, 5, 30), False),
        (HalfYearEnd(1, startingMonth=1), datetime(2008, 5, 31), False),
        (HalfYearEnd(1, startingMonth=1), datetime(2007, 6, 29), False),
        (HalfYearEnd(1, startingMonth=1), datetime(2007, 6, 30), False),
        (HalfYearEnd(1, startingMonth=3), datetime(2008, 1, 31), False),
        (HalfYearEnd(1, startingMonth=3), datetime(2007, 12, 31), False),
        (HalfYearEnd(1, startingMonth=3), datetime(2008, 2, 29), False),
        (HalfYearEnd(1, startingMonth=3), datetime(2007, 3, 30), False),
        (HalfYearEnd(1, startingMonth=3), datetime(2007, 3, 31), True),
        (HalfYearEnd(1, startingMonth=3), datetime(2008, 4, 30), False),
        (HalfYearEnd(1, startingMonth=3), datetime(2008, 5, 30), False),
        (HalfYearEnd(1, startingMonth=3), datetime(2008, 5, 31), False),
        (HalfYearEnd(1, startingMonth=3), datetime(2007, 6, 29), False),
        (HalfYearEnd(1, startingMonth=3), datetime(2007, 6, 30), False),
        (HalfYearEnd(1, startingMonth=6), datetime(2008, 1, 31), False),
        (HalfYearEnd(1, startingMonth=6), datetime(2007, 12, 31), True),
        (HalfYearEnd(1, startingMonth=6), datetime(2008, 2, 29), False),
        (HalfYearEnd(1, startingMonth=6), datetime(2007, 3, 30), False),
        (HalfYearEnd(1, startingMonth=6), datetime(2007, 3, 31), False),
        (HalfYearEnd(1, startingMonth=6), datetime(2008, 4, 30), False),
        (HalfYearEnd(1, startingMonth=6), datetime(2008, 5, 30), False),
        (HalfYearEnd(1, startingMonth=6), datetime(2008, 5, 31), False),
        (HalfYearEnd(1, startingMonth=6), datetime(2007, 6, 29), False),
        (HalfYearEnd(1, startingMonth=6), datetime(2007, 6, 30), True),
    ]

    @pytest.mark.parametrize("case", on_offset_cases)
    def test_is_on_offset(self, case):
        offset, dt, expected = case
        assert_is_on_offset(offset, dt, expected)
