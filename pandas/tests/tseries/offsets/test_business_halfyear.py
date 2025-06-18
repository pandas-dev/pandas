"""
Tests for the following offsets:
- BHalfYearBegin
- BHalfYearEnd
"""

from __future__ import annotations

from datetime import datetime

import pytest

from pandas.tests.tseries.offsets.common import (
    assert_is_on_offset,
    assert_offset_equal,
)

from pandas.tseries.offsets import (
    BHalfYearBegin,
    BHalfYearEnd,
)


@pytest.mark.parametrize("klass", (BHalfYearBegin, BHalfYearEnd))
def test_halfyearly_dont_normalize(klass):
    date = datetime(2012, 3, 31, 5, 30)
    result = date + klass()
    assert result.time() == date.time()


@pytest.mark.parametrize("offset", [BHalfYearBegin(), BHalfYearEnd()])
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


class TestBHalfYearBegin:
    def test_repr(self):
        expected = "<BusinessHalfYearBegin: startingMonth=1>"
        assert repr(BHalfYearBegin()) == expected
        expected = "<BusinessHalfYearBegin: startingMonth=3>"
        assert repr(BHalfYearBegin(startingMonth=3)) == expected
        expected = "<BusinessHalfYearBegin: startingMonth=1>"
        assert repr(BHalfYearBegin(startingMonth=1)) == expected

    def test_offset_corner_case(self):
        # corner
        offset = BHalfYearBegin(n=-1, startingMonth=1)
        assert datetime(2010, 2, 1) + offset == datetime(2010, 1, 1)

    offset_cases = []
    offset_cases.append(
        (
            BHalfYearBegin(startingMonth=1),
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
            BHalfYearBegin(startingMonth=2),
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
            BHalfYearBegin(startingMonth=1, n=0),
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
            BHalfYearBegin(startingMonth=1, n=-1),
            {
                datetime(2008, 1, 1): datetime(2007, 7, 2),
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
            BHalfYearBegin(startingMonth=1, n=2),
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
        (BHalfYearBegin(1, startingMonth=1), datetime(2008, 1, 1), True),
        (BHalfYearBegin(1, startingMonth=1), datetime(2007, 12, 1), False),
        (BHalfYearBegin(1, startingMonth=1), datetime(2008, 2, 1), False),
        (BHalfYearBegin(1, startingMonth=1), datetime(2007, 3, 1), False),
        (BHalfYearBegin(1, startingMonth=1), datetime(2008, 4, 1), False),
        (BHalfYearBegin(1, startingMonth=1), datetime(2008, 5, 1), False),
        (BHalfYearBegin(1, startingMonth=1), datetime(2007, 6, 1), False),
        (BHalfYearBegin(1, startingMonth=3), datetime(2008, 1, 1), False),
        (BHalfYearBegin(1, startingMonth=3), datetime(2007, 12, 1), False),
        (BHalfYearBegin(1, startingMonth=3), datetime(2008, 2, 1), False),
        (BHalfYearBegin(1, startingMonth=3), datetime(2007, 3, 1), True),
        (BHalfYearBegin(1, startingMonth=3), datetime(2008, 4, 1), False),
        (BHalfYearBegin(1, startingMonth=3), datetime(2008, 5, 1), False),
        (BHalfYearBegin(1, startingMonth=3), datetime(2008, 5, 2), False),
        (BHalfYearBegin(1, startingMonth=3), datetime(2007, 6, 1), False),
        (BHalfYearBegin(1, startingMonth=3), datetime(2007, 6, 2), False),
        (BHalfYearBegin(1, startingMonth=6), datetime(2008, 1, 1), False),
        (BHalfYearBegin(1, startingMonth=6), datetime(2007, 12, 3), True),
        (BHalfYearBegin(1, startingMonth=6), datetime(2008, 2, 1), False),
        (BHalfYearBegin(1, startingMonth=6), datetime(2007, 3, 1), False),
        (BHalfYearBegin(1, startingMonth=6), datetime(2007, 3, 2), False),
        (BHalfYearBegin(1, startingMonth=6), datetime(2008, 4, 1), False),
        (BHalfYearBegin(1, startingMonth=6), datetime(2008, 5, 1), False),
        (BHalfYearBegin(1, startingMonth=6), datetime(2008, 5, 2), False),
        (BHalfYearBegin(1, startingMonth=6), datetime(2007, 6, 1), True),
    ]

    @pytest.mark.parametrize("case", on_offset_cases)
    def test_is_on_offset(self, case):
        offset, dt, expected = case
        assert_is_on_offset(offset, dt, expected)


class TestBHalfYearEnd:
    def test_repr(self):
        expected = "<BusinessHalfYearEnd: startingMonth=6>"
        assert repr(BHalfYearEnd()) == expected
        expected = "<BusinessHalfYearEnd: startingMonth=3>"
        assert repr(BHalfYearEnd(startingMonth=3)) == expected
        expected = "<BusinessHalfYearEnd: startingMonth=1>"
        assert repr(BHalfYearEnd(startingMonth=1)) == expected

    def test_offset_corner_case(self):
        # corner
        offset = BHalfYearEnd(n=-1, startingMonth=1)
        assert datetime(2010, 1, 30) + offset == datetime(2010, 1, 29)

    offset_cases = []
    offset_cases.append(
        (
            BHalfYearEnd(startingMonth=1),
            {
                datetime(2008, 1, 1): datetime(2008, 1, 31),
                datetime(2008, 1, 31): datetime(2008, 7, 31),
                datetime(2008, 2, 15): datetime(2008, 7, 31),
                datetime(2008, 2, 29): datetime(2008, 7, 31),
                datetime(2008, 3, 15): datetime(2008, 7, 31),
                datetime(2008, 3, 31): datetime(2008, 7, 31),
                datetime(2008, 4, 15): datetime(2008, 7, 31),
                datetime(2008, 7, 31): datetime(2009, 1, 30),
            },
        )
    )

    offset_cases.append(
        (
            BHalfYearEnd(startingMonth=2),
            {
                datetime(2008, 1, 1): datetime(2008, 2, 29),
                datetime(2008, 1, 31): datetime(2008, 2, 29),
                datetime(2008, 2, 15): datetime(2008, 2, 29),
                datetime(2008, 2, 29): datetime(2008, 8, 29),
                datetime(2008, 3, 15): datetime(2008, 8, 29),
                datetime(2008, 3, 31): datetime(2008, 8, 29),
                datetime(2008, 4, 15): datetime(2008, 8, 29),
                datetime(2008, 8, 28): datetime(2008, 8, 29),
                datetime(2008, 8, 29): datetime(2009, 2, 27),
            },
        )
    )

    offset_cases.append(
        (
            BHalfYearEnd(startingMonth=1, n=0),
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
            BHalfYearEnd(startingMonth=1, n=-1),
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
            BHalfYearEnd(startingMonth=6, n=2),
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
        (BHalfYearEnd(1, startingMonth=1), datetime(2008, 1, 31), True),
        (BHalfYearEnd(1, startingMonth=1), datetime(2007, 12, 31), False),
        (BHalfYearEnd(1, startingMonth=1), datetime(2008, 2, 29), False),
        (BHalfYearEnd(1, startingMonth=1), datetime(2007, 3, 30), False),
        (BHalfYearEnd(1, startingMonth=1), datetime(2007, 3, 31), False),
        (BHalfYearEnd(1, startingMonth=1), datetime(2008, 4, 30), False),
        (BHalfYearEnd(1, startingMonth=1), datetime(2008, 5, 30), False),
        (BHalfYearEnd(1, startingMonth=1), datetime(2008, 5, 31), False),
        (BHalfYearEnd(1, startingMonth=1), datetime(2007, 6, 29), False),
        (BHalfYearEnd(1, startingMonth=1), datetime(2007, 6, 30), False),
        (BHalfYearEnd(1, startingMonth=3), datetime(2008, 1, 31), False),
        (BHalfYearEnd(1, startingMonth=3), datetime(2007, 12, 31), False),
        (BHalfYearEnd(1, startingMonth=3), datetime(2008, 2, 29), False),
        (BHalfYearEnd(1, startingMonth=3), datetime(2007, 3, 30), True),
        (BHalfYearEnd(1, startingMonth=3), datetime(2007, 3, 31), False),
        (BHalfYearEnd(1, startingMonth=3), datetime(2008, 4, 30), False),
        (BHalfYearEnd(1, startingMonth=3), datetime(2008, 5, 30), False),
        (BHalfYearEnd(1, startingMonth=3), datetime(2008, 5, 31), False),
        (BHalfYearEnd(1, startingMonth=3), datetime(2007, 6, 29), False),
        (BHalfYearEnd(1, startingMonth=3), datetime(2007, 6, 30), False),
        (BHalfYearEnd(1, startingMonth=6), datetime(2008, 1, 31), False),
        (BHalfYearEnd(1, startingMonth=6), datetime(2007, 12, 31), True),
        (BHalfYearEnd(1, startingMonth=6), datetime(2008, 2, 29), False),
        (BHalfYearEnd(1, startingMonth=6), datetime(2007, 3, 30), False),
        (BHalfYearEnd(1, startingMonth=6), datetime(2007, 3, 31), False),
        (BHalfYearEnd(1, startingMonth=6), datetime(2008, 4, 30), False),
        (BHalfYearEnd(1, startingMonth=6), datetime(2008, 5, 30), False),
        (BHalfYearEnd(1, startingMonth=6), datetime(2008, 5, 31), False),
        (BHalfYearEnd(1, startingMonth=6), datetime(2007, 6, 29), True),
        (BHalfYearEnd(1, startingMonth=6), datetime(2007, 6, 30), False),
    ]

    @pytest.mark.parametrize("case", on_offset_cases)
    def test_is_on_offset(self, case):
        offset, dt, expected = case
        assert_is_on_offset(offset, dt, expected)
