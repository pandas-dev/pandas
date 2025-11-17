"""
Tests for the following offsets:
- Easter
"""

from __future__ import annotations

from datetime import datetime

from dateutil.easter import (
    EASTER_ORTHODOX,
    EASTER_WESTERN,
)
import pytest

from pandas.tests.tseries.offsets.common import assert_offset_equal

from pandas.tseries.offsets import Easter


class TestEaster:
    @pytest.mark.parametrize(
        "offset,date,expected",
        [
            (Easter(), datetime(2010, 1, 1), datetime(2010, 4, 4)),
            (Easter(), datetime(2010, 4, 5), datetime(2011, 4, 24)),
            (Easter(2), datetime(2010, 1, 1), datetime(2011, 4, 24)),
            (Easter(), datetime(2010, 4, 4), datetime(2011, 4, 24)),
            (Easter(2), datetime(2010, 4, 4), datetime(2012, 4, 8)),
            (-Easter(), datetime(2011, 1, 1), datetime(2010, 4, 4)),
            (-Easter(), datetime(2010, 4, 5), datetime(2010, 4, 4)),
            (-Easter(2), datetime(2011, 1, 1), datetime(2009, 4, 12)),
            (-Easter(), datetime(2010, 4, 4), datetime(2009, 4, 12)),
            (-Easter(2), datetime(2010, 4, 4), datetime(2008, 3, 23)),
        ],
    )
    def test_offset(self, offset, date, expected):
        assert_offset_equal(offset, date, expected)

    @pytest.mark.parametrize(
        "offset,date,expected",
        [
            (Easter(method=EASTER_WESTERN), datetime(2010, 1, 1), datetime(2010, 4, 4)),
            (
                Easter(method=EASTER_WESTERN),
                datetime(2010, 4, 5),
                datetime(2011, 4, 24),
            ),
            (
                Easter(2, method=EASTER_WESTERN),
                datetime(2010, 1, 1),
                datetime(2011, 4, 24),
            ),
            (
                Easter(method=EASTER_WESTERN),
                datetime(2010, 4, 4),
                datetime(2011, 4, 24),
            ),
            (
                Easter(2, method=EASTER_WESTERN),
                datetime(2010, 4, 4),
                datetime(2012, 4, 8),
            ),
            (
                -Easter(method=EASTER_WESTERN),
                datetime(2011, 1, 1),
                datetime(2010, 4, 4),
            ),
            (
                -Easter(method=EASTER_WESTERN),
                datetime(2010, 4, 5),
                datetime(2010, 4, 4),
            ),
            (
                -Easter(2, method=EASTER_WESTERN),
                datetime(2011, 1, 1),
                datetime(2009, 4, 12),
            ),
            (
                -Easter(method=EASTER_WESTERN),
                datetime(2010, 4, 4),
                datetime(2009, 4, 12),
            ),
            (
                -Easter(2, method=EASTER_WESTERN),
                datetime(2010, 4, 4),
                datetime(2008, 3, 23),
            ),
        ],
    )
    def test_western_easter_offset(self, offset, date, expected):
        assert_offset_equal(offset, date, expected)

    @pytest.mark.parametrize(
        "offset,date,expected",
        [
            (
                Easter(method=EASTER_ORTHODOX),
                datetime(2010, 1, 1),
                datetime(2010, 4, 4),
            ),
            (
                Easter(method=EASTER_ORTHODOX),
                datetime(2010, 4, 5),
                datetime(2011, 4, 24),
            ),
            (
                Easter(2, method=EASTER_ORTHODOX),
                datetime(2010, 1, 1),
                datetime(2011, 4, 24),
            ),
            (
                Easter(method=EASTER_ORTHODOX),
                datetime(2010, 4, 4),
                datetime(2011, 4, 24),
            ),
            (
                Easter(2, method=EASTER_ORTHODOX),
                datetime(2010, 4, 4),
                datetime(2012, 4, 15),
            ),
            (
                -Easter(method=EASTER_ORTHODOX),
                datetime(2011, 1, 1),
                datetime(2010, 4, 4),
            ),
            (
                -Easter(method=EASTER_ORTHODOX),
                datetime(2010, 4, 5),
                datetime(2010, 4, 4),
            ),
            (
                -Easter(2, method=EASTER_ORTHODOX),
                datetime(2011, 1, 1),
                datetime(2009, 4, 19),
            ),
            (
                -Easter(method=EASTER_ORTHODOX),
                datetime(2010, 4, 4),
                datetime(2009, 4, 19),
            ),
            (
                -Easter(2, method=EASTER_ORTHODOX),
                datetime(2010, 4, 4),
                datetime(2008, 4, 27),
            ),
        ],
    )
    def test_orthodox_easter_offset(self, offset, date, expected):
        assert_offset_equal(offset, date, expected)
