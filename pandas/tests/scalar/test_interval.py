from __future__ import division

from pandas import Interval

import pytest
import pandas.util.testing as tm


@pytest.fixture
def interval():
    return Interval(0, 1)


class TestInterval(object):

    def test_properties(self, interval):
        assert interval.closed == 'right'
        assert interval.left == 0
        assert interval.right == 1
        assert interval.mid == 0.5

    def test_repr(self, interval):
        assert repr(interval) == "Interval(0, 1, closed='right')"
        assert str(interval) == "(0, 1]"

        interval_left = Interval(0, 1, closed='left')
        assert repr(interval_left) == "Interval(0, 1, closed='left')"
        assert str(interval_left) == "[0, 1)"

    def test_contains(self, interval):
        assert 0.5 in interval
        assert 1 in interval
        assert 0 not in interval

        msg = "__contains__ not defined for two intervals"
        with tm.assert_raises_regex(TypeError, msg):
            interval in interval

        interval_both = Interval(0, 1, closed='both')
        assert 0 in interval_both
        assert 1 in interval_both

        interval_neither = Interval(0, 1, closed='neither')
        assert 0 not in interval_neither
        assert 0.5 in interval_neither
        assert 1 not in interval_neither

    def test_equal(self):
        assert Interval(0, 1) == Interval(0, 1, closed='right')
        assert Interval(0, 1) != Interval(0, 1, closed='left')
        assert Interval(0, 1) != 0

    def test_comparison(self):
        with tm.assert_raises_regex(TypeError, 'unorderable types'):
            Interval(0, 1) < 2

        assert Interval(0, 1) < Interval(1, 2)
        assert Interval(0, 1) < Interval(0, 2)
        assert Interval(0, 1) < Interval(0.5, 1.5)
        assert Interval(0, 1) <= Interval(0, 1)
        assert Interval(0, 1) > Interval(-1, 2)
        assert Interval(0, 1) >= Interval(0, 1)

    def test_hash(self, interval):
        # should not raise
        hash(interval)

    def test_math_add(self, interval):
        expected = Interval(1, 2)
        actual = interval + 1
        assert expected == actual

        expected = Interval(1, 2)
        actual = 1 + interval
        assert expected == actual

        actual = interval
        actual += 1
        assert expected == actual

        msg = "unsupported operand type\(s\) for \+"
        with tm.assert_raises_regex(TypeError, msg):
            interval + Interval(1, 2)

        with tm.assert_raises_regex(TypeError, msg):
            interval + 'foo'

    def test_math_sub(self, interval):
        expected = Interval(-1, 0)
        actual = interval - 1
        assert expected == actual

        actual = interval
        actual -= 1
        assert expected == actual

        msg = "unsupported operand type\(s\) for -"
        with tm.assert_raises_regex(TypeError, msg):
            interval - Interval(1, 2)

        with tm.assert_raises_regex(TypeError, msg):
            interval - 'foo'

    def test_math_mult(self, interval):
        expected = Interval(0, 2)
        actual = interval * 2
        assert expected == actual

        expected = Interval(0, 2)
        actual = 2 * interval
        assert expected == actual

        actual = interval
        actual *= 2
        assert expected == actual

        msg = "unsupported operand type\(s\) for \*"
        with tm.assert_raises_regex(TypeError, msg):
            interval * Interval(1, 2)

        msg = "can\'t multiply sequence by non-int"
        with tm.assert_raises_regex(TypeError, msg):
            interval * 'foo'

    def test_math_div(self, interval):
        expected = Interval(0, 0.5)
        actual = interval / 2.0
        assert expected == actual

        actual = interval
        actual /= 2.0
        assert expected == actual

        msg = "unsupported operand type\(s\) for /"
        with tm.assert_raises_regex(TypeError, msg):
            interval / Interval(1, 2)

        with tm.assert_raises_regex(TypeError, msg):
            interval / 'foo'
