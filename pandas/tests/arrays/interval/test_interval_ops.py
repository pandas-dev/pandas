"""Tests for Interval-Interval operations, such as overlaps, contains, etc."""
import numpy as np
import pytest

import pandas.util.testing as tm
from pandas import Interval, IntervalIndex, Timedelta, Timestamp
from pandas.core.arrays import IntervalArray


@pytest.fixture(params=[IntervalArray, IntervalIndex])
def constructor(request):
    """
    Fixture for IntervalArray and IntervalIndex class constructors
    """
    return request.param


@pytest.fixture(params=[
    (Timedelta('0 days'), Timedelta('1 day')),
    (Timestamp('2018-01-01'), Timedelta('1 day')),
    (0, 1)], ids=lambda x: type(x[0]).__name__)
def start_shift(request):
    """
    Fixture for generating intervals of types from a start value and a shift
    value that can be added to start to generate an endpoint
    """
    return request.param


class TestOverlaps(object):

    def test_overlaps_interval(
            self, constructor, start_shift, closed, other_closed):
        start, shift = start_shift
        interval = Interval(start, start + 3 * shift, other_closed)

        # container: identical, nested, spanning, partial, adjacent, disjoint
        tuples = [(start, start + 3 * shift),
                  (start + shift, start + 2 * shift),
                  (start - shift, start + 4 * shift),
                  (start + 2 * shift, start + 4 * shift),
                  (start + 3 * shift, start + 4 * shift),
                  (start + 4 * shift, start + 5 * shift)]
        interval_container = constructor.from_tuples(tuples, closed)

        adjacent = (interval.closed_right and interval_container.closed_left)
        expected = np.array([True, True, True, True, adjacent, False])
        result = interval_container.overlaps(interval)
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize('other_constructor', [
        IntervalArray, IntervalIndex])
    def test_overlaps_interval_container(self, constructor, other_constructor):
        # TODO: modify this test when implemented
        interval_container = constructor.from_breaks(range(5))
        other_container = other_constructor.from_breaks(range(5))
        with pytest.raises(NotImplementedError):
            interval_container.overlaps(other_container)

    def test_overlaps_na(self, constructor, start_shift):
        """NA values are marked as False"""
        start, shift = start_shift
        interval = Interval(start, start + shift)

        tuples = [(start, start + shift),
                  np.nan,
                  (start + 2 * shift, start + 3 * shift)]
        interval_container = constructor.from_tuples(tuples)

        expected = np.array([True, False, False])
        result = interval_container.overlaps(interval)
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize('other', [
        10, True, 'foo', Timedelta('1 day'), Timestamp('2018-01-01')],
        ids=lambda x: type(x).__name__)
    def test_overlaps_invalid_type(self, constructor, other):
        interval_container = constructor.from_breaks(range(5))
        msg = '`other` must be Interval-like, got {other}'.format(
            other=type(other).__name__)
        with tm.assert_raises_regex(TypeError, msg):
            interval_container.overlaps(other)
