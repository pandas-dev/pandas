import numpy as np
import pytest

from pandas import (
    Interval,
    Timedelta,
    Timestamp,
)
import pandas._testing as tm


@pytest.fixture(
    params=[
        (Timedelta("0 days"), Timedelta("1 day")),
        (Timestamp("2018-01-01"), Timedelta("1 day")),
        (0, 1),
    ],
    ids=lambda x: type(x[0]).__name__,
)
def start_shift(request):
    """
    Fixture for generating intervals of types from a start value and a shift
    value that can be added to start to generate an endpoint
    """
    return request.param


class TestOverlaps:
    def test_overlaps_self(self, start_shift, closed):
        start, shift = start_shift
        interval = Interval(start, start + shift, closed)
        assert interval.overlaps(interval)

    def test_overlaps_nested(self, start_shift, closed, other_closed):
        start, shift = start_shift
        interval1 = Interval(start, start + 3 * shift, other_closed)
        interval2 = Interval(start + shift, start + 2 * shift, closed)

        # nested intervals should always overlap
        assert interval1.overlaps(interval2)

    def test_overlaps_disjoint(self, start_shift, closed, other_closed):
        start, shift = start_shift
        interval1 = Interval(start, start + shift, other_closed)
        interval2 = Interval(start + 2 * shift, start + 3 * shift, closed)

        # disjoint intervals should never overlap
        assert not interval1.overlaps(interval2)

    def test_overlaps_endpoint(self, start_shift, closed, other_closed):
        start, shift = start_shift
        interval1 = Interval(start, start + shift, other_closed)
        interval2 = Interval(start + shift, start + 2 * shift, closed)

        # overlap if shared endpoint is closed for both (overlap at a point)
        result = interval1.overlaps(interval2)
        expected = interval1.closed_right and interval2.closed_left
        assert result == expected

    @pytest.mark.parametrize(
        "other",
        [10, True, "foo", Timedelta("1 day"), Timestamp("2018-01-01")],
        ids=lambda x: type(x).__name__,
    )
    def test_overlaps_invalid_type(self, other):
        interval = Interval(0, 1)
        msg = f"`other` must be an Interval, got {type(other).__name__}"
        with pytest.raises(TypeError, match=msg):
            interval.overlaps(other)


class TestIntersection:
    def test_intersection_self(self):
        interval = Interval(1, 8, "left")
        assert interval.intersection(interval) == interval

    def test_intersection_include_limits(self):
        other = Interval(1, 8, "left")

        intervals = np.array(
            [
                Interval(7, 9, "left"),  # include left
                Interval(0, 2, "right"),  # include right
                Interval(1, 8, "right"),  # open limit
            ]
        )

        expected = np.array(
            [
                Interval(7, 8, "left"),
                Interval(1, 2, "both"),
                Interval(1, 8, "neither"),
            ]
        )

        result = np.array([interval.intersection(other) for interval in intervals])
        tm.assert_numpy_array_equal(result, expected)

    def test_intersection_overlapping(self):
        other = Interval(1, 8, "left")

        intervals = np.array(
            [
                Interval(2, 4, "both"),  # nested
                Interval(0, 9, "both"),  # spanning
                Interval(4, 10, "both"),  # partial
            ]
        )

        expected = np.array(
            [
                Interval(2, 4, "both"),
                Interval(1, 8, "left"),
                Interval(4, 8, "left"),
            ]
        )

        result = np.array([interval.intersection(other) for interval in intervals])
        tm.assert_numpy_array_equal(result, expected)

    def test_intersection_adjacent(self):
        other = Interval(1, 8, "left")

        intervals = np.array(
            [
                Interval(-5, 1, "both"),  # adjacent closed
                Interval(8, 10, "both"),  # adjacent open
                Interval(10, 15, "both"),  # disjoint
            ]
        )

        expected = np.array(
            [
                Interval(1, 1, "both"),
                None,
                None,
            ]
        )

        result = np.array([interval.intersection(other) for interval in intervals])
        tm.assert_numpy_array_equal(result, expected)

    def test_intersection_timestamps(self):
        year_2020 = Interval(
            Timestamp("2020-01-01 00:00:00"),
            Timestamp("2021-01-01 00:00:00"),
            closed="left",
        )

        march_2020 = Interval(
            Timestamp("2020-03-01 00:00:00"),
            Timestamp("2020-04-01 00:00:00"),
            closed="left",
        )

        result = year_2020.intersection(march_2020)
        assert result == march_2020

    @pytest.mark.parametrize(
        "other",
        [10, True, "foo", Timedelta("1 day"), Timestamp("2018-01-01")],
        ids=lambda x: type(x).__name__,
    )
    def test_intersection_invalid_type(self, other):
        interval = Interval(0, 1)
        msg = f"`other` must be an Interval, got {type(other).__name__}"
        with pytest.raises(TypeError, match=msg):
            interval.intersection(other)


class TestUnion:
    def test_union_self(self):
        interval = Interval(1, 8, "left")

        result = interval.union(interval)

        expected = np.array([interval], dtype=object)
        tm.assert_numpy_array_equal(result, expected)

    def test_union_include_limits(self):
        other = Interval(1, 8, "left")

        intervals = np.array(
            [
                Interval(7, 9, "left"),  # include left
                Interval(0, 2, "right"),  # include right
                Interval(1, 8, "right"),  # open limit
            ]
        )

        expected = np.array(
            [
                np.array([Interval(1, 9, "left")], dtype=object),
                np.array([Interval(0, 8, "neither")], dtype=object),
                np.array([Interval(1, 8, "both")], dtype=object),
            ],
            dtype=object,
        )

        result = np.array([interval.union(other) for interval in intervals])
        tm.assert_numpy_array_equal(result, expected)

    def test_union_overlapping(self):
        other = Interval(1, 8, "left")

        intervals = np.array(
            [
                Interval(2, 4, "both"),  # nested
                Interval(0, 9, "both"),  # spanning
                Interval(4, 10, "both"),  # partial
            ]
        )

        expected = np.array(
            [
                np.array([Interval(1, 8, "left")], dtype=object),
                np.array([Interval(0, 9, "both")], dtype=object),
                np.array([Interval(1, 10, "both")], dtype=object),
            ],
            dtype=object,
        )

        result = np.array([interval.union(other) for interval in intervals])
        tm.assert_numpy_array_equal(result, expected)

    def test_union_adjacent(self):
        other = Interval(1, 8, "left")

        intervals = np.array(
            [
                Interval(-5, 1, "both"),  # adjacent closed
                Interval(8, 10, "both"),  # adjacent open
                Interval(10, 15, "both"),  # disjoint
            ]
        )

        expected = np.array(
            [
                np.array([Interval(-5, 8, "left")], dtype=object),
                np.array([Interval(1, 10, "both")], dtype=object),
                np.array([other, Interval(10, 15, "both")], dtype=object),
            ],
            dtype=object,
        )

        result = np.array(
            [interval.union(other) for interval in intervals], dtype=object
        )
        tm.assert_numpy_array_equal(result, expected)

    def test_union_timestamps(self):
        year_2020 = Interval(
            Timestamp("2020-01-01 00:00:00"),
            Timestamp("2021-01-01 00:00:00"),
            closed="left",
        )

        year_2021 = Interval(
            Timestamp("2021-01-01 00:00:00"),
            Timestamp("2022-01-01 00:00:00"),
            closed="left",
        )

        expected = np.array(
            [
                Interval(
                    Timestamp("2020-01-01 00:00:00"),
                    Timestamp("2022-01-01 00:00:00"),
                    closed="left",
                )
            ],
            dtype=object,
        )

        result = year_2020.union(year_2021)
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize(
        "other",
        [10, True, "foo", Timedelta("1 day"), Timestamp("2018-01-01")],
        ids=lambda x: type(x).__name__,
    )
    def test_union_invalid_type(self, other):
        interval = Interval(0, 1)
        msg = f"`other` must be an Interval, got {type(other).__name__}"
        with pytest.raises(TypeError, match=msg):
            interval.union(other)
