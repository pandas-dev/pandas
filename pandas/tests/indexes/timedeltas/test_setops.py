import numpy as np
import pytest

import pandas as pd
import pandas.util.testing as tm
from pandas import Int64Index, TimedeltaIndex, timedelta_range


class TestTimedeltaIndex(object):

    def test_union(self):

        i1 = timedelta_range('1day', periods=5)
        i2 = timedelta_range('3day', periods=5)
        result = i1.union(i2)
        expected = timedelta_range('1day', periods=7)
        tm.assert_index_equal(result, expected)

        i1 = Int64Index(np.arange(0, 20, 2))
        i2 = timedelta_range(start='1 day', periods=10, freq='D')
        i1.union(i2)  # Works
        i2.union(i1)  # Fails with "AttributeError: can't set attribute"

    def test_union_coverage(self):

        idx = TimedeltaIndex(['3d', '1d', '2d'])
        ordered = TimedeltaIndex(idx.sort_values(), freq='infer')
        result = ordered.union(idx)
        tm.assert_index_equal(result, ordered)

        result = ordered[:0].union(ordered)
        tm.assert_index_equal(result, ordered)
        assert result.freq == ordered.freq

    def test_union_bug_1730(self):

        rng_a = timedelta_range('1 day', periods=4, freq='3H')
        rng_b = timedelta_range('1 day', periods=4, freq='4H')

        result = rng_a.union(rng_b)
        exp = TimedeltaIndex(sorted(set(list(rng_a)) | set(list(rng_b))))
        tm.assert_index_equal(result, exp)

    def test_union_bug_1745(self):

        left = TimedeltaIndex(['1 day 15:19:49.695000'])
        right = TimedeltaIndex(['2 day 13:04:21.322000',
                                '1 day 15:27:24.873000',
                                '1 day 15:31:05.350000'])

        result = left.union(right)
        exp = TimedeltaIndex(sorted(set(list(left)) | set(list(right))))
        tm.assert_index_equal(result, exp)

    def test_union_bug_4564(self):

        left = timedelta_range("1 day", "30d")
        right = left + pd.offsets.Minute(15)

        result = left.union(right)
        exp = TimedeltaIndex(sorted(set(list(left)) | set(list(right))))
        tm.assert_index_equal(result, exp)

    def test_intersection_bug_1708(self):
        index_1 = timedelta_range('1 day', periods=4, freq='h')
        index_2 = index_1 + pd.offsets.Hour(5)

        result = index_1 & index_2
        assert len(result) == 0

        index_1 = timedelta_range('1 day', periods=4, freq='h')
        index_2 = index_1 + pd.offsets.Hour(1)

        result = index_1 & index_2
        expected = timedelta_range('1 day 01:00:00', periods=3, freq='h')
        tm.assert_index_equal(result, expected)


@pytest.mark.parametrize('idx1,idx2,expected', [
    (pd.to_timedelta(range(2, 6), unit='s'),
     pd.to_timedelta(range(3), unit='s'),
     TimedeltaIndex(['00:00:002'])),
    (pd.to_timedelta(range(3), unit='s'),
     pd.to_timedelta(range(2, 6), unit='s'),
     TimedeltaIndex(['00:00:002'])),
])
def test_intersection_intersects_ascending(idx1, idx2, expected):
    result = idx1.intersection(idx2)
    assert result.equals(expected)


@pytest.mark.parametrize('idx1,idx2,expected', [
    (pd.to_timedelta(range(6, 3, -1), unit='s'),
     pd.to_timedelta(range(5, 1, -1), unit='s'),
     TimedeltaIndex(['00:00:05', '00:00:04'])),
    (pd.to_timedelta(range(5, 1, -1), unit='s'),
     pd.to_timedelta(range(6, 3, -1), unit='s'),
     TimedeltaIndex(['00:00:05', '00:00:04'])),
])
def test_intersection_intersects_descending(idx1, idx2, expected):
    # GH 17391
    result = idx1.intersection(idx2)
    assert result.equals(expected)


def test_intersection_intersects_descending_no_intersect():
    idx1 = pd.to_timedelta(range(6, 4, -1), unit='s')
    idx2 = pd.to_timedelta(range(4, 1, -1), unit='s')
    result = idx1.intersection(idx2)
    assert len(result) == 0


def test_intersection_intersects_len_1():
    idx1 = pd.to_timedelta(range(1, 2), unit='s')
    idx2 = pd.to_timedelta(range(1, 0, -1), unit='s')
    intersection = idx1.intersection(idx2)
    expected = TimedeltaIndex(['00:00:01'],
                              dtype='timedelta64[ns]')
    tm.assert_index_equal(intersection, expected)


def test_intersection_can_intersect_self():
    idx = pd.to_timedelta(range(1, 2), unit='s')
    result = idx.intersection(idx)
    tm.assert_index_equal(idx, result)


def test_intersection_not_sorted():
    idx1 = pd.to_timedelta((1, 3, 2, 5, 4), unit='s')
    idx2 = pd.to_timedelta((1, 2, 3, 5, 4), unit='s')
    result = idx1.intersection(idx2)
    expected = idx1
    tm.assert_index_equal(result, expected)


def test_intersection_not_unique():
    idx1 = pd.to_timedelta((1, 2, 2, 3, 3, 5), unit='s')
    idx2 = pd.to_timedelta((1, 2, 3, 4), unit='s')
    result = idx1.intersection(idx2)
    expected = pd.to_timedelta((1, 2, 2, 3, 3), unit='s')
    tm.assert_index_equal(result, expected)

    result = idx2.intersection(idx1)
    expected = pd.to_timedelta((1, 2, 2, 3, 3), unit='s')
    tm.assert_index_equal(result, expected)


@pytest.mark.parametrize("index1, index2, expected", [
    (pd.to_timedelta((1, 2, 3, 4, 5, 6, 7, 8), unit='s'),
     pd.to_timedelta((2, 3, 4, 8), unit='s'),
     pd.to_timedelta((2, 3, 4, 8), unit='s')),
    (pd.to_timedelta((1, 2, 3, 4, 5), unit='s'),
     pd.to_timedelta((2, 3, 4), unit='s'),
     pd.to_timedelta((2, 3, 4), unit='s')),
    (pd.to_timedelta((2, 4, 5, 6), unit='s'),
     pd.to_timedelta((2, 3, 4), unit='s'),
     pd.to_timedelta((2, 4), unit='s')),
])
def test_intersection_different_lengths(index1, index2, expected):
    def intersect(idx1, idx2, expected):
        result = idx1.intersection(idx2)
        tm.assert_index_equal(result, expected)
        result = idx2.intersection(idx1)
        tm.assert_index_equal(result, expected)

    intersect(index1, index2, expected)
    intersect(index1.sort_values(ascending=False),
              index2.sort_values(ascending=False),
              expected.sort_values(ascending=False)
              )


@pytest.mark.parametrize("index1, index2, expected", [
    (pd.to_timedelta((2, 4, 5, 6), unit='s'),
     pd.to_timedelta((2, 3, 4, 6), unit='s'),
     pd.to_timedelta((2, 4, 6), unit='s')),
    (pd.to_timedelta((2, 4, 5), unit='s'),
     pd.to_timedelta((3, 4, 5, 6), unit='s'),
     pd.to_timedelta((4, 5), unit='s')),
])
def test_intersection_not_a_subset(index1, index2, expected):
    def intersect(idx1, idx2, expected):
        result = idx1.intersection(idx2)
        tm.assert_index_equal(result, expected)
        result = idx2.intersection(idx1)
        tm.assert_index_equal(result, expected)

    intersect(index1, index2, expected)
    intersect(index1.sort_values(ascending=False),
              index2.sort_values(ascending=False),
              expected.sort_values(ascending=False)
              )
