import numpy as np
import pytest

import pandas as pd
from pandas import Int64Index, TimedeltaIndex, timedelta_range
import pandas.util.testing as tm


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

    @pytest.mark.parametrize("sort", [None, False])
    def test_intersection_equal(self, sort):
        # for equal indicies intersection should return the original index
        first = timedelta_range('1 day', periods=4, freq='h')
        second = timedelta_range('1 day', periods=4, freq='h')
        intersect = first.intersection(second, sort=sort)
        if sort is None:
            tm.assert_index_equal(intersect, second.sort_values())
        assert tm.equalContents(intersect, second)

        # Corner cases
        inter = first.intersection(first, sort=sort)
        assert inter is first

    @pytest.mark.parametrize("sort", [None, False])
    def test_intersection_zero_length(self, sort):
        index_1 = timedelta_range('1 day', periods=4, freq='h')
        index_2 = timedelta_range('1 day', periods=0, freq='h')
        inter = index_1.intersection(index_2, sort=sort)
        tm.assert_index_equal(index_2, inter)
        inter_2 = index_2.intersection(index_1, sort=sort)
        tm.assert_index_equal(index_2, inter_2)

    @pytest.mark.parametrize("sort", [None, False])
    def test_intersection(self, sort):
        # GH 4690 (with tz)
        base = timedelta_range('1 day', periods=4, freq='h', name='idx')

        # if target has the same name, it is preserved
        rng2 = timedelta_range('1 day', periods=5, freq='h', name='idx')
        expected2 = timedelta_range('1 day', periods=4, freq='h', name='idx')

        # if target name is different, it will be reset
        rng3 = timedelta_range('1 day', periods=5, freq='h', name='other')
        expected3 = timedelta_range('1 day', periods=4, freq='h', name=None)

        rng4 = timedelta_range('1 day', periods=10, freq='h', name='idx')[5:]
        expected4 = TimedeltaIndex([], name='idx')

        for (rng, expected) in [(rng2, expected2), (rng3, expected3),
                                (rng4, expected4)]:
            result = base.intersection(rng)
            tm.assert_index_equal(result, expected)
            assert result.name == expected.name
            assert result.freq == expected.freq

    @pytest.mark.parametrize("sort", [None, False])
    def intersection_non_monotonic(self, sort):
        # non-monotonic
        base = TimedeltaIndex(['1 hour', '2 hour',
                              '4 hour', '3 hour'],
                              name='idx')

        rng2 = TimedeltaIndex(['5 hour', '2 hour',
                              '4 hour', '9 hour'],
                              name='idx')
        expected2 = TimedeltaIndex(['2 hour', '4 hour'],
                                   name='idx')

        rng3 = TimedeltaIndex(['2 hour', '5 hour',
                              '5 hour', '1 hour'],
                              name='other')
        expected3 = TimedeltaIndex(['1 hour', '2 hour'],
                                   name=None)

        rng4 = base[::-1]
        expected4 = base

        for (rng, expected) in [(rng2, expected2), (rng3, expected3),
                                (rng4, expected4)]:
            result = base.intersection(rng, sort=sort)
            if sort is None:
                expected = expected.sort_values()
            tm.assert_index_equal(result, expected)
            assert result.name == expected.name
            assert result.freq is None
