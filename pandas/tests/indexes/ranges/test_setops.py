from datetime import datetime, timedelta

import numpy as np
import pytest

from pandas import Index, Int64Index, RangeIndex
import pandas._testing as tm


class TestRangeIndexSetOps:
    @pytest.mark.parametrize("sort", [None, False])
    def test_intersection(self, sort):
        # intersect with Int64Index
        index = RangeIndex(start=0, stop=20, step=2)
        other = Index(np.arange(1, 6))
        result = index.intersection(other, sort=sort)
        expected = Index(np.sort(np.intersect1d(index.values, other.values)))
        tm.assert_index_equal(result, expected)

        result = other.intersection(index, sort=sort)
        expected = Index(
            np.sort(np.asarray(np.intersect1d(index.values, other.values)))
        )
        tm.assert_index_equal(result, expected)

        # intersect with increasing RangeIndex
        other = RangeIndex(1, 6)
        result = index.intersection(other, sort=sort)
        expected = Index(np.sort(np.intersect1d(index.values, other.values)))
        tm.assert_index_equal(result, expected)

        # intersect with decreasing RangeIndex
        other = RangeIndex(5, 0, -1)
        result = index.intersection(other, sort=sort)
        expected = Index(np.sort(np.intersect1d(index.values, other.values)))
        tm.assert_index_equal(result, expected)

        # reversed (GH 17296)
        result = other.intersection(index, sort=sort)
        tm.assert_index_equal(result, expected)

        # GH 17296: intersect two decreasing RangeIndexes
        first = RangeIndex(10, -2, -2)
        other = RangeIndex(5, -4, -1)
        expected = first.astype(int).intersection(other.astype(int), sort=sort)
        result = first.intersection(other, sort=sort).astype(int)
        tm.assert_index_equal(result, expected)

        # reversed
        result = other.intersection(first, sort=sort).astype(int)
        tm.assert_index_equal(result, expected)

        index = RangeIndex(5)

        # intersect of non-overlapping indices
        other = RangeIndex(5, 10, 1)
        result = index.intersection(other, sort=sort)
        expected = RangeIndex(0, 0, 1)
        tm.assert_index_equal(result, expected)

        other = RangeIndex(-1, -5, -1)
        result = index.intersection(other, sort=sort)
        expected = RangeIndex(0, 0, 1)
        tm.assert_index_equal(result, expected)

        # intersection of empty indices
        other = RangeIndex(0, 0, 1)
        result = index.intersection(other, sort=sort)
        expected = RangeIndex(0, 0, 1)
        tm.assert_index_equal(result, expected)

        result = other.intersection(index, sort=sort)
        tm.assert_index_equal(result, expected)

        # intersection of non-overlapping values based on start value and gcd
        index = RangeIndex(1, 10, 2)
        other = RangeIndex(0, 10, 4)
        result = index.intersection(other, sort=sort)
        expected = RangeIndex(0, 0, 1)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("sort", [False, None])
    def test_union_noncomparable(self, sort):
        # corner case, non-Int64Index
        index = RangeIndex(start=0, stop=20, step=2)
        other = Index([datetime.now() + timedelta(i) for i in range(4)], dtype=object)
        result = index.union(other, sort=sort)
        expected = Index(np.concatenate((index, other)))
        tm.assert_index_equal(result, expected)

        result = other.union(index, sort=sort)
        expected = Index(np.concatenate((other, index)))
        tm.assert_index_equal(result, expected)

    @pytest.fixture(
        params=[
            (
                RangeIndex(0, 10, 1),
                RangeIndex(0, 10, 1),
                RangeIndex(0, 10, 1),
                RangeIndex(0, 10, 1),
            ),
            (
                RangeIndex(0, 10, 1),
                RangeIndex(5, 20, 1),
                RangeIndex(0, 20, 1),
                Int64Index(range(20)),
            ),
            (
                RangeIndex(0, 10, 1),
                RangeIndex(10, 20, 1),
                RangeIndex(0, 20, 1),
                Int64Index(range(20)),
            ),
            (
                RangeIndex(0, -10, -1),
                RangeIndex(0, -10, -1),
                RangeIndex(0, -10, -1),
                RangeIndex(0, -10, -1),
            ),
            (
                RangeIndex(0, -10, -1),
                RangeIndex(-10, -20, -1),
                RangeIndex(-19, 1, 1),
                Int64Index(range(0, -20, -1)),
            ),
            (
                RangeIndex(0, 10, 2),
                RangeIndex(1, 10, 2),
                RangeIndex(0, 10, 1),
                Int64Index(list(range(0, 10, 2)) + list(range(1, 10, 2))),
            ),
            (
                RangeIndex(0, 11, 2),
                RangeIndex(1, 12, 2),
                RangeIndex(0, 12, 1),
                Int64Index(list(range(0, 11, 2)) + list(range(1, 12, 2))),
            ),
            (
                RangeIndex(0, 21, 4),
                RangeIndex(-2, 24, 4),
                RangeIndex(-2, 24, 2),
                Int64Index(list(range(0, 21, 4)) + list(range(-2, 24, 4))),
            ),
            (
                RangeIndex(0, -20, -2),
                RangeIndex(-1, -21, -2),
                RangeIndex(-19, 1, 1),
                Int64Index(list(range(0, -20, -2)) + list(range(-1, -21, -2))),
            ),
            (
                RangeIndex(0, 100, 5),
                RangeIndex(0, 100, 20),
                RangeIndex(0, 100, 5),
                Int64Index(range(0, 100, 5)),
            ),
            (
                RangeIndex(0, -100, -5),
                RangeIndex(5, -100, -20),
                RangeIndex(-95, 10, 5),
                Int64Index(list(range(0, -100, -5)) + [5]),
            ),
            (
                RangeIndex(0, -11, -1),
                RangeIndex(1, -12, -4),
                RangeIndex(-11, 2, 1),
                Int64Index(list(range(0, -11, -1)) + [1, -11]),
            ),
            (RangeIndex(0), RangeIndex(0), RangeIndex(0), RangeIndex(0)),
            (
                RangeIndex(0, -10, -2),
                RangeIndex(0),
                RangeIndex(0, -10, -2),
                RangeIndex(0, -10, -2),
            ),
            (
                RangeIndex(0, 100, 2),
                RangeIndex(100, 150, 200),
                RangeIndex(0, 102, 2),
                Int64Index(range(0, 102, 2)),
            ),
            (
                RangeIndex(0, -100, -2),
                RangeIndex(-100, 50, 102),
                RangeIndex(-100, 4, 2),
                Int64Index(list(range(0, -100, -2)) + [-100, 2]),
            ),
            (
                RangeIndex(0, -100, -1),
                RangeIndex(0, -50, -3),
                RangeIndex(-99, 1, 1),
                Int64Index(list(range(0, -100, -1))),
            ),
            (
                RangeIndex(0, 1, 1),
                RangeIndex(5, 6, 10),
                RangeIndex(0, 6, 5),
                Int64Index([0, 5]),
            ),
            (
                RangeIndex(0, 10, 5),
                RangeIndex(-5, -6, -20),
                RangeIndex(-5, 10, 5),
                Int64Index([0, 5, -5]),
            ),
            (
                RangeIndex(0, 3, 1),
                RangeIndex(4, 5, 1),
                Int64Index([0, 1, 2, 4]),
                Int64Index([0, 1, 2, 4]),
            ),
            (
                RangeIndex(0, 10, 1),
                Int64Index([]),
                RangeIndex(0, 10, 1),
                RangeIndex(0, 10, 1),
            ),
            (
                RangeIndex(0),
                Int64Index([1, 5, 6]),
                Int64Index([1, 5, 6]),
                Int64Index([1, 5, 6]),
            ),
        ]
    )
    def unions(self, request):
        """Inputs and expected outputs for RangeIndex.union tests"""
        return request.param

    def test_union_sorted(self, unions):

        idx1, idx2, expected_sorted, expected_notsorted = unions

        res1 = idx1.union(idx2, sort=None)
        tm.assert_index_equal(res1, expected_sorted, exact=True)

        res1 = idx1.union(idx2, sort=False)
        tm.assert_index_equal(res1, expected_notsorted, exact=True)

        res2 = idx2.union(idx1, sort=None)
        res3 = idx1._int64index.union(idx2, sort=None)
        tm.assert_index_equal(res2, expected_sorted, exact=True)
        tm.assert_index_equal(res3, expected_sorted)
