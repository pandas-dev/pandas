from __future__ import division

import pytest
import numpy as np

from pandas import Interval, IntervalIndex, Int64Index
import pandas.util.testing as tm


pytestmark = pytest.mark.skip(reason="new indexing tests for issue 16316")


class TestIntervalIndex(object):

    @pytest.mark.parametrize("side", ['right', 'left', 'both', 'neither'])
    def test_get_loc_interval(self, closed, side):

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3)], closed=closed)

        for bound in [[0, 1], [1, 2], [2, 3], [3, 4],
                      [0, 2], [2.5, 3], [-1, 4]]:
            # if get_loc is supplied an interval, it should only search
            # for exact matches, not overlaps or covers, else KeyError.
            if closed == side:
                if bound == [0, 1]:
                    assert idx.get_loc(Interval(0, 1, closed=side)) == 0
                elif bound == [2, 3]:
                    assert idx.get_loc(Interval(2, 3, closed=side)) == 1
                else:
                    with pytest.raises(KeyError):
                        idx.get_loc(Interval(*bound, closed=side))
            else:
                with pytest.raises(KeyError):
                    idx.get_loc(Interval(*bound, closed=side))

    @pytest.mark.parametrize("scalar", [-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5])
    def test_get_loc_scalar(self, closed, scalar):

        # correct = {side: {query: answer}}.
        # If query is not in the dict, that query should raise a KeyError
        correct = {'right': {0.5: 0, 1: 0, 2.5: 1, 3: 1},
                   'left': {0: 0, 0.5: 0, 2: 1, 2.5: 1},
                   'both': {0: 0, 0.5: 0, 1: 0, 2: 1, 2.5: 1, 3: 1},
                   'neither': {0.5: 0, 2.5: 1}}

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3)], closed=closed)

        # if get_loc is supplied a scalar, it should return the index of
        # the interval which contains the scalar, or KeyError.
        if scalar in correct[closed].keys():
            assert idx.get_loc(scalar) == correct[closed][scalar]
        else:
            pytest.raises(KeyError, idx.get_loc, scalar)

    def test_slice_locs_with_interval(self):

        # increasing monotonically
        index = IntervalIndex.from_tuples([(0, 2), (1, 3), (2, 4)])

        assert index.slice_locs(
            start=Interval(0, 2), end=Interval(2, 4)) == (0, 3)
        assert index.slice_locs(start=Interval(0, 2)) == (0, 3)
        assert index.slice_locs(end=Interval(2, 4)) == (0, 3)
        assert index.slice_locs(end=Interval(0, 2)) == (0, 1)
        assert index.slice_locs(
            start=Interval(2, 4), end=Interval(0, 2)) == (2, 1)

        # decreasing monotonically
        index = IntervalIndex.from_tuples([(2, 4), (1, 3), (0, 2)])

        assert index.slice_locs(
            start=Interval(0, 2), end=Interval(2, 4)) == (2, 1)
        assert index.slice_locs(start=Interval(0, 2)) == (2, 3)
        assert index.slice_locs(end=Interval(2, 4)) == (0, 1)
        assert index.slice_locs(end=Interval(0, 2)) == (0, 3)
        assert index.slice_locs(
            start=Interval(2, 4), end=Interval(0, 2)) == (0, 3)

        # sorted duplicates
        index = IntervalIndex.from_tuples([(0, 2), (0, 2), (2, 4)])

        assert index.slice_locs(
            start=Interval(0, 2), end=Interval(2, 4)) == (0, 3)
        assert index.slice_locs(start=Interval(0, 2)) == (0, 3)
        assert index.slice_locs(end=Interval(2, 4)) == (0, 3)
        assert index.slice_locs(end=Interval(0, 2)) == (0, 2)
        assert index.slice_locs(
            start=Interval(2, 4), end=Interval(0, 2)) == (2, 2)

        # unsorted duplicates
        index = IntervalIndex.from_tuples([(0, 2), (2, 4), (0, 2)])

        pytest.raises(KeyError, index.slice_locs(
            start=Interval(0, 2), end=Interval(2, 4)))
        pytest.raises(KeyError, index.slice_locs(start=Interval(0, 2)))
        assert index.slice_locs(end=Interval(2, 4)) == (0, 2)
        pytest.raises(KeyError, index.slice_locs(end=Interval(0, 2)))
        pytest.raises(KeyError, index.slice_locs(
            start=Interval(2, 4), end=Interval(0, 2)))

        # another unsorted duplicates
        index = IntervalIndex.from_tuples([(0, 2), (0, 2), (2, 4), (1, 3)])

        assert index.slice_locs(
            start=Interval(0, 2), end=Interval(2, 4)) == (0, 3)
        assert index.slice_locs(start=Interval(0, 2)) == (0, 4)
        assert index.slice_locs(end=Interval(2, 4)) == (0, 3)
        assert index.slice_locs(end=Interval(0, 2)) == (0, 2)
        assert index.slice_locs(
            start=Interval(2, 4), end=Interval(0, 2)) == (2, 2)

    def test_slice_locs_with_ints_and_floats_succeeds(self):

        # increasing non-overlapping
        index = IntervalIndex.from_tuples([(0, 1), (1, 2), (3, 4)])

        assert index.slice_locs(0, 1) == (0, 1)
        assert index.slice_locs(0, 2) == (0, 2)
        assert index.slice_locs(0, 3) == (0, 2)
        assert index.slice_locs(3, 1) == (2, 1)
        assert index.slice_locs(3, 4) == (2, 3)
        assert index.slice_locs(0, 4) == (0, 3)

        # decreasing non-overlapping
        index = IntervalIndex.from_tuples([(3, 4), (1, 2), (0, 1)])
        assert index.slice_locs(0, 1) == (3, 2)
        assert index.slice_locs(0, 2) == (3, 1)
        assert index.slice_locs(0, 3) == (3, 1)
        assert index.slice_locs(3, 1) == (1, 2)
        assert index.slice_locs(3, 4) == (1, 0)
        assert index.slice_locs(0, 4) == (3, 0)

    @pytest.mark.parametrize("query", [
        [0, 1], [0, 2], [0, 3], [3, 1], [3, 4], [0, 4]])
    @pytest.mark.parametrize("tuples", [
        [(0, 2), (1, 3), (2, 4)], [(2, 4), (1, 3), (0, 2)],
        [(0, 2), (0, 2), (2, 4)], [(0, 2), (2, 4), (0, 2)],
        [(0, 2), (0, 2), (2, 4), (1, 3)]])
    def test_slice_locs_with_ints_and_floats_errors(self, tuples, query):
        index = IntervalIndex.from_tuples(tuples)
        with pytest.raises(KeyError):
            index.slice_locs(query)

    @pytest.mark.parametrize('query, expected', [
        ([Interval(1, 3, closed='right')], [1]),
        ([Interval(1, 3, closed='left')], [-1]),
        ([Interval(1, 3, closed='both')], [-1]),
        ([Interval(1, 3, closed='neither')], [-1]),
        ([Interval(1, 4, closed='right')], [-1]),
        ([Interval(0, 4, closed='right')], [-1]),
        ([Interval(1, 2, closed='right')], [-1]),
        ([Interval(2, 4, closed='right'), Interval(1, 3, closed='right')],
         [2, 1]),
        ([Interval(1, 3, closed='right'), Interval(0, 2, closed='right')],
         [1, -1]),
        ([Interval(1, 3, closed='right'), Interval(1, 3, closed='left')],
         [1, -1])])
    def test_get_indexer_with_interval(self, query, expected):

        tuples = [(0, 2.5), (1, 3), (2, 4)]
        index = IntervalIndex.from_tuples(tuples, closed='right')

        result = index.get_indexer(query)
        expected = np.array(expected, dtype='intp')
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize('query, expected', [
        ([-0.5], [-1]),
        ([0], [-1]),
        ([0.5], [0]),
        ([1], [0]),
        ([1.5], [1]),
        ([2], [1]),
        ([2.5], [-1]),
        ([3], [-1]),
        ([3.5], [2]),
        ([4], [2]),
        ([4.5], [-1]),
        ([1, 2], [0, 1]),
        ([1, 2, 3], [0, 1, -1]),
        ([1, 2, 3, 4], [0, 1, -1, 2]),
        ([1, 2, 3, 4, 2], [0, 1, -1, 2, 1])])
    def test_get_indexer_with_int_and_float(self, query, expected):

        tuples = [(0, 1), (1, 2), (3, 4)]
        index = IntervalIndex.from_tuples(tuples, closed='right')

        result = index.get_indexer(query)
        expected = np.array(expected, dtype='intp')
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize('tuples, closed', [
        ([(0, 2), (1, 3), (3, 4)], 'neither'),
        ([(0, 5), (1, 4), (6, 7)], 'left'),
        ([(0, 1), (0, 1), (1, 2)], 'right'),
        ([(0, 1), (2, 3), (3, 4)], 'both')])
    def test_get_indexer_errors(self, tuples, closed):
        # IntervalIndex needs non-overlapping for uniqueness when querying
        index = IntervalIndex.from_tuples(tuples, closed=closed)

        msg = ('cannot handle overlapping indices; use '
               'IntervalIndex.get_indexer_non_unique')
        with tm.assert_raises_regex(ValueError, msg):
            index.get_indexer([0, 2])

    @pytest.mark.parametrize('query, expected', [
        ([-0.5], ([-1], [0])),
        ([0], ([0], [])),
        ([0.5], ([0], [])),
        ([1], ([0, 1], [])),
        ([1.5], ([0, 1], [])),
        ([2], ([0, 1, 2], [])),
        ([2.5], ([1, 2], [])),
        ([3], ([2], [])),
        ([3.5], ([2], [])),
        ([4], ([-1], [0])),
        ([4.5], ([-1], [0])),
        ([1, 2], ([0, 1, 0, 1, 2], [])),
        ([1, 2, 3], ([0, 1, 0, 1, 2, 2], [])),
        ([1, 2, 3, 4], ([0, 1, 0, 1, 2, 2, -1], [3])),
        ([1, 2, 3, 4, 2], ([0, 1, 0, 1, 2, 2, -1, 0, 1, 2], [3]))])
    def test_get_indexer_non_unique_with_int_and_float(self, query, expected):

        tuples = [(0, 2.5), (1, 3), (2, 4)]
        index = IntervalIndex.from_tuples(tuples, closed='left')

        result_indexer, result_missing = index.get_indexer_non_unique(query)
        expected_indexer = Int64Index(expected[0])
        expected_missing = np.array(expected[1], dtype='intp')

        tm.assert_index_equal(result_indexer, expected_indexer)
        tm.assert_numpy_array_equal(result_missing, expected_missing)

        # TODO we may also want to test get_indexer for the case when
        # the intervals are duplicated, decreasing, non-monotonic, etc..

    def test_contains(self):

        index = IntervalIndex.from_arrays([0, 1], [1, 2], closed='right')

        # __contains__ requires perfect matches to intervals.
        assert 0 not in index
        assert 1 not in index
        assert 2 not in index

        assert Interval(0, 1, closed='right') in index
        assert Interval(0, 2, closed='right') not in index
        assert Interval(0, 0.5, closed='right') not in index
        assert Interval(3, 5, closed='right') not in index
        assert Interval(-1, 0, closed='left') not in index
        assert Interval(0, 1, closed='left') not in index
        assert Interval(0, 1, closed='both') not in index

    def test_contains_method(self):

        index = IntervalIndex.from_arrays([0, 1], [1, 2], closed='right')

        assert not index.contains(0)
        assert index.contains(0.1)
        assert index.contains(0.5)
        assert index.contains(1)

        assert index.contains(Interval(0, 1, closed='right'))
        assert not index.contains(Interval(0, 1, closed='left'))
        assert not index.contains(Interval(0, 1, closed='both'))
        assert not index.contains(Interval(0, 2, closed='right'))

        assert not index.contains(Interval(0, 3, closed='right'))
        assert not index.contains(Interval(1, 3, closed='right'))

        assert not index.contains(20)
        assert not index.contains(-20)
