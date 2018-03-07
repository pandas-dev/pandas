from __future__ import division

import pytest
import numpy as np

from pandas import Interval, IntervalIndex, Int64Index
import pandas.util.testing as tm


pytestmark = pytest.mark.skip(reason="new indexing tests for issue 16316")


class TestIntervalIndex(object):

    def _compare_tuple_of_numpy_array(self, result, expected):
        lidx, ridx = result
        lidx_expected, ridx_expected = expected

        tm.assert_numpy_array_equal(lidx, lidx_expected)
        tm.assert_numpy_array_equal(ridx, ridx_expected)

    def check_get_loc_result(self, idx, key, correct):
        """
        helper function to check the result for get_loc related tests

        Parameters
        ----------

        idx: IntervalIndex
            IntervalIndex over which get_loc checks will be run

        key: scalar or Interval

        correct: dict
            dictionary of the form {key: expected} denoting the expected loc's
            to be returned for each given key; keys not in the dictionary are
            assumed to raise a KeyError
        """
        expected = correct.get(key, None)
        if expected is None:
            # no expected in dict means KeyError
            with pytest.raises(KeyError):
                idx.get_loc(key)
        elif isinstance(expected, list):
            # multiple loc's returned -> numpy array
            expected = np.array(expected, dtype='int64')
            result = idx.get_loc(key)
            tm.assert_numpy_array_equal(result, expected)
        else:
            # otherwise single loc returned as an integer
            result = idx.get_loc(key)
            assert result == expected

    @pytest.mark.parametrize('ii_closed', ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize('iv_closed', ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize('iv_tuple', [
        (-2, -1), (-1, 0.5), (0, 1), (0.25, 0.75), (0, 2), (0, 3), (1, 3),
        (2, 3), (2.5, 4), (4, 5)], ids=str)
    def test_get_loc_interval(self, ii_closed, iv_closed, iv_tuple):
        """
        if get_loc is supplied an interval, it should only search for exact
        matches, not overlaps or covers, else KeyError
        """
        # construct interval to test against
        iv = Interval(iv_tuple[0], iv_tuple[1], iv_closed)

        # non-overlapping monotonic
        idx = IntervalIndex.from_tuples([(0, 1), (2, 3)], closed=ii_closed)
        correct = {Interval(0, 1, ii_closed): 0, Interval(2, 3, ii_closed): 1}
        self.check_get_loc_result(idx, iv, correct)

        # non-monotonic with dupes
        idx = IntervalIndex.from_tuples([(0, 1), (2, 3), (0, 1)],
                                        closed=ii_closed)
        correct = {Interval(0, 1, ii_closed): [0, 2],
                   Interval(2, 3, ii_closed): 1}
        self.check_get_loc_result(idx, iv, correct)

        # overlapping
        idx = IntervalIndex.from_tuples([(0, 2), (1, 3)], closed=ii_closed)
        correct = {Interval(0, 2, ii_closed): 0, Interval(1, 3, ii_closed): 1}
        self.check_get_loc_result(idx, iv, correct)

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
        Interval(1, 3, closed='right'),
        Interval(1, 3, closed='left'),
        Interval(1, 3, closed='both'),
        Interval(1, 3, closed='neither'),
        Interval(1, 4, closed='right'),
        Interval(0, 4, closed='right'),
        Interval(1, 2, closed='right')])
    @pytest.mark.parametrize("expected_result", [1, -1, -1, -1, -1, -1, -1])
    def test_get_indexer_with_interval_single_queries(
            self, query, expected_result):

        index = IntervalIndex.from_tuples(
            [(0, 2.5), (1, 3), (2, 4)], closed='right')

        result = index.get_indexer([query])
        expect = np.array([expected_result], dtype='intp')
        tm.assert_numpy_array_equal(result, expect)

    @pytest.mark.parametrize("query", [
        [Interval(2, 4, closed='right'), Interval(1, 3, closed='right')],
        [Interval(1, 3, closed='right'), Interval(0, 2, closed='right')],
        [Interval(1, 3, closed='right'), Interval(1, 3, closed='left')]])
    @pytest.mark.parametrize("expected_result", [[2, 1], [1, -1], [1, -1]])
    def test_get_indexer_with_interval_multiple_queries(
            self, query, expected_result):

        index = IntervalIndex.from_tuples(
            [(0, 2.5), (1, 3), (2, 4)], closed='right')

        result = index.get_indexer(query)
        expect = np.array(expected_result, dtype='intp')
        tm.assert_numpy_array_equal(result, expect)

    @pytest.mark.parametrize(
        "query",
        [-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5])
    @pytest.mark.parametrize(
        "expected_result",
        [-1, -1, 0, 0, 1, 1, -1, -1, 2, 2, -1])
    def test_get_indexer_with_ints_and_floats_single_queries(
            self, query, expected_result):

        index = IntervalIndex.from_tuples(
            [(0, 1), (1, 2), (3, 4)], closed='right')

        result = index.get_indexer([query])
        expect = np.array([expected_result], dtype='intp')
        tm.assert_numpy_array_equal(result, expect)

    @pytest.mark.parametrize(
        "query",
        [[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 3, 4, 2]])
    @pytest.mark.parametrize(
        "expected_result",
        [[0, 1], [0, 1, -1], [0, 1, -1, 2], [0, 1, -1, 2, 1]])
    def test_get_indexer_with_ints_and_floats_multiple_queries(
            self, query, expected_result):

        index = IntervalIndex.from_tuples(
            [(0, 1), (1, 2), (3, 4)], closed='right')

        result = index.get_indexer(query)
        expect = np.array(expected_result, dtype='intp')
        tm.assert_numpy_array_equal(result, expect)

        index = IntervalIndex.from_tuples([(0, 2), (1, 3), (2, 4)])
        # TODO: @shoyer believes this should raise, master branch doesn't

    @pytest.mark.parametrize(
        "query",
        [-0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5])
    @pytest.mark.parametrize("expected_result", [
        (Int64Index([], dtype='int64'), np.array([0])),
        (Int64Index([0], dtype='int64'), np.array([])),
        (Int64Index([0], dtype='int64'), np.array([])),
        (Int64Index([0, 1], dtype='int64'), np.array([])),
        (Int64Index([0, 1], dtype='int64'), np.array([])),
        (Int64Index([0, 1, 2], dtype='int64'), np.array([])),
        (Int64Index([1, 2], dtype='int64'), np.array([])),
        (Int64Index([2], dtype='int64'), np.array([])),
        (Int64Index([2], dtype='int64'), np.array([])),
        (Int64Index([], dtype='int64'), np.array([0])),
        (Int64Index([], dtype='int64'), np.array([0]))])
    def test_get_indexer_non_unique_with_ints_and_floats_single_queries(
            self, query, expected_result):

        index = IntervalIndex.from_tuples(
            [(0, 2.5), (1, 3), (2, 4)], closed='left')

        result = index.get_indexer_non_unique([query])
        tm.assert_numpy_array_equal(result, expected_result)

    @pytest.mark.parametrize(
        "query",
        [[1, 2], [1, 2, 3], [1, 2, 3, 4], [1, 2, 3, 4, 2]])
    @pytest.mark.parametrize("expected_result", [
        (Int64Index([0, 1, 0, 1, 2], dtype='int64'), np.array([])),
        (Int64Index([0, 1, 0, 1, 2, 2], dtype='int64'), np.array([])),
        (Int64Index([0, 1, 0, 1, 2, 2, -1], dtype='int64'), np.array([3])),
        (Int64Index([0, 1, 0, 1, 2, 2, -1, 0, 1, 2], dtype='int64'),
            np.array([3]))])
    def test_get_indexer_non_unique_with_ints_and_floats_multiple_queries(
            self, query, expected_result):

        index = IntervalIndex.from_tuples(
            [(0, 2.5), (1, 3), (2, 4)], closed='left')

        result = index.get_indexer_non_unique(query)
        tm.assert_numpy_array_equal(result, expected_result)

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
