from __future__ import division

import pytest
import numpy as np

from pandas import Interval, IntervalIndex, Int64Index
from pandas.tests.indexes.common import Base
import pandas.util.testing as tm


pytestmark = pytest.mark.skip(reason="new indexing tests for issue 16316")


class TestIntervalIndex(Base):

    def _compare_tuple_of_numpy_array(self, result, expected):
        lidx, ridx = result
        lidx_expected, ridx_expected = expected

        tm.assert_numpy_array_equal(lidx, lidx_expected)
        tm.assert_numpy_array_equal(ridx, ridx_expected)

    @pytest.mark.parametrize("ivl_side", ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize("other_side", ['right', 'left', 'both', 'neither'])
    def test_interval_covers_interval(self, ivl_side, other_side):

        # class Interval:
        #     def covers(self, other: Interval) -> bool

        assert Interval(1, 3).covers(Interval(1.5, 2.5))
        assert Interval(1, 3).covers(Interval(1, 2))
        assert Interval(1, 3).covers(Interval(2, 3))
        assert not Interval(1, 3).covers(Interval(0.5, 2.5))
        assert not Interval(1, 3).covers(Interval(1.5, 3.5))

        ivl = Interval(1, 3, closed=ivl_side)
        other = Interval(1, 3, closed=other_side)

        should_cover = {
            'right': {
                'right': True, 'left': False, 'both': False, 'neither': True},
            'left': {
                'right': False, 'left': True, 'both': False, 'neither': True},
            'both': {
                'right': True, 'left': True, 'both': True, 'neither': True},
            'neither': {
                'right': False, 'left': False, 'both': False, 'neither': True}
        }

        assert ivl.covers(other) == should_cover[ivl_side][other_side]

    @pytest.mark.parametrize("ivl_side", ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize("other_side", ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize("ivl_range", [(1, 3), (-1, 1), (3, 5)])
    def test_interval_overlaps_interval(self, ivl_side, other_side, ivl_range):

        # class Interval:
        #     def overlaps(self, other: Interval) -> bool

        assert Interval(1, 3).overlaps(Interval(1.5, 2.5))
        assert Interval(1, 3).overlaps(Interval(1, 2))
        assert Interval(1, 3).overlaps(Interval(2, 3))
        assert Interval(1, 3).overlaps(Interval(0.5, 2.5))
        assert Interval(1, 3).overlaps(Interval(1.5, 3.5))

        assert not Interval(1, 3).overlaps(Interval(-1, 1))
        assert not Interval(1, 3).overlaps(Interval(3, 5))

        ivl = Interval(*ivl_range, closed=ivl_side)
        other = Interval(1, 3, closed=other_side)

        should_overlap = {
            # idx_side:
            #   ivl_side: {ivl_range: expected_result}
            'right': {
                'right': {(-1, 1): False, (1, 3): True, (3, 5): False},
                'left':{(-1, 1): False, (1, 3): True, (3, 5): True},
                'both': {(-1, 1): False, (1, 3): True, (3, 5): True},
                'neither': {(-1, 1): False, (1, 3): True, (3, 5): False},
            },
            'left': {
                'right': {(-1, 1): True, (1, 3): True, (3, 5): False},
                'left':{(-1, 1): False, (1, 3): True, (3, 5): False},
                'both': {(-1, 1): True, (1, 3): True, (3, 5): False},
                'neither': {(-1, 1): False, (1, 3): True, (3, 5): False},
            },
            'both': {
                'right': {(-1, 1): True, (1, 3): True, (3, 5): False},
                'left':{(-1, 1): False, (1, 3): True, (3, 5): True},
                'both': {(-1, 1): True, (1, 3): True, (3, 5): True},
                'neither': {(-1, 1): False, (1, 3): True, (3, 5): False},
            },
            'neither': {
                'right': {(-1, 1): False, (1, 3): True, (3, 5): False},
                'left':{(-1, 1): False, (1, 3): True, (3, 5): False},
                'both': {(-1, 1): False, (1, 3): True, (3, 5): False},
                'neither': {(-1, 1): False, (1, 3): True, (3, 5): False},
            }
        }

        assert ivl.overlaps(other) == should_overlap[other_side][ivl_side][ivl_range] == other.overlaps(ivl)

    @pytest.mark.parametrize("idx_side", ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize("ivl_side", ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize("ivl_range", [(1, 3), (0, 3), (0, 2), (2, 4)])
    def test_interval_covers_intervalIndex(self, idx_side, ivl_side, ivl_range):

        # class Interval:
        #     def covers(self, other: IntervalIndex) -> IntegerArray1D

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)], closed=idx_side)
        ivl = Interval(*ivl_range, closed=ivl_side)

        should_cover = {
            # idx_side:
            #   ivl_side: {ivl_range: expected_result}
            'right': {
                'right': {
                (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0], (2, 4): [1]},
                'left': {
                (1, 3): [], (0, 3): [0], (0, 2): [0], (2, 4): [1]},
                'both': {
                (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0], (2, 4): [1]},
                'neither': {
                (1, 3): [], (0, 3): [0], (0, 2): [0], (2, 4): [1]}
            },
            'left': {
                'right': {
                (1, 3): [1], (0, 3): [1, 2], (0, 2): [], (2, 4): []},
                'left': {
                (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0], (2, 4): [1]},
                'both': {
                (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0], (2, 4): [1]},
                'neither': {
                (1, 3): [1], (0, 3): [1, 2], (0, 2): [], (2, 4): []}
            },
            'both': {
                'right': {
                (1, 3): [1], (0, 3): [1, 2], (0, 2): [], (2, 4): []},
                'left': {
                (1, 3): [], (0, 3): [0], (0, 2): [0], (2, 4): [1]},
                'both': {
                (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0], (2, 4): [1]},
                'neither': {
                (1, 3): [], (0, 3): [], (0, 2): [], (2, 4): []}
            },
            'neither': {
                'right': {
                (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0], (2, 4): [1]},
                'left': {
                (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0], (2, 4): [1]},
                'both': {
                (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0], (2, 4): [1]},
                'neither': {
                (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0], (2, 4): [1]}
            }
        }

        tm.assert_numpy_array_equal(ivl.covers(idx),
                np.array(should_cover[idx_side][ivl_side][ivl_range]))
        tm.assert_numpy_array_equal(idx.covers(ivl),
                np.array(should_cover[idx_side][ivl_side][ivl_range]))

    @pytest.mark.parametrize("idx_side", ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize("ivl_side", ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize("ivl_range", [(1, 3), (1, 2), (0, 2), (3, 4)])
    def test_interval_overlaps_intervalIndex(self, idx_side, ivl_side, ivl_range):

        # class Interval:
        #     def overlaps(self, other: IntervalIndex) -> IntegerArray1D

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)], closed=idx_side)
        ivl = Interval(*ivl_range, closed=ivl_side)

        should_overlap = {
            # idx_side:
            #   ivl_side: {ivl_range: expected_result}
            'right': {
                'right': {
                (1, 3): [1, 2], (1, 2): [2], (0, 2): [0, 2], (3, 4): []},
                'left': {
                (1, 3): [0, 1, 2], (1, 2): [0, 2], (0, 2): [0, 2], (3, 4): [1, 2]},
                'both': {
                (1, 3): [0, 1, 2], (1, 2): [0, 2], (0, 2): [0, 2], (3, 4): [1, 2]},
                'neither': {
                (1, 3): [1, 2], (1, 2): [2], (0, 2): [0, 2], (3, 4): []},
            },
            'left': {
                'right': {
                (1, 3): [1, 2], (1, 2): [1, 2], (0, 2): [0, 1, 2], (3, 4): []},
                'left': {
                (1, 3): [1, 2], (1, 2): [2], (0, 2): [0, 2], (3, 4): []},
                'both': {
                (1, 3): [1, 2], (1, 2): [1, 2], (0, 2): [0, 1, 2], (3, 4): []},
                'neither': {
                (1, 3): [1, 2], (1, 2): [2], (0, 2): [0, 2], (3, 4): []},
            },
            'both': {
                'right': {
                (1, 3): [1, 2], (1, 2): [1, 2], (0, 2): [0, 1, 2], (3, 4): []},
                'left': {
                (1, 3): [0, 1, 2], (1, 2): [0, 2], (0, 2): [0, 2], (3, 4): [1, 2]},
                'both': {
                (1, 3): [0, 1, 2], (1, 2): [0, 1, 2], (0, 2): [0, 1, 2], (3, 4): [1, 2]},
                'neither': {
                (1, 3): [1, 2], (1, 2): [2], (0, 2): [0, 2], (3, 4): []},
            },
            'neither': {
                'right': {
                (1, 3): [1, 2], (1, 2): [2], (0, 2): [0, 2], (3, 4): []},
                'left': {
                (1, 3): [1, 2], (1, 2): [2], (0, 2): [0, 2], (3, 4): []},
                'both': {
                (1, 3): [1, 2], (1, 2): [2], (0, 2): [0, 2], (3, 4): []},
                'neither': {
                (1, 3): [1, 2], (1, 2): [2], (0, 2): [0, 2], (3, 4): []},
            }
        }

        tm.assert_numpy_array_equal(ivl.overlaps(idx),
                np.array(should_overlap[idx_side][ivl_side][ivl_range]))
        tm.assert_numpy_array_equal(idx.overlaps(ivl),
                np.array(should_overlap[idx_side][ivl_side][ivl_range]))

    def test_intervalIndex_covers_intervalIndex(self):

        # class IntervalIndex:
        #     def covers(self, other: IntervalIndex) -> Tuple[IntegerArray1D, IntegerArray1D]

        idx1 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)], closed='right')
        idx2 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)], closed='left')
        idx3 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)], closed='both')

        self._compare_tuple_of_numpy_array(idx.covers(idx1), (np.array([0,1,2,2]), np.array([0,1,1,2])))
        self._compare_tuple_of_numpy_array(idx.covers(idx2), (np.array([2]), np.array([1])))
        self._compare_tuple_of_numpy_array(idx.covers(idx3), (np.array([0,1,2,2]), np.array([0,1,1,2])))

    def test_intervalIndex_overlaps_intervalIndex(self):

        # class IntervalIndex:
        #     def overlaps(self, other: IntervalIndex) -> Tuple[IntegerArray1D, IntegerArray1D]

        idx1 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)], closed='right')
        idx2 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)], closed='left')
        idx3 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)], closed='both')

        self._compare_tuple_of_numpy_array(idx.overlaps(idx1), (np.array([0,1,2,2]), np.array([0,1,1,2])))
        self._compare_tuple_of_numpy_array(idx.overlaps(idx2), (np.array([0,0,1,1,2,2]), np.array([0,2,1,2,1,2])))
        self._compare_tuple_of_numpy_array(idx.overlaps(idx3), (np.array([0,0,1,1,2,2]), np.array([0,2,1,2,1,2])))
