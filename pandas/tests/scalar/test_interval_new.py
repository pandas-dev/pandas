from __future__ import division

import pytest
import numpy as np

from pandas import Interval, IntervalIndex
from pandas.tests.indexes.common import Base
import pandas.util.testing as tm

pytestmark = pytest.mark.skip(reason="new indexing tests for issue 16316")


class TestIntervalIndex(Base):

    @pytest.mark.parametrize("ivl_side", ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize("oth_side", ['right', 'left', 'both', 'neither'])
    def test_interval_covers_interval(self, ivl_side, oth_side):

        # class Interval:
        #     def covers(self, other: Interval) -> bool

        assert Interval(1, 3).covers(Interval(1.5, 2.5))
        assert Interval(1, 3).covers(Interval(1, 2))
        assert Interval(1, 3).covers(Interval(2, 3))
        assert not Interval(1, 3).covers(Interval(0.5, 2.5))
        assert not Interval(1, 3).covers(Interval(1.5, 3.5))

        ivl = Interval(1, 3, closed=ivl_side)
        other = Interval(1, 3, closed=oth_side)

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

        result = ivl.covers(other)
        expected = should_cover[ivl_side][oth_side]
        assert result == expected

    @pytest.mark.parametrize("ivl_side", ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize("oth_side", ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize("ivl_range", [(1, 3), (-1, 1), (3, 5)])
    def test_interval_overlaps_interval(self, ivl_side, oth_side, ivl_range):

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
        other = Interval(1, 3, closed=oth_side)

        should_overlap = {
            # idx_side:
            #   ivl_side: {ivl_range: expected_result}
            'right': {
                'right': {(-1, 1): False, (1, 3): True, (3, 5): False},
                'left': {(-1, 1): False, (1, 3): True, (3, 5): True},
                'both': {(-1, 1): False, (1, 3): True, (3, 5): True},
                'neither': {(-1, 1): False, (1, 3): True, (3, 5): False},
            },
            'left': {
                'right': {(-1, 1): True, (1, 3): True, (3, 5): False},
                'left': {(-1, 1): False, (1, 3): True, (3, 5): False},
                'both': {(-1, 1): True, (1, 3): True, (3, 5): False},
                'neither': {(-1, 1): False, (1, 3): True, (3, 5): False},
            },
            'both': {
                'right': {(-1, 1): True, (1, 3): True, (3, 5): False},
                'left': {(-1, 1): False, (1, 3): True, (3, 5): True},
                'both': {(-1, 1): True, (1, 3): True, (3, 5): True},
                'neither': {(-1, 1): False, (1, 3): True, (3, 5): False},
            },
            'neither': {
                'right': {(-1, 1): False, (1, 3): True, (3, 5): False},
                'left': {(-1, 1): False, (1, 3): True, (3, 5): False},
                'both': {(-1, 1): False, (1, 3): True, (3, 5): False},
                'neither': {(-1, 1): False, (1, 3): True, (3, 5): False},
            }
        }

        result = ivl.overlaps(other)
        expected = should_overlap[oth_side][ivl_side][ivl_range]
        other_result = other.overlaps(ivl)

        assert result == expected == other_result

    @pytest.mark.parametrize("idx_side", ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize("ivl_side", ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize("ivl_range", [(1, 3), (0, 3), (0, 2), (2, 4)])
    def test_interval_covers_intervalIndex(self, idx_side, ivl_side,
                                           ivl_range):

        # class Interval:
        #     def covers(self, other: IntervalIndex) -> IntegerArray1D

        # class IntervalIndex:
        #     def covers(self, other: Interval) -> IntegerArray1D

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)],
                                        closed=idx_side)
        ivl = Interval(*ivl_range, closed=ivl_side)

        should_cover = {
            # idx_side:
            #   ivl_side: {ivl_range: expected_result}
            'right': {
                'right': {
                    (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0],
                    (2, 4): [1]},
                'left': {
                    (1, 3): [], (0, 3): [0], (0, 2): [0], (2, 4): [1]},
                'both': {
                    (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0],
                    (2, 4): [1]},
                'neither': {
                    (1, 3): [], (0, 3): [0], (0, 2): [0], (2, 4): [1]}
            },
            'left': {
                'right': {
                    (1, 3): [1], (0, 3): [1, 2], (0, 2): [], (2, 4): []},
                'left': {
                    (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0],
                    (2, 4): [1]},
                'both': {
                    (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0],
                    (2, 4): [1]},
                'neither': {
                    (1, 3): [1], (0, 3): [1, 2], (0, 2): [], (2, 4): []}
            },
            'both': {
                'right': {
                    (1, 3): [1], (0, 3): [1, 2], (0, 2): [], (2, 4): []},
                'left': {
                    (1, 3): [], (0, 3): [0], (0, 2): [0], (2, 4): [1]},
                'both': {
                    (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0],
                    (2, 4): [1]},
                'neither': {
                    (1, 3): [], (0, 3): [], (0, 2): [], (2, 4): []}
            },
            'neither': {
                'right': {
                    (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0],
                    (2, 4): [1]},
                'left': {
                    (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0],
                    (2, 4): [1]},
                'both': {
                    (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0],
                    (2, 4): [1]},
                'neither': {
                    (1, 3): [1, 2], (0, 3): [0, 1, 2], (0, 2): [0],
                    (2, 4): [1]}
            }
        }

        result = ivl.covers(idx)
        expected = np.array(should_cover[idx_side][ivl_side][ivl_range])
        other_result = idx.covers(ivl)

        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(other_result, expected)

    @pytest.mark.parametrize("idx_side", ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize("ivl_side", ['right', 'left', 'both', 'neither'])
    @pytest.mark.parametrize("ivl_range", [(1, 3), (1, 2), (0, 2), (3, 4)])
    def test_interval_overlaps_intervalIndex(self, idx_side, ivl_side,
                                             ivl_range):

        # class Interval:
        #     def overlaps(self, other: IntervalIndex) -> IntegerArray1D

        # class IntervalIndex:
        #     def overlaps(self, other: Interval) -> IntegerArray1D

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)],
                                        closed=idx_side)
        ivl = Interval(*ivl_range, closed=ivl_side)

        should_overlap = {
            # idx_side:
            #   ivl_side: {ivl_range: expected_result}
            'right': {
                'right': {
                    (1, 3): [1, 2], (1, 2): [2], (0, 2): [0, 2], (3, 4): []},
                'left': {
                    (1, 3): [0, 1, 2], (1, 2): [0, 2], (0, 2): [0, 2],
                    (3, 4): [1, 2]},
                'both': {
                    (1, 3): [0, 1, 2], (1, 2): [0, 2], (0, 2): [0, 2],
                    (3, 4): [1, 2]},
                'neither': {
                    (1, 3): [1, 2], (1, 2): [2], (0, 2): [0, 2], (3, 4): []},
            },
            'left': {
                'right': {
                    (1, 3): [1, 2], (1, 2): [1, 2], (0, 2): [0, 1, 2],
                    (3, 4): []},
                'left': {
                    (1, 3): [1, 2], (1, 2): [2], (0, 2): [0, 2], (3, 4): []},
                'both': {
                    (1, 3): [1, 2], (1, 2): [1, 2], (0, 2): [0, 1, 2],
                    (3, 4): []},
                'neither': {
                    (1, 3): [1, 2], (1, 2): [2], (0, 2): [0, 2], (3, 4): []},
            },
            'both': {
                'right': {
                    (1, 3): [1, 2], (1, 2): [1, 2], (0, 2): [0, 1, 2],
                    (3, 4): []},
                'left': {
                    (1, 3): [0, 1, 2], (1, 2): [0, 2], (0, 2): [0, 2],
                    (3, 4): [1, 2]},
                'both': {
                    (1, 3): [0, 1, 2], (1, 2): [0, 1, 2], (0, 2): [0, 1, 2],
                    (3, 4): [1, 2]},
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

        result = ivl.overlaps(idx)
        expected = np.array(should_overlap[idx_side][ivl_side][ivl_range])
        other_result = idx.overlaps(ivl)

        tm.assert_numpy_array_equal(result, expected)
        tm.assert_numpy_array_equal(other_result, expected)
