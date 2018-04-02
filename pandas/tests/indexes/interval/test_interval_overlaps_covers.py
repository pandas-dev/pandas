from __future__ import division

import pytest
import numpy as np

from pandas import IntervalIndex
from pandas.tests.indexes.common import Base
import pandas.util.testing as tm

pytestmark = pytest.mark.skip(reason="new indexing tests for issue 16316")


class TestIntervalIndex(Base):

    def _compare_tuple_of_numpy_array(self, result, expected):
        lidx, ridx = result
        lidx_expected, ridx_expected = expected

        tm.assert_numpy_array_equal(lidx, lidx_expected)
        tm.assert_numpy_array_equal(ridx, ridx_expected)

    @pytest.mark.parametrize("oth_side", ['right', 'left', 'both'])
    def test_intervalIndex_covers_intervalIndex(self, idx):

        # class IntervalIndex:
        #     def covers(self, other: IntervalIndex) -> Tuple[IntegerArray1D,
        #                                                       IntegerArray1D]

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)],
                                        closed="right")
        other = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)],
                                        closed=oth_side)

        result = idx.covers(idx)
        expected = {
            "right": (np.array([0, 1, 2, 2]), np.array([0, 1, 1, 2])),
            "left": (np.array([2]), np.array([1])),
            "both": (np.array([0, 1, 2, 2]), np.array([0, 1, 1, 2]))
        }

        self._compare_tuple_of_numpy_array(result, expected[oth_side])

    @pytest.mark.parametrize("oth_side", ['right', 'left', 'both'])
    def test_intervalIndex_overlaps_intervalIndex(self):

        # class IntervalIndex:
        #     def overlaps(self, other: IntervalIndex) -> Tuple[IntegerArray1D,
        #                                                       IntegerArray1D]

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)],
                                        closed="right")
        other = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)],
                                        closed=oth_side)

        result = idx.overlaps(idx)
        expected = {
            "right": (np.array([0, 1, 2, 2]), np.array([0, 1, 1, 2])),
            "left": (np.array([0, 0, 1, 1, 2, 2]),
                     np.array([0, 2, 1, 2, 1, 2])),
            "both": (np.array([0, 0, 1, 1, 2, 2]),
                     np.array([0, 2, 1, 2, 1, 2]))
        }
        self._compare_tuple_of_numpy_array(result, expected[oth_side])
