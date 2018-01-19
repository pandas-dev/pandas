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

    def test_intervalIndex_covers_intervalIndex(self):

        # class IntervalIndex:
        #     def covers(self, other: IntervalIndex) -> Tuple[IntegerArray1D,
        #                                                       IntegerArray1D]

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)],
                                        closed="right")

        result = idx.covers(idx)
        expected = (np.array([0, 1, 2, 2]), np.array([0, 1, 1, 2]))
        self._compare_tuple_of_numpy_array(result, expected)

        idx2 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)],
                                        closed="left")

        result = idx.covers(idx2),
        expected = (np.array([2]), np.array([1]))
        self._compare_tuple_of_numpy_array(result, expected)

        idx3 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)],
                                        closed="both")

        result = idx.covers(idx3),
        expected = (np.array([0, 1, 2, 2]), np.array([0, 1, 1, 2]))
        self._compare_tuple_of_numpy_array(result, expected)

    def test_intervalIndex_overlaps_intervalIndex(self):

        # class IntervalIndex:
        #     def overlaps(self, other: IntervalIndex) -> Tuple[IntegerArray1D,
        #                                                       IntegerArray1D]

        idx = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)],
                                        closed="right")

        result = idx.overlaps(idx)
        expected = (np.array([0, 1, 2, 2]), np.array([0, 1, 1, 2]))
        self._compare_tuple_of_numpy_array(result, expected)

        idx2 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)],
                                        closed="left")

        result = idx.overlaps(idx2)
        expected = (np.array([0, 0, 1, 1, 2, 2]), np.array([0, 2, 1, 2, 1, 2]))
        self._compare_tuple_of_numpy_array(result, expected)

        idx3 = IntervalIndex.from_tuples([(0, 1), (2, 3), (1, 3)],
                                        closed="both")

        result = idx.overlaps(idx3)
        expected = (np.array([0, 0, 1, 1, 2, 2]), np.array([0, 2, 1, 2, 1, 2]))
        self._compare_tuple_of_numpy_array(result, expected)
