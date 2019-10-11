import numpy as np
import pytest

from pandas import IntervalIndex, Series, date_range
from pandas.tests.indexes.common import Base
import pandas.util.testing as tm


class TestBase(Base):
    """
    Tests specific to the shared common index tests; unrelated tests should be placed
    in test_interval.py or the specific test file (e.g. test_astype.py)
    """

    _holder = IntervalIndex

    @pytest.fixture
    def indices(self):
        return tm.makeIntervalIndex(10)

    def create_index(self, closed="right"):
        return IntervalIndex.from_breaks(range(11), closed=closed)

    def test_equals(self, closed):
        expected = IntervalIndex.from_breaks(np.arange(5), closed=closed)
        assert expected.equals(expected)
        assert expected.equals(expected.copy())

        assert not expected.equals(expected.astype(object))
        assert not expected.equals(np.array(expected))
        assert not expected.equals(list(expected))

        assert not expected.equals([1, 2])
        assert not expected.equals(np.array([1, 2]))
        assert not expected.equals(date_range("20130101", periods=2))

        expected_name1 = IntervalIndex.from_breaks(
            np.arange(5), closed=closed, name="foo"
        )
        expected_name2 = IntervalIndex.from_breaks(
            np.arange(5), closed=closed, name="bar"
        )
        assert expected.equals(expected_name1)
        assert expected_name1.equals(expected_name2)

        for other_closed in {"left", "right", "both", "neither"} - {closed}:
            expected_other_closed = IntervalIndex.from_breaks(
                np.arange(5), closed=other_closed
            )
            assert not expected.equals(expected_other_closed)

    def test_repr_max_seq_item_setting(self):
        # override base test: not a valid repr as we use interval notation
        pass

    def test_repr_roundtrip(self):
        # override base test: not a valid repr as we use interval notation
        pass

    def test_take(self, closed):
        index = self.create_index(closed=closed)

        result = index.take(range(10))
        tm.assert_index_equal(result, index)

        result = index.take([0, 0, 1])
        expected = IntervalIndex.from_arrays([0, 0, 1], [1, 1, 2], closed=closed)
        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("klass", [list, tuple, np.array, Series])
    def test_where(self, closed, klass):
        idx = self.create_index(closed=closed)
        cond = [True] * len(idx)
        expected = idx
        result = expected.where(klass(cond))
        tm.assert_index_equal(result, expected)

        cond = [False] + [True] * len(idx[1:])
        expected = IntervalIndex([np.nan] + idx[1:].tolist())
        result = idx.where(klass(cond))
        tm.assert_index_equal(result, expected)
