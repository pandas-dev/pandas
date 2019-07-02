"""
Tests for reshaping utilities.
"""
import pytest

from pandas.core.arrays._reshaping import tuplify_shape


class TestTuplify:
    def test_tuplify_single_arg(self):
        # Single-tuple cases, i.e.
        #  arr.reshape((x, y))
        shape = tuplify_shape(3, ((3,),))
        assert shape == (3,)

        shape = tuplify_shape(3, ((1, 3),))
        assert shape == (1, 3)

        shape = tuplify_shape(3, ((3, 1),))
        assert shape == (3, 1)

    def test_tuplify_multi_arg(self):
        # Multi-arg cases, i.e.
        #  arr.reshape(x, y)
        shape = tuplify_shape(3, (3,))
        assert shape == (3,)

        shape = tuplify_shape(3, (3, 1))
        assert shape == (3, 1)

        shape = tuplify_shape(3, (1, 3))
        assert shape == (1, 3)

    def test_tuplify_minus_one(self):
        shape = tuplify_shape(4, (1, -1))
        assert shape == (1, 4)

        shape = tuplify_shape(4, (-1, 1))
        assert shape == (4, 1)

    def test_tuplify_minus_one_factors(self):
        shape = tuplify_shape(4, (1, -1, 2), restrict=False)
        assert shape == (1, 2, 2)

    def test_tuplify_multiple_minus_ones(self):
        # No more than 1 "-1"
        with pytest.raises(ValueError, match="Invalid shape"):
            tuplify_shape(99, (-1, -1))

    def test_tuplify_negative(self):
        # Nothing less than -1 in a shape
        with pytest.raises(ValueError, match="Invalid shape"):
            tuplify_shape(99, (-2, 3))

    def test_tuplify_size_match(self):
        # must match original size
        with pytest.raises(ValueError, match="Product of shape"):
            tuplify_shape(3, (2, 2))
