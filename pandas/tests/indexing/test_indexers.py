# Tests aimed at pandas.core.indexers
import numpy as np
import pytest

from pandas.core.indexers import (
    is_scalar_indexer,
    length_of_indexer,
    validate_indices,
)


def test_length_of_indexer():
    arr = np.zeros(4, dtype=bool)
    arr[0] = 1
    result = length_of_indexer(arr)
    assert result == 1


@pytest.mark.parametrize(
    "start, stop, step, target_len, expected",
    [
        # GH#66100 negative step with None bounds previously returned
        # negative lengths
        (None, None, -1, 7, 7),
        (None, None, -2, 7, 4),
        (None, None, -3, 7, 3),
        (4, None, -2, 7, 3),
        (None, 0, -2, 7, 3),
        (None, 1, -1, 7, 5),
        (3, None, -1, 7, 4),
        # explicit (non-None) bounds with negative step already worked
        (4, 0, -2, 7, 2),
        (6, 0, -1, 7, 6),
        (5, 1, -2, 7, 2),
        (-1, -5, -1, 7, 4),
        # positive step, including GH#9995 non-divisible-step case
        (1, None, 2, 6, 3),
        (None, None, 1, 7, 7),
        (1, 5, 2, 7, 2),
        (None, None, None, 7, 7),
        (2, None, None, 7, 5),
    ],
)
def test_length_of_indexer_slice(start, stop, step, target_len, expected):
    target = np.arange(target_len)
    indexer = slice(start, stop, step)
    result = length_of_indexer(indexer, target)
    assert result == expected
    # cross-check against the actual sliced length, the real source of truth
    assert result == len(target[indexer])


@pytest.mark.parametrize(
    "start, stop, step, expected",
    [
        (10, 0, -3, 4),
        (0, 1, 5, 1),
    ],
)
def test_length_of_indexer_range_overflow(start, stop, step, expected):
    # https://github.com/pandas-dev/pandas/pull/63872
    indexer = range(start, stop, step)
    result = length_of_indexer(indexer)
    assert result == expected


def test_is_scalar_indexer():
    indexer = (0, 1)
    assert is_scalar_indexer(indexer, 2)
    assert not is_scalar_indexer(indexer[0], 2)

    indexer = (np.array([2]), 1)
    assert not is_scalar_indexer(indexer, 2)

    indexer = (np.array([2]), np.array([3]))
    assert not is_scalar_indexer(indexer, 2)

    indexer = (np.array([2]), np.array([3, 4]))
    assert not is_scalar_indexer(indexer, 2)

    assert not is_scalar_indexer(slice(None), 1)

    indexer = 0
    assert is_scalar_indexer(indexer, 1)

    indexer = (0,)
    assert is_scalar_indexer(indexer, 1)


class TestValidateIndices:
    def test_validate_indices_ok(self):
        indices = np.asarray([0, 1])
        validate_indices(indices, 2)
        validate_indices(indices[:0], 0)
        validate_indices(np.array([-1, -1]), 0)

    def test_validate_indices_low(self):
        indices = np.asarray([0, -2])
        with pytest.raises(ValueError, match="'indices' contains"):
            validate_indices(indices, 2)

    def test_validate_indices_high(self):
        indices = np.asarray([0, 1, 2])
        with pytest.raises(IndexError, match="indices are out"):
            validate_indices(indices, 2)

    def test_validate_indices_empty(self):
        with pytest.raises(IndexError, match="indices are out"):
            validate_indices(np.array([0, 1]), 0)
