# Tests aimed at pandas.core.indexers
import numpy as np

from pandas.core.indexers import is_scalar_indexer, length_of_indexer


def test_length_of_indexer():
    arr = np.zeros(4, dtype=bool)
    arr[0] = 1
    result = length_of_indexer(arr)
    assert result == 1


def test_is_scalar_indexer():
    indexer = (0, 1)
    assert is_scalar_indexer(indexer, 2)
    assert not is_scalar_indexer(indexer[0], 2)

    indexer = (np.array([2]), 1)
    assert is_scalar_indexer(indexer, 2)

    indexer = (np.array([2]), np.array([3]))
    assert is_scalar_indexer(indexer, 2)

    indexer = (np.array([2]), np.array([3, 4]))
    assert not is_scalar_indexer(indexer, 2)

    assert not is_scalar_indexer(slice(None), 1)
