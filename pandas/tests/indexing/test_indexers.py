# Tests aimed at pandas.core.indexers
import numpy as np

from pandas.core.indexers import length_of_indexer


def test_length_of_indexer():
    arr = np.zeros(4, dtype=bool)
    arr[0] = 1
    result = length_of_indexer(arr)
    assert result == 1
