import pytest

from pandas import Series
from pandas.api.indexers import BaseIndexer


def test_bad_get_window_bounds_signature():
    class BadIndexer(BaseIndexer):
        def get_window_bounds(self):
            return None

    indexer = BadIndexer()
    with pytest.raises(
        ValueError,
        match="BadIndexer does not implement the correct signature "
              "forget_window_bounds.",
    ):
        Series(range(5)).rolling(indexer)
