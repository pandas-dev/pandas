import numpy as np
import pytest

from pandas import DataFrame, Index, MultiIndex


@pytest.fixture
def multiindex_dataframe_random_data():
    """DataFrame with 2 level MultiIndex with random data"""
    index = MultiIndex(
        levels=[["foo", "bar", "baz", "qux"], ["one", "two", "three"]],
        codes=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3], [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
        names=["first", "second"],
    )
    return DataFrame(
        np.random.randn(10, 3), index=index, columns=Index(["A", "B", "C"], name="exp")
    )
