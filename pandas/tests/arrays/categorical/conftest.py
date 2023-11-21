import pytest

from pandas import Categorical


@pytest.fixture
def factor():
    """Fixture returning  a Categorical object"""
    return Categorical(["a", "b", "b", "a", "a", "c", "c", "c"], ordered=True)
