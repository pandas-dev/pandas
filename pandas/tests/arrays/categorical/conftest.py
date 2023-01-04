import pytest

from pandas import Categorical


@pytest.fixture(name="allow_fill", params=[True, False])
def fixture_allow_fill(request):
    """Boolean 'allow_fill' parameter for Categorical.take"""
    return request.param


@pytest.fixture(name="factor")
def fixture_factor():
    """Fixture returning  a Categorical object"""
    return Categorical(["a", "b", "b", "a", "a", "c", "c", "c"], ordered=True)
