import pytest

import pandas._testing as tm


@pytest.fixture
def string_series():
    """
    Fixture for Series of floats with Index of unique strings
    """
    s = tm.makeStringSeries()
    s.name = "series"
    return s


@pytest.fixture
def object_series():
    """
    Fixture for Series of dtype object with Index of unique strings
    """
    s = tm.makeObjectSeries()
    s.name = "objects"
    return s
