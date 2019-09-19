import pytest

import pandas.util.testing as tm


@pytest.fixture
def datetime_series():
    """
    Fixture for Series of floats with DatetimeIndex
    """
    s = tm.makeTimeSeries()
    s.name = "ts"
    return s


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
