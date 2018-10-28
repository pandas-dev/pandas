import pytest

import pandas.util.testing as tm
from pandas import Series


@pytest.fixture
def datetime_series():
    """
    Fixture for Series of floats with DatetimeIndex
    """
    s = tm.makeTimeSeries()
    s.name = 'ts'
    return s


@pytest.fixture
def string_series():
    """
    Fixture for Series of floats with Index of unique strings
    """
    s = tm.makeStringSeries()
    s.name = 'series'
    return s


@pytest.fixture
def object_series():
    """
    Fixture for Series of dtype datetime64[ns] with Index of unique strings
    """
    s = tm.makeObjectSeries()
    s.name = 'objects'
    return s


@pytest.fixture
def empty_series():
    """
    Fixture for empty Series
    """
    return Series([], index=[])
