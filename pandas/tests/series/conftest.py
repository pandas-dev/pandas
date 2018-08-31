import pytest

import pandas.util.testing as tm

from pandas import Series


@pytest.fixture
def datetime_series():
    """
    Fixture for Series of floats with DatetimeIndex

    See pandas.util.testing.makeTimeSeries
    """
    s = tm.makeTimeSeries()
    s.name = 'ts'
    return s


@pytest.fixture
def string_series():
    """
    Fixture for Series of floats with Index of unique strings

    See pandas.util.testing.makeStringSeries
    """
    s = tm.makeStringSeries()
    s.name = 'series'
    return s


@pytest.fixture
def object_series():
    """
    Fixture for Series of dtype datetime64[ns] with Index of unique strings

    See pandas.util.testing.makeObjectSeries
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
