import pytest

import pandas.util.testing as tm

from pandas import Series


@pytest.fixture
def ts():
    ts = tm.makeTimeSeries()
    ts.name = 'ts'
    return ts


@pytest.fixture
def series():
    series = tm.makeStringSeries()
    series.name = 'series'
    return series


@pytest.fixture
def objSeries():
    objSeries = tm.makeObjectSeries()
    objSeries.name = 'objects'
    return objSeries


@pytest.fixture
def empty():
    return Series([], index=[])
