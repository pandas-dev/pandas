import numpy as np
import pytest

from pandas import Series
import pandas.util.testing as tm

_all_indices = [tm.makeBoolIndex(10, name='a'),
                tm.makeIntIndex(10, name='a'),
                tm.makeFloatIndex(10, name='a'),
                tm.makeDateIndex(10, name='a'),
                tm.makeDateIndex(10, name='a').tz_localize(tz='US/Eastern'),
                tm.makePeriodIndex(10, name='a'),
                tm.makeStringIndex(10, name='a'),
                tm.makeUnicodeIndex(10, name='a')]

_all_indices_ids = [
    'testing boolean index',
    'testing integer index',
    'testing float index',
    'testing date index',
    'testing date index with timezone',
    'testing period index',
    'testing string index',
    'testing unicode index',
]

_arr = np.arange(10)

_all_series = [Series(_arr, index=tm.makeBoolIndex(10, name='a')),
               Series(_arr, index=tm.makeIntIndex(10, name='a')),
               Series(_arr, index=tm.makeFloatIndex(10, name='a')),
               Series(_arr, index=tm.makeDateIndex(10, name='a')),
               tm.makeDateIndex(10, name='a').tz_localize(
                   tz='US/Eastern').to_series(keep_tz=True),
               Series(_arr, index=tm.makePeriodIndex(10, name='a')),
               Series(_arr, index=tm.makeStringIndex(10, name='a')),
               Series(_arr, index=tm.makeUnicodeIndex(10, name='a'))]

_all_series_ids = [
    'testing pandas series with boolean index',
    'testing pandas series with integer index',
    'testing pandas series with float index',
    'testing pandas series with date index',
    'testing pandas series with date index with timezone',
    'testing pandas series with period index',
    'testing pandas series with string index',
    'testing pandas series with unicode index',
]


@pytest.fixture(params=_all_indices, ids=_all_indices_ids)
def all_indices_fixture(request):
    """Fixture for testing index operations"""
    return request.param


@pytest.fixture(params=_all_series, ids=_all_series_ids)
def all_series_fixture(request):
    """Fixture for testing pandas Series operations"""
    return request.param


@pytest.fixture(params=_all_indices + _all_series,
                ids=_all_indices_ids + _all_series_ids)
def all_indices_and_series_fixture(request):
    """Fixture for testing index and pandas Series operations"""
    return request.param
