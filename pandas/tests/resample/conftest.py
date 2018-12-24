from datetime import datetime

import numpy as np
import pytest

from pandas import DataFrame, Series
from pandas.core.indexes.datetimes import date_range
from pandas.core.indexes.period import period_range

# The various methods we support
downsample_methods = ['min', 'max', 'first', 'last', 'sum', 'mean', 'sem',
                      'median', 'prod', 'var', 'std', 'ohlc', 'quantile']
upsample_methods = ['count', 'size']
series_methods = ['nunique']
resample_methods = downsample_methods + upsample_methods + series_methods


@pytest.fixture(params=downsample_methods)
def downsample_method(request):
    """Fixture for parametrization of Grouper downsample methods."""
    return request.param


@pytest.fixture(params=upsample_methods)
def upsample_method(request):
    """Fixture for parametrization of Grouper upsample methods."""
    return request.param


@pytest.fixture(params=resample_methods)
def resample_method(request):
    """Fixture for parametrization of Grouper resample methods."""
    return request.param


@pytest.fixture
def simple_date_range_series():
    """
    Series with date range index and random data for test purposes.
    """
    def _simple_date_range_series(start, end, freq='D'):
        rng = date_range(start, end, freq=freq)
        return Series(np.random.randn(len(rng)), index=rng)
    return _simple_date_range_series


@pytest.fixture
def simple_period_range_series():
    """
    Series with period range index and random data for test purposes.
    """
    def _simple_period_range_series(start, end, freq='D'):
        rng = period_range(start, end, freq=freq)
        return Series(np.random.randn(len(rng)), index=rng)
    return _simple_period_range_series


@pytest.fixture
def _index_start():
    return datetime(2005, 1, 1)


@pytest.fixture
def _index_end():
    return datetime(2005, 1, 10)


@pytest.fixture
def _index_freq():
    return 'D'


@pytest.fixture
def index(_index_factory, _index_start, _index_end, _index_freq):
    return _index_factory(_index_start, _index_end, freq=_index_freq)


@pytest.fixture
def _static_values(index):
    return np.arange(len(index))


@pytest.fixture
def series(index, _series_name, _static_values):
    return Series(_static_values, index=index, name=_series_name)


@pytest.fixture
def empty_series(series):
    return series[:0]


@pytest.fixture
def frame(index, _series_name, _static_values):
    # _series_name is intentionally unused
    return DataFrame({'value': _static_values}, index=index)


@pytest.fixture
def empty_frame(series):
    index = series.index[:0]
    return DataFrame(index=index)


@pytest.fixture(params=[Series, DataFrame])
def series_and_frame(request, series, frame):
    if request.param == Series:
        return series
    if request.param == DataFrame:
        return frame
