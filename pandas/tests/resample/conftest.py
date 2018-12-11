import numpy as np
import pytest

from pandas import Series
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


@pytest.fixture()
def simple_date_range_series():
    """
    Series with date range index and random data for test purposes.
    """
    def _simple_date_range_series(start, end, freq='D'):
        rng = date_range(start, end, freq=freq)
        return Series(np.random.randn(len(rng)), index=rng)
    return _simple_date_range_series


@pytest.fixture()
def simple_period_range_series():
    """
    Series with period range index and random data for test purposes.
    """
    def _simple_period_range_series(start, end, freq='D'):
        rng = period_range(start, end, freq=freq)
        return Series(np.random.randn(len(rng)), index=rng)
    return _simple_period_range_series
