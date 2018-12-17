from datetime import datetime

import numpy as np
import pytest

from pandas import DataFrame, Series
from pandas.core.indexes.datetimes import date_range
from pandas.core.indexes.period import period_range
from pandas.core.indexes.timedeltas import timedelta_range

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


class ResampleFixture(object):
    def __init__(self, _index_factory, _series_name, _index_start=None,
                 _index_end=None, _index_freq=None):
        self._index_factory = _index_factory
        self._series_name = _series_name
        self._index_start = _index_start or datetime(2005, 1, 1)
        self._index_end = _index_end or datetime(2005, 1, 10)
        self._index_freq = _index_freq or 'D'

    @property
    def index(self):
        return self._index_factory(
            self._index_start, self._index_end, freq=self._index_freq)

    @property
    def _static_values(self):
        return np.arange(len(self.index))

    @property
    def series(self):
        return Series(
            self._static_values, index=self.index, name=self._series_name)

    @property
    def frame(self):
        return DataFrame({'value': self._static_values}, index=self.index)

    def create_index(self, *args, **kwargs):
        """ return the _index_factory created using the args, kwargs """
        return self._index_factory(*args, **kwargs)


@pytest.fixture
def date_range_fixture():
    return ResampleFixture(date_range, 'dti')


@pytest.fixture
def period_range_fixture():
    return ResampleFixture(period_range, 'pi')


@pytest.fixture
def timedelta_range_fixture():
    return ResampleFixture(timedelta_range, 'tdi', '1 day', '10 day')


@pytest.fixture
def index(resample_fixture):
    return resample_fixture.index


@pytest.fixture
def series(resample_fixture):
    return resample_fixture.series


@pytest.fixture
def frame(resample_fixture):
    return resample_fixture.frame


@pytest.fixture(params=['series', 'frame'])
def series_and_frame(request, resample_fixture):
    return getattr(resample_fixture, request.param)
