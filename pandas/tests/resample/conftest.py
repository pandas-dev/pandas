from datetime import datetime

import numpy as np
import pytest

from pandas import Series

# The various methods we support
downsample_methods = [
    "min",
    "max",
    "first",
    "last",
    "sum",
    "mean",
    "sem",
    "median",
    "prod",
    "var",
    "std",
    "ohlc",
    "quantile",
]
upsample_methods = ["count", "size"]
series_methods = ["nunique"]
resample_methods = downsample_methods + upsample_methods + series_methods


@pytest.fixture(params=downsample_methods)
def downsample_method(request):
    """Fixture for parametrization of Grouper downsample methods."""
    return request.param


@pytest.fixture(params=resample_methods)
def resample_method(request):
    """Fixture for parametrization of Grouper resample methods."""
    return request.param


@pytest.fixture
def _index_start():
    """Fixture for parametrization of index, series and frame."""
    return datetime(2005, 1, 1)


@pytest.fixture
def _index_end():
    """Fixture for parametrization of index, series and frame."""
    return datetime(2005, 1, 10)


@pytest.fixture
def _index_freq():
    """Fixture for parametrization of index, series and frame."""
    return "D"


@pytest.fixture
def index(_index_factory, _index_start, _index_end):
    """
    Fixture for parametrization of date_range, period_range and
    timedelta_range indexes
    """
    return _index_factory(_index_start, _index_end, freq="D", name=None)


@pytest.fixture
def _static_values(index):
    """
    Fixture for parametrization of values used in parametrization of
    Series and DataFrames with date_range, period_range and
    timedelta_range indexes
    """
    return np.arange(len(index))


@pytest.fixture
def _series_name():
    """
    Fixture for parametrization of Series name for Series used with
    date_range, period_range and timedelta_range indexes
    """
    return None


@pytest.fixture
def series(index, _series_name, _static_values):
    """
    Fixture for parametrization of Series with date_range, period_range and
    timedelta_range indexes
    """
    return Series(_static_values, index=index, name=_series_name)
