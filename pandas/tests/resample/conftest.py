import pytest

from pandas.tests.resample.test_base import (
    downsample_methods, resample_methods, upsample_methods)


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
