import pytest

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
