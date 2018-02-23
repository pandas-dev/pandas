import pytest
import pandas.tseries.offsets as offsets


@pytest.fixture(params=[getattr(offsets, o) for o in offsets.__all__])
def offset_types(request):
    """
    Fixture for all the datetime offsets available for a time series.
    """
    return request.param


@pytest.fixture(params=[getattr(offsets, o) for o in offsets.__all__ if
                        issubclass(getattr(offsets, o), offsets.MonthOffset)
                        and o != 'MonthOffset'])
def month_classes(request):
    """
    Fixture for all the month based datetime offsets available for a time series.
    """
    return request.param


@pytest.fixture(params=[getattr(offsets, o) for o in offsets.__all__ if
                        issubclass(getattr(offsets, o), offsets.Tick)])
def tick_classes(request):
    """
    Fixture for all the Tick (day, hour, ..., nano) based datetime offsets available for a time series.
    """
    return request.param


@pytest.fixture(params=[None, 'UTC', 'Asia/Tokyo', 'US/Eastern',
                        'dateutil/Asia/Tokyo', 'dateutil/US/Pacific'])
def tz(request):
    """
    Fixture for different timezones strings.
    """
    return request.param
