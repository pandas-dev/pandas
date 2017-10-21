import pytest
import pandas.tseries.offsets as offsets


@pytest.fixture(params=[getattr(offsets, o) for o in offsets.__all__])
def offset_types(request):
    return request.param


@pytest.fixture(params=[None, 'UTC', 'Asia/Tokyo', 'US/Eastern',
                        'dateutil/Asia/Tokyo', 'dateutil/US/Pacific'])
def tz(request):
    return request.param
