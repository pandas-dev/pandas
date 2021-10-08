import numpy as np
import pytest

from pandas import (
    Series,
    array,
)


@pytest.fixture(params=[None, False])
def sort(request):
    """
    Valid values for the 'sort' parameter used in the Index
    setops methods (intersection, union, etc.)

    Caution:
        Don't confuse this one with the "sort" fixture used
        for DataFrame.append or concat. That one has
        parameters [True, False].

        We can't combine them as sort=True is not permitted
        in the Index setops methods.
    """
    return request.param


@pytest.fixture(params=["D", "3D", "-3D", "H", "2H", "-2H", "T", "2T", "S", "-3S"])
def freq_sample(request):
    """
    Valid values for 'freq' parameter used to create date_range and
    timedelta_range..
    """
    return request.param


@pytest.fixture(params=[list, np.array, array, Series])
def listlike_box(request):
    """
    Types that may be passed as the indexer to searchsorted.
    """
    return request.param


# TODO: not clear if this _needs_ to be different from listlike_box or
#  if that is just a historical artifact
@pytest.fixture(params=[list, tuple, np.array, Series])
def listlike_box_with_tuple(request):
    """
    Types that may be passed as the indexer to searchsorted.
    """
    return request.param
