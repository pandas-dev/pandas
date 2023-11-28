import datetime

import pytest

from pandas._libs.tslibs import Timestamp


@pytest.fixture
def dt():
    """
    Fixture for common Timestamp.
    """
    return Timestamp(datetime.datetime(2008, 1, 2))
