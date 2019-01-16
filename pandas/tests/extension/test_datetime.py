"""
Extension Tests for *tz-naive* arrays.DatetimeArray.

Currently, we only run the Dtype tests, as we do not allow a
tz-naive DatetimeArray inside internals.
"""
import numpy as np
import pytest

import pandas as pd
from pandas.tests.extension import base


@pytest.fixture
def dtype():
    return pd.DatetimeDtype()


@pytest.fixture
def data():
    return pd.arrays.DatetimeArray(np.arange(0, 100, dtype='M8[ns]'))


class BaseTimedeltaTests(object):
    pass


class TestDtype(BaseTimedeltaTests, base.BaseDtypeTests):
    pass
