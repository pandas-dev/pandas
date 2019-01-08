"""
Extension Tests for arrays.TimedeltaArray.

Currently, we only run the Dtype tests, as we do not allow a
TimedeltaArray inside internals.
"""
import numpy as np
import pytest

import pandas as pd
from pandas.tests.extension import base


@pytest.fixture
def dtype():
    return pd.TimedeltaDtype()


@pytest.fixture
def data():
    return pd.arrays.TimedeltaArray(np.arange(100))


class BaseTimedeltaTests(object):
    pass


class TestDtype(BaseTimedeltaTests, base.BaseDtypeTests):
    pass
