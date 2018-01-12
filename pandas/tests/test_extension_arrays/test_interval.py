import pytest

import pandas as pd
from pandas.core.interval import IntervalArray

from .base import BaseArrayTests


@pytest.fixture
def test_data():
    """Length-100 PeriodArray for semantics test."""
    return IntervalArray(pd.interval_range(0, periods=100))


class TestPeriod(BaseArrayTests):
    pass
