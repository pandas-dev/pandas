import numpy as np

import pandas as pd


class TestIsMonotonic:
    def test_is_monotonic_numeric(self):
        ser = pd.Series(np.random.default_rng(2).integers(0, 10, size=1000))
        assert not ser.is_monotonic_increasing
        ser = pd.Series(np.arange(1000))
        assert ser.is_monotonic_increasing is True
        assert ser.is_monotonic_increasing is True
        ser = pd.Series(np.arange(1000, 0, -1))
        assert ser.is_monotonic_decreasing is True

    def test_is_monotonic_dt64(self):
        ser = pd.Series(pd.date_range("20130101", periods=10))
        assert ser.is_monotonic_increasing is True
        assert ser.is_monotonic_increasing is True

        ser = pd.Series(list(reversed(ser)))
        assert ser.is_monotonic_increasing is False
        assert ser.is_monotonic_decreasing is True
