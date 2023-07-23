import pytest

import pandas as pd


class TestResolution:
    @pytest.mark.parametrize(
        "freq,expected",
        [
            ("A", "year"),
            ("Q", "quarter"),
            ("M", "month"),
            ("D", "day"),
            ("H", "hour"),
            ("min", "minute"),
            ("S", "second"),
            ("ms", "millisecond"),
            ("U", "microsecond"),
        ],
    )
    def test_resolution(self, freq, expected):
        idx = pd.period_range(start="2013-04-01", periods=30, freq=freq)
        assert idx.resolution == expected
