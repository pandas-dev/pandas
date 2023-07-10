import pytest

import pandas as pd
import pandas._testing as tm


class TestResolution:
    @pytest.mark.parametrize(
        "freq,expected",
        [
            ("A", "year"),
            ("Q", "quarter"),
            ("M", "month"),
            ("D", "day"),
            ("H", "hour"),
            ("T", "minute"),
            ("S", "second"),
            ("L", "millisecond"),
            ("U", "microsecond"),
        ],
    )
    def test_resolution(self, freq, expected):
        msg = f"Code freq={freq} is deprecated and will be removed in a future version."

        if freq in {"T", "L"}:
            with tm.assert_produces_warning(FutureWarning, match=msg):
                idx = pd.period_range(start="2013-04-01", periods=30, freq=freq)
                assert idx.resolution == expected
        else:
            idx = pd.period_range(start="2013-04-01", periods=30, freq=freq)
            assert idx.resolution == expected
