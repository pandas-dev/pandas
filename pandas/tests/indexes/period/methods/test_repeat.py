import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestRepeat:
    @pytest.mark.parametrize("use_numpy", [True, False])
    @pytest.mark.parametrize(
        "index",
        [
            pd.period_range("2000-01-01", periods=3, freq="D"),
            pd.period_range("2001-01-01", periods=3, freq="2D"),
            pd.PeriodIndex(["2001-01", "NaT", "2003-01"], freq="M"),
        ],
    )
    def test_repeat_freqstr(self, index, use_numpy):
        # GH#10183
        expected = pd.PeriodIndex([per for per in index for _ in range(3)])
        result = np.repeat(index, 3) if use_numpy else index.repeat(3)
        tm.assert_index_equal(result, expected)
        assert result.freqstr == index.freqstr
