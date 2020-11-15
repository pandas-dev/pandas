from datetime import timedelta

import pytest

from pandas import NA, DataFrame, Index, date_range
import pandas._testing as tm


class TestReindex:
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"method": "pad", "tolerance": timedelta(seconds=9)},
            {"method": "backfill", "tolerance": timedelta(seconds=9)},
            {"method": "nearest"},
            {"method": None},
        ],
    )
    def test_reindex_empty_frame(self, kwargs):
        # GH#27315
        idx = date_range(start="2020", freq="30s", periods=3)
        df = DataFrame([], index=Index([], name="time"), columns=["a"])
        result = df.reindex(idx, **kwargs)
        expected = DataFrame({"a": [NA] * 3}, index=idx)
        tm.assert_frame_equal(result, expected)
