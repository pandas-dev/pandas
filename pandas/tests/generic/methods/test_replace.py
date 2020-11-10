import re

import pytest

import pandas as pd
import pandas._testing as tm


class TestReplace:
    @pytest.mark.parametrize("value", [pd.Period("2020-01"), pd.Interval(0, 5)])
    def test_replace_ea_ignore_float(self, frame_or_series, value):
        # GH#34871
        df = frame_or_series([value] * 3)
        result = df.replace(1.0, 0.0)
        expected = frame_or_series([value] * 3)
        tm.assert_equal(expected, result)

    def test_replace_with_compiled_regex(self):
        # https://github.com/pandas-dev/pandas/issues/35680
        s = pd.Series(["a", "b", "c"])
        regex = re.compile("^a$")
        result = s.replace({regex: "z"}, regex=True)
        expected = pd.Series(["z", "b", "c"])
        tm.assert_series_equal(result, expected)