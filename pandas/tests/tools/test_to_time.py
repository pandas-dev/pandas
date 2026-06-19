from datetime import time

import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import Series
import pandas._testing as tm
from pandas.core.tools.times import to_time


class TestToTime:
    @pytest.mark.parametrize(
        "time_string",
        [
            "14:15",
            "1415",
            pytest.param("2:15pm", marks=td.skip_if_not_us_locale),
            pytest.param("0215pm", marks=td.skip_if_not_us_locale),
            "14:15:00",
            "141500",
            pytest.param("2:15:00pm", marks=td.skip_if_not_us_locale),
            pytest.param("021500pm", marks=td.skip_if_not_us_locale),
            time(14, 15),
        ],
    )
    def test_parsers_time(self, time_string):
        # GH#11818
        assert to_time(time_string) == time(14, 15)

    def test_odd_format(self):
        new_string = "14.15"
        assert to_time(new_string, format="%H.%M") == time(14, 15)

    def test_arraylike(self):
        arg = ["14:15", "20:20"]
        expected_arr = [time(14, 15), time(20, 20)]
        assert to_time(arg) == expected_arr
        assert to_time(arg, format="%H:%M") == expected_arr
        assert to_time(arg, infer_time_format=True) == expected_arr
        assert to_time(arg, format="%I:%M%p", errors="coerce") == [None, None]

        with pytest.raises(ValueError, match="errors must be"):
            to_time(arg, format="%I:%M%p", errors="ignore")

        msg = "Cannot convert.+to a time with given format"
        with pytest.raises(ValueError, match=msg):
            to_time(arg, format="%I:%M%p", errors="raise")

        tm.assert_series_equal(
            to_time(Series(arg, name="test")), Series(expected_arr, name="test")
        )

        res = to_time(np.array(arg))
        assert isinstance(res, list)
        assert res == expected_arr
