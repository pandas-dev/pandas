import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas import DataFrame, Series, date_range


class TestBetweenTime:
    @td.skip_if_has_locale
    def test_between_time_formats(self):
        # GH11818
        rng = date_range("1/1/2000", "1/5/2000", freq="5min")
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)

        strings = [
            ("2:00", "2:30"),
            ("0200", "0230"),
            ("2:00am", "2:30am"),
            ("0200am", "0230am"),
            ("2:00:00", "2:30:00"),
            ("020000", "023000"),
            ("2:00:00am", "2:30:00am"),
            ("020000am", "023000am"),
        ]
        expected_length = 28

        for time_string in strings:
            assert len(ts.between_time(*time_string)) == expected_length

    def test_between_time_axis(self):
        # issue 8839
        rng = date_range("1/1/2000", periods=100, freq="10min")
        ts = Series(np.random.randn(len(rng)), index=rng)
        stime, etime = ("08:00:00", "09:00:00")
        expected_length = 7

        assert len(ts.between_time(stime, etime)) == expected_length
        assert len(ts.between_time(stime, etime, axis=0)) == expected_length
        msg = "No axis named 1 for object type Series"
        with pytest.raises(ValueError, match=msg):
            ts.between_time(stime, etime, axis=1)
