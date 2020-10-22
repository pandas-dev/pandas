from datetime import timedelta

import numpy as np
import pytest

import pandas as pd
from pandas import Index, Timedelta, TimedeltaIndex, timedelta_range
import pandas._testing as tm


class TestTimedeltaIndexInsert:
    def test_insert(self):

        idx = TimedeltaIndex(["4day", "1day", "2day"], name="idx")

        result = idx.insert(2, timedelta(days=5))
        exp = TimedeltaIndex(["4day", "1day", "5day", "2day"], name="idx")
        tm.assert_index_equal(result, exp)

        # insertion of non-datetime should coerce to object index
        result = idx.insert(1, "inserted")
        expected = Index(
            [Timedelta("4day"), "inserted", Timedelta("1day"), Timedelta("2day")],
            name="idx",
        )
        assert not isinstance(result, TimedeltaIndex)
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name

        idx = timedelta_range("1day 00:00:01", periods=3, freq="s", name="idx")

        # preserve freq
        expected_0 = TimedeltaIndex(
            ["1day", "1day 00:00:01", "1day 00:00:02", "1day 00:00:03"],
            name="idx",
            freq="s",
        )
        expected_3 = TimedeltaIndex(
            ["1day 00:00:01", "1day 00:00:02", "1day 00:00:03", "1day 00:00:04"],
            name="idx",
            freq="s",
        )

        # reset freq to None
        expected_1_nofreq = TimedeltaIndex(
            ["1day 00:00:01", "1day 00:00:01", "1day 00:00:02", "1day 00:00:03"],
            name="idx",
            freq=None,
        )
        expected_3_nofreq = TimedeltaIndex(
            ["1day 00:00:01", "1day 00:00:02", "1day 00:00:03", "1day 00:00:05"],
            name="idx",
            freq=None,
        )

        cases = [
            (0, Timedelta("1day"), expected_0),
            (-3, Timedelta("1day"), expected_0),
            (3, Timedelta("1day 00:00:04"), expected_3),
            (1, Timedelta("1day 00:00:01"), expected_1_nofreq),
            (3, Timedelta("1day 00:00:05"), expected_3_nofreq),
        ]

        for n, d, expected in cases:
            result = idx.insert(n, d)
            tm.assert_index_equal(result, expected)
            assert result.name == expected.name
            assert result.freq == expected.freq

    @pytest.mark.parametrize(
        "null", [None, np.nan, np.timedelta64("NaT"), pd.NaT, pd.NA]
    )
    def test_insert_nat(self, null):
        # GH 18295 (test missing)
        idx = timedelta_range("1day", "3day")
        result = idx.insert(1, null)
        expected = TimedeltaIndex(["1day", pd.NaT, "2day", "3day"])
        tm.assert_index_equal(result, expected)

    def test_insert_invalid_na(self):
        idx = TimedeltaIndex(["4day", "1day", "2day"], name="idx")
        msg = r"value should be a 'Timedelta' or 'NaT'\. Got 'datetime64' instead\."
        with pytest.raises(TypeError, match=msg):
            idx.insert(0, np.datetime64("NaT"))

    @pytest.mark.parametrize(
        "item", [0, np.int64(0), np.float64(0), np.array(0), np.datetime64(456, "us")]
    )
    def test_insert_mismatched_types_raises(self, item):
        # GH#33703 dont cast these to td64
        tdi = TimedeltaIndex(["4day", "1day", "2day"], name="idx")

        msg = r"value should be a 'Timedelta' or 'NaT'\. Got '.*' instead\."
        with pytest.raises(TypeError, match=msg):
            tdi.insert(1, item)

    def test_insert_dont_cast_strings(self):
        # To match DatetimeIndex and PeriodIndex behavior, dont try to
        #  parse strings to Timedelta
        idx = timedelta_range("1day", "3day")

        result = idx.insert(0, "1 Day")
        assert result.dtype == object
        assert result[0] == "1 Day"

    def test_insert_empty(self):
        # Corner case inserting with length zero doesnt raise IndexError
        # GH#33573 for freq preservation
        idx = timedelta_range("1 Day", periods=3)
        td = idx[0]

        result = idx[:0].insert(0, td)
        assert result.freq == "D"

        result = idx[:0].insert(1, td)
        assert result.freq == "D"

        result = idx[:0].insert(-1, td)
        assert result.freq == "D"
