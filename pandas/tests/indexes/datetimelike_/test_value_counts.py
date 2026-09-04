import numpy as np

import pandas as pd
import pandas._testing as tm


class TestValueCounts:
    # GH#7735

    def test_value_counts_unique_datetimeindex(self, tz_naive_fixture):
        tz = tz_naive_fixture
        orig = pd.date_range("2011-01-01 09:00", freq="h", periods=10, tz=tz)
        self._check_value_counts_with_repeats(orig)

    def test_value_counts_unique_timedeltaindex(self):
        orig = pd.timedelta_range("1 days 09:00:00", freq="h", periods=10)
        self._check_value_counts_with_repeats(orig)

    def test_value_counts_unique_periodindex(self):
        orig = pd.period_range("2011-01-01 09:00", freq="h", periods=10)
        self._check_value_counts_with_repeats(orig)

    def _check_value_counts_with_repeats(self, orig):
        # create repeated values, 'n'th element is repeated by n+1 times
        idx = type(orig)(
            np.repeat(orig._values, range(1, len(orig) + 1)), dtype=orig.dtype
        )

        exp_idx = orig[::-1]
        if not isinstance(exp_idx, pd.PeriodIndex):
            exp_idx = exp_idx._with_freq(None)
        expected = pd.Series(
            range(10, 0, -1), index=exp_idx, dtype="int64", name="count"
        )

        for obj in [idx, pd.Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        tm.assert_index_equal(idx.unique(), orig, check_freq=False)

    def test_value_counts_unique_datetimeindex2(self, tz_naive_fixture):
        tz = tz_naive_fixture
        idx = pd.DatetimeIndex(
            [
                "2013-01-01 09:00",
                "2013-01-01 09:00",
                "2013-01-01 09:00",
                "2013-01-01 08:00",
                "2013-01-01 08:00",
                pd.NaT,
            ],
            tz=tz,
        )
        self._check_value_counts_dropna(idx)

    def test_value_counts_unique_timedeltaindex2(self):
        idx = pd.TimedeltaIndex(
            [
                "1 days 09:00:00",
                "1 days 09:00:00",
                "1 days 09:00:00",
                "1 days 08:00:00",
                "1 days 08:00:00",
                pd.NaT,
            ]
        )
        self._check_value_counts_dropna(idx)

    def test_value_counts_unique_periodindex2(self):
        idx = pd.PeriodIndex(
            [
                "2013-01-01 09:00",
                "2013-01-01 09:00",
                "2013-01-01 09:00",
                "2013-01-01 08:00",
                "2013-01-01 08:00",
                pd.NaT,
            ],
            freq="h",
        )
        self._check_value_counts_dropna(idx)

    def _check_value_counts_dropna(self, idx):
        exp_idx = idx[[2, 3]]
        expected = pd.Series([3, 2], index=exp_idx, name="count")

        for obj in [idx, pd.Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        exp_idx = idx[[2, 3, -1]]
        expected = pd.Series([3, 2, 1], index=exp_idx, name="count")

        for obj in [idx, pd.Series(idx)]:
            tm.assert_series_equal(obj.value_counts(dropna=False), expected)

        tm.assert_index_equal(idx.unique(), exp_idx)
