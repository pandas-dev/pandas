import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestTimedeltaIndex:
    def test_misc_coverage(self):
        rng = pd.timedelta_range("1 day", periods=5)
        result = rng.groupby(rng.days)
        assert isinstance(next(iter(result.values()))[0], pd.Timedelta)

    def test_map(self):
        # test_map_dictlike generally tests

        rng = pd.timedelta_range("1 day", periods=10)

        f = lambda x: x.days
        result = rng.map(f)
        exp = pd.Index([f(x) for x in rng], dtype=np.int64)
        tm.assert_index_equal(result, exp)

    def test_fields(self):
        rng = pd.timedelta_range("1 days, 10:11:12.100123456", periods=2, freq="s")
        tm.assert_index_equal(rng.days, pd.Index([1, 1], dtype=np.int64))
        tm.assert_index_equal(
            rng.seconds,
            pd.Index(
                [10 * 3600 + 11 * 60 + 12, 10 * 3600 + 11 * 60 + 13], dtype=np.int32
            ),
        )
        tm.assert_index_equal(
            rng.microseconds,
            pd.Index([100 * 1000 + 123, 100 * 1000 + 123], dtype=np.int32),
        )
        tm.assert_index_equal(rng.nanoseconds, pd.Index([456, 456], dtype=np.int32))

        msg = "'TimedeltaIndex' object has no attribute '{}'"
        with pytest.raises(AttributeError, match=msg.format("hours")):
            rng.hours
        with pytest.raises(AttributeError, match=msg.format("minutes")):
            rng.minutes
        with pytest.raises(AttributeError, match=msg.format("milliseconds")):
            rng.milliseconds

        # with nat
        s = pd.Series(rng)
        s[1] = np.nan

        tm.assert_series_equal(s.dt.days, pd.Series([1, np.nan], index=range(2)))
        tm.assert_series_equal(
            s.dt.seconds, pd.Series([10 * 3600 + 11 * 60 + 12, np.nan], index=range(2))
        )

        # preserve name (GH15589)
        rng.name = "name"
        assert rng.days.name == "name"
