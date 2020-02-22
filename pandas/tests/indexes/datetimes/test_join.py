from datetime import datetime

import pytest

from pandas import DatetimeIndex, Index, date_range, to_datetime
import pandas._testing as tm

from pandas.tseries.offsets import BMonthEnd


class TestJoin:
    def test_join_utc_convert(self, join_type):
        rng = date_range("1/1/2011", periods=100, freq="H", tz="utc")

        left = rng.tz_convert("US/Eastern")
        right = rng.tz_convert("Europe/Berlin")

        result = left.join(left[:-5], how=join_type)
        assert isinstance(result, DatetimeIndex)
        assert result.tz == left.tz

        result = left.join(right[:-5], how=join_type)
        assert isinstance(result, DatetimeIndex)
        assert result.tz.zone == "UTC"

    @pytest.mark.parametrize("sort", [None, False])
    def test_datetimeindex_union_join_empty(self, sort):
        dti = date_range(start="1/1/2001", end="2/1/2001", freq="D")
        empty = Index([])

        result = dti.union(empty, sort=sort)
        expected = dti.astype("O")
        tm.assert_index_equal(result, expected)

        result = dti.join(empty)
        assert isinstance(result, DatetimeIndex)
        tm.assert_index_equal(result, dti)

    def test_join_nonunique(self):
        idx1 = to_datetime(["2012-11-06 16:00:11.477563", "2012-11-06 16:00:11.477563"])
        idx2 = to_datetime(["2012-11-06 15:11:09.006507", "2012-11-06 15:11:09.006507"])
        rs = idx1.join(idx2, how="outer")
        assert rs.is_monotonic

    @pytest.mark.parametrize("freq", ["B", "C"])
    def test_outer_join(self, freq):
        # should just behave as union
        start, end = datetime(2009, 1, 1), datetime(2010, 1, 1)
        rng = date_range(start=start, end=end, freq=freq)

        # overlapping
        left = rng[:10]
        right = rng[5:10]

        the_join = left.join(right, how="outer")
        assert isinstance(the_join, DatetimeIndex)

        # non-overlapping, gap in middle
        left = rng[:5]
        right = rng[10:]

        the_join = left.join(right, how="outer")
        assert isinstance(the_join, DatetimeIndex)
        assert the_join.freq is None

        # non-overlapping, no gap
        left = rng[:5]
        right = rng[5:10]

        the_join = left.join(right, how="outer")
        assert isinstance(the_join, DatetimeIndex)

        # overlapping, but different offset
        rng = date_range(start, end, freq=BMonthEnd())

        the_join = rng.join(rng, how="outer")
        assert isinstance(the_join, DatetimeIndex)
        assert the_join.freq is None
