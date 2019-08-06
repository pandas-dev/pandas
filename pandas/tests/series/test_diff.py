from numpy import nan
from pandas import Series, date_range
from pandas.tests.series.common import TestData
from pandas.core.indexes.timedeltas import TimedeltaIndex
from pandas.util.testing import assert_series_equal


class TestDiff(TestData):
    def test_ts_diff(self):
        # Just run the function
        self.ts.diff()

        # int dtype
        a = 10000000000000000
        b = a + 1
        s = Series([a, b])

        rs = s.diff()
        assert rs[1] == 1

        # neg n
        rs = self.ts.diff(-1)
        xp = self.ts - self.ts.shift(-1)
        assert_series_equal(rs, xp)

        # 0
        rs = self.ts.diff(0)
        xp = self.ts - self.ts
        assert_series_equal(rs, xp)

        # datetime diff (GH3100)
        s = Series(date_range("20130102", periods=5))
        rs = s - s.shift(1)
        xp = s.diff()
        assert_series_equal(rs, xp)

        # timedelta diff
        nrs = rs - rs.shift(1)
        nxp = xp.diff()
        assert_series_equal(nrs, nxp)

        # with tz
        s = Series(
            date_range("2000-01-01 09:00:00", periods=5, tz="US/Eastern"), name="foo"
        )
        result = s.diff()
        assert_series_equal(
            result, Series(TimedeltaIndex(["NaT"] + ["1 days"] * 4), name="foo")
        )
        
    def test_boolean_diff(self):
        # boolean series
        s = Series([False, True, True, False, False])
        result = s.diff()
        assert_series_equal(result, Series[nan, True, False, True, False])
