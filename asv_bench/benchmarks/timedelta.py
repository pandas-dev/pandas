"""
Timedelta benchmarks with non-tslibs dependencies.  See
benchmarks.tslibs.timedelta for benchmarks that rely only on tslibs.
"""

from pandas import (
    DataFrame,
    Series,
    timedelta_range,
)


class DatetimeAccessor:
    def setup_cache(self):
        N = 100000
        series = Series(timedelta_range("1 days", periods=N, freq="h"))
        return series

    def time_dt_accessor(self, series):
        series.dt

    def time_timedelta_days(self, series):
        series.dt.days

    def time_timedelta_seconds(self, series):
        series.dt.seconds

    def time_timedelta_microseconds(self, series):
        series.dt.microseconds

    def time_timedelta_nanoseconds(self, series):
        series.dt.nanoseconds


class TimedeltaIndexing:
    def setup(self):
        self.index = timedelta_range(start="1985", periods=1000, freq="D")
        self.index2 = timedelta_range(start="1986", periods=1000, freq="D")
        self.series = Series(range(1000), index=self.index)
        self.timedelta = self.index[500]

    def time_get_loc(self):
        self.index.get_loc(self.timedelta)

    def time_shallow_copy(self):
        self.index._view()

    def time_series_loc(self):
        self.series.loc[self.timedelta]

    def time_align(self):
        DataFrame({"a": self.series, "b": self.series[:500]})

    def time_intersection(self):
        self.index.intersection(self.index2)

    def time_union(self):
        self.index.union(self.index2)

    def time_unique(self):
        self.index.unique()
