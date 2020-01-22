"""
Timedelta benchmarks with non-tslibs dependencies.  See
benchmarks.tslibs.timedelta for benchmarks that rely only on tslibs.
"""

import numpy as np

from pandas import DataFrame, Series, Timestamp, timedelta_range, to_timedelta


class ToTimedelta:
    def setup(self):
        self.ints = np.random.randint(0, 60, size=10000)
        self.str_days = []
        self.str_seconds = []
        for i in self.ints:
            self.str_days.append(f"{i} days")
            self.str_seconds.append(f"00:00:{i:02d}")

    def time_convert_int(self):
        to_timedelta(self.ints, unit="s")

    def time_convert_string_days(self):
        to_timedelta(self.str_days)

    def time_convert_string_seconds(self):
        to_timedelta(self.str_seconds)


class ToTimedeltaErrors:

    params = ["coerce", "ignore"]
    param_names = ["errors"]

    def setup(self, errors):
        ints = np.random.randint(0, 60, size=10000)
        self.arr = [f"{i} days" for i in ints]
        self.arr[-1] = "apple"

    def time_convert(self, errors):
        to_timedelta(self.arr, errors=errors)


class TimedeltaOps:
    def setup(self):
        self.td = to_timedelta(np.arange(1000000))
        self.ts = Timestamp("2000")

    def time_add_td_ts(self):
        self.td + self.ts


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

    def time_shape(self):
        self.index.shape

    def time_shallow_copy(self):
        self.index._shallow_copy()

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
