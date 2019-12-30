"""
Period benchmarks with non-tslibs dependencies.  See
benchmarks.tslibs.period for benchmarks that rely only on tslibs.
"""
from pandas import DataFrame, Period, PeriodIndex, Series, date_range, period_range

from pandas.tseries.frequencies import to_offset


class PeriodIndexConstructor:

    params = [["D"], [True, False]]
    param_names = ["freq", "is_offset"]

    def setup(self, freq, is_offset):
        self.rng = date_range("1985", periods=1000)
        self.rng2 = date_range("1985", periods=1000).to_pydatetime()
        self.ints = list(range(2000, 3000))
        self.daily_ints = (
            date_range("1/1/2000", periods=1000, freq=freq).strftime("%Y%m%d").map(int)
        )
        if is_offset:
            self.freq = to_offset(freq)
        else:
            self.freq = freq

    def time_from_date_range(self, freq, is_offset):
        PeriodIndex(self.rng, freq=freq)

    def time_from_pydatetime(self, freq, is_offset):
        PeriodIndex(self.rng2, freq=freq)

    def time_from_ints(self, freq, is_offset):
        PeriodIndex(self.ints, freq=freq)

    def time_from_ints_daily(self, freq, is_offset):
        PeriodIndex(self.daily_ints, freq=freq)


class DataFramePeriodColumn:
    def setup(self):
        self.rng = period_range(start="1/1/1990", freq="S", periods=20000)
        self.df = DataFrame(index=range(len(self.rng)))

    def time_setitem_period_column(self):
        self.df["col"] = self.rng

    def time_set_index(self):
        # GH#21582 limited by comparisons of Period objects
        self.df["col2"] = self.rng
        self.df.set_index("col2", append=True)


class Algorithms:

    params = ["index", "series"]
    param_names = ["typ"]

    def setup(self, typ):
        data = [
            Period("2011-01", freq="M"),
            Period("2011-02", freq="M"),
            Period("2011-03", freq="M"),
            Period("2011-04", freq="M"),
        ]

        if typ == "index":
            self.vector = PeriodIndex(data * 1000, freq="M")
        elif typ == "series":
            self.vector = Series(data * 1000)

    def time_drop_duplicates(self, typ):
        self.vector.drop_duplicates()

    def time_value_counts(self, typ):
        self.vector.value_counts()


class Indexing:
    def setup(self):
        self.index = period_range(start="1985", periods=1000, freq="D")
        self.series = Series(range(1000), index=self.index)
        self.period = self.index[500]

    def time_get_loc(self):
        self.index.get_loc(self.period)

    def time_shape(self):
        self.index.shape

    def time_shallow_copy(self):
        self.index._shallow_copy()

    def time_series_loc(self):
        self.series.loc[self.period]

    def time_align(self):
        DataFrame({"a": self.series, "b": self.series[:500]})

    def time_intersection(self):
        self.index[:750].intersection(self.index[250:])

    def time_unique(self):
        self.index.unique()
