"""
Period benchmarks that rely only on tslibs.  See benchmarks.period for
Period benchmarks that rely on other parts fo pandas.
"""
from pandas import Period

from pandas.tseries.frequencies import to_offset


class PeriodProperties:

    params = (
        ["M", "min"],
        [
            "year",
            "month",
            "day",
            "hour",
            "minute",
            "second",
            "is_leap_year",
            "quarter",
            "qyear",
            "week",
            "daysinmonth",
            "dayofweek",
            "dayofyear",
            "start_time",
            "end_time",
        ],
    )
    param_names = ["freq", "attr"]

    def setup(self, freq, attr):
        self.per = Period("2012-06-01", freq=freq)

    def time_property(self, freq, attr):
        getattr(self.per, attr)


class PeriodUnaryMethods:

    params = ["M", "min"]
    param_names = ["freq"]

    def setup(self, freq):
        self.per = Period("2012-06-01", freq=freq)

    def time_to_timestamp(self, freq):
        self.per.to_timestamp()

    def time_now(self, freq):
        self.per.now(freq)

    def time_asfreq(self, freq):
        self.per.asfreq("A")


class PeriodConstructor:
    params = [["D"], [True, False]]
    param_names = ["freq", "is_offset"]

    def setup(self, freq, is_offset):
        if is_offset:
            self.freq = to_offset(freq)
        else:
            self.freq = freq

    def time_period_constructor(self, freq, is_offset):
        Period("2012-06-01", freq=freq)
