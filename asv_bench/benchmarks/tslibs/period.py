"""
Period benchmarks that rely only on tslibs.  See benchmarks.period for
Period benchmarks that rely on other parts fo pandas.
"""

import numpy as np

from pandas._libs.tslibs.period import (
    Period,
    periodarr_to_dt64arr,
)

from pandas.tseries.frequencies import to_offset

from .tslib import (
    _sizes,
    _tzs,
    tzlocal_obj,
)

try:
    from pandas._libs.tslibs.vectorized import dt64arr_to_periodarr
except ImportError:
    from pandas._libs.tslibs.period import dt64arr_to_periodarr


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


_freq_ints = [
    1000,
    1011,  # Annual - November End
    2000,
    2011,  # Quarterly - November End
    3000,
    4000,
    4006,  # Weekly - Saturday End
    5000,
    6000,
    7000,
    8000,
    9000,
    10000,
    11000,
    12000,
]


class TimePeriodArrToDT64Arr:
    params = [
        _sizes,
        _freq_ints,
    ]
    param_names = ["size", "freq"]

    def setup(self, size, freq):
        arr = np.arange(10, dtype="i8").repeat(size // 10)
        self.i8values = arr

    def time_periodarray_to_dt64arr(self, size, freq):
        periodarr_to_dt64arr(self.i8values, freq)


class TimeDT64ArrToPeriodArr:
    params = [
        _sizes,
        _freq_ints,
        _tzs,
    ]
    param_names = ["size", "freq", "tz"]

    def setup(self, size, freq, tz):
        if size == 10 ** 6 and tz is tzlocal_obj:
            # tzlocal is cumbersomely slow, so skip to keep runtime in check
            raise NotImplementedError

        arr = np.arange(10, dtype="i8").repeat(size // 10)
        self.i8values = arr

    def time_dt64arr_to_periodarr(self, size, freq, tz):
        dt64arr_to_periodarr(self.i8values, freq, tz)
