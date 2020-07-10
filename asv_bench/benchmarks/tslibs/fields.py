import numpy as np

from pandas._libs.tslibs.fields import (
    get_date_field,
    get_start_end_field,
    get_timedelta_field,
)

from .tslib import _sizes


class TimeGetTimedeltaField:
    params = [
        _sizes,
        ["days", "h", "s", "seconds", "ms", "microseconds", "us", "ns", "nanoseconds"],
    ]
    param_names = ["size", "field"]

    def setup(self, size, field):
        arr = np.random.randint(0, 10, size=size, dtype="i8")
        self.i8data = arr

    def time_get_timedelta_field(self, size, field):
        get_timedelta_field(self.i8data, field)


class TimeGetDateField:
    params = [
        _sizes,
        [
            "Y",
            "M",
            "D",
            "h",
            "m",
            "s",
            "us",
            "ns",
            "doy",
            "dow",
            "woy",
            "q",
            "dim",
            "is_leap_year",
        ],
    ]
    param_names = ["size", "field"]

    def setup(self, size, field):
        arr = np.random.randint(0, 10, size=size, dtype="i8")
        self.i8data = arr

    def time_get_date_field(self, size, field):
        get_date_field(self.i8data, field)


class TimeGetStartEndField:
    params = [
        _sizes,
        ["start", "end"],
        ["month", "quarter", "year"],
        ["B", None, "QS"],
        [12, 3, 5],
    ]
    param_names = ["size", "side", "period", "freqstr", "month_kw"]

    def setup(self, size, side, period, freqstr, month_kw):
        arr = np.random.randint(0, 10, size=size, dtype="i8")
        self.i8data = arr

        self.attrname = f"is_{period}_{side}"

    def time_get_start_end_field(self, size, side, period, freqstr, month_kw):
        get_start_end_field(self.i8data, self.attrname, freqstr, month_kw=month_kw)
