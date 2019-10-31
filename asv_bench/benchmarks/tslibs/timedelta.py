"""
Timedelta benchmarks that rely only on tslibs.  See benchmarks.timedeltas for
Timedelta benchmarks that rely on other parts fo pandas.
"""
import datetime

import numpy as np

from pandas import Timedelta


class TimedeltaConstructor:
    def time_from_int(self):
        Timedelta(123456789)

    def time_from_unit(self):
        Timedelta(1, unit="d")

    def time_from_components(self):
        Timedelta(
            days=1,
            hours=2,
            minutes=3,
            seconds=4,
            milliseconds=5,
            microseconds=6,
            nanoseconds=7,
        )

    def time_from_datetime_timedelta(self):
        Timedelta(datetime.timedelta(days=1, seconds=1))

    def time_from_np_timedelta(self):
        Timedelta(np.timedelta64(1, "ms"))

    def time_from_string(self):
        Timedelta("1 days")

    def time_from_iso_format(self):
        Timedelta("P4DT12H30M5S")

    def time_from_missing(self):
        Timedelta("nat")


class TimedeltaProperties:
    def setup_cache(self):
        td = Timedelta(days=365, minutes=35, seconds=25, milliseconds=35)
        return td

    def time_timedelta_days(self, td):
        td.days

    def time_timedelta_seconds(self, td):
        td.seconds

    def time_timedelta_microseconds(self, td):
        td.microseconds

    def time_timedelta_nanoseconds(self, td):
        td.nanoseconds
