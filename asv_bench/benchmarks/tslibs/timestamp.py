from datetime import datetime

import numpy as np
import pytz

from pandas import Timestamp

from .tslib import _tzs


class TimestampConstruction:
    def setup(self):
        self.npdatetime64 = np.datetime64("2020-01-01 00:00:00")
        self.dttime_unaware = datetime(2020, 1, 1, 0, 0, 0)
        self.dttime_aware = datetime(2020, 1, 1, 0, 0, 0, 0, pytz.UTC)
        self.ts = Timestamp("2020-01-01 00:00:00")

    def time_parse_iso8601_no_tz(self):
        Timestamp("2017-08-25 08:16:14")

    def time_parse_iso8601_tz(self):
        Timestamp("2017-08-25 08:16:14-0500")

    def time_parse_dateutil(self):
        Timestamp("2017/08/25 08:16:14 AM")

    def time_parse_today(self):
        Timestamp("today")

    def time_parse_now(self):
        Timestamp("now")

    def time_fromordinal(self):
        Timestamp.fromordinal(730120)

    def time_fromtimestamp(self):
        Timestamp.fromtimestamp(1515448538)

    def time_from_npdatetime64(self):
        Timestamp(self.npdatetime64)

    def time_from_datetime_unaware(self):
        Timestamp(self.dttime_unaware)

    def time_from_datetime_aware(self):
        Timestamp(self.dttime_aware)

    def time_from_pd_timestamp(self):
        Timestamp(self.ts)


class TimestampProperties:
    params = [_tzs]
    param_names = ["tz"]

    def setup(self, tz):
        self.ts = Timestamp("2017-08-25 08:16:14", tzinfo=tz)

    def time_tz(self, tz):
        self.ts.tz

    def time_dayofweek(self, tz):
        self.ts.dayofweek

    def time_dayofyear(self, tz):
        self.ts.dayofyear

    def time_week(self, tz):
        self.ts.week

    def time_quarter(self, tz):
        self.ts.quarter

    def time_days_in_month(self, tz):
        self.ts.days_in_month

    def time_is_month_start(self, tz):
        self.ts.is_month_start

    def time_is_month_end(self, tz):
        self.ts.is_month_end

    def time_is_quarter_start(self, tz):
        self.ts.is_quarter_start

    def time_is_quarter_end(self, tz):
        self.ts.is_quarter_end

    def time_is_year_start(self, tz):
        self.ts.is_year_start

    def time_is_year_end(self, tz):
        self.ts.is_year_end

    def time_is_leap_year(self, tz):
        self.ts.is_leap_year

    def time_microsecond(self, tz):
        self.ts.microsecond

    def time_month_name(self, tz):
        self.ts.month_name()

    def time_weekday_name(self, tz):
        self.ts.day_name()


class TimestampOps:
    params = _tzs
    param_names = ["tz"]

    def setup(self, tz):
        self.ts = Timestamp("2017-08-25 08:16:14", tz=tz)

    def time_replace_tz(self, tz):
        self.ts.replace(tzinfo=pytz.timezone("US/Eastern"))

    def time_replace_None(self, tz):
        self.ts.replace(tzinfo=None)

    def time_to_pydatetime(self, tz):
        self.ts.to_pydatetime()

    def time_normalize(self, tz):
        self.ts.normalize()

    def time_tz_convert(self, tz):
        if self.ts.tz is not None:
            self.ts.tz_convert(tz)

    def time_tz_localize(self, tz):
        if self.ts.tz is None:
            self.ts.tz_localize(tz)

    def time_to_julian_date(self, tz):
        self.ts.to_julian_date()

    def time_floor(self, tz):
        self.ts.floor("5T")

    def time_ceil(self, tz):
        self.ts.ceil("5T")


class TimestampAcrossDst:
    def setup(self):
        dt = datetime(2016, 3, 27, 1)
        self.tzinfo = pytz.timezone("CET").localize(dt, is_dst=False).tzinfo
        self.ts2 = Timestamp(dt)

    def time_replace_across_dst(self):
        self.ts2.replace(tzinfo=self.tzinfo)
