import datetime

from pandas import Timestamp
import pytz


class TimestampConstruction(object):

    def time_parse_iso8601_no_tz(self):
        Timestamp('2017-08-25 08:16:14')

    def time_parse_iso8601_tz(self):
        Timestamp('2017-08-25 08:16:14-0500')

    def time_parse_dateutil(self):
        Timestamp('2017/08/25 08:16:14 AM')

    def time_parse_today(self):
        Timestamp('today')

    def time_parse_now(self):
        Timestamp('now')

    def time_fromordinal(self):
        Timestamp.fromordinal(730120)

    def time_fromtimestamp(self):
        Timestamp.fromtimestamp(1515448538)


class TimestampProperties(object):
    goal_time = 0.2

    _tzs = [None, pytz.timezone('Europe/Amsterdam')]
    _freqs = [None, 'B']
    params = [_tzs, _freqs]
    param_names = ['tz', 'freq']

    def setup(self, tz, freq):
        self.ts = Timestamp('2017-08-25 08:16:14', tzinfo=tz, freq=freq)

    def time_tz(self, tz, freq):
        self.ts.tz

    def time_dayofweek(self, tz, freq):
        self.ts.dayofweek

    def time_weekday_name(self, tz, freq):
        self.ts.weekday_name

    def time_dayofyear(self, tz, freq):
        self.ts.dayofyear

    def time_week(self, tz, freq):
        self.ts.week

    def time_quarter(self, tz, freq):
        self.ts.quarter

    def time_days_in_month(self, tz, freq):
        self.ts.days_in_month

    def time_freqstr(self, tz, freq):
        self.ts.freqstr

    def time_is_month_start(self, tz, freq):
        self.ts.is_month_start

    def time_is_month_end(self, tz, freq):
        self.ts.is_month_end

    def time_is_quarter_start(self, tz, freq):
        self.ts.is_quarter_start

    def time_is_quarter_end(self, tz, freq):
        self.ts.is_quarter_end

    def time_is_year_start(self, tz, freq):
        self.ts.is_quarter_end

    def time_is_year_end(self, tz, freq):
        self.ts.is_quarter_end

    def time_is_leap_year(self, tz, freq):
        self.ts.is_quarter_end

    def time_microsecond(self, tz, freq):
        self.ts.microsecond


class TimestampOps(object):
    goal_time = 0.2

    params = [None, 'US/Eastern']
    param_names = ['tz']

    def setup(self, tz):
        self.ts = Timestamp('2017-08-25 08:16:14', tz=tz)

    def time_replace_tz(self, tz):
        self.ts.replace(tzinfo=pytz.timezone('US/Eastern'))

    def time_replace_None(self, tz):
        self.ts.replace(tzinfo=None)

    def time_to_pydatetime(self, tz):
        self.ts.to_pydatetime()


class TimestampAcrossDst(object):
    goal_time = 0.2

    def setup(self):
        dt = datetime.datetime(2016, 3, 27, 1)
        self.tzinfo = pytz.timezone('CET').localize(dt, is_dst=False).tzinfo
        self.ts2 = Timestamp(dt)

    def time_replace_across_dst(self):
        self.ts2.replace(tzinfo=self.tzinfo)
