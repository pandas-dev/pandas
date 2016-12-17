from pandas.tseries.converter import DatetimeConverter
from .pandas_vb_common import *
import pandas as pd
from datetime import timedelta
import datetime as dt
try:
    import pandas.tseries.holiday
except ImportError:
    pass
from pandas.tseries.frequencies import infer_freq
import numpy as np

if hasattr(Series, 'convert'):
    Series.resample = Series.convert


class DatetimeIndex(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.rng = date_range(start='1/1/2000', periods=self.N, freq='T')
        self.delta_offset = pd.offsets.Day()
        self.fast_offset = pd.offsets.DateOffset(months=2, days=2)
        self.slow_offset = pd.offsets.BusinessDay()

        self.rng2 = date_range(start='1/1/2000 9:30', periods=10000, freq='S', tz='US/Eastern')

        self.index_repeated = date_range(start='1/1/2000', periods=1000, freq='T').repeat(10)

        self.rng3 = date_range(start='1/1/2000', periods=1000, freq='H')
        self.df = DataFrame(np.random.randn(len(self.rng3), 2), self.rng3)

        self.rng4 = date_range(start='1/1/2000', periods=1000, freq='H', tz='US/Eastern')
        self.df2 = DataFrame(np.random.randn(len(self.rng4), 2), index=self.rng4)

        N = 100000
        self.dti = pd.date_range('2011-01-01', freq='H', periods=N).repeat(5)
        self.dti_tz = pd.date_range('2011-01-01', freq='H', periods=N,
                                    tz='Asia/Tokyo').repeat(5)

        self.rng5 = date_range(start='1/1/2000', end='3/1/2000', tz='US/Eastern')

        self.dst_rng = date_range(start='10/29/2000 1:00:00', end='10/29/2000 1:59:59', freq='S')
        self.index = date_range(start='10/29/2000', end='10/29/2000 00:59:59', freq='S')
        self.index = self.index.append(self.dst_rng)
        self.index = self.index.append(self.dst_rng)
        self.index = self.index.append(date_range(start='10/29/2000 2:00:00', end='10/29/2000 3:00:00', freq='S'))

        self.N = 10000
        self.rng6 = date_range(start='1/1/1', periods=self.N, freq='B')

        self.rng7 = date_range(start='1/1/1700', freq='D', periods=100000)
        self.a = self.rng7[:50000].append(self.rng7[50002:])

    def time_add_timedelta(self):
        (self.rng + timedelta(minutes=2))

    def time_add_offset_delta(self):
        (self.rng + self.delta_offset)

    def time_add_offset_fast(self):
        (self.rng + self.fast_offset)

    def time_add_offset_slow(self):
        (self.rng + self.slow_offset)

    def time_normalize(self):
        self.rng2.normalize()

    def time_unique(self):
        self.index_repeated.unique()

    def time_reset_index(self):
        self.df.reset_index()

    def time_reset_index_tz(self):
        self.df2.reset_index()

    def time_dti_factorize(self):
        self.dti.factorize()

    def time_dti_tz_factorize(self):
        self.dti_tz.factorize()

    def time_timestamp_tzinfo_cons(self):
        self.rng5[0]

    def time_infer_dst(self):
        self.index.tz_localize('US/Eastern', infer_dst=True)

    def time_timeseries_is_month_start(self):
        self.rng6.is_month_start

    def time_infer_freq(self):
        infer_freq(self.a)


class TimeDatetimeConverter(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.rng = date_range(start='1/1/2000', periods=self.N, freq='T')

    def time_convert(self):
        DatetimeConverter.convert(self.rng, None, None)


class Iteration(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.M = 10000
        self.idx1 = date_range(start='20140101', freq='T', periods=self.N)
        self.idx2 = period_range(start='20140101', freq='T', periods=self.N)

    def iter_n(self, iterable, n=None):
        self.i = 0
        for _ in iterable:
            self.i += 1
            if ((n is not None) and (self.i > n)):
                break

    def time_iter_datetimeindex(self):
        self.iter_n(self.idx1)

    def time_iter_datetimeindex_preexit(self):
        self.iter_n(self.idx1, self.M)

    def time_iter_periodindex(self):
        self.iter_n(self.idx2)

    def time_iter_periodindex_preexit(self):
        self.iter_n(self.idx2, self.M)


#----------------------------------------------------------------------
# Resampling

class ResampleDataFrame(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range(start='20130101', periods=100000, freq='50L')
        self.df = DataFrame(np.random.randn(100000, 2), index=self.rng)

    def time_max_numpy(self):
        self.df.resample('1s', how=np.max)

    def time_max_string(self):
        self.df.resample('1s', how='max')

    def time_mean_numpy(self):
        self.df.resample('1s', how=np.mean)

    def time_mean_string(self):
        self.df.resample('1s', how='mean')

    def time_min_numpy(self):
        self.df.resample('1s', how=np.min)

    def time_min_string(self):
        self.df.resample('1s', how='min')


class ResampleSeries(object):
    goal_time = 0.2

    def setup(self):
        self.rng1 = period_range(start='1/1/2000', end='1/1/2001', freq='T')
        self.ts1 = Series(np.random.randn(len(self.rng1)), index=self.rng1)

        self.rng2 = date_range(start='1/1/2000', end='1/1/2001', freq='T')
        self.ts2 = Series(np.random.randn(len(self.rng2)), index=self.rng2)

        self.rng3 = date_range(start='2000-01-01 00:00:00', end='2000-01-01 10:00:00', freq='555000U')
        self.int_ts = Series(5, self.rng3, dtype='int64')
        self.dt_ts = self.int_ts.astype('datetime64[ns]')

    def time_period_downsample_mean(self):
        self.ts1.resample('D', how='mean')

    def time_timestamp_downsample_mean(self):
        self.ts2.resample('D', how='mean')

    def time_resample_datetime64(self):
        # GH 7754
        self.dt_ts.resample('1S', how='last')

    def time_1min_5min_mean(self):
        self.ts2[:10000].resample('5min', how='mean')

    def time_1min_5min_ohlc(self):
        self.ts2[:10000].resample('5min', how='ohlc')


class AsOf(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.rng = date_range(start='1/1/1990', periods=self.N, freq='53s')
        self.ts = Series(np.random.randn(self.N), index=self.rng)
        self.dates = date_range(start='1/1/1990', periods=(self.N * 10), freq='5s')
        self.ts2 = self.ts.copy()
        self.ts2[250:5000] = np.nan
        self.ts3 = self.ts.copy()
        self.ts3[-5000:] = np.nan

    # test speed of pre-computing NAs.
    def time_asof(self):
        self.ts.asof(self.dates)

    # should be roughly the same as above.
    def time_asof_nan(self):
        self.ts2.asof(self.dates)

    # test speed of the code path for a scalar index
    # without *while* loop
    def time_asof_single(self):
        self.ts.asof(self.dates[0])

    # test speed of the code path for a scalar index
    # before the start. should be the same as above.
    def time_asof_single_early(self):
        self.ts.asof(self.dates[0] - dt.timedelta(10))

    # test the speed of the code path for a scalar index
    # with a long *while* loop. should still be much
    # faster than pre-computing all the NAs.
    def time_asof_nan_single(self):
        self.ts3.asof(self.dates[-1])


class AsOfDataFrame(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.M = 100
        self.rng = date_range(start='1/1/1990', periods=self.N, freq='53s')
        self.dates = date_range(start='1/1/1990', periods=(self.N * 10), freq='5s')
        self.ts = DataFrame(np.random.randn(self.N, self.M), index=self.rng)
        self.ts2 = self.ts.copy()
        self.ts2.iloc[250:5000] = np.nan
        self.ts3 = self.ts.copy()
        self.ts3.iloc[-5000:] = np.nan

    # test speed of pre-computing NAs.
    def time_asof(self):
        self.ts.asof(self.dates)

    # should be roughly the same as above.
    def time_asof_nan(self):
        self.ts2.asof(self.dates)

    # test speed of the code path for a scalar index
    # with pre-computing all NAs.
    def time_asof_single(self):
        self.ts.asof(self.dates[0])

    # should be roughly the same as above.
    def time_asof_nan_single(self):
        self.ts3.asof(self.dates[-1])

    # test speed of the code path for a scalar index
    # before the start. should be without the cost of
    # pre-computing all the NAs.
    def time_asof_single_early(self):
        self.ts.asof(self.dates[0] - dt.timedelta(10))


class TimeSeries(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.rng = date_range(start='1/1/2000', periods=self.N, freq='s')
        self.rng = self.rng.take(np.random.permutation(self.N))
        self.ts = Series(np.random.randn(self.N), index=self.rng)

        self.rng2 = date_range(start='1/1/2000', periods=self.N, freq='T')
        self.ts2 = Series(np.random.randn(self.N), index=self.rng2)

        self.lindex = np.random.permutation(self.N)[:(self.N // 2)]
        self.rindex = np.random.permutation(self.N)[:(self.N // 2)]
        self.left = Series(self.ts2.values.take(self.lindex), index=self.ts2.index.take(self.lindex))
        self.right = Series(self.ts2.values.take(self.rindex), index=self.ts2.index.take(self.rindex))

        self.rng3 = date_range(start='1/1/2000', periods=1500000, freq='S')
        self.ts3 = Series(1, index=self.rng3)

    def time_sort_index(self):
        self.ts.sort_index()

    def time_timeseries_slice_minutely(self):
        self.ts2[:10000]

    def time_add_irregular(self):
        (self.left + self.right)

    def time_large_lookup_value(self):
        self.ts3[self.ts3.index[(len(self.ts3) // 2)]]
        self.ts3.index._cleanup()


class SeriesArithmetic(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.s = Series(date_range(start='20140101', freq='T', periods=self.N))
        self.delta_offset = pd.offsets.Day()
        self.fast_offset = pd.offsets.DateOffset(months=2, days=2)
        self.slow_offset = pd.offsets.BusinessDay()

    def time_add_offset_delta(self):
        (self.s + self.delta_offset)

    def time_add_offset_fast(self):
        (self.s + self.fast_offset)

    def time_add_offset_slow(self):
        (self.s + self.slow_offset)


class ToDatetime(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range(start='1/1/2000', periods=10000, freq='D')
        self.stringsD = Series((((self.rng.year * 10000) + (self.rng.month * 100)) + self.rng.day), dtype=np.int64).apply(str)

        self.rng = date_range(start='1/1/2000', periods=20000, freq='H')
        self.strings = [x.strftime('%Y-%m-%d %H:%M:%S') for x in self.rng]
        self.strings_nosep = [x.strftime('%Y%m%d %H:%M:%S') for x in self.rng]
        self.strings_tz_space = [x.strftime('%Y-%m-%d %H:%M:%S') + ' -0800'
                                 for x in self.rng]

        self.s = Series((['19MAY11', '19MAY11:00:00:00'] * 100000))
        self.s2 = self.s.str.replace(':\\S+$', '')

    def time_format_YYYYMMDD(self):
        to_datetime(self.stringsD, format='%Y%m%d')

    def time_iso8601(self):
        to_datetime(self.strings)

    def time_iso8601_nosep(self):
        to_datetime(self.strings_nosep)

    def time_iso8601_format(self):
        to_datetime(self.strings, format='%Y-%m-%d %H:%M:%S')

    def time_iso8601_format_no_sep(self):
        to_datetime(self.strings_nosep, format='%Y%m%d %H:%M:%S')

    def time_iso8601_tz_spaceformat(self):
        to_datetime(self.strings_tz_space)

    def time_format_exact(self):
        to_datetime(self.s2, format='%d%b%y')

    def time_format_no_exact(self):
        to_datetime(self.s, format='%d%b%y', exact=False)


class Offsets(object):
    goal_time = 0.2

    def setup(self):
        self.date = dt.datetime(2011, 1, 1)
        self.dt64 = np.datetime64('2011-01-01 09:00Z')
        self.hcal = pd.tseries.holiday.USFederalHolidayCalendar()
        self.day = pd.offsets.Day()
        self.year = pd.offsets.YearBegin()
        self.cday = pd.offsets.CustomBusinessDay()
        self.cmb = pd.offsets.CustomBusinessMonthBegin(calendar=self.hcal)
        self.cme = pd.offsets.CustomBusinessMonthEnd(calendar=self.hcal)
        self.cdayh = pd.offsets.CustomBusinessDay(calendar=self.hcal)

    def time_timeseries_day_apply(self):
        self.day.apply(self.date)

    def time_timeseries_day_incr(self):
        (self.date + self.day)

    def time_timeseries_year_apply(self):
        self.year.apply(self.date)

    def time_timeseries_year_incr(self):
        (self.date + self.year)

    # custom business offsets

    def time_custom_bday_decr(self):
        (self.date - self.cday)

    def time_custom_bday_incr(self):
        (self.date + self.cday)

    def time_custom_bday_apply(self):
        self.cday.apply(self.date)

    def time_custom_bday_apply_dt64(self):
        self.cday.apply(self.dt64)

    def time_custom_bday_cal_incr(self):
        self.date + 1 * self.cdayh

    def time_custom_bday_cal_decr(self):
        self.date - 1 * self.cdayh

    def time_custom_bday_cal_incr_n(self):
        self.date + 10 * self.cdayh

    def time_custom_bday_cal_incr_neg_n(self):
        self.date - 10 * self.cdayh

    # Increment custom business month

    def time_custom_bmonthend_incr(self):
        (self.date + self.cme)

    def time_custom_bmonthend_incr_n(self):
        (self.date + (10 * self.cme))

    def time_custom_bmonthend_decr_n(self):
        (self.date - (10 * self.cme))

    def time_custom_bmonthbegin_decr_n(self):
        (self.date - (10 * self.cmb))

    def time_custom_bmonthbegin_incr_n(self):
        (self.date + (10 * self.cmb))


class SemiMonthOffset(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.rng = date_range(start='1/1/2000', periods=self.N, freq='T')
        # date is not on an offset which will be slowest case
        self.date = dt.datetime(2011, 1, 2)
        self.semi_month_end = pd.offsets.SemiMonthEnd()
        self.semi_month_begin = pd.offsets.SemiMonthBegin()

    def time_end_apply(self):
        self.semi_month_end.apply(self.date)

    def time_end_incr(self):
        self.date + self.semi_month_end

    def time_end_incr_n(self):
        self.date + 10 * self.semi_month_end

    def time_end_decr(self):
        self.date - self.semi_month_end

    def time_end_decr_n(self):
        self.date - 10 * self.semi_month_end

    def time_end_apply_index(self):
        self.semi_month_end.apply_index(self.rng)

    def time_end_incr_rng(self):
        self.rng + self.semi_month_end

    def time_end_decr_rng(self):
        self.rng - self.semi_month_end

    def time_begin_apply(self):
        self.semi_month_begin.apply(self.date)

    def time_begin_incr(self):
        self.date + self.semi_month_begin

    def time_begin_incr_n(self):
        self.date + 10 * self.semi_month_begin

    def time_begin_decr(self):
        self.date - self.semi_month_begin

    def time_begin_decr_n(self):
        self.date - 10 * self.semi_month_begin

    def time_begin_apply_index(self):
        self.semi_month_begin.apply_index(self.rng)

    def time_begin_incr_rng(self):
        self.rng + self.semi_month_begin

    def time_begin_decr_rng(self):
        self.rng - self.semi_month_begin
