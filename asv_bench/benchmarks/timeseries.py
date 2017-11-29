try:
    from pandas.plotting._converter import DatetimeConverter
except ImportError:
    from pandas.tseries.converter import DatetimeConverter

import pandas as pd
from pandas import to_datetime, date_range, Series, DataFrame, period_range

import datetime as dt
from pandas.tseries.frequencies import infer_freq
import numpy as np

if hasattr(Series, 'convert'):
    Series.resample = Series.convert


class DatetimeIndex(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.rng = date_range(start='1/1/2000', periods=self.N, freq='T')

        self.rng2 = date_range(start='1/1/2000 9:30', periods=10000,
                               freq='S', tz='US/Eastern')

        self.index_repeated = date_range(start='1/1/2000',
                                         periods=1000, freq='T').repeat(10)

        self.rng3 = date_range(start='1/1/2000', periods=1000, freq='H')
        self.df = DataFrame(np.random.randn(len(self.rng3), 2), self.rng3)

        self.rng4 = date_range(start='1/1/2000', periods=1000,
                               freq='H', tz='US/Eastern')
        self.df2 = DataFrame(np.random.randn(len(self.rng4), 2),
                             index=self.rng4)

        N = 100000
        self.dti = pd.date_range('2011-01-01', freq='H', periods=N).repeat(5)
        self.dti_tz = pd.date_range('2011-01-01', freq='H', periods=N,
                                    tz='Asia/Tokyo').repeat(5)

        self.rng5 = date_range(start='1/1/2000',
                               end='3/1/2000', tz='US/Eastern')

        self.dst_rng = date_range(start='10/29/2000 1:00:00',
                                  end='10/29/2000 1:59:59', freq='S')
        self.index = date_range(start='10/29/2000',
                                end='10/29/2000 00:59:59', freq='S')
        self.index = self.index.append(self.dst_rng)
        self.index = self.index.append(self.dst_rng)
        self.index = self.index.append(date_range(start='10/29/2000 2:00:00',
                                                  end='10/29/2000 3:00:00',
                                                  freq='S'))

        self.N = 10000
        self.rng6 = date_range(start='1/1/1', periods=self.N, freq='B')

        self.rng7 = date_range(start='1/1/1700', freq='D', periods=100000)
        self.no_freq = self.rng7[:50000].append(self.rng7[50002:])
        self.d_freq = self.rng7[:50000].append(self.rng7[50000:])

        self.rng8 = date_range(start='1/1/1700', freq='B', periods=75000)
        self.b_freq = self.rng8[:50000].append(self.rng8[50000:])

    def time_add_timedelta(self):
        (self.rng + dt.timedelta(minutes=2))

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

    def time_dti_time(self):
        self.rng.time

    def time_timestamp_tzinfo_cons(self):
        self.rng5[0]

    def time_infer_dst(self):
        self.index.tz_localize('US/Eastern', infer_dst=True)

    def time_timeseries_is_month_start(self):
        self.rng6.is_month_start

    def time_infer_freq_none(self):
        infer_freq(self.no_freq)

    def time_infer_freq_daily(self):
        infer_freq(self.d_freq)

    def time_infer_freq_business(self):
        infer_freq(self.b_freq)

    def time_to_date(self):
        self.rng.date

    def time_to_pydatetime(self):
        self.rng.to_pydatetime()


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


# ----------------------------------------------------------------------
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

        self.rng3 = date_range(start='2000-01-01 00:00:00',
                               end='2000-01-01 10:00:00', freq='555000U')
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
        self.dates = date_range(start='1/1/1990',
                                periods=(self.N * 10), freq='5s')
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
        self.dates = date_range(start='1/1/1990',
                                periods=(self.N * 10), freq='5s')
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
        self.left = Series(self.ts2.values.take(self.lindex),
                           index=self.ts2.index.take(self.lindex))
        self.right = Series(self.ts2.values.take(self.rindex),
                            index=self.ts2.index.take(self.rindex))

        self.rng3 = date_range(start='1/1/2000', periods=1500000, freq='S')
        self.ts3 = Series(1, index=self.rng3)

    def time_sort_index_monotonic(self):
        self.ts2.sort_index()

    def time_sort_index_non_monotonic(self):
        self.ts.sort_index()

    def time_timeseries_slice_minutely(self):
        self.ts2[:10000]

    def time_add_irregular(self):
        (self.left + self.right)

    def time_large_lookup_value(self):
        self.ts3[self.ts3.index[(len(self.ts3) // 2)]]
        self.ts3.index._cleanup()


class ToDatetime(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range(start='1/1/2000', periods=10000, freq='D')
        self.stringsD = Series(self.rng.strftime('%Y%m%d'))

        self.rng = date_range(start='1/1/2000', periods=20000, freq='H')
        self.strings = self.rng.strftime('%Y-%m-%d %H:%M:%S').tolist()
        self.strings_nosep = self.rng.strftime('%Y%m%d %H:%M:%S').tolist()
        self.strings_tz_space = [x.strftime('%Y-%m-%d %H:%M:%S') + ' -0800'
                                 for x in self.rng]

        self.s = Series((['19MAY11', '19MAY11:00:00:00'] * 100000))
        self.s2 = self.s.str.replace(':\\S+$', '')

        self.unique_numeric_seconds = range(10000)
        self.dup_numeric_seconds = [1000] * 10000
        self.dup_string_dates = ['2000-02-11'] * 10000
        self.dup_string_with_tz = ['2000-02-11 15:00:00-0800'] * 10000

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

    def time_cache_true_with_unique_seconds_and_unit(self):
        to_datetime(self.unique_numeric_seconds, unit='s', cache=True)

    def time_cache_false_with_unique_seconds_and_unit(self):
        to_datetime(self.unique_numeric_seconds, unit='s', cache=False)

    def time_cache_true_with_dup_seconds_and_unit(self):
        to_datetime(self.dup_numeric_seconds, unit='s', cache=True)

    def time_cache_false_with_dup_seconds_and_unit(self):
        to_datetime(self.dup_numeric_seconds, unit='s', cache=False)

    def time_cache_true_with_dup_string_dates(self):
        to_datetime(self.dup_string_dates, cache=True)

    def time_cache_false_with_dup_string_dates(self):
        to_datetime(self.dup_string_dates, cache=False)

    def time_cache_true_with_dup_string_dates_and_format(self):
        to_datetime(self.dup_string_dates, format='%Y-%m-%d', cache=True)

    def time_cache_false_with_dup_string_dates_and_format(self):
        to_datetime(self.dup_string_dates, format='%Y-%m-%d', cache=False)

    def time_cache_true_with_dup_string_tzoffset_dates(self):
        to_datetime(self.dup_string_with_tz, cache=True)

    def time_cache_false_with_dup_string_tzoffset_dates(self):
        to_datetime(self.dup_string_with_tz, cache=False)


class DatetimeAccessor(object):
    def setup(self):
        self.N = 100000
        self.series = pd.Series(
            pd.date_range(start='1/1/2000', periods=self.N, freq='T')
        )

    def time_dt_accessor(self):
        self.series.dt

    def time_dt_accessor_normalize(self):
        self.series.dt.normalize()
