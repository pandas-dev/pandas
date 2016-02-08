from vbench.api import Benchmark
from datetime import datetime
from pandas import *

N = 100000
try:
    rng = date_range(start='1/1/2000', periods=N, freq='min')
except NameError:
    rng = DatetimeIndex(start='1/1/2000', periods=N, freq='T')
    def date_range(start=None, end=None, periods=None, freq=None):
        return DatetimeIndex(start=start, end=end, periods=periods, offset=freq)


common_setup = """from .pandas_vb_common import *
from datetime import timedelta
N = 100000

rng = date_range(start='1/1/2000', periods=N, freq='T')

if hasattr(Series, 'convert'):
    Series.resample = Series.convert

ts = Series(np.random.randn(N), index=rng)
"""

#----------------------------------------------------------------------
# Lookup value in large time series, hash map population

setup = common_setup + """
rng = date_range(start='1/1/2000', periods=1500000, freq='S')
ts = Series(1, index=rng)
"""

stmt = "ts[ts.index[len(ts) // 2]]; ts.index._cleanup()"
timeseries_large_lookup_value = Benchmark(stmt, setup,
                                          start_date=datetime(2012, 1, 1))

#----------------------------------------------------------------------
# Test slice minutely series

timeseries_slice_minutely = Benchmark('ts[:10000]', common_setup)

#----------------------------------------------------------------------
# Test conversion

setup = common_setup + """

"""

timeseries_1min_5min_ohlc = Benchmark(
    "ts[:10000].resample('5min', how='ohlc')",
    common_setup,
    start_date=datetime(2012, 5, 1))

timeseries_1min_5min_mean = Benchmark(
    "ts[:10000].resample('5min', how='mean')",
    common_setup,
    start_date=datetime(2012, 5, 1))

#----------------------------------------------------------------------
# Irregular alignment

setup = common_setup + """
lindex = np.random.permutation(N)[:N // 2]
rindex = np.random.permutation(N)[:N // 2]
left = Series(ts.values.take(lindex), index=ts.index.take(lindex))
right = Series(ts.values.take(rindex), index=ts.index.take(rindex))
"""

timeseries_add_irregular = Benchmark('left + right', setup)

#----------------------------------------------------------------------
# Sort large irregular time series

setup = common_setup + """
N = 100000
rng = date_range(start='1/1/2000', periods=N, freq='s')
rng = rng.take(np.random.permutation(N))
ts = Series(np.random.randn(N), index=rng)
"""

timeseries_sort_index = Benchmark('ts.sort_index()', setup,
                                  start_date=datetime(2012, 4, 1))

#----------------------------------------------------------------------
# Shifting, add offset

setup = common_setup + """
rng = date_range(start='1/1/2000', periods=10000, freq='T')
"""

datetimeindex_add_offset = Benchmark('rng + timedelta(minutes=2)', setup,
                                     start_date=datetime(2012, 4, 1))

setup = common_setup + """
N = 10000
rng = date_range(start='1/1/1990', periods=N, freq='53s')
ts = Series(np.random.randn(N), index=rng)
dates = date_range(start='1/1/1990', periods=N * 10, freq='5s')
"""
timeseries_asof_single = Benchmark('ts.asof(dates[0])', setup,
                                   start_date=datetime(2012, 4, 27))

timeseries_asof = Benchmark('ts.asof(dates)', setup,
                            start_date=datetime(2012, 4, 27))

setup = setup + 'ts[250:5000] = np.nan'

timeseries_asof_nan = Benchmark('ts.asof(dates)', setup,
                                start_date=datetime(2012, 4, 27))

#----------------------------------------------------------------------
# Time zone

setup = common_setup + """
rng = date_range(start='1/1/2000', end='3/1/2000', tz='US/Eastern')
"""

timeseries_timestamp_tzinfo_cons = \
    Benchmark('rng[0]', setup, start_date=datetime(2012, 5, 5))

#----------------------------------------------------------------------
# Resampling period

setup = common_setup + """
rng = period_range(start='1/1/2000', end='1/1/2001', freq='T')
ts = Series(np.random.randn(len(rng)), index=rng)
"""

timeseries_period_downsample_mean = \
    Benchmark("ts.resample('D', how='mean')", setup,
              start_date=datetime(2012, 4, 25))

setup = common_setup + """
rng = date_range(start='1/1/2000', end='1/1/2001', freq='T')
ts = Series(np.random.randn(len(rng)), index=rng)
"""

timeseries_timestamp_downsample_mean = \
    Benchmark("ts.resample('D', how='mean')", setup,
              start_date=datetime(2012, 4, 25))

# GH 7754
setup = common_setup + """
rng = date_range(start='2000-01-01 00:00:00',
                    end='2000-01-01 10:00:00', freq='555000U')
int_ts = Series(5, rng, dtype='int64')
ts = int_ts.astype('datetime64[ns]')
"""

timeseries_resample_datetime64 = Benchmark("ts.resample('1S', how='last')", setup)

#----------------------------------------------------------------------
# to_datetime

setup = common_setup + """
rng = date_range(start='1/1/2000', periods=20000, freq='H')
strings = [x.strftime('%Y-%m-%d %H:%M:%S') for x in rng]
"""

timeseries_to_datetime_iso8601 = \
    Benchmark('to_datetime(strings)', setup,
              start_date=datetime(2012, 7, 11))

timeseries_to_datetime_iso8601_format = \
    Benchmark("to_datetime(strings, format='%Y-%m-%d %H:%M:%S')", setup,
              start_date=datetime(2012, 7, 11))

setup = common_setup + """
rng = date_range(start='1/1/2000', periods=10000, freq='D')
strings = Series(rng.year*10000+rng.month*100+rng.day,dtype=np.int64).apply(str)
"""

timeseries_to_datetime_YYYYMMDD = \
    Benchmark('to_datetime(strings,format="%Y%m%d")', setup,
              start_date=datetime(2012, 7, 1))

setup = common_setup + """
s = Series(['19MAY11','19MAY11:00:00:00']*100000)
"""
timeseries_with_format_no_exact = Benchmark("to_datetime(s,format='%d%b%y',exact=False)", \
     setup, start_date=datetime(2014, 11, 26))
timeseries_with_format_replace = Benchmark("to_datetime(s.str.replace(':\S+$',''),format='%d%b%y')", \
     setup, start_date=datetime(2014, 11, 26))

# ---- infer_freq
# infer_freq

setup = common_setup + """
from pandas.tseries.frequencies import infer_freq
rng = date_range(start='1/1/1700', freq='D', periods=100000)
a = rng[:50000].append(rng[50002:])
"""

timeseries_infer_freq = \
    Benchmark('infer_freq(a)', setup, start_date=datetime(2012, 7, 1))

# setitem PeriodIndex

setup = common_setup + """
rng = period_range(start='1/1/1990', freq='S', periods=20000)
df = DataFrame(index=range(len(rng)))
"""

period_setitem = \
    Benchmark("df['col'] = rng", setup,
              start_date=datetime(2012, 8, 1))

setup = common_setup + """
rng = date_range(start='1/1/2000 9:30', periods=10000, freq='S', tz='US/Eastern')
"""

datetimeindex_normalize = \
    Benchmark('rng.normalize()', setup,
              start_date=datetime(2012, 9, 1))

setup = common_setup + """
from pandas.tseries.offsets import Second
s1 = date_range(start='1/1/2000', periods=100, freq='S')
curr = s1[-1]
slst = []
for i in range(100):
    slst.append(curr + Second()), periods=100, freq='S')
    curr = slst[-1][-1]
"""

# dti_append_tz = \
#     Benchmark('s1.append(slst)', setup, start_date=datetime(2012, 9, 1))


setup = common_setup + """
rng = date_range(start='1/1/2000', periods=1000, freq='H')
df = DataFrame(np.random.randn(len(rng), 2), rng)
"""

dti_reset_index = \
    Benchmark('df.reset_index()', setup, start_date=datetime(2012, 9, 1))

setup = common_setup + """
rng = date_range(start='1/1/2000', periods=1000, freq='H',
                 tz='US/Eastern')
df = DataFrame(np.random.randn(len(rng), 2), index=rng)
"""

dti_reset_index_tz = \
    Benchmark('df.reset_index()', setup, start_date=datetime(2012, 9, 1))

setup = common_setup + """
rng = date_range(start='1/1/2000', periods=1000, freq='T')
index = rng.repeat(10)
"""

datetimeindex_unique = Benchmark('index.unique()', setup,
                                 start_date=datetime(2012, 7, 1))

# tz_localize with infer argument.  This is an attempt to emulate the results
# of read_csv with duplicated data.  Not passing infer_dst will fail
setup = common_setup + """
dst_rng = date_range(start='10/29/2000 1:00:00',
                     end='10/29/2000 1:59:59', freq='S')
index = date_range(start='10/29/2000', end='10/29/2000 00:59:59', freq='S')
index = index.append(dst_rng)
index = index.append(dst_rng)
index = index.append(date_range(start='10/29/2000 2:00:00',
                                end='10/29/2000 3:00:00', freq='S'))
"""

datetimeindex_infer_dst = \
Benchmark('index.tz_localize("US/Eastern", infer_dst=True)',
          setup, start_date=datetime(2013, 9, 30))


#----------------------------------------------------------------------
# Resampling: fast-path various functions

setup = common_setup + """
rng = date_range(start='20130101',periods=100000,freq='50L')
df = DataFrame(np.random.randn(100000,2),index=rng)
"""

dataframe_resample_mean_string = \
    Benchmark("df.resample('1s', how='mean')", setup)

dataframe_resample_mean_numpy = \
    Benchmark("df.resample('1s', how=np.mean)", setup)

dataframe_resample_min_string = \
    Benchmark("df.resample('1s', how='min')", setup)

dataframe_resample_min_numpy = \
    Benchmark("df.resample('1s', how=np.min)", setup)

dataframe_resample_max_string = \
    Benchmark("df.resample('1s', how='max')", setup)

dataframe_resample_max_numpy = \
    Benchmark("df.resample('1s', how=np.max)", setup)


#----------------------------------------------------------------------
# DatetimeConverter

setup = common_setup + """
from pandas.tseries.converter import DatetimeConverter
"""

datetimeindex_converter = \
    Benchmark('DatetimeConverter.convert(rng, None, None)',
              setup, start_date=datetime(2013, 1, 1))

# Adding custom business day
setup = common_setup + """
import datetime as dt
import pandas as pd
try:
    import pandas.tseries.holiday
except ImportError:
    pass
import numpy as np

date = dt.datetime(2011,1,1)
dt64 = np.datetime64('2011-01-01 09:00Z')
hcal = pd.tseries.holiday.USFederalHolidayCalendar()

day = pd.offsets.Day()
year = pd.offsets.YearBegin()
cday = pd.offsets.CustomBusinessDay()
cmb = pd.offsets.CustomBusinessMonthBegin(calendar=hcal)
cme = pd.offsets.CustomBusinessMonthEnd(calendar=hcal)

cdayh = pd.offsets.CustomBusinessDay(calendar=hcal)
"""
timeseries_day_incr = Benchmark("date + day",setup)

timeseries_day_apply = Benchmark("day.apply(date)",setup)

timeseries_year_incr = Benchmark("date + year",setup)

timeseries_year_apply = Benchmark("year.apply(date)",setup)

timeseries_custom_bday_incr = \
    Benchmark("date + cday",setup)

timeseries_custom_bday_decr = \
    Benchmark("date - cday",setup)

timeseries_custom_bday_apply = \
    Benchmark("cday.apply(date)",setup)

timeseries_custom_bday_apply_dt64 = \
    Benchmark("cday.apply(dt64)",setup)

timeseries_custom_bday_cal_incr = \
    Benchmark("date + 1 * cdayh",setup)

timeseries_custom_bday_cal_decr = \
    Benchmark("date - 1 * cdayh",setup)

timeseries_custom_bday_cal_incr_n = \
    Benchmark("date + 10 * cdayh",setup)

timeseries_custom_bday_cal_incr_neg_n = \
    Benchmark("date - 10 * cdayh",setup)

# Increment custom business month
timeseries_custom_bmonthend_incr = \
    Benchmark("date + cme",setup)

timeseries_custom_bmonthend_incr_n = \
    Benchmark("date + 10 * cme",setup)

timeseries_custom_bmonthend_decr_n = \
    Benchmark("date - 10 * cme",setup)

timeseries_custom_bmonthbegin_incr_n = \
    Benchmark("date + 10 * cmb",setup)

timeseries_custom_bmonthbegin_decr_n = \
    Benchmark("date - 10 * cmb",setup)


#----------------------------------------------------------------------
# month/quarter/year start/end accessors

setup = common_setup + """
N = 10000
rng = date_range(start='1/1/1', periods=N, freq='B')
"""

timeseries_is_month_start = Benchmark('rng.is_month_start', setup,
                                  start_date=datetime(2014, 4, 1))

#----------------------------------------------------------------------
# iterate over DatetimeIndex/PeriodIndex
setup = common_setup + """
N = 1000000
M = 10000
idx1 = date_range(start='20140101', freq='T', periods=N)
idx2 = period_range(start='20140101', freq='T', periods=N)

def iter_n(iterable, n=None):
  i = 0
  for _ in iterable:
    i += 1
    if n is not None and i > n:
      break
"""

timeseries_iter_datetimeindex = Benchmark('iter_n(idx1)', setup)

timeseries_iter_periodindex = Benchmark('iter_n(idx2)', setup)

timeseries_iter_datetimeindex_preexit = Benchmark('iter_n(idx1, M)', setup)

timeseries_iter_periodindex_preexit = Benchmark('iter_n(idx2, M)', setup)


#----------------------------------------------------------------------
# apply an Offset to a  DatetimeIndex
setup = common_setup + """
N = 100000
idx1 = date_range(start='20140101', freq='T', periods=N)
delta_offset = pd.offsets.Day()
fast_offset = pd.offsets.DateOffset(months=2, days=2)
slow_offset = pd.offsets.BusinessDay()

"""

timeseries_datetimeindex_offset_delta = Benchmark('idx1 + delta_offset', setup)
timeseries_datetimeindex_offset_fast = Benchmark('idx1 + fast_offset', setup)
timeseries_datetimeindex_offset_slow = Benchmark('idx1 + slow_offset', setup)

# apply an Offset to a Series containing datetime64 values
setup = common_setup + """
N = 100000
s = Series(date_range(start='20140101', freq='T', periods=N))
delta_offset = pd.offsets.Day()
fast_offset = pd.offsets.DateOffset(months=2, days=2)
slow_offset = pd.offsets.BusinessDay()

"""

timeseries_series_offset_delta = Benchmark('s + delta_offset', setup)
timeseries_series_offset_fast = Benchmark('s + fast_offset', setup)
timeseries_series_offset_slow = Benchmark('s + slow_offset', setup)
