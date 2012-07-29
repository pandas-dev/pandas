from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
from datetime import timedelta
N = 100000

try:
    rng = date_range('1/1/2000', periods=N, freq='min')
except NameError:
    rng = DateRange('1/1/2000', periods=N, offset=datetools.Minute())
    def date_range(start=None, end=None, periods=None, freq=None):
        return DateRange(start, end, periods=periods, offset=freq)

if hasattr(Series, 'convert'):
    Series.resample = Series.convert

ts = Series(np.random.randn(N), index=rng)
"""

#----------------------------------------------------------------------
# Lookup value in large time series, hash map population

setup = common_setup + """
rng = date_range('1/1/2000', periods=1500000, freq='s')
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

timeseries_1min_5min_ohlc = Benchmark("ts[:10000].resample('5min', how='ohlc')",
                                      common_setup,
                                      start_date=datetime(2012, 5, 1))

timeseries_1min_5min_mean = Benchmark("ts[:10000].resample('5min', how='mean')",
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
rng = date_range('1/1/2000', periods=N, freq='s')
rng = rng.take(np.random.permutation(N))
ts = Series(np.random.randn(N), index=rng)
"""

timeseries_sort_index = Benchmark('ts.sort_index()', setup,
                                  start_date=datetime(2012, 4, 1))

#----------------------------------------------------------------------
# Shifting, add offset

setup = common_setup + """
rng = date_range('1/1/2000', periods=10000, freq='T')
"""

datetimeindex_add_offset = Benchmark('rng + timedelta(minutes=2)', setup,
                                     start_date=datetime(2012, 4, 1))

setup = common_setup + """
N = 10000
rng = date_range('1/1/1990', periods=N, freq='53s')
ts = Series(np.random.randn(N), index=rng)
dates = date_range('1/1/1990', periods=N * 10, freq='5s')
"""
timeseries_asof_single = Benchmark('ts.asof(dates[0])', setup,
                                 start_date=datetime(2012, 4, 27))

timeseries_asof = Benchmark('ts.asof(dates)', setup,
                            start_date=datetime(2012, 4, 27))

setup = setup + 'ts[250:5000] = np.nan'

timeseries_asof_nan = Benchmark('ts.asof(dates)', setup,
                                start_date=datetime(2012, 4, 27))

#----------------------------------------------------------------------
# Time zone stuff

setup = common_setup + """
rng = date_range('1/1/2000', '3/1/2000', tz='US/Eastern')
"""

timeseries_timestamp_tzinfo_cons = \
    Benchmark('rng[0]', setup, start_date=datetime(2012, 5, 5))

#----------------------------------------------------------------------
# Resampling period

setup = common_setup + """
rng = period_range('1/1/2000', '1/1/2001', freq='T')
ts = Series(np.random.randn(len(rng)), index=rng)
"""

timeseries_period_downsample_mean = \
    Benchmark("ts.resample('D', how='mean')", setup,
              start_date=datetime(2012, 4, 25))

setup = common_setup + """
rng = date_range('1/1/2000', '1/1/2001', freq='T')
ts = Series(np.random.randn(len(rng)), index=rng)
"""

timeseries_timestamp_downsample_mean = \
    Benchmark("ts.resample('D', how='mean')", setup,
              start_date=datetime(2012, 4, 25))

#----------------------------------------------------------------------
# to_datetime

setup = common_setup + """
rng = date_range('1/1/2000', periods=20000, freq='h')
strings = [x.strftime('%Y-%m-%d %H:%M:%S') for x in rng]
"""

timeseries_to_datetime_iso8601 = \
    Benchmark('to_datetime(strings)', setup,
              start_date=datetime(2012, 7, 11))

# ---- infer_freq
# infer_freq

setup = common_setup + """
from pandas.tseries.frequencies import infer_freq
rng = date_range('1/1/1700', freq='D', periods=100000)
a = rng[:50000].append(rng[50002:])
"""

timeseries_infer_freq = \
    Benchmark('infer_freq(a)', setup, start_date=datetime(2012, 7, 1))
