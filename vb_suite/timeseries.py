from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
from datetime import timedelta
N = 100000

try:
    rng = date_range('1/1/2000', periods=N, freq='min')
except NameError:
    rng = DateRange('1/1/2000', periods=N, offset=datetools.Minute())
    date_range = DateRange

ts = Series(np.random.randn(N), index=rng)
"""

#----------------------------------------------------------------------
# Test slice minutely series

timeseries_slice_minutely = Benchmark('ts[:10000]', common_setup)

#----------------------------------------------------------------------
# Test conversion

setup = common_setup + """

"""

timeseries_1min_5min_ohlc = Benchmark("ts[:10000].convert('5min', how='ohlc')",
                                      common_setup)

timeseries_1min_5min_mean = Benchmark("ts[:10000].convert('5min', how='mean')",
                                      common_setup)

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
N = 1000
rng = date_range('1/1/1990', periods=N, freq='53s')
ts = Series(np.random.randn(N), index=rng)
dates = date_range('1/1/1990', periods=N * 10, freq='5s')
"""
timeseries_asof_loop = Benchmark('[ts.asof(d) for d in dates]', setup,
                                 start_date=datetime(2012, 4, 27))

timeseries_asof = Benchmark('ts.asof(dates)', setup,
                            start_date=datetime(2012, 4, 27))
