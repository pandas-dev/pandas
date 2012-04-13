from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
N = 100000

try:
    rng = date_range('1/1/2000', periods=N, freq='min')
except NameError:
    rng = DateRange('1/1/2000', periods=N, offset=datetools.Minute())

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
