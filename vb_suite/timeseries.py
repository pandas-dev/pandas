from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
try:
    rng = date_range('1/1/2000', periods=100000, freq='min')
except NameError:
    rng = DateRange('1/1/2000', periods=100000,
                    offset=datetools.Minute())

ts = Series(np.random.randn(100000), index=rng)
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
