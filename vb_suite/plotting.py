from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from .pandas_vb_common import *

try:
    from pandas import date_range
except ImportError:
    def date_range(start=None, end=None, periods=None, freq=None):
        return DatetimeIndex(start, end, periods=periods, offset=freq)

"""

#-----------------------------------------------------------------------------
# Timeseries plotting

setup = common_setup + """
N = 2000
M = 5
df = DataFrame(np.random.randn(N,M), index=date_range('1/1/1975', periods=N))
"""

plot_timeseries_period = Benchmark("df.plot()", setup=setup, 
                                   name='plot_timeseries_period')

