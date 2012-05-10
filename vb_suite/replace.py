from vbench.api import Benchmark

common_setup = """from pandas_vb_common import *
from datetime import timedelta
import pandas._tseries as lib
N = 1000000

try:
    rng = date_range('1/1/2000', periods=N, freq='min')
except NameError:
    rng = DateRange('1/1/2000', periods=N, offset=datetools.Minute())
    date_range = DateRange

ts = Series(np.random.randn(N), index=rng)

def replace_slow(ser, old, new):
    lib.slow_replace(ser.values, old, new)
    return ser
"""

replace_fillna = Benchmark('ts.fillna(0., inplace=True)', common_setup)
replace_replacena = Benchmark('ts.replace(np.nan, 0., inplace=True)',
                              common_setup)
replace_putmask = Benchmark('replace_slow(ts, np.nan, 0.)', common_setup)
