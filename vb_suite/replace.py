from vbench.api import Benchmark
from datetime import datetime

common_setup = """from .pandas_vb_common import *
from datetime import timedelta

N = 1000000

try:
    rng = date_range('1/1/2000', periods=N, freq='min')
except NameError:
    rng = DatetimeIndex('1/1/2000', periods=N, offset=datetools.Minute())
    date_range = DateRange

ts = Series(np.random.randn(N), index=rng)
"""

large_dict_setup = """from .pandas_vb_common import *
from pandas.compat import range
n = 10 ** 6
start_value = 10 ** 5
to_rep = dict((i, start_value + i) for i in range(n))
s = Series(np.random.randint(n, size=10 ** 3))
"""

replace_fillna = Benchmark('ts.fillna(0., inplace=True)', common_setup,
                           name='replace_fillna',
                           start_date=datetime(2012, 4, 4))
replace_replacena = Benchmark('ts.replace(np.nan, 0., inplace=True)',
                              common_setup,
                              name='replace_replacena',
                              start_date=datetime(2012, 5, 15))
replace_large_dict = Benchmark('s.replace(to_rep, inplace=True)',
                               large_dict_setup,
                               name='replace_large_dict',
                               start_date=datetime(2014, 4, 6))
