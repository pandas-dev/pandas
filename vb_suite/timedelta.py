from vbench.api import Benchmark
from datetime import datetime

common_setup = """from .pandas_vb_common import *
from pandas import to_timedelta
"""

#----------------------------------------------------------------------
# conversion

setup = common_setup + """
arr = np.random.randint(0,1000,size=10000)
"""

stmt = "to_timedelta(arr,unit='s')"
timedelta_convert_int = Benchmark(stmt, setup, start_date=datetime(2014, 1, 1))

setup = common_setup + """
arr = np.random.randint(0,1000,size=10000)
arr = [ '{0} days'.format(i) for i in arr ]
"""

stmt = "to_timedelta(arr)"
timedelta_convert_string = Benchmark(stmt, setup, start_date=datetime(2014, 1, 1))

setup = common_setup + """
arr = np.random.randint(0,60,size=10000)
arr = [ '00:00:{0:02d}'.format(i) for i in arr ]
"""

stmt = "to_timedelta(arr)"
timedelta_convert_string_seconds = Benchmark(stmt, setup, start_date=datetime(2014, 1, 1))
