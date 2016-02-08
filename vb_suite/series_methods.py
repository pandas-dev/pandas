from vbench.api import Benchmark
from datetime import datetime

common_setup = """from .pandas_vb_common import *
"""

setup = common_setup + """
s1 = Series(np.random.randn(10000))
s2 = Series(np.random.randint(1, 10, 10000))
s3 = Series(np.random.randint(1, 10, 100000)).astype('int64')
values = [1,2]
s4 = s3.astype('object')
"""

series_nlargest1 = Benchmark('s1.nlargest(3, take_last=True);'
                             's1.nlargest(3, take_last=False)',
                             setup,
                             start_date=datetime(2014, 1, 25))
series_nlargest2 = Benchmark('s2.nlargest(3, take_last=True);'
                             's2.nlargest(3, take_last=False)',
                             setup,
                             start_date=datetime(2014, 1, 25))

series_nsmallest2 = Benchmark('s1.nsmallest(3, take_last=True);'
                              's1.nsmallest(3, take_last=False)',
                              setup,
                              start_date=datetime(2014, 1, 25))

series_nsmallest2 = Benchmark('s2.nsmallest(3, take_last=True);'
                              's2.nsmallest(3, take_last=False)',
                              setup,
                              start_date=datetime(2014, 1, 25))

series_isin_int64 = Benchmark('s3.isin(values)',
                              setup,
                              start_date=datetime(2014, 1, 25))
series_isin_object = Benchmark('s4.isin(values)',
                               setup,
                               start_date=datetime(2014, 1, 25))
