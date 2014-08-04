from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

setup = common_setup + """
s1 = Series(np.random.randn(10000))
s2 = Series(np.random.randint(1, 10, 10000))
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

series_nunique1 = Benchmark('s1.nunique()',
                            setup,
                            start_date=datetime(2014, 1, 25))

series_nunique2 = Benchmark('s2.nunique()',
                            setup,
                            start_date=datetime(2014, 1, 25))
