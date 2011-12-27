from vbench.benchmark import Benchmark
from datetime import datetime

SECTION = "Index / MultiIndex objects"


common_setup = """from pandas_vb_common import *
"""

#----------------------------------------------------------------------
# intersection, union

setup = common_setup + """
rng = DateRange('1/1/2000', periods=10000, offset=datetools.Minute())
rng = rng.view(Index)
rng2 = rng[:-1]
"""

index_datetime_intersection = Benchmark("rng.intersection(rng2)", setup,
                                        name='index_datetime_intersection')
index_datetime_union = Benchmark("rng.union(rng2)", setup,
                                 name='index_datetime_union')
