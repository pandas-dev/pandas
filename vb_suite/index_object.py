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

index_datetime_intersection = Benchmark("rng.intersection(rng2)", setup)
index_datetime_union = Benchmark("rng.union(rng2)", setup)

# integers
setup = common_setup + """
N = 1000000
options = np.arange(N)

left = Index(options.take(np.random.permutation(N)[:N // 2]))
right = Index(options.take(np.random.permutation(N)[:N // 2]))
"""

index_int64_union = Benchmark('left.union(right)', setup,
                              start_date=datetime(2011, 1, 1))

index_int64_intersection = Benchmark('left.intersection(right)', setup,
                                     start_date=datetime(2011, 1, 1))
