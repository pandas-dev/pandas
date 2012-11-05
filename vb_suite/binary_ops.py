from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

#----------------------------------------------------------------------
# data alignment

setup = common_setup + """n = 1000000
# indices = Index([rands(10) for _ in xrange(n)])
def sample(values, k):
    sampler = np.random.permutation(len(values))
    return values.take(sampler[:k])
sz = 500000
rng = np.arange(0, 10000000000000, 10000000)
stamps = np.datetime64(datetime.now()).view('i8') + rng
idx1 = np.sort(sample(stamps, sz))
idx2 = np.sort(sample(stamps, sz))
ts1 = Series(np.random.randn(sz), idx1)
ts2 = Series(np.random.randn(sz), idx2)
"""
stmt = "ts1 + ts2"
series_align_int64_index_binop = Benchmark(stmt, setup,
                                           start_date=datetime(2010, 6, 1),
                                           logy=True)
