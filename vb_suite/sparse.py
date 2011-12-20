from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

#----------------------------------------------------------------------

setup = common_setup + """
from pandas.core.sparse import SparseSeries, SparseDataFrame

K = 50
N = 50000
rng = np.asarray(DateRange('1/1/2000', periods=N,
                           offset=datetools.Minute()))

# rng2 = np.asarray(rng).astype('M8').astype('i8')

series = {}
for i in range(1, K + 1):
    data = np.random.randn(N)[:-i]
    this_rng = rng[:-i]
    data[100:] = np.nan
    series[i] = SparseSeries(data, index=this_rng)
"""
stmt = "SparseDataFrame(series)"

bm_sparse1 = Benchmark(stmt, setup, name="sparse_series_to_frame",
                      start_date=datetime(2011, 6, 1))
