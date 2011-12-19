from vbench.benchmark import Benchmark
from vbench.db import BenchmarkDB
from vbench.runner import BenchmarkRunner
from vbench.git import GitRepo

from datetime import datetime


common_setup = """
from pandas import *
import pandas.util.testing as tm
import random
import numpy as np
"""

#----------------------------------------------------------------------
# Series.__getitem__, get_value

setup = common_setup + """
tm.N = 1000
ts = tm.makeTimeSeries()
dt = ts.index[500]
"""
statement = "ts[dt]"

bm_getitem = Benchmark(statement, setup, ncalls=100000,
                       name='Series.__getitem__ scalar')

#----------------------------------------------------------------------
# DataFrame reindex columns

setup = common_setup + """
index = [tm.rands(10) for _ in xrange(1000)]
columns = [tm.rands(10) for _ in xrange(30)]
df = DataFrame(np.random.rand(1000, 30), index=index,
               columns=columns)
idx = index[100]
col = columns[10]
"""
statement = "df[col][idx]"
bm_df_getitem = Benchmark(statement, setup,
                        name='DataFrame get scalar value')

setup = common_setup + """
try:
    klass = DataMatrix
except:
    klass = DataFrame

index = [tm.rands(10) for _ in xrange(1000)]
columns = [tm.rands(10) for _ in xrange(30)]
df = klass(np.random.rand(1000, 30), index=index,
               columns=columns)
idx = index[100]
col = columns[10]
"""
statement = "df[col][idx]"
bm_df_getitem2 = Benchmark(statement, setup,
                        name='DataMatrix get scalar value')

setup = common_setup + """
try:
    klass = DataMatrix
except:
    klass = DataFrame

index = [tm.rands(10) for _ in xrange(1000)]
columns = [tm.rands(10) for _ in xrange(30)]
df = klass(np.random.rand(1000, 30), index=index,
               columns=columns)
idx = index[100]
col = columns[10]
"""
statement = "df.get_value(idx, col)"
bm_df_getitem3 = Benchmark(statement, setup, name='DataFrame get_value',
                           start_date=datetime(2011, 11, 12))

#----------------------------------------------------------------------
# DataFrame reindex columns

setup = common_setup + """
df = DataFrame(index=range(10000), data=np.random.rand(10000,30),
               columns=range(30))
"""
statement = "df.reindex(columns=df.columns[1:5])"

bm_reindex1 = Benchmark(statement, setup, name='DataFrame.reindex columns')

#----------------------------------------------------------------------

setup = common_setup + """
rng = DateRange('1/1/1970', periods=10000, offset=datetools.Minute())
df = DataFrame(np.random.rand(10000, 10), index=rng,
               columns=range(10))
df['foo'] = 'bar'
rng2 = Index(rng[::2])
"""
statement = "df.reindex(rng2)"
bm_reindex2 = Benchmark(statement, setup,
                        name='DataFrame.reindex index daterange')

#----------------------------------------------------------------------
setup = common_setup + """

N = 100000
ngroups = 100

def get_test_data(ngroups=100, n=N):
    unique_groups = range(ngroups)
    arr = np.asarray(np.tile(unique_groups, n / ngroups), dtype=object)

    if len(arr) < n:
        arr = np.asarray(list(arr) + unique_groups[:n - len(arr)],
                         dtype=object)

    random.shuffle(arr)
    return arr

df = DataFrame({'key1' : get_test_data(ngroups=ngroups),
                'key2' : get_test_data(ngroups=ngroups),
                'data' : np.random.randn(N)})
def f():
    df.groupby(['key1', 'key2']).agg(lambda x: x.values.sum())
"""
stmt = "df.groupby(['key1', 'key2'])['data'].agg(lambda x: x.values.sum())"

bm_groupby1 = Benchmark(stmt, setup, name="GroupBy test foo",
                        start_date=datetime(2011, 7, 1))

stmt2 = "df.groupby(['key1', 'key2']).sum()"
bm_groupby2 = Benchmark(stmt2, setup, name="GroupBy test 2",
                        start_date=datetime(2011, 7, 1))

#----------------------------------------------------------------------
# data alignment

setup = common_setup + """
from pandas import *
from pandas.util.testing import rands

n = 1000000
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
bm_align1 = Benchmark(stmt, setup, name="int64 alignment 1",
                      start_date=datetime(2011, 3, 1))

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

bm_sparse1 = Benchmark(stmt, setup, name="SparseSeries to SparseDataFrame",
                      start_date=datetime(2011, 6, 1))

#----------------------------------------------------------------------
# Actually running the benchmarks

benchmarks = [v for v in locals().values() if isinstance(v, Benchmark)]

REPO_PATH = '/home/wesm/code/pandas'
REPO_URL = 'git@github.com:wesm/pandas.git'
DB_PATH = '/home/wesm/code/pandas/gb_suite/benchmarks.db'
TMP_DIR = '/home/wesm/tmp/gb_pandas'
PREPARE = """
python setup.py clean
"""
BUILD = """
python setup.py build_ext --inplace
"""
START_DATE = datetime(2011, 3, 1)

repo = GitRepo(REPO_PATH)

to_consider = repo.shas.truncate(START_DATE)

runner = BenchmarkRunner(benchmarks, REPO_PATH, REPO_URL,
                         BUILD, DB_PATH, TMP_DIR, PREPARE,
                         run_option='eod', start_date=START_DATE)

runner.run()
