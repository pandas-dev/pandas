from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

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

# aggregate multiple columns
df = DataFrame({'key1' : get_test_data(ngroups=ngroups),
                'key2' : get_test_data(ngroups=ngroups),
                'data1' : np.random.randn(N),
                'data2' : np.random.randn(N)})
def f():
    df.groupby(['key1', 'key2']).agg(lambda x: x.values.sum())

simple_series = Series(np.random.randn(N))
key1 = df['key1']
"""

stmt1 = "df.groupby(['key1', 'key2'])['data1'].agg(lambda x: x.values.sum())"
groupby_multi_python = Benchmark(stmt1, setup,
                                 start_date=datetime(2011, 7, 1))

stmt3 = "df.groupby(['key1', 'key2']).sum()"
groupby_multi_cython = Benchmark(stmt3, setup,
                                 start_date=datetime(2011, 7, 1))

stmt = "df.groupby(['key1', 'key2'])['data1'].agg(np.std)"
groupby_multi_series_op = Benchmark(stmt, setup,
                                    start_date=datetime(2011, 8, 1))

groupby_series_simple_cython = \
    Benchmark('simple_series.groupby(key1).sum()', setup,
              start_date=datetime(2011, 3, 1))

#----------------------------------------------------------------------
# 2d grouping, aggregate many columns

setup = common_setup + """
labels = np.random.randint(0, 100, size=1000)
df = DataFrame(randn(1000, 1000))
"""

groupby_frame_cython_many_columns = Benchmark('df.groupby(labels).sum()', setup,
                                              start_date=datetime(2011, 8, 1),
                                              logy=True)

#----------------------------------------------------------------------
# single key, long, integer key

setup = common_setup + """
data = np.random.randn(100000, 1)
labels = np.random.randint(0, 1000, size=100000)
df = DataFrame(data)
"""

groupby_frame_singlekey_integer = \
    Benchmark('df.groupby(labels).sum()', setup,
              start_date=datetime(2011, 8, 1), logy=True)

