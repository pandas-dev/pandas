from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

setup = common_setup + """

N = 100000
ngroups = 5

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

stmt1 = "df.groupby(['key1', 'key2'])['data'].agg(lambda x: x.values.sum())"
bm_groupby1 = Benchmark(stmt1, setup,
                        name="groupby_multi_python",
                        start_date=datetime(2011, 7, 1))

stmt3 = "df.groupby(['key1', 'key2']).sum()"
bm_groupby3 = Benchmark(stmt3, setup,
                        name="groupby_multi_cython",
                        start_date=datetime(2011, 7, 1))

