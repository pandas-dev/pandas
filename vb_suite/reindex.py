from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

#----------------------------------------------------------------------
# DataFrame reindex columns

setup = common_setup + """
df = DataFrame(index=range(10000), data=np.random.rand(10000,30),
               columns=range(30))
"""
statement = "df.reindex(columns=df.columns[1:5])"

reindex_frame_columns = Benchmark(statement, setup,
                                  name='dataframe_reindex_columns')

#----------------------------------------------------------------------

setup = common_setup + """
rng = DateRange('1/1/1970', periods=10000, offset=datetools.Minute())
df = DataFrame(np.random.rand(10000, 10), index=rng,
               columns=range(10))
df['foo'] = 'bar'
rng2 = Index(rng[::2])
"""
statement = "df.reindex(rng2)"
reindex_frame_daterange = Benchmark(statement, setup,
                                    name='dataframe_reindex_daterange')

#----------------------------------------------------------------------
# multiindex reindexing

setup = common_setup + """
N = 1000
K = 20

level1 = np.array([tm.rands(10) for _ in xrange(N)], dtype='O').repeat(K)
level2 = np.tile(np.array([tm.rands(10) for _ in xrange(K)], dtype='O'),
                 N)
index = MultiIndex.from_arrays([level1, level2])

s1 = Series(np.random.randn(N * K), index=index)
s2 = s1[::2]
"""
statement = "s1.reindex(s2.index)"
reindex_multi = Benchmark(statement, setup,
                          name='reindex_multiindex',
                          start_date=datetime(2011, 8, 1))
