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

bm_reindex1 = Benchmark(statement, setup,
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
bm_reindex2 = Benchmark(statement, setup,
                        name='dataframe_reindex_daterange')
