from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

#----------------------------------------------------------------------
# nanops

setup = common_setup + """
s = Series(np.random.randn(100000))
s[::2] = np.nan
"""

stat_ops_series_std = Benchmark("s.std()", setup)

#----------------------------------------------------------------------
# ops by level

setup = common_setup + """
index = MultiIndex(levels=[np.arange(10), np.arange(100), np.arange(100)],
                   labels=[np.arange(10).repeat(10000),
                           np.tile(np.arange(100).repeat(100), 10),
                           np.tile(np.tile(np.arange(100), 100), 10)])
random.shuffle(index.values)
df = DataFrame(np.random.randn(len(index), 4), index=index)
df_level = DataFrame(np.random.randn(100, 4), index=index.levels[1])
"""

stat_ops_level_frame_sum = \
    Benchmark("df.sum(level=1)", setup,
              start_date=datetime(2011, 11, 15))

stat_ops_level_frame_sum_multiple = \
    Benchmark("df.sum(level=[0, 1])", setup, repeat=1,
              start_date=datetime(2011, 11, 15))

stat_ops_level_series_sum = \
    Benchmark("df[1].sum(level=1)", setup,
              start_date=datetime(2011, 11, 15))

stat_ops_level_series_sum_multiple = \
    Benchmark("df[1].sum(level=[0, 1])", setup, repeat=1,
              start_date=datetime(2011, 11, 15))

