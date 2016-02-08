from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from .pandas_vb_common import *
"""

#----------------------------------------------------------------------
# nanops

setup = common_setup + """
s = Series(np.random.randn(100000), index=np.arange(100000))
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

sum_setup = common_setup + """
df = DataFrame(np.random.randn(100000, 4))
dfi = DataFrame(np.random.randint(1000, size=df.shape))
"""

stat_ops_frame_sum_int_axis_0 = \
    Benchmark("dfi.sum()", sum_setup, start_date=datetime(2013, 7, 25))

stat_ops_frame_sum_float_axis_0 = \
    Benchmark("df.sum()", sum_setup, start_date=datetime(2013, 7, 25))

stat_ops_frame_mean_int_axis_0 = \
    Benchmark("dfi.mean()", sum_setup, start_date=datetime(2013, 7, 25))

stat_ops_frame_mean_float_axis_0 = \
    Benchmark("df.mean()", sum_setup, start_date=datetime(2013, 7, 25))

stat_ops_frame_sum_int_axis_1 = \
    Benchmark("dfi.sum(1)", sum_setup, start_date=datetime(2013, 7, 25))

stat_ops_frame_sum_float_axis_1 = \
    Benchmark("df.sum(1)", sum_setup, start_date=datetime(2013, 7, 25))

stat_ops_frame_mean_int_axis_1 = \
    Benchmark("dfi.mean(1)", sum_setup, start_date=datetime(2013, 7, 25))

stat_ops_frame_mean_float_axis_1 = \
    Benchmark("df.mean(1)", sum_setup, start_date=datetime(2013, 7, 25))

#----------------------------------------------------------------------
# rank

setup = common_setup + """
values = np.concatenate([np.arange(100000),
                         np.random.randn(100000),
                         np.arange(100000)])
s = Series(values)
"""

stats_rank_average = Benchmark('s.rank()', setup,
                               start_date=datetime(2011, 12, 12))

stats_rank_pct_average = Benchmark('s.rank(pct=True)', setup,
                                   start_date=datetime(2014, 1, 16))
stats_rank_pct_average_old = Benchmark('s.rank() / len(s)', setup,
                                       start_date=datetime(2014, 1, 16))
setup = common_setup + """
values = np.random.randint(0, 100000, size=200000)
s = Series(values)
"""

stats_rank_average_int = Benchmark('s.rank()', setup,
                                   start_date=datetime(2011, 12, 12))

setup = common_setup + """
df = DataFrame(np.random.randn(5000, 50))
"""

stats_rank2d_axis1_average = Benchmark('df.rank(1)', setup,
                                       start_date=datetime(2011, 12, 12))

stats_rank2d_axis0_average = Benchmark('df.rank()', setup,
                                       start_date=datetime(2011, 12, 12))

# rolling functions

setup = common_setup + """
arr = np.random.randn(100000)
"""

stats_rolling_mean = Benchmark('rolling_mean(arr, 100)', setup,
                               start_date=datetime(2011, 6, 1))

# spearman correlation

setup = common_setup + """
df = DataFrame(np.random.randn(1000, 30))
"""

stats_corr_spearman = Benchmark("df.corr(method='spearman')", setup,
                                start_date=datetime(2011, 12, 4))
