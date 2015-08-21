from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from .pandas_vb_common import *
import pandas as pd
df  = DataFrame(np.random.randn(20000, 100))
df2 = DataFrame(np.random.randn(20000, 100))
df3 = DataFrame(np.random.randn(20000, 100))
df4 = DataFrame(np.random.randn(20000, 100))
"""

setup = common_setup + """
import pandas.computation.expressions as expr
expr.set_numexpr_threads(1)
"""

SECTION = 'Eval'

#----------------------------------------------------------------------
# binary ops

#----------------------------------------------------------------------
# add
eval_frame_add_all_threads = \
    Benchmark("pd.eval('df + df2 + df3 + df4')", common_setup,
              name='eval_frame_add_all_threads',
              start_date=datetime(2013, 7, 21))



eval_frame_add_one_thread = \
    Benchmark("pd.eval('df + df2 + df3 + df4')", setup,
              name='eval_frame_add_one_thread',
              start_date=datetime(2013, 7, 26))

eval_frame_add_python = \
    Benchmark("pd.eval('df + df2 + df3 + df4', engine='python')", common_setup,
              name='eval_frame_add_python', start_date=datetime(2013, 7, 21))

eval_frame_add_python_one_thread = \
    Benchmark("pd.eval('df + df2 + df3 + df4', engine='python')", setup,
              name='eval_frame_add_python_one_thread',
              start_date=datetime(2013, 7, 26))
#----------------------------------------------------------------------
# mult

eval_frame_mult_all_threads = \
    Benchmark("pd.eval('df * df2 * df3 * df4')", common_setup,
              name='eval_frame_mult_all_threads',
              start_date=datetime(2013, 7, 21))

eval_frame_mult_one_thread = \
    Benchmark("pd.eval('df * df2 * df3 * df4')", setup,
              name='eval_frame_mult_one_thread',
              start_date=datetime(2013, 7, 26))

eval_frame_mult_python = \
    Benchmark("pd.eval('df * df2 * df3 * df4', engine='python')",
              common_setup,
              name='eval_frame_mult_python', start_date=datetime(2013, 7, 21))

eval_frame_mult_python_one_thread = \
    Benchmark("pd.eval('df * df2 * df3 * df4', engine='python')", setup,
              name='eval_frame_mult_python_one_thread',
              start_date=datetime(2013, 7, 26))

#----------------------------------------------------------------------
# multi and

eval_frame_and_all_threads = \
    Benchmark("pd.eval('(df > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)')",
              common_setup,
              name='eval_frame_and_all_threads',
              start_date=datetime(2013, 7, 21))

eval_frame_and_one_thread = \
    Benchmark("pd.eval('(df > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)')", setup,
              name='eval_frame_and_one_thread',
              start_date=datetime(2013, 7, 26))

eval_frame_and_python = \
    Benchmark("pd.eval('(df > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)', engine='python')",
              common_setup, name='eval_frame_and_python',
              start_date=datetime(2013, 7, 21))

eval_frame_and_one_thread = \
    Benchmark("pd.eval('(df > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)', engine='python')",
              setup,
              name='eval_frame_and_python_one_thread',
              start_date=datetime(2013, 7, 26))

#--------------------------------------------------------------------
# chained comp
eval_frame_chained_cmp_all_threads = \
    Benchmark("pd.eval('df < df2 < df3 < df4')", common_setup,
              name='eval_frame_chained_cmp_all_threads',
              start_date=datetime(2013, 7, 21))

eval_frame_chained_cmp_one_thread = \
    Benchmark("pd.eval('df < df2 < df3 < df4')", setup,
              name='eval_frame_chained_cmp_one_thread',
              start_date=datetime(2013, 7, 26))

eval_frame_chained_cmp_python = \
    Benchmark("pd.eval('df < df2 < df3 < df4', engine='python')",
              common_setup, name='eval_frame_chained_cmp_python',
              start_date=datetime(2013, 7, 26))

eval_frame_chained_cmp_one_thread = \
    Benchmark("pd.eval('df < df2 < df3 < df4', engine='python')", setup,
              name='eval_frame_chained_cmp_python_one_thread',
              start_date=datetime(2013, 7, 26))


common_setup = """from .pandas_vb_common import *
"""

setup = common_setup + """
N = 1000000
halfway = N // 2 - 1
index = date_range('20010101', periods=N, freq='T')
s = Series(index)
ts = s.iloc[halfway]
"""

series_setup = setup + """
df = DataFrame({'dates': s.values})
"""

query_datetime_series = Benchmark("df.query('dates < @ts')",
                                  series_setup,
                                  start_date=datetime(2013, 9, 27))

index_setup = setup + """
df = DataFrame({'a': np.random.randn(N)}, index=index)
"""

query_datetime_index = Benchmark("df.query('index < @ts')",
                                 index_setup, start_date=datetime(2013, 9, 27))

setup = setup + """
N = 1000000
df = DataFrame({'a': np.random.randn(N)})
min_val = df['a'].min()
max_val = df['a'].max()
"""

query_with_boolean_selection = Benchmark("df.query('(a >= @min_val) & (a <= @max_val)')",
                                         setup, start_date=datetime(2013, 9, 27))

