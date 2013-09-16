from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
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
              start_date=datetime(2012, 7, 21))

eval_frame_mult_one_thread = \
    Benchmark("pd.eval('df * df2 * df3 * df4')", setup,
              name='eval_frame_mult_one_thread',
              start_date=datetime(2012, 7, 26))

eval_frame_mult_python = \
    Benchmark("pdl.eval('df * df2 * df3 * df4', engine='python')",
              common_setup,
              name='eval_frame_mult_python', start_date=datetime(2013, 7, 21))

eval_frame_mult_python_one_thread = \
    Benchmark("pd.eval('df * df2 * df3 * df4', engine='python')", setup,
              name='eval_frame_mult_python_one_thread',
              start_date=datetime(2012, 7, 26))

#----------------------------------------------------------------------
# multi and

eval_frame_and_all_threads = \
    Benchmark("pd.eval('(df > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)')",
              common_setup,
              name='eval_frame_and_all_threads',
              start_date=datetime(2012, 7, 21))

eval_frame_and_one_thread = \
    Benchmark("pd.eval('(df > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)')", setup,
              name='eval_frame_and_one_thread',
              start_date=datetime(2012, 7, 26))

setup = common_setup
eval_frame_and_python = \
    Benchmark("pd.eval('(df > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)', engine='python')",
              common_setup, name='eval_frame_and_python',
              start_date=datetime(2013, 7, 21))

eval_frame_and_one_thread = \
    Benchmark("pd.eval('(df > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)', engine='python')",
              setup,
              name='eval_frame_and_python_one_thread',
              start_date=datetime(2012, 7, 26))

#--------------------------------------------------------------------
# chained comp
eval_frame_chained_cmp_all_threads = \
    Benchmark("pd.eval('df < df2 < df3 < df4')", common_setup,
              name='eval_frame_chained_cmp_all_threads',
              start_date=datetime(2012, 7, 21))

eval_frame_chained_cmp_one_thread = \
    Benchmark("pd.eval('df < df2 < df3 < df4')", setup,
              name='eval_frame_chained_cmp_one_thread',
              start_date=datetime(2012, 7, 26))

setup = common_setup
eval_frame_chained_cmp_python = \
    Benchmark("pd.eval('df < df2 < df3 < df4', engine='python')",
              common_setup, name='eval_frame_chained_cmp_python',
              start_date=datetime(2013, 7, 26))

eval_frame_chained_cmp_one_thread = \
    Benchmark("pd.eval('df < df2 < df3 < df4', engine='python')", setup,
              name='eval_frame_chained_cmp_python_one_thread',
              start_date=datetime(2012, 7, 26))
