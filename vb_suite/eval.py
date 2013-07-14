from vbench.benchmark import Benchmark
from datetime import datetime

setup = """from pandas_vb_common import *
import pandas as pd
df  = DataFrame(np.random.randn(20000, 100))
df2 = DataFrame(np.random.randn(20000, 100))
df3 = DataFrame(np.random.randn(20000, 100))
df4 = DataFrame(np.random.randn(20000, 100))
"""

SECTION = 'Eval'

#----------------------------------------------------------------------
# binary ops

#----------------------------------------------------------------------
# add

frame_add_eval = \
    Benchmark("pd.eval('df + df2 + df3 + df4')", setup, name='frame_add_eval',
              start_date=datetime(2013, 7, 21))

frame_add_python = \
    Benchmark("pd.eval('df + df2 + df3 + df4', engine='python')", setup,
              name='frame_add_python', start_date=datetime(2013, 7, 21))

#----------------------------------------------------------------------
# mult

frame_mult_eval = \
    Benchmark("pd.eval('df * df2 * df3 * df4')", setup, name='frame_mult_eval',
              start_date=datetime(2012, 7, 21))

frame_mult_python = \
    Benchmark("pdl.eval('df * df2 * df3 * df4', engine='python')", setup,
              name='frame_mult_python', start_date=datetime(2013, 7, 21))

#----------------------------------------------------------------------
# multi and

frame_and_eval = \
    Benchmark("pd.eval('(df > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)')", setup,
              name='frame_and_eval', start_date=datetime(2012, 7, 21))

frame_and_python = \
    Benchmark("pd.eval('(df > 0) & (df2 > 0) & (df3 > 0) & (df4 > 0)', "
              "engine='python')", setup, name='frame_and_python',
              start_date=datetime(2013, 7, 21))
