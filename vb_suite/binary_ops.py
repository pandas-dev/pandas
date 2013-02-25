from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

SECTION = 'Binary ops'

#----------------------------------------------------------------------
# binary ops

setup = common_setup + """
df  = DataFrame(np.random.randn(100000, 100))
df2 = DataFrame(np.random.randn(100000, 100))
"""
frame_add = \
    Benchmark("df + df2", setup, name='frame_add',
              start_date=datetime(2012, 1, 1))

setup = common_setup + """
df  = DataFrame(np.random.randn(100000, 100))
df2 = DataFrame(np.random.randn(100000, 100))
"""
frame_mult = \
    Benchmark("df * df2", setup, name='frame_mult',
              start_date=datetime(2012, 1, 1))

setup = common_setup + """
df  = DataFrame(np.random.randn(100000, 100))
df2 = DataFrame(np.random.randn(100000, 100))
"""
frame_multi_and = \
    Benchmark("df[(df>0) & (df2>0)]", setup, name='frame_multi_and',
              start_date=datetime(2012, 1, 1))

