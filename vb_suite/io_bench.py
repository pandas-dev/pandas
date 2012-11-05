from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

#----------------------------------------------------------------------
# read_csv

setup1 = common_setup + """
index = [rands(10) for _ in xrange(10000)]
df = DataFrame({'float1' : randn(10000),
                'float2' : randn(10000),
                'string1' : ['foo'] * 10000,
                'bool1' : [True] * 10000,
                'int1' : np.random.randint(0, 100000, size=10000)},
               index=index)
df.to_csv('__test__.csv')
"""

read_csv_standard = Benchmark("read_csv('__test__.csv')", setup1,
                              start_date=datetime(2011, 9, 15))


#----------------------------------------------------------------------
# write_csv

setup2 = common_setup + """
index = [rands(10) for _ in xrange(10000)]
df = DataFrame({'float1' : randn(10000),
                'float2' : randn(10000),
                'string1' : ['foo'] * 10000,
                'bool1' : [True] * 10000,
                'int1' : np.random.randint(0, 100000, size=10000)},
               index=index)
"""

write_csv_standard = Benchmark("df.to_csv('__test__.csv')", setup2,
                               start_date=datetime(2011, 9, 15))

#----------------------------------
setup = common_setup + """
df = DataFrame(np.random.randn(3000, 30))
"""
frame_to_csv = Benchmark("df.to_csv('__test__.csv')", setup,
                         start_date=datetime(2011, 1, 1))
