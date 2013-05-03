from vbench.api import Benchmark
from datetime import datetime

start_date = datetime(2013, 5, 1)

common_setup = """from pandas_vb_common import *
import os
from pandas.io import packers
from pandas.core import common as com

f = '__test__.msg'
def remove(f):
   try:
       os.remove(f)
   except:
       pass

"""

#----------------------------------------------------------------------
# read a pack

setup1 = common_setup + """
index = date_range('20000101',periods=25000,freq='H')
df = DataFrame({'float1' : randn(25000),
                'float2' : randn(25000)},
               index=index)
remove(f)
packers.save(f,df)
"""

read_pack = Benchmark("packers.load(f)", setup1, 
                       start_date=start_date)


#----------------------------------------------------------------------
# write to a pack

setup2 = common_setup + """
index = date_range('20000101',periods=25000,freq='H')
df = DataFrame({'float1' : randn(25000),
                'float2' : randn(25000)},
               index=index)
remove(f)
"""

write_pack = Benchmark(
    "packers.save(f,df)", setup2, cleanup="remove(f)",
    start_date=start_date)

#----------------------------------------------------------------------
# read a pickle

setup1 = common_setup + """
index = date_range('20000101',periods=25000,freq='H')
df = DataFrame({'float1' : randn(25000),
                'float2' : randn(25000)},
               index=index)
remove(f)
df.save(f)
"""

read_pickle = Benchmark("com.load(f)", setup1, 
                       start_date=start_date)


#----------------------------------------------------------------------
# write to a pickle

setup2 = common_setup + """
index = date_range('20000101',periods=25000,freq='H')
df = DataFrame({'float1' : randn(25000),
                'float2' : randn(25000)},
               index=index)
remove(f)
"""

write_pickle = Benchmark(
    "df.save(f)", setup2, cleanup="remove(f)",
    start_date=start_date)
