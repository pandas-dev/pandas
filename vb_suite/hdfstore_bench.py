from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
import os

f = '__test__.h5'
def remove(f):
   try:
       os.remove(f)
   except:
       pass

"""

#----------------------------------------------------------------------
# get from a store

setup1 = common_setup + """
index = [rands(10) for _ in xrange(100000)]
df = DataFrame({'float1' : randn(100000),
                'float2' : randn(100000)},
               index=index)
remove(f)
store = HDFStore(f)
store.put('df1',df)
"""

read_store = Benchmark("store.get('df1')", setup1, cleanup = "store.close()",
                       start_date=datetime(2011, 9, 15))


#----------------------------------------------------------------------
# write to a store

setup2 = common_setup + """
index = [rands(10) for _ in xrange(100000)]
df = DataFrame({'float1' : randn(100000),
                'float2' : randn(100000)},
               index=index)
remove(f)
store = HDFStore(f)
"""

write_store = Benchmark("store.put('df2',df)", setup2, cleanup = "store.close()",
                        start_date=datetime(2011, 9, 15))

#----------------------------------------------------------------------
# get from a store (mixed)

setup3 = common_setup + """
index = [rands(10) for _ in xrange(100000)]
df = DataFrame({'float1' : randn(100000),
                'float2' : randn(100000),
                'string1' : ['foo'] * 100000,
                'bool1' : [True] * 100000,
                'int1' : np.random.randint(0, 1000000, size=100000)},
               index=index)
remove(f)
store = HDFStore(f)
store.put('df3',df)
"""

read_store_mixed = Benchmark("store.get('df3')", setup3, cleanup = "store.close()",
                             start_date=datetime(2011, 9, 15))


#----------------------------------------------------------------------
# write to a store (mixed)

setup4 = common_setup + """
index = [rands(10) for _ in xrange(100000)]
df = DataFrame({'float1' : randn(100000),
                'float2' : randn(100000),
                'string1' : ['foo'] * 100000,
                'bool1' : [True] * 100000,
                'int1' : np.random.randint(0, 1000000, size=100000)},
               index=index)
remove(f)
store = HDFStore(f)
"""

write_store_mixed = Benchmark("store.put('df4',df)", setup4, cleanup = "store.close()",
                               start_date=datetime(2011, 9, 15))

#----------------------------------------------------------------------
# get from a table (mixed)

setup5 = common_setup + """
index = [rands(10) for _ in xrange(100000)]
df = DataFrame({'float1' : randn(100000),
                'float2' : randn(100000),
                'string1' : ['foo'] * 100000,
                'bool1' : [True] * 100000,
                'int1' : np.random.randint(0, 1000000, size=100000)},
               index=index)

remove(f)
store = HDFStore(f)
store.append('df5',df)
"""

read_store_table_mixed = Benchmark("store.select('df5')", setup5, cleanup = "store.close()",
                                   start_date=datetime(2011, 9, 15))


#----------------------------------------------------------------------
# write to a table (mixed)

setup6 = common_setup + """
index = [rands(10) for _ in xrange(100000)]
df = DataFrame({'float1' : randn(100000),
                'float2' : randn(100000),
                'string1' : ['foo'] * 100000,
                'bool1' : [True] * 100000,
                'int1' : np.random.randint(0, 100000, size=100000)},
               index=index)
remove(f)
store = HDFStore(f)
"""

write_store_table_mixed = Benchmark("store.append('df6',df)", setup6, cleanup = "store.close()",
                                    start_date=datetime(2011, 9, 15))

#----------------------------------------------------------------------
# select from a table

setup7 = common_setup + """
index = [rands(10) for _ in xrange(100000)]
df = DataFrame({'float1' : randn(100000),
                'float2' : randn(100000) },
               index=index)

remove(f)
store = HDFStore(f)
store.append('df7',df)
"""

read_store_table = Benchmark("store.select('df7')", setup7, cleanup = "store.close()",
                             start_date=datetime(2011, 9, 15))


#----------------------------------------------------------------------
# write to a table

setup8 = common_setup + """
index = [rands(10) for _ in xrange(100000)]
df = DataFrame({'float1' : randn(100000),
                'float2' : randn(100000) },
               index=index)
remove(f)
store = HDFStore(f)
"""

write_store_table = Benchmark("store.append('df8',df)", setup8, cleanup = "store.close()",
                              start_date=datetime(2011, 9, 15))

#----------------------------------------------------------------------
# get from a table (wide)

setup9 = common_setup + """
df = DataFrame(np.random.randn(100000,200))

remove(f)
store = HDFStore(f)
store.append('df9',df)
"""

read_store_table_wide = Benchmark("store.select('df9')", setup9, cleanup = "store.close()",
                             start_date=datetime(2011, 9, 15))


#----------------------------------------------------------------------
# write to a table (wide)

setup10 = common_setup + """
df = DataFrame(np.random.randn(100000,200))

remove(f)
store = HDFStore(f)
"""

write_store_table_wide = Benchmark("store.append('df10',df)", setup10, cleanup = "store.close()",
                              start_date=datetime(2011, 9, 15))

#----------------------------------------------------------------------
# get from a table (wide) (indexed)

setup11 = common_setup + """
index = date_range('1/1/2000', periods = 100000)
df = DataFrame(np.random.randn(100000,200), index = index)

remove(f)
store = HDFStore(f)
store.append('df11',df)
store.create_table_index('df11')
"""

query_store_table_wide = Benchmark("store.select('df11', [ ('index', '>', df.index[10000]), ('index', '<', df.index[15000]) ])", setup11, cleanup = "store.close()",
                                   start_date=datetime(2011, 9, 15))


#----------------------------------------------------------------------
# query from a table (indexed)

setup12 = common_setup + """
index = date_range('1/1/2000', periods = 100000)
df = DataFrame({'float1' : randn(100000),
                'float2' : randn(100000) },
               index=index)

remove(f)
store = HDFStore(f)
store.append('df12',df)
store.create_table_index('df12')
"""

query_store_table = Benchmark("store.select('df12', [ ('index', '>', df.index[10000]), ('index', '<', df.index[15000]) ])", setup12, cleanup = "store.close()",
                              start_date=datetime(2011, 9, 15))

