from vbench.api import Benchmark
from datetime import datetime

start_date = datetime(2012, 7, 1)

common_setup = """from .pandas_vb_common import *
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
index = tm.makeStringIndex(25000)
df = DataFrame({'float1' : randn(25000),
                'float2' : randn(25000)},
               index=index)
remove(f)
store = HDFStore(f)
store.put('df1',df)
"""

read_store = Benchmark("store.get('df1')", setup1, cleanup="store.close()",
                       start_date=start_date)


#----------------------------------------------------------------------
# write to a store

setup2 = common_setup + """
index = tm.makeStringIndex(25000)
df = DataFrame({'float1' : randn(25000),
                'float2' : randn(25000)},
               index=index)
remove(f)
store = HDFStore(f)
"""

write_store = Benchmark(
    "store.put('df2',df)", setup2, cleanup="store.close()",
    start_date=start_date)

#----------------------------------------------------------------------
# get from a store (mixed)

setup3 = common_setup + """
index = tm.makeStringIndex(25000)
df = DataFrame({'float1' : randn(25000),
                'float2' : randn(25000),
                'string1' : ['foo'] * 25000,
                'bool1' : [True] * 25000,
                'int1' : np.random.randint(0, 250000, size=25000)},
               index=index)
remove(f)
store = HDFStore(f)
store.put('df3',df)
"""

read_store_mixed = Benchmark(
    "store.get('df3')", setup3, cleanup="store.close()",
    start_date=start_date)


#----------------------------------------------------------------------
# write to a store (mixed)

setup4 = common_setup + """
index = tm.makeStringIndex(25000)
df = DataFrame({'float1' : randn(25000),
                'float2' : randn(25000),
                'string1' : ['foo'] * 25000,
                'bool1' : [True] * 25000,
                'int1' : np.random.randint(0, 250000, size=25000)},
               index=index)
remove(f)
store = HDFStore(f)
"""

write_store_mixed = Benchmark(
    "store.put('df4',df)", setup4, cleanup="store.close()",
    start_date=start_date)

#----------------------------------------------------------------------
# get from a table (mixed)

setup5 = common_setup + """
N=10000
index = tm.makeStringIndex(N)
df = DataFrame({'float1' : randn(N),
                'float2' : randn(N),
                'string1' : ['foo'] * N,
                'bool1' : [True] * N,
                'int1' : np.random.randint(0, N, size=N)},
               index=index)

remove(f)
store = HDFStore(f)
store.append('df5',df)
"""

read_store_table_mixed = Benchmark(
    "store.select('df5')", setup5, cleanup="store.close()",
    start_date=start_date)


#----------------------------------------------------------------------
# write to a table (mixed)

setup6 = common_setup + """
index = tm.makeStringIndex(25000)
df = DataFrame({'float1' : randn(25000),
                'float2' : randn(25000),
                'string1' : ['foo'] * 25000,
                'bool1' : [True] * 25000,
                'int1' : np.random.randint(0, 25000, size=25000)},
               index=index)
remove(f)
store = HDFStore(f)
"""

write_store_table_mixed = Benchmark(
    "store.append('df6',df)", setup6, cleanup="store.close()",
    start_date=start_date)

#----------------------------------------------------------------------
# select from a table

setup7 = common_setup + """
index = tm.makeStringIndex(25000)
df = DataFrame({'float1' : randn(25000),
                'float2' : randn(25000) },
               index=index)

remove(f)
store = HDFStore(f)
store.append('df7',df)
"""

read_store_table = Benchmark(
    "store.select('df7')", setup7, cleanup="store.close()",
    start_date=start_date)


#----------------------------------------------------------------------
# write to a table

setup8 = common_setup + """
index = tm.makeStringIndex(25000)
df = DataFrame({'float1' : randn(25000),
                'float2' : randn(25000) },
               index=index)
remove(f)
store = HDFStore(f)
"""

write_store_table = Benchmark(
    "store.append('df8',df)", setup8, cleanup="store.close()",
    start_date=start_date)

#----------------------------------------------------------------------
# get from a table (wide)

setup9 = common_setup + """
df = DataFrame(np.random.randn(25000,100))

remove(f)
store = HDFStore(f)
store.append('df9',df)
"""

read_store_table_wide = Benchmark(
    "store.select('df9')", setup9, cleanup="store.close()",
    start_date=start_date)


#----------------------------------------------------------------------
# write to a table (wide)

setup10 = common_setup + """
df = DataFrame(np.random.randn(25000,100))

remove(f)
store = HDFStore(f)
"""

write_store_table_wide = Benchmark(
    "store.append('df10',df)", setup10, cleanup="store.close()",
    start_date=start_date)

#----------------------------------------------------------------------
# get from a table (wide)

setup11 = common_setup + """
index = date_range('1/1/2000', periods = 25000)
df = DataFrame(np.random.randn(25000,100), index = index)

remove(f)
store = HDFStore(f)
store.append('df11',df)
"""

query_store_table_wide = Benchmark(
    "store.select('df11', [ ('index', '>', df.index[10000]), ('index', '<', df.index[15000]) ])", setup11, cleanup="store.close()",
    start_date=start_date)


#----------------------------------------------------------------------
# query from a table

setup12 = common_setup + """
index = date_range('1/1/2000', periods = 25000)
df = DataFrame({'float1' : randn(25000),
                'float2' : randn(25000) },
               index=index)

remove(f)
store = HDFStore(f)
store.append('df12',df)
"""

query_store_table = Benchmark(
    "store.select('df12', [ ('index', '>', df.index[10000]), ('index', '<', df.index[15000]) ])", setup12, cleanup="store.close()",
    start_date=start_date)

#----------------------------------------------------------------------
# select from a panel table

setup13 = common_setup + """
p = Panel(randn(20, 1000, 25), items= [ 'Item%03d' % i for i in range(20) ],
                   major_axis=date_range('1/1/2000', periods=1000), minor_axis = [ 'E%03d' % i for i in range(25) ])

remove(f)
store = HDFStore(f)
store.append('p1',p)
"""

read_store_table_panel = Benchmark(
    "store.select('p1')", setup13, cleanup="store.close()",
    start_date=start_date)


#----------------------------------------------------------------------
# write to a panel table

setup14 = common_setup + """
p = Panel(randn(20, 1000, 25), items= [ 'Item%03d' % i for i in range(20) ],
                   major_axis=date_range('1/1/2000', periods=1000), minor_axis = [ 'E%03d' % i for i in range(25) ])

remove(f)
store = HDFStore(f)
"""

write_store_table_panel = Benchmark(
    "store.append('p2',p)", setup14, cleanup="store.close()",
    start_date=start_date)

#----------------------------------------------------------------------
# write to a table (data_columns)

setup15 = common_setup + """
df = DataFrame(np.random.randn(10000,10),columns = [ 'C%03d' % i for i in range(10) ])

remove(f)
store = HDFStore(f)
"""

write_store_table_dc = Benchmark(
    "store.append('df15',df,data_columns=True)", setup15, cleanup="store.close()",
    start_date=start_date)

