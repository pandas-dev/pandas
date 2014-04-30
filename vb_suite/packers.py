from vbench.api import Benchmark
from datetime import datetime

start_date = datetime(2013, 5, 1)

common_setup = """from pandas_vb_common import *
import os
import pandas as pd
from pandas.core import common as com
from random import randrange

f = '__test__.msg'
def remove(f):
   try:
       os.remove(f)
   except:
       pass

N=100000
C=5
index = date_range('20000101',periods=N,freq='H')
df = DataFrame(dict([ ("float{0}".format(i),randn(N)) for i in range(C) ]),
               index=index)

N=100000
C=5
index = date_range('20000101',periods=N,freq='H')
df2 = DataFrame(dict([ ("float{0}".format(i),randn(N)) for i in range(C) ]),
                index=index)
df2['object'] = ['%08x'%randrange(16**8) for _ in range(N)]
remove(f)
"""

#----------------------------------------------------------------------
# msgpack

setup = common_setup + """
df2.to_msgpack(f)
"""

packers_read_pack = Benchmark("pd.read_msgpack(f)", setup, start_date=start_date)

setup = common_setup + """
"""

packers_write_pack = Benchmark("df2.to_msgpack(f)", setup, cleanup="remove(f)", start_date=start_date)

#----------------------------------------------------------------------
# pickle

setup = common_setup + """
df2.to_pickle(f)
"""

packers_read_pickle = Benchmark("pd.read_pickle(f)", setup, start_date=start_date)

setup = common_setup + """
"""

packers_write_pickle = Benchmark("df2.to_pickle(f)", setup, cleanup="remove(f)", start_date=start_date)

#----------------------------------------------------------------------
# csv

setup = common_setup + """
df.to_csv(f)
"""

packers_read_csv = Benchmark("pd.read_csv(f)", setup, start_date=start_date)

setup = common_setup + """
"""

packers_write_csv = Benchmark("df.to_csv(f)", setup, cleanup="remove(f)", start_date=start_date)

#----------------------------------------------------------------------
# hdf store

setup = common_setup + """
df2.to_hdf(f,'df')
"""

packers_read_hdf_store = Benchmark("pd.read_hdf(f,'df')", setup, start_date=start_date)

setup = common_setup + """
"""

packers_write_hdf_store = Benchmark("df2.to_hdf(f,'df')", setup, cleanup="remove(f)", start_date=start_date)

#----------------------------------------------------------------------
# hdf table

setup = common_setup + """
df2.to_hdf(f,'df',table=True)
"""

packers_read_hdf_table = Benchmark("pd.read_hdf(f,'df')", setup, start_date=start_date)

setup = common_setup + """
"""

packers_write_hdf_table = Benchmark("df2.to_hdf(f,'df',table=True)", setup, cleanup="remove(f)", start_date=start_date)

#----------------------------------------------------------------------
# json

setup_int_index = """
import numpy as np
df.index = np.arange(N)
"""

setup = common_setup + """
df.to_json(f,orient='split')
"""
packers_read_json_date_index = Benchmark("pd.read_json(f, orient='split')", setup, start_date=start_date)
setup = setup + setup_int_index
packers_read_json = Benchmark("pd.read_json(f, orient='split')", setup, start_date=start_date)

setup = common_setup + """
"""
packers_write_json_date_index = Benchmark("df.to_json(f,orient='split')", setup, cleanup="remove(f)", start_date=start_date)
setup = setup + setup_int_index
packers_write_json = Benchmark("df.to_json(f,orient='split')", setup, cleanup="remove(f)", start_date=start_date)
