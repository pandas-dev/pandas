from vbench.api import Benchmark
from datetime import datetime

start_date = datetime(2013, 5, 1)

common_setup = """from pandas_vb_common import *
import os
import pandas as pd
from pandas.core import common as com
from pandas.compat import BytesIO
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
df2.to_hdf(f,'df',format='table')
"""

packers_read_hdf_table = Benchmark("pd.read_hdf(f,'df')", setup, start_date=start_date)

setup = common_setup + """
"""

packers_write_hdf_table = Benchmark("df2.to_hdf(f,'df',table=True)", setup, cleanup="remove(f)", start_date=start_date)

#----------------------------------------------------------------------
# sql

setup = common_setup + """
import sqlite3
from sqlalchemy import create_engine
engine = create_engine('sqlite:///:memory:')

df2.to_sql('table', engine, if_exists='replace')
"""

packers_read_sql= Benchmark("pd.read_sql_table('table', engine)", setup, start_date=start_date)

setup = common_setup + """
import sqlite3
from sqlalchemy import create_engine
engine = create_engine('sqlite:///:memory:')
"""

packers_write_sql = Benchmark("df2.to_sql('table', engine, if_exists='replace')", setup, start_date=start_date)

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
packers_write_json_T = Benchmark("df.to_json(f,orient='columns')", setup, cleanup="remove(f)", start_date=start_date)

setup = common_setup + """
from numpy.random import randint
from collections import OrderedDict

cols = [
  lambda i: ("{0}_timedelta".format(i), [pd.Timedelta('%d seconds' % randrange(1e6)) for _ in range(N)]),
  lambda i: ("{0}_int".format(i), randint(1e8, size=N)),
  lambda i: ("{0}_timestamp".format(i), [pd.Timestamp( 1418842918083256000 + randrange(1e9, 1e18, 200)) for _ in range(N)])
  ]
df_mixed = DataFrame(OrderedDict([cols[i % len(cols)](i) for i in range(C)]),
                     index=index)
"""
packers_write_json_mixed_delta_int_tstamp = Benchmark("df_mixed.to_json(f,orient='split')", setup, cleanup="remove(f)", start_date=start_date)

setup = common_setup + """
from numpy.random import randint
from collections import OrderedDict
cols = [
  lambda i: ("{0}_float".format(i), randn(N)),
  lambda i: ("{0}_int".format(i), randint(1e8, size=N))
  ]
df_mixed = DataFrame(OrderedDict([cols[i % len(cols)](i) for i in range(C)]),
                     index=index)
"""
packers_write_json_mixed_float_int = Benchmark("df_mixed.to_json(f,orient='index')", setup, cleanup="remove(f)", start_date=start_date)
packers_write_json_mixed_float_int_T = Benchmark("df_mixed.to_json(f,orient='columns')", setup, cleanup="remove(f)", start_date=start_date)

setup = common_setup + """
from numpy.random import randint
from collections import OrderedDict
cols = [
  lambda i: ("{0}_float".format(i), randn(N)),
  lambda i: ("{0}_int".format(i), randint(1e8, size=N)),
  lambda i: ("{0}_str".format(i), ['%08x'%randrange(16**8) for _ in range(N)])
  ]
df_mixed = DataFrame(OrderedDict([cols[i % len(cols)](i) for i in range(C)]),
                     index=index)
"""
packers_write_json_mixed_float_int_str = Benchmark("df_mixed.to_json(f,orient='split')", setup, cleanup="remove(f)", start_date=start_date)

#----------------------------------------------------------------------
# stata

setup = common_setup + """
df.to_stata(f, {'index': 'tc'})
"""
packers_read_stata = Benchmark("pd.read_stata(f)", setup, start_date=start_date)

packers_write_stata = Benchmark("df.to_stata(f, {'index': 'tc'})", setup, cleanup="remove(f)", start_date=start_date)

setup = common_setup + """
df['int8_'] = [randint(np.iinfo(np.int8).min, np.iinfo(np.int8).max - 27) for _ in range(N)]
df['int16_'] = [randint(np.iinfo(np.int16).min, np.iinfo(np.int16).max - 27) for _ in range(N)]
df['int32_'] = [randint(np.iinfo(np.int32).min, np.iinfo(np.int32).max - 27) for _ in range(N)]
df['float32_'] = np.array(randn(N), dtype=np.float32)
df.to_stata(f, {'index': 'tc'})
"""

packers_read_stata_with_validation = Benchmark("pd.read_stata(f)", setup, start_date=start_date)

packers_write_stata_with_validation = Benchmark("df.to_stata(f, {'index': 'tc'})", setup, cleanup="remove(f)", start_date=start_date)

#----------------------------------------------------------------------
# Excel - alternative writers
setup = common_setup + """
bio = BytesIO()
"""

excel_writer_bench = """
bio.seek(0)
writer = pd.io.excel.ExcelWriter(bio, engine='{engine}')
df[:2000].to_excel(writer)
writer.save()
"""

benchmark_xlsxwriter = excel_writer_bench.format(engine='xlsxwriter')

packers_write_excel_xlsxwriter = Benchmark(benchmark_xlsxwriter, setup)

benchmark_openpyxl = excel_writer_bench.format(engine='openpyxl')

packers_write_excel_openpyxl = Benchmark(benchmark_openpyxl, setup)

benchmark_xlwt = excel_writer_bench.format(engine='xlwt')

packers_write_excel_xlwt = Benchmark(benchmark_xlwt, setup)


#----------------------------------------------------------------------
# Excel - reader

setup = common_setup + """
bio = BytesIO()
writer = pd.io.excel.ExcelWriter(bio, engine='xlsxwriter')
df[:2000].to_excel(writer)
writer.save()
"""

benchmark_read_excel="""
bio.seek(0)
pd.read_excel(bio)
"""

packers_read_excel = Benchmark(benchmark_read_excel, setup)
