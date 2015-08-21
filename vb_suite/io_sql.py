from vbench.api import Benchmark
from datetime import datetime

common_setup = """from .pandas_vb_common import *
import sqlite3
import sqlalchemy
from sqlalchemy import create_engine

engine = create_engine('sqlite:///:memory:')
con = sqlite3.connect(':memory:')
"""

sdate = datetime(2014, 6, 1)


#-------------------------------------------------------------------------------
# to_sql

setup = common_setup + """
index = tm.makeStringIndex(10000)
df = DataFrame({'float1' : randn(10000),
                'float2' : randn(10000),
                'string1' : ['foo'] * 10000,
                'bool1' : [True] * 10000,
                'int1' : np.random.randint(0, 100000, size=10000)},
               index=index)
"""

sql_write_sqlalchemy = Benchmark("df.to_sql('test1', engine, if_exists='replace')",
                                 setup, start_date=sdate)

sql_write_fallback = Benchmark("df.to_sql('test1', con, if_exists='replace')",
                               setup, start_date=sdate)


#-------------------------------------------------------------------------------
# read_sql

setup = common_setup + """
index = tm.makeStringIndex(10000)
df = DataFrame({'float1' : randn(10000),
                'float2' : randn(10000),
                'string1' : ['foo'] * 10000,
                'bool1' : [True] * 10000,
                'int1' : np.random.randint(0, 100000, size=10000)},
            index=index)
df.to_sql('test2', engine, if_exists='replace')
df.to_sql('test2', con, if_exists='replace')
"""

sql_read_query_sqlalchemy = Benchmark("read_sql_query('SELECT * FROM test2', engine)",
                                      setup, start_date=sdate)

sql_read_query_fallback = Benchmark("read_sql_query('SELECT * FROM test2', con)",
                                    setup, start_date=sdate)

sql_read_table_sqlalchemy = Benchmark("read_sql_table('test2', engine)",
                                      setup, start_date=sdate)


#-------------------------------------------------------------------------------
# type specific write

setup = common_setup + """
df = DataFrame({'float' : randn(10000),
                'string' : ['foo'] * 10000,
                'bool' : [True] * 10000,
                'datetime' : date_range('2000-01-01', periods=10000, freq='s')})
df.loc[1000:3000, 'float'] = np.nan
"""

sql_float_write_sqlalchemy = \
    Benchmark("df[['float']].to_sql('test_float', engine, if_exists='replace')",
              setup, start_date=sdate)

sql_float_write_fallback = \
    Benchmark("df[['float']].to_sql('test_float', con, if_exists='replace')",
              setup, start_date=sdate)

sql_string_write_sqlalchemy = \
    Benchmark("df[['string']].to_sql('test_string', engine, if_exists='replace')",
              setup, start_date=sdate)

sql_string_write_fallback = \
    Benchmark("df[['string']].to_sql('test_string', con, if_exists='replace')",
              setup, start_date=sdate)

sql_datetime_write_sqlalchemy = \
    Benchmark("df[['datetime']].to_sql('test_datetime', engine, if_exists='replace')",
            setup, start_date=sdate)

#sql_datetime_write_fallback = \
#   Benchmark("df[['datetime']].to_sql('test_datetime', con, if_exists='replace')",
#           setup3, start_date=sdate)

#-------------------------------------------------------------------------------
# type specific read

setup = common_setup + """
df = DataFrame({'float' : randn(10000),
                'datetime' : date_range('2000-01-01', periods=10000, freq='s')})
df['datetime_string'] = df['datetime'].map(str)

df.to_sql('test_type', engine, if_exists='replace')
df[['float', 'datetime_string']].to_sql('test_type', con, if_exists='replace')
"""

sql_float_read_query_sqlalchemy = \
    Benchmark("read_sql_query('SELECT float FROM test_type', engine)",
            setup, start_date=sdate)

sql_float_read_table_sqlalchemy = \
    Benchmark("read_sql_table('test_type', engine, columns=['float'])",
            setup, start_date=sdate)

sql_float_read_query_fallback = \
    Benchmark("read_sql_query('SELECT float FROM test_type', con)",
            setup, start_date=sdate)

sql_datetime_read_as_native_sqlalchemy = \
    Benchmark("read_sql_table('test_type', engine, columns=['datetime'])",
            setup, start_date=sdate)

sql_datetime_read_and_parse_sqlalchemy = \
    Benchmark("read_sql_table('test_type', engine, columns=['datetime_string'], parse_dates=['datetime_string'])",
            setup, start_date=sdate)
