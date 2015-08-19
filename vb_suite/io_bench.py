from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
from StringIO import StringIO
"""

#----------------------------------------------------------------------
# read_csv

setup1 = common_setup + """
index = tm.makeStringIndex(10000)
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

#----------------------------------
# skiprows

setup1 = common_setup + """
index = tm.makeStringIndex(20000)
df = DataFrame({'float1' : randn(20000),
                'float2' : randn(20000),
                'string1' : ['foo'] * 20000,
                'bool1' : [True] * 20000,
                'int1' : np.random.randint(0, 200000, size=20000)},
               index=index)
df.to_csv('__test__.csv')
"""

read_csv_skiprows = Benchmark("read_csv('__test__.csv', skiprows=10000)", setup1,
                              start_date=datetime(2011, 9, 15))

#----------------------------------------------------------------------
# write_csv

setup2 = common_setup + """
index = tm.makeStringIndex(10000)
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
#----------------------------------

setup = common_setup + """
df=DataFrame({'A':range(50000)})
df['B'] = df.A + 1.0
df['C'] = df.A + 2.0
df['D'] = df.A + 3.0
"""
frame_to_csv2 = Benchmark("df.to_csv('__test__.csv')", setup,
                         start_date=datetime(2011, 1, 1))

#----------------------------------
setup = common_setup + """
from pandas import concat, Timestamp

def create_cols(name):
    return [ "%s%03d" % (name,i) for i in xrange(5) ]
df_float  = DataFrame(np.random.randn(5000, 5),dtype='float64',columns=create_cols('float'))
df_int    = DataFrame(np.random.randn(5000, 5),dtype='int64',columns=create_cols('int'))
df_bool   = DataFrame(True,index=df_float.index,columns=create_cols('bool'))
df_object = DataFrame('foo',index=df_float.index,columns=create_cols('object'))
df_dt     = DataFrame(Timestamp('20010101'),index=df_float.index,columns=create_cols('date'))

# add in some nans
df_float.ix[30:500,1:3] = np.nan

df        = concat([ df_float, df_int, df_bool, df_object, df_dt ], axis=1)

"""
frame_to_csv_mixed = Benchmark("df.to_csv('__test__.csv')", setup,
                               start_date=datetime(2012, 6, 1))

#----------------------------------------------------------------------
# parse dates, ISO8601 format

setup = common_setup + """
rng = date_range('1/1/2000', periods=1000)
data = '\\n'.join(rng.map(lambda x: x.strftime("%Y-%m-%d %H:%M:%S")))
"""

stmt = ("read_csv(StringIO(data), header=None, names=['foo'], "
        "         parse_dates=['foo'])")
read_parse_dates_iso8601 = Benchmark(stmt, setup,
                                     start_date=datetime(2012, 3, 1))

setup = common_setup + """
rng = date_range('1/1/2000', periods=1000)
data = DataFrame(rng, index=rng)
"""

stmt = ("data.to_csv('__test__.csv', date_format='%Y%m%d')")

frame_to_csv_date_formatting = Benchmark(stmt, setup,
                                     start_date=datetime(2013, 9, 1))

#----------------------------------------------------------------------
# infer datetime format

setup = common_setup + """
rng = date_range('1/1/2000', periods=1000)
data = '\\n'.join(rng.map(lambda x: x.strftime("%Y-%m-%d %H:%M:%S")))
"""

stmt = ("read_csv(StringIO(data), header=None, names=['foo'], "
        "         parse_dates=['foo'], infer_datetime_format=True)")

read_csv_infer_datetime_format_iso8601 = Benchmark(stmt, setup)

setup = common_setup + """
rng = date_range('1/1/2000', periods=1000)
data = '\\n'.join(rng.map(lambda x: x.strftime("%Y%m%d")))
"""

stmt = ("read_csv(StringIO(data), header=None, names=['foo'], "
        "         parse_dates=['foo'], infer_datetime_format=True)")

read_csv_infer_datetime_format_ymd = Benchmark(stmt, setup)

setup = common_setup + """
rng = date_range('1/1/2000', periods=1000)
data = '\\n'.join(rng.map(lambda x: x.strftime("%m/%d/%Y %H:%M:%S.%f")))
"""

stmt = ("read_csv(StringIO(data), header=None, names=['foo'], "
        "         parse_dates=['foo'], infer_datetime_format=True)")

read_csv_infer_datetime_format_custom = Benchmark(stmt, setup)
