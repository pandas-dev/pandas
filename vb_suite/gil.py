from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

basic = common_setup + """
from pandas.util.testing import test_parallel

N = 1000000
ngroups = 1000
np.random.seed(1234)

df = DataFrame({'key' : np.random.randint(0,ngroups,size=N),
                'data' : np.random.randn(N) })
"""

setup = basic + """

def f():
    df.groupby('key')['data'].sum()

# run consecutivily
def g2():
    for i in range(2):
        f()
def g4():
    for i in range(4):
        f()
def g8():
    for i in range(8):
        f()

# run in parallel
@test_parallel(num_threads=2)
def pg2():
    f()

@test_parallel(num_threads=4)
def pg4():
    f()

@test_parallel(num_threads=8)
def pg8():
    f()

"""

nogil_groupby_sum_4 = Benchmark(
    'pg4()', setup,
    start_date=datetime(2015, 1, 1))

nogil_groupby_sum_8 = Benchmark(
    'pg8()', setup,
    start_date=datetime(2015, 1, 1))


#### test all groupby funcs ####

setup = basic + """

@test_parallel(num_threads=2)
def pg2():
    df.groupby('key')['data'].func()

"""

for f in ['sum','prod','var','count','min','max','mean','last']:

    name = "nogil_groupby_{f}_2".format(f=f)
    bmark = Benchmark('pg2()', setup.replace('func',f), start_date=datetime(2015, 1, 1))
    bmark.name = name
    globals()[name] = bmark

del bmark


#### test take_1d ####
setup = basic + """
from pandas.core import common as com

N = 1e7
df = DataFrame({'int64' : np.arange(N,dtype='int64'),
                'float64' : np.arange(N,dtype='float64')})
indexer = np.arange(100,len(df)-100)

@test_parallel(num_threads=2)
def take_1d_pg2_int64():
    com.take_1d(df.int64.values,indexer)

@test_parallel(num_threads=2)
def take_1d_pg2_float64():
    com.take_1d(df.float64.values,indexer)

"""

nogil_take1d_float64 = Benchmark('take_1d_pg2_int64()', setup, start_date=datetime(2015, 1, 1))
nogil_take1d_int64 = Benchmark('take_1d_pg2_float64()', setup, start_date=datetime(2015, 1, 1))
