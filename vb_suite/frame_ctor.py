from vbench.benchmark import Benchmark
from datetime import datetime
try:
    import pandas.tseries.offsets as offsets
except:
    import pandas.core.datetools as offsets

common_setup = """from pandas_vb_common import *
try:
    from pandas.tseries.offsets import *
except:
    from pandas.core.datetools import *
"""

#----------------------------------------------------------------------
# Creation from nested dict

setup = common_setup + """
N, K = 5000, 50
index = [rands(10) for _ in xrange(N)]
columns = [rands(10) for _ in xrange(K)]
frame = DataFrame(np.random.randn(N, K), index=index, columns=columns)

try:
    data = frame.to_dict()
except:
    data = frame.toDict()

some_dict = data.values()[0]
dict_list = [dict(zip(columns, row)) for row in frame.values]
"""

frame_ctor_nested_dict = Benchmark("DataFrame(data)", setup)

# From JSON-like stuff
frame_ctor_list_of_dict = Benchmark("DataFrame(dict_list)", setup,
                                    start_date=datetime(2011, 12, 20))

series_ctor_from_dict = Benchmark("Series(some_dict)", setup)

# nested dict, integer indexes, regression described in #621
setup = common_setup + """
data = dict((i,dict((j,float(j)) for j in xrange(100))) for i in xrange(2000))
"""
frame_ctor_nested_dict_int64 = Benchmark("DataFrame(data)", setup)

# dynamically generate benchmarks for every offset
dynamic_benchmarks = {}
n_steps = [1, 2]
for offset in offsets.__all__:
    for n in n_steps:
        setup = common_setup + """
df = DataFrame(np.random.randn(1000,10),index=date_range('1/1/1900',periods=1000,freq={}({})))
d = dict([ (col,df[col]) for col in df.columns ])
""".format(offset, n)
        key = 'frame_ctor_dtindex_{}({})'.format(offset, n)
        dynamic_benchmarks[key] = Benchmark("DataFrame(d)", setup, name=key)

# Have to stuff them in globals() so vbench detects them
globals().update(dynamic_benchmarks)

# from a mi-series
setup = common_setup + """
mi = MultiIndex.from_tuples([(x,y) for x in range(100) for y in range(100)])
s = Series(randn(10000), index=mi)
"""
frame_from_series = Benchmark("DataFrame(s)", setup)

#----------------------------------------------------------------------
# get_numeric_data

setup = common_setup + """
df = DataFrame(randn(10000, 25))
df['foo'] = 'bar'
df['bar'] = 'baz'
df = df.consolidate()
"""

frame_get_numeric_data = Benchmark('df._get_numeric_data()', setup,
                                   start_date=datetime(2011, 11, 1))
