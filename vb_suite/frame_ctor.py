from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
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
