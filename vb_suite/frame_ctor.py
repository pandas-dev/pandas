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
data = frame.to_dict()
some_dict = data.values()[0]
dict_list = [dict(zip(columns, row)) for row in frame.values]
"""

frame_ctor_nested_dict = Benchmark("DataFrame(data)", setup)

# From JSON-like stuff

frame_ctor_list_of_dict = Benchmark("DataFrame(dict_list)", setup,
                                    start_date=datetime(2011, 12, 20))

series_ctor_from_dict = Benchmark("Series(some_dict)", setup)
