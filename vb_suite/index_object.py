from vbench.benchmark import Benchmark
from datetime import datetime

SECTION = "Index / MultiIndex objects"


common_setup = """from pandas_vb_common import *
"""

#----------------------------------------------------------------------
# intersection, union

setup = common_setup + """
rng = DatetimeIndex(start='1/1/2000', periods=10000, freq=datetools.Minute())
if rng.dtype == object:
    rng = rng.view(Index)
else:
    rng = rng.asobject
rng2 = rng[:-1]
"""

index_datetime_intersection = Benchmark("rng.intersection(rng2)", setup)
index_datetime_union = Benchmark("rng.union(rng2)", setup)

setup = common_setup + """
rng = date_range('1/1/2000', periods=10000, freq='T')
rng2 = rng[:-1]
"""

datetime_index_intersection = Benchmark("rng.intersection(rng2)", setup,
                                        start_date=datetime(2013, 9, 27))
datetime_index_union = Benchmark("rng.union(rng2)", setup,
                                 start_date=datetime(2013, 9, 27))

# integers
setup = common_setup + """
N = 1000000
options = np.arange(N)

left = Index(options.take(np.random.permutation(N)[:N // 2]))
right = Index(options.take(np.random.permutation(N)[:N // 2]))
"""

index_int64_union = Benchmark('left.union(right)', setup,
                              start_date=datetime(2011, 1, 1))

index_int64_intersection = Benchmark('left.intersection(right)', setup,
                                     start_date=datetime(2011, 1, 1))

#----------------------------------------------------------------------
# string index slicing
setup = common_setup + """
idx = tm.makeStringIndex(1000000)

mask = np.arange(1000000) % 3 == 0
series_mask = Series(mask)
"""
index_str_slice_indexer_basic = Benchmark('idx[:-1]', setup)
index_str_slice_indexer_even = Benchmark('idx[::2]', setup)
index_str_boolean_indexer = Benchmark('idx[mask]', setup)
index_str_boolean_series_indexer = Benchmark('idx[series_mask]', setup)

#----------------------------------------------------------------------
# float64 index
#----------------------------------------------------------------------
# construction
setup = common_setup + """
baseidx = np.arange(1e6)
"""

index_float64_construct = Benchmark('Index(baseidx)', setup,
                                    name='index_float64_construct',
                                    start_date=datetime(2014, 4, 13))

setup = common_setup + """
idx = tm.makeFloatIndex(1000000)

mask = np.arange(idx.size) % 3 == 0
series_mask = Series(mask)
"""
#----------------------------------------------------------------------
# getting
index_float64_get = Benchmark('idx[1]', setup, name='index_float64_get',
                              start_date=datetime(2014, 4, 13))


#----------------------------------------------------------------------
# slicing
index_float64_slice_indexer_basic = Benchmark('idx[:-1]', setup,
                                              name='index_float64_slice_indexer_basic',
                                              start_date=datetime(2014, 4, 13))
index_float64_slice_indexer_even = Benchmark('idx[::2]', setup,
                                             name='index_float64_slice_indexer_even',
                                             start_date=datetime(2014, 4, 13))
index_float64_boolean_indexer = Benchmark('idx[mask]', setup,
                                          name='index_float64_boolean_indexer',
                                          start_date=datetime(2014, 4, 13))
index_float64_boolean_series_indexer = Benchmark('idx[series_mask]', setup,
                                                 name='index_float64_boolean_series_indexer',
                                                 start_date=datetime(2014, 4, 13))

#----------------------------------------------------------------------
# arith ops
index_float64_mul = Benchmark('idx * 2', setup, name='index_float64_mul',
                              start_date=datetime(2014, 4, 13))
index_float64_div = Benchmark('idx / 2', setup, name='index_float64_div',
                              start_date=datetime(2014, 4, 13))
