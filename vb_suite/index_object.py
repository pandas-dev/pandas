from vbench.benchmark import Benchmark
from datetime import datetime

SECTION = "Index / MultiIndex objects"


common_setup = """from .pandas_vb_common import *
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


# Constructing MultiIndex from cartesian product of iterables
#

setup = common_setup + """
iterables = [tm.makeStringIndex(10000), range(20)]
"""

multiindex_from_product = Benchmark('MultiIndex.from_product(iterables)',
                                    setup, name='multiindex_from_product',
                                    start_date=datetime(2014, 6, 30))

#----------------------------------------------------------------------
# MultiIndex with DatetimeIndex level

setup = common_setup + """
level1 = range(1000)
level2 = date_range(start='1/1/2012', periods=100)
mi = MultiIndex.from_product([level1, level2])
"""

multiindex_with_datetime_level_full = \
    Benchmark("mi.copy().values", setup,
              name='multiindex_with_datetime_level_full',
              start_date=datetime(2014, 10, 11))


multiindex_with_datetime_level_sliced = \
    Benchmark("mi[:10].values", setup,
              name='multiindex_with_datetime_level_sliced',
              start_date=datetime(2014, 10, 11))

# multi-index duplicated
setup = common_setup + """
n, k = 200, 5000
levels = [np.arange(n), tm.makeStringIndex(n).values, 1000 + np.arange(n)]
labels = [np.random.choice(n, k * n) for lev in levels]
mi = MultiIndex(levels=levels, labels=labels)
"""

multiindex_duplicated = Benchmark('mi.duplicated()', setup,
                                  name='multiindex_duplicated')

#----------------------------------------------------------------------
# repr

setup = common_setup + """
dr = pd.date_range('20000101', freq='D', periods=100000)
"""

datetime_index_repr = \
    Benchmark("dr._is_dates_only", setup,
              start_date=datetime(2012, 1, 11))

setup = common_setup + """
n = 3 * 5 * 7 * 11 * (1 << 10)
low, high = - 1 << 12, 1 << 12
f = lambda k: np.repeat(np.random.randint(low, high, n // k), k)

i = np.random.permutation(n)
mi = MultiIndex.from_arrays([f(11), f(7), f(5), f(3), f(1)])[i]
"""

multiindex_sortlevel_int64 = Benchmark('mi.sortlevel()', setup,
                                       name='multiindex_sortlevel_int64')
