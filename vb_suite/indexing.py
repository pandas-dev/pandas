from vbench.benchmark import Benchmark
from datetime import datetime

SECTION = 'Indexing and scalar value access'

common_setup = """from pandas_vb_common import *
"""

#----------------------------------------------------------------------
# Series.__getitem__, get_value

setup = common_setup + """
tm.N = 1000
ts = tm.makeTimeSeries()
dt = ts.index[500]
"""
statement = "ts[dt]"

bm_getitem = Benchmark(statement, setup, ncalls=100000,
                       name='series_getitem_scalar')

setup = common_setup + """
index = [tm.rands(10) for _ in xrange(1000)]
s = Series(np.random.rand(1000), index=index)
idx = index[100]
"""
statement = "s.get_value(idx)"
bm_df_getitem3 = Benchmark(statement, setup,
                           name='series_get_value',
                           start_date=datetime(2011, 11, 12))

#----------------------------------------------------------------------
# DataFrame __getitem__

setup = common_setup + """
index = [tm.rands(10) for _ in xrange(1000)]
columns = [tm.rands(10) for _ in xrange(30)]
df = DataFrame(np.random.rand(1000, 30), index=index,
               columns=columns)
idx = index[100]
col = columns[10]
"""
statement = "df[col][idx]"
bm_df_getitem = Benchmark(statement, setup,
                        name='dataframe_getitem_scalar')

setup = common_setup + """
try:
    klass = DataMatrix
except:
    klass = DataFrame

index = [tm.rands(10) for _ in xrange(1000)]
columns = [tm.rands(10) for _ in xrange(30)]
df = klass(np.random.rand(1000, 30), index=index,
               columns=columns)
idx = index[100]
col = columns[10]
"""
statement = "df[col][idx]"
bm_df_getitem2 = Benchmark(statement, setup,
                        name='datamatrix_getitem_scalar')

setup = common_setup + """
try:
    klass = DataMatrix
except:
    klass = DataFrame

index = [tm.rands(10) for _ in xrange(1000)]
columns = [tm.rands(10) for _ in xrange(30)]
df = klass(np.random.rand(1000, 30), index=index,
               columns=columns)
idx = index[100]
col = columns[10]
"""
statement = "df.get_value(idx, col)"
bm_df_getitem3 = Benchmark(statement, setup,
                           name='dataframe_get_value',
                           start_date=datetime(2011, 11, 12))

#----------------------------------------------------------------------
# Boolean DataFrame row selection

setup = common_setup + """
df = DataFrame(np.random.randn(10000, 4), columns=['A', 'B', 'C', 'D'])
indexer = df['B'] > 0
obj_indexer = indexer.astype('O')
"""
indexing_dataframe_boolean_rows = \
    Benchmark("df[indexer]", setup, name='indexing_dataframe_boolean_rows')

indexing_dataframe_boolean_rows_object = \
    Benchmark("df[obj_indexer]", setup,
              name='indexing_dataframe_boolean_rows_object')

#----------------------------------------------------------------------
# MultiIndex sortlevel

setup = common_setup + """
a = np.repeat(np.arange(100), 1000)
b = np.tile(np.arange(1000), 100)
midx = MultiIndex.from_arrays([a, b])
midx = midx.take(np.random.permutation(np.arange(100000)))
"""
sort_level_zero = Benchmark("midx.sortlevel(0)", setup,
                            start_date=datetime(2012,1,1))
sort_level_one = Benchmark("midx.sortlevel(1)", setup,
                           start_date=datetime(2012,1,1))

#----------------------------------------------------------------------
# Panel subset selection

setup = common_setup + """
p = Panel(np.random.randn(500, 500, 500))
inds = range(0, 500, 10)
"""

indexing_panel_subset = Benchmark('p.ix[inds, inds, inds]', setup,
                                  start_date=datetime(2012, 1, 1))
