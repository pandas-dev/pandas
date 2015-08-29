from vbench.benchmark import Benchmark
from datetime import datetime

SECTION = 'Indexing and scalar value access'

common_setup = """from .pandas_vb_common import *
"""

#----------------------------------------------------------------------
# Series.__getitem__, get_value, __getitem__(slice)

setup = common_setup + """
tm.N = 1000
ts = tm.makeTimeSeries()
dt = ts.index[500]
"""
statement = "ts[dt]"
bm_getitem = Benchmark(statement, setup, ncalls=100000,
                       name='time_series_getitem_scalar')

setup = common_setup + """
index = tm.makeStringIndex(1000)
s = Series(np.random.rand(1000), index=index)
idx = index[100]
"""
statement = "s.get_value(idx)"
bm_get_value = Benchmark(statement, setup,
                         name='series_get_value',
                         start_date=datetime(2011, 11, 12))


setup = common_setup + """
index = tm.makeStringIndex(1000000)
s = Series(np.random.rand(1000000), index=index)
"""
series_getitem_pos_slice = Benchmark("s[:800000]", setup,
                                     name="series_getitem_pos_slice")


setup = common_setup + """
index = tm.makeStringIndex(1000000)
s = Series(np.random.rand(1000000), index=index)
lbl = s.index[800000]
"""
series_getitem_label_slice = Benchmark("s[:lbl]", setup,
                                       name="series_getitem_label_slice")


#----------------------------------------------------------------------
# DataFrame __getitem__

setup = common_setup + """
index = tm.makeStringIndex(1000)
columns = tm.makeStringIndex(30)
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

index = tm.makeStringIndex(1000)
columns = tm.makeStringIndex(30)
df = klass(np.random.rand(1000, 30), index=index, columns=columns)
idx = index[100]
col = columns[10]
"""
statement = "df[col][idx]"
bm_df_getitem2 = Benchmark(statement, setup,
                           name='datamatrix_getitem_scalar')


#----------------------------------------------------------------------
# ix get scalar

setup = common_setup + """
index = tm.makeStringIndex(1000)
columns = tm.makeStringIndex(30)
df = DataFrame(np.random.randn(1000, 30), index=index, columns=columns)
idx = index[100]
col = columns[10]
"""

indexing_frame_get_value_ix = Benchmark("df.ix[idx,col]", setup,
                                        name='indexing_frame_get_value_ix',
                                        start_date=datetime(2011, 11, 12))

indexing_frame_get_value = Benchmark("df.get_value(idx,col)", setup,
                                     name='indexing_frame_get_value',
                                     start_date=datetime(2011, 11, 12))

setup = common_setup + """
mi = MultiIndex.from_tuples([(x,y) for x in range(1000) for y in range(1000)])
s =  Series(np.random.randn(1000000), index=mi)
"""

series_xs_mi_ix = Benchmark("s.ix[999]", setup,
                            name='series_xs_mi_ix',
                            start_date=datetime(2013, 1, 1))

setup = common_setup + """
mi = MultiIndex.from_tuples([(x,y) for x in range(1000) for y in range(1000)])
s =  Series(np.random.randn(1000000), index=mi)
df = DataFrame(s)
"""

frame_xs_mi_ix = Benchmark("df.ix[999]", setup,
                           name='frame_xs_mi_ix',
                           start_date=datetime(2013, 1, 1))

#----------------------------------------------------------------------
# Boolean DataFrame row selection

setup = common_setup + """
df  = DataFrame(np.random.randn(10000, 4), columns=['A', 'B', 'C', 'D'])
indexer = df['B'] > 0
obj_indexer = indexer.astype('O')
"""
indexing_dataframe_boolean_rows = \
    Benchmark("df[indexer]", setup, name='indexing_dataframe_boolean_rows')

indexing_dataframe_boolean_rows_object = \
    Benchmark("df[obj_indexer]", setup,
              name='indexing_dataframe_boolean_rows_object')

setup = common_setup + """
df  = DataFrame(np.random.randn(50000, 100))
df2 = DataFrame(np.random.randn(50000, 100))
"""
indexing_dataframe_boolean = \
    Benchmark("df > df2", setup, name='indexing_dataframe_boolean',
              start_date=datetime(2012, 1, 1))

setup = common_setup + """
try:
    import pandas.computation.expressions as expr
except:
    expr = None

if expr is None:
    raise NotImplementedError
df  = DataFrame(np.random.randn(50000, 100))
df2 = DataFrame(np.random.randn(50000, 100))
expr.set_numexpr_threads(1)
"""

indexing_dataframe_boolean_st = \
    Benchmark("df > df2", setup, name='indexing_dataframe_boolean_st',cleanup="expr.set_numexpr_threads()",
              start_date=datetime(2013, 2, 26))


setup = common_setup + """
try:
    import pandas.computation.expressions as expr
except:
    expr = None

if expr is None:
    raise NotImplementedError
df  = DataFrame(np.random.randn(50000, 100))
df2 = DataFrame(np.random.randn(50000, 100))
expr.set_use_numexpr(False)
"""

indexing_dataframe_boolean_no_ne = \
    Benchmark("df > df2", setup, name='indexing_dataframe_boolean_no_ne',cleanup="expr.set_use_numexpr(True)",
              start_date=datetime(2013, 2, 26))
#----------------------------------------------------------------------
# MultiIndex sortlevel

setup = common_setup + """
a = np.repeat(np.arange(100), 1000)
b = np.tile(np.arange(1000), 100)
midx = MultiIndex.from_arrays([a, b])
midx = midx.take(np.random.permutation(np.arange(100000)))
"""
sort_level_zero = Benchmark("midx.sortlevel(0)", setup,
                            start_date=datetime(2012, 1, 1))
sort_level_one = Benchmark("midx.sortlevel(1)", setup,
                           start_date=datetime(2012, 1, 1))

#----------------------------------------------------------------------
# Panel subset selection

setup = common_setup + """
p = Panel(np.random.randn(100, 100, 100))
inds = range(0, 100, 10)
"""

indexing_panel_subset = Benchmark('p.ix[inds, inds, inds]', setup,
                                  start_date=datetime(2012, 1, 1))

#----------------------------------------------------------------------
# Iloc

setup = common_setup + """
df = DataFrame({'A' : [0.1] * 3000, 'B' : [1] * 3000})
idx = np.array(range(30)) * 99
df2 = DataFrame({'A' : [0.1] * 1000, 'B' : [1] * 1000})
df2 = concat([df2, 2*df2, 3*df2])
"""

frame_iloc_dups = Benchmark('df2.iloc[idx]', setup,
                            start_date=datetime(2013, 1, 1))

frame_loc_dups = Benchmark('df2.loc[idx]', setup,
                            start_date=datetime(2013, 1, 1))

setup = common_setup + """
df = DataFrame(dict( A = [ 'foo'] * 1000000))
"""

frame_iloc_big = Benchmark('df.iloc[:100,0]', setup,
                            start_date=datetime(2013, 1, 1))

#----------------------------------------------------------------------
# basic tests for [], .loc[], .iloc[] and .ix[]

setup = common_setup + """
s = Series(np.random.rand(1000000))
"""

series_getitem_scalar = Benchmark("s[800000]", setup)
series_getitem_slice = Benchmark("s[:800000]", setup)
series_getitem_list_like = Benchmark("s[[800000]]", setup)
series_getitem_array = Benchmark("s[np.arange(10000)]", setup)

series_loc_scalar = Benchmark("s.loc[800000]", setup)
series_loc_slice = Benchmark("s.loc[:800000]", setup)
series_loc_list_like = Benchmark("s.loc[[800000]]", setup)
series_loc_array = Benchmark("s.loc[np.arange(10000)]", setup)

series_iloc_scalar = Benchmark("s.iloc[800000]", setup)
series_iloc_slice = Benchmark("s.iloc[:800000]", setup)
series_iloc_list_like = Benchmark("s.iloc[[800000]]", setup)
series_iloc_array = Benchmark("s.iloc[np.arange(10000)]", setup)

series_ix_scalar = Benchmark("s.ix[800000]", setup)
series_ix_slice = Benchmark("s.ix[:800000]", setup)
series_ix_list_like = Benchmark("s.ix[[800000]]", setup)
series_ix_array = Benchmark("s.ix[np.arange(10000)]", setup)


# multi-index slicing
setup = common_setup + """
np.random.seed(1234)
idx=pd.IndexSlice
n=100000
mdt = pandas.DataFrame()
mdt['A'] = np.random.choice(range(10000,45000,1000), n)
mdt['B'] = np.random.choice(range(10,400), n)
mdt['C'] = np.random.choice(range(1,150), n)
mdt['D'] = np.random.choice(range(10000,45000), n)
mdt['x'] = np.random.choice(range(400), n)
mdt['y'] = np.random.choice(range(25), n)


test_A = 25000
test_B = 25
test_C = 40
test_D = 35000

eps_A = 5000
eps_B = 5
eps_C = 5
eps_D = 5000
mdt2 = mdt.set_index(['A','B','C','D']).sortlevel()
"""

multiindex_slicers = Benchmark('mdt2.loc[idx[test_A-eps_A:test_A+eps_A,test_B-eps_B:test_B+eps_B,test_C-eps_C:test_C+eps_C,test_D-eps_D:test_D+eps_D],:]', setup,
                               start_date=datetime(2015, 1, 1))

#----------------------------------------------------------------------
# take

setup = common_setup + """
s = Series(np.random.rand(100000))
ts = Series(np.random.rand(100000),
            index=date_range('2011-01-01', freq='S', periods=100000))
indexer = [True, False, True, True, False] * 20000
"""

series_take_intindex = Benchmark("s.take(indexer)", setup)
series_take_dtindex = Benchmark("ts.take(indexer)", setup)
