from vbench.api import Benchmark
from datetime import datetime

common_setup = """from .pandas_vb_common import *
"""

#----------------------------------------------------------------------
# lookup

setup = common_setup + """
df = DataFrame(np.random.randn(10000, 8), columns=list('abcdefgh'))
df['foo'] = 'bar'

row_labels = list(df.index[::10])[:900]
col_labels = list(df.columns) * 100
row_labels_all = np.array(list(df.index) * len(df.columns), dtype='object')
col_labels_all = np.array(list(df.columns) * len(df.index), dtype='object')
"""

frame_fancy_lookup = Benchmark('df.lookup(row_labels, col_labels)', setup,
                               start_date=datetime(2012, 1, 12))

frame_fancy_lookup_all = Benchmark('df.lookup(row_labels_all, col_labels_all)',
                                   setup,
                                   start_date=datetime(2012, 1, 12))

#----------------------------------------------------------------------
# fillna in place

setup = common_setup + """
df = DataFrame(randn(10000, 100))
df.values[::2] = np.nan
"""

frame_fillna_inplace = Benchmark('df.fillna(0, inplace=True)', setup,
                                 start_date=datetime(2012, 4, 4))


#----------------------------------------------------------------------
# reindex both axes

setup = common_setup + """
df = DataFrame(randn(10000, 10000))
idx = np.arange(4000, 7000)
"""

frame_reindex_axis0 = Benchmark('df.reindex(idx)', setup)

frame_reindex_axis1 = Benchmark('df.reindex(columns=idx)', setup)

frame_reindex_both_axes = Benchmark('df.reindex(index=idx, columns=idx)',
                                    setup, start_date=datetime(2011, 1, 1))

frame_reindex_both_axes_ix = Benchmark('df.ix[idx, idx]', setup,
                                       start_date=datetime(2011, 1, 1))

#----------------------------------------------------------------------
# reindex with upcasts
setup = common_setup + """
df=DataFrame(dict([(c, {
        0: randint(0, 2, 1000).astype(np.bool_),
        1: randint(0, 1000, 1000).astype(np.int16),
        2: randint(0, 1000, 1000).astype(np.int32),
        3: randint(0, 1000, 1000).astype(np.int64)
    }[randint(0, 4)]) for c in range(1000)]))
"""

frame_reindex_upcast = Benchmark('df.reindex(permutation(range(1200)))', setup)

#----------------------------------------------------------------------
# boolean indexing

setup = common_setup + """
df = DataFrame(randn(10000, 100))
bool_arr = np.zeros(10000, dtype=bool)
bool_arr[:1000] = True
"""

frame_boolean_row_select = Benchmark('df[bool_arr]', setup,
                                     start_date=datetime(2011, 1, 1))

#----------------------------------------------------------------------
# iteritems (monitor no-copying behaviour)

setup = common_setup + """
df = DataFrame(randn(10000, 1000))
df2 = DataFrame(randn(3000,1),columns=['A'])
df3 = DataFrame(randn(3000,1))

def f():
    if hasattr(df, '_item_cache'):
        df._item_cache.clear()
    for name, col in df.iteritems():
        pass

def g():
    for name, col in df.iteritems():
        pass

def h():
    for i in range(10000):
        df2['A']

def j():
    for i in range(10000):
        df3[0]

"""

# as far back as the earliest test currently in the suite
frame_iteritems = Benchmark('f()', setup,
                            start_date=datetime(2010, 6, 1))

frame_iteritems_cached = Benchmark('g()', setup,
                                   start_date=datetime(2010, 6, 1))

frame_getitem_single_column = Benchmark('h()', setup,
                                        start_date=datetime(2010, 6, 1))

frame_getitem_single_column2 = Benchmark('j()', setup,
                                         start_date=datetime(2010, 6, 1))

#----------------------------------------------------------------------
# assignment

setup = common_setup + """
idx = date_range('1/1/2000', periods=100000, freq='D')
df = DataFrame(randn(100000, 1),columns=['A'],index=idx)
def f(df):
    x = df.copy()
    x['date'] = x.index
"""

frame_assign_timeseries_index = Benchmark('f(df)', setup,
                                          start_date=datetime(2013, 10, 1))


#----------------------------------------------------------------------
# to_string

setup = common_setup + """
df = DataFrame(randn(100, 10))
"""

frame_to_string_floats = Benchmark('df.to_string()', setup,
                                   start_date=datetime(2010, 6, 1))

#----------------------------------------------------------------------
# to_html

setup = common_setup + """
nrows=500
df = DataFrame(randn(nrows, 10))
df[0]=period_range("2000","2010",nrows)
df[1]=range(nrows)

"""

frame_to_html_mixed = Benchmark('df.to_html()', setup,
                                   start_date=datetime(2011, 11, 18))


# truncated repr_html, single index

setup = common_setup + """
nrows=10000
data=randn(nrows,10)
idx=MultiIndex.from_arrays(np.tile(randn(3,nrows/100),100))
df=DataFrame(data,index=idx)

"""

frame_html_repr_trunc_mi = Benchmark('df._repr_html_()', setup,
                                   start_date=datetime(2013, 11, 25))

# truncated repr_html, MultiIndex

setup = common_setup + """
nrows=10000
data=randn(nrows,10)
idx=randn(nrows)
df=DataFrame(data,index=idx)

"""

frame_html_repr_trunc_si = Benchmark('df._repr_html_()', setup,
                                   start_date=datetime(2013, 11, 25))


# insert many columns

setup = common_setup + """
N = 1000

def f(K=500):
    df = DataFrame(index=range(N))
    new_col = np.random.randn(N)
    for i in range(K):
        df[i] = new_col
"""

frame_insert_500_columns_end = Benchmark('f()', setup, start_date=datetime(2011, 1, 1))

setup = common_setup + """
N = 1000

def f(K=100):
    df = DataFrame(index=range(N))
    new_col = np.random.randn(N)
    for i in range(K):
        df.insert(0,i,new_col)
"""

frame_insert_100_columns_begin = Benchmark('f()', setup, start_date=datetime(2011, 1, 1))

#----------------------------------------------------------------------
# strings methods, #2602

setup = common_setup + """
s = Series(['abcdefg', np.nan]*500000)
"""

series_string_vector_slice = Benchmark('s.str[:5]', setup,
                                       start_date=datetime(2012, 8, 1))

#----------------------------------------------------------------------
# df.info() and get_dtype_counts() # 2807

setup = common_setup + """
df = pandas.DataFrame(np.random.randn(10,10000))
"""

frame_get_dtype_counts = Benchmark('df.get_dtype_counts()', setup,
                                       start_date=datetime(2012, 8, 1))

##
setup = common_setup + """
df = pandas.DataFrame(np.random.randn(10,10000))
"""

frame_repr_wide = Benchmark('repr(df)', setup,
                            start_date=datetime(2012, 8, 1))

##
setup = common_setup + """
df = pandas.DataFrame(np.random.randn(10000, 10))
"""

frame_repr_tall = Benchmark('repr(df)', setup,
                            start_date=datetime(2012, 8, 1))

##
setup = common_setup + """
df = DataFrame(randn(100000, 1))
"""

frame_xs_row = Benchmark('df.xs(50000)', setup)

##
setup = common_setup + """
df = DataFrame(randn(1,100000))
"""

frame_xs_col = Benchmark('df.xs(50000,axis = 1)', setup)

#----------------------------------------------------------------------
# nulls/masking

## masking
setup = common_setup + """
data = np.random.randn(1000, 500)
df = DataFrame(data)
df = df.where(df > 0) # create nans
bools = df > 0
mask = isnull(df)
"""

frame_mask_bools = Benchmark('bools.mask(mask)', setup,
                             start_date=datetime(2013,1,1))

frame_mask_floats  = Benchmark('bools.astype(float).mask(mask)', setup,
                             start_date=datetime(2013,1,1))

## isnull
setup = common_setup + """
data = np.random.randn(1000, 1000)
df = DataFrame(data)
"""
frame_isnull  = Benchmark('isnull(df)', setup,
                           start_date=datetime(2012,1,1))

## dropna
dropna_setup = common_setup + """
data = np.random.randn(10000, 1000)
df = DataFrame(data)
df.ix[50:1000,20:50] = np.nan
df.ix[2000:3000] = np.nan
df.ix[:,60:70] = np.nan
"""
frame_dropna_axis0_any  = Benchmark('df.dropna(how="any",axis=0)', dropna_setup,
                                     start_date=datetime(2012,1,1))
frame_dropna_axis0_all  = Benchmark('df.dropna(how="all",axis=0)', dropna_setup,
                                     start_date=datetime(2012,1,1))

frame_dropna_axis1_any  = Benchmark('df.dropna(how="any",axis=1)', dropna_setup,
                                    start_date=datetime(2012,1,1))

frame_dropna_axis1_all  = Benchmark('df.dropna(how="all",axis=1)', dropna_setup,
                                    start_date=datetime(2012,1,1))

# dropna on mixed dtypes
dropna_mixed_setup = common_setup + """
data = np.random.randn(10000, 1000)
df = DataFrame(data)
df.ix[50:1000,20:50] = np.nan
df.ix[2000:3000] = np.nan
df.ix[:,60:70] = np.nan
df['foo'] = 'bar'
"""
frame_dropna_axis0_any_mixed_dtypes  = Benchmark('df.dropna(how="any",axis=0)', dropna_mixed_setup,
                                                 start_date=datetime(2012,1,1))
frame_dropna_axis0_all_mixed_dtypes  = Benchmark('df.dropna(how="all",axis=0)', dropna_mixed_setup,
                                                 start_date=datetime(2012,1,1))

frame_dropna_axis1_any_mixed_dtypes  = Benchmark('df.dropna(how="any",axis=1)', dropna_mixed_setup,
                                                 start_date=datetime(2012,1,1))

frame_dropna_axis1_all_mixed_dtypes  = Benchmark('df.dropna(how="all",axis=1)', dropna_mixed_setup,
                                                 start_date=datetime(2012,1,1))

## dropna multi
dropna_setup = common_setup + """
data = np.random.randn(10000, 1000)
df = DataFrame(data)
df.ix[50:1000,20:50] = np.nan
df.ix[2000:3000] = np.nan
df.ix[:,60:70] = np.nan
df.index = MultiIndex.from_tuples(df.index.map(lambda x: (x, x)))
df.columns = MultiIndex.from_tuples(df.columns.map(lambda x: (x, x)))
"""
frame_count_level_axis0_multi = Benchmark('df.count(axis=0, level=1)', dropna_setup,
                                          start_date=datetime(2012,1,1))

frame_count_level_axis1_multi = Benchmark('df.count(axis=1, level=1)', dropna_setup,
                                          start_date=datetime(2012,1,1))

# dropna on mixed dtypes
dropna_mixed_setup = common_setup + """
data = np.random.randn(10000, 1000)
df = DataFrame(data)
df.ix[50:1000,20:50] = np.nan
df.ix[2000:3000] = np.nan
df.ix[:,60:70] = np.nan
df['foo'] = 'bar'
df.index = MultiIndex.from_tuples(df.index.map(lambda x: (x, x)))
df.columns = MultiIndex.from_tuples(df.columns.map(lambda x: (x, x)))
"""
frame_count_level_axis0_mixed_dtypes_multi  = Benchmark('df.count(axis=0, level=1)', dropna_mixed_setup,
                                                        start_date=datetime(2012,1,1))

frame_count_level_axis1_mixed_dtypes_multi  = Benchmark('df.count(axis=1, level=1)', dropna_mixed_setup,
                                                        start_date=datetime(2012,1,1))

#----------------------------------------------------------------------
# apply

setup = common_setup + """
s = Series(np.arange(1028.))
df = DataFrame({ i:s for i in range(1028) })
"""
frame_apply_user_func = Benchmark('df.apply(lambda x: np.corrcoef(x,s)[0,1])', setup,
                           name = 'frame_apply_user_func',
                           start_date=datetime(2012,1,1))

setup = common_setup + """
df = DataFrame(np.random.randn(1000,100))
"""
frame_apply_lambda_mean = Benchmark('df.apply(lambda x: x.sum())', setup,
                                    name = 'frame_apply_lambda_mean',
                                    start_date=datetime(2012,1,1))
setup = common_setup + """
df = DataFrame(np.random.randn(1000,100))
"""
frame_apply_np_mean = Benchmark('df.apply(np.mean)', setup,
                               name = 'frame_apply_np_mean',
                               start_date=datetime(2012,1,1))

setup = common_setup + """
df = DataFrame(np.random.randn(1000,100))
"""
frame_apply_pass_thru = Benchmark('df.apply(lambda x: x)', setup,
                                  name = 'frame_apply_pass_thru',
                                  start_date=datetime(2012,1,1))

setup = common_setup + """
df = DataFrame(np.random.randn(1000,100))
"""
frame_apply_axis_1 = Benchmark('df.apply(lambda x: x+1,axis=1)', setup,
                               name = 'frame_apply_axis_1',
                               start_date=datetime(2012,1,1))

setup = common_setup + """
df = DataFrame(np.random.randn(1000,3),columns=list('ABC'))
"""
frame_apply_ref_by_name = Benchmark('df.apply(lambda x: x["A"] + x["B"],axis=1)', setup,
                                     name = 'frame_apply_ref_by_name',
                                     start_date=datetime(2012,1,1))

#----------------------------------------------------------------------
# dtypes

setup = common_setup + """
df = DataFrame(np.random.randn(1000,1000))
"""
frame_dtypes = Benchmark('df.dtypes', setup,
                         start_date=datetime(2012,1,1))

#----------------------------------------------------------------------
# equals
setup = common_setup + """
def make_pair(frame):
    df = frame
    df2 = df.copy()
    df2.ix[-1,-1] = np.nan
    return df, df2

def test_equal(name):
    df, df2 = pairs[name]
    return df.equals(df)

def test_unequal(name):
    df, df2 = pairs[name]
    return df.equals(df2)

float_df = DataFrame(np.random.randn(1000, 1000))
object_df = DataFrame([['foo']*1000]*1000)
nonunique_cols = object_df.copy()
nonunique_cols.columns = ['A']*len(nonunique_cols.columns)

pairs = dict([(name, make_pair(frame))
         for name, frame in (('float_df', float_df), ('object_df', object_df), ('nonunique_cols', nonunique_cols))])
"""
frame_float_equal = Benchmark('test_equal("float_df")', setup)
frame_object_equal = Benchmark('test_equal("object_df")', setup)
frame_nonunique_equal = Benchmark('test_equal("nonunique_cols")', setup)

frame_float_unequal = Benchmark('test_unequal("float_df")', setup)
frame_object_unequal = Benchmark('test_unequal("object_df")', setup)
frame_nonunique_unequal = Benchmark('test_unequal("nonunique_cols")', setup)

#-----------------------------------------------------------------------------
# interpolate
# this is the worst case, where every column has NaNs.
setup = common_setup + """
df = DataFrame(randn(10000, 100))
df.values[::2] = np.nan
"""

frame_interpolate = Benchmark('df.interpolate()', setup,
                               start_date=datetime(2014, 2, 7))

setup = common_setup + """
df = DataFrame({'A': np.arange(0, 10000),
                'B': np.random.randint(0, 100, 10000),
                'C': randn(10000),
                'D': randn(10000)})
df.loc[1::5, 'A'] = np.nan
df.loc[1::5, 'C'] = np.nan
"""

frame_interpolate_some_good = Benchmark('df.interpolate()', setup,
                                        start_date=datetime(2014, 2, 7))
frame_interpolate_some_good_infer = Benchmark('df.interpolate(downcast="infer")',
                                              setup,
                                              start_date=datetime(2014, 2, 7))


#-------------------------------------------------------------------------
# frame shift speedup issue-5609

setup = common_setup + """
df = DataFrame(np.random.rand(10000,500))
# note: df._data.blocks are f_contigous
"""
frame_shift_axis0 = Benchmark('df.shift(1,axis=0)', setup,
                    start_date=datetime(2014,1,1))
frame_shift_axis1 = Benchmark('df.shift(1,axis=1)', setup,
                    name = 'frame_shift_axis_1',
                    start_date=datetime(2014,1,1))


#-----------------------------------------------------------------------------
# from_records issue-6700

setup = common_setup + """
def get_data(n=100000):
    return ((x, x*20, x*100) for x in range(n))
"""

frame_from_records_generator = Benchmark('df = DataFrame.from_records(get_data())',
                                setup,
                                name='frame_from_records_generator',
                                start_date=datetime(2013,10,4))  # issue-4911

frame_from_records_generator_nrows = Benchmark('df = DataFrame.from_records(get_data(), nrows=1000)',
                                setup,
                                name='frame_from_records_generator_nrows',
                                start_date=datetime(2013,10,04))  # issue-4911

#-----------------------------------------------------------------------------
# duplicated

setup = common_setup + '''
n = 1 << 20

t = date_range('2015-01-01', freq='S', periods=n // 64)
xs = np.random.randn(n // 64).round(2)

df = DataFrame({'a':np.random.randint(- 1 << 8, 1 << 8, n),
                'b':np.random.choice(t, n),
                'c':np.random.choice(xs, n)})
'''

frame_duplicated = Benchmark('df.duplicated()', setup,
                             name='frame_duplicated')
