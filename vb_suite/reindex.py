from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from .pandas_vb_common import *
"""

#----------------------------------------------------------------------
# DataFrame reindex columns

setup = common_setup + """
df = DataFrame(index=range(10000), data=np.random.rand(10000,30),
               columns=range(30))
"""
statement = "df.reindex(columns=df.columns[1:5])"

frame_reindex_columns = Benchmark(statement, setup)

#----------------------------------------------------------------------

setup = common_setup + """
rng = DatetimeIndex(start='1/1/1970', periods=10000, freq=datetools.Minute())
df = DataFrame(np.random.rand(10000, 10), index=rng,
               columns=range(10))
df['foo'] = 'bar'
rng2 = Index(rng[::2])
"""
statement = "df.reindex(rng2)"
dataframe_reindex = Benchmark(statement, setup)

#----------------------------------------------------------------------
# multiindex reindexing

setup = common_setup + """
N = 1000
K = 20

level1 = tm.makeStringIndex(N).values.repeat(K)
level2 = np.tile(tm.makeStringIndex(K).values, N)
index = MultiIndex.from_arrays([level1, level2])

s1 = Series(np.random.randn(N * K), index=index)
s2 = s1[::2]
"""
statement = "s1.reindex(s2.index)"
reindex_multi = Benchmark(statement, setup,
                          name='reindex_multiindex',
                          start_date=datetime(2011, 9, 1))

#----------------------------------------------------------------------
# Pad / backfill

def pad(source_series, target_index):
    try:
        source_series.reindex(target_index, method='pad')
    except:
        source_series.reindex(target_index, fillMethod='pad')

def backfill(source_series, target_index):
    try:
        source_series.reindex(target_index, method='backfill')
    except:
        source_series.reindex(target_index, fillMethod='backfill')

setup = common_setup + """
rng = date_range('1/1/2000', periods=100000, freq=datetools.Minute())

ts = Series(np.random.randn(len(rng)), index=rng)
ts2 = ts[::2]
ts3 = ts2.reindex(ts.index)
ts4 = ts3.astype('float32')

def pad(source_series, target_index):
    try:
        source_series.reindex(target_index, method='pad')
    except:
        source_series.reindex(target_index, fillMethod='pad')
def backfill(source_series, target_index):
    try:
        source_series.reindex(target_index, method='backfill')
    except:
        source_series.reindex(target_index, fillMethod='backfill')
"""

statement = "pad(ts2, ts.index)"
reindex_daterange_pad = Benchmark(statement, setup,
                                  name="reindex_daterange_pad")

statement = "backfill(ts2, ts.index)"
reindex_daterange_backfill = Benchmark(statement, setup,
                                       name="reindex_daterange_backfill")

reindex_fillna_pad = Benchmark("ts3.fillna(method='pad')", setup,
                               name="reindex_fillna_pad",
                               start_date=datetime(2011, 3, 1))

reindex_fillna_pad_float32 = Benchmark("ts4.fillna(method='pad')", setup,
                                       name="reindex_fillna_pad_float32",
                                       start_date=datetime(2013, 1, 1))

reindex_fillna_backfill = Benchmark("ts3.fillna(method='backfill')", setup,
                                    name="reindex_fillna_backfill",
                                    start_date=datetime(2011, 3, 1))
reindex_fillna_backfill_float32 = Benchmark("ts4.fillna(method='backfill')", setup,
                                            name="reindex_fillna_backfill_float32",
                                            start_date=datetime(2013, 1, 1))

#----------------------------------------------------------------------
# align on level

setup = common_setup + """
index = MultiIndex(levels=[np.arange(10), np.arange(100), np.arange(100)],
                   labels=[np.arange(10).repeat(10000),
                           np.tile(np.arange(100).repeat(100), 10),
                           np.tile(np.tile(np.arange(100), 100), 10)])
random.shuffle(index.values)
df = DataFrame(np.random.randn(len(index), 4), index=index)
df_level = DataFrame(np.random.randn(100, 4), index=index.levels[1])
"""

reindex_frame_level_align = \
    Benchmark("df.align(df_level, level=1, copy=False)", setup,
              name='reindex_frame_level_align',
              start_date=datetime(2011, 12, 27))

reindex_frame_level_reindex = \
    Benchmark("df_level.reindex(df.index, level=1)", setup,
              name='reindex_frame_level_reindex',
              start_date=datetime(2011, 12, 27))


#----------------------------------------------------------------------
# sort_index, drop_duplicates

# pathological, but realistic
setup = common_setup + """
N = 10000
K = 10

key1 = tm.makeStringIndex(N).values.repeat(K)
key2 = tm.makeStringIndex(N).values.repeat(K)

df = DataFrame({'key1' : key1, 'key2' : key2,
                'value' : np.random.randn(N * K)})
col_array_list = list(df.values.T)
"""
statement = "df.sort_index(by=['key1', 'key2'])"
frame_sort_index_by_columns = Benchmark(statement, setup,
                                        start_date=datetime(2011, 11, 1))

# drop_duplicates

statement = "df.drop_duplicates(['key1', 'key2'])"
frame_drop_duplicates = Benchmark(statement, setup,
                                  start_date=datetime(2011, 11, 15))

statement = "df.drop_duplicates(['key1', 'key2'], inplace=True)"
frame_drop_dup_inplace = Benchmark(statement, setup,
                                   start_date=datetime(2012, 5, 16))

lib_fast_zip = Benchmark('lib.fast_zip(col_array_list)', setup,
                         name='lib_fast_zip',
                         start_date=datetime(2012, 1, 1))

setup = setup + """
df.ix[:10000, :] = np.nan
"""
statement2 = "df.drop_duplicates(['key1', 'key2'])"
frame_drop_duplicates_na = Benchmark(statement2, setup,
                                     start_date=datetime(2012, 5, 15))

lib_fast_zip_fillna = Benchmark('lib.fast_zip_fillna(col_array_list)', setup,
                                start_date=datetime(2012, 5, 15))

statement2 = "df.drop_duplicates(['key1', 'key2'], inplace=True)"
frame_drop_dup_na_inplace = Benchmark(statement2, setup,
                                      start_date=datetime(2012, 5, 16))

setup = common_setup + """
s = Series(np.random.randint(0, 1000, size=10000))
s2 = Series(np.tile(tm.makeStringIndex(1000).values, 10))
"""

series_drop_duplicates_int = Benchmark('s.drop_duplicates()', setup,
                                       start_date=datetime(2012, 11, 27))

series_drop_duplicates_string = \
    Benchmark('s2.drop_duplicates()', setup,
              start_date=datetime(2012, 11, 27))

#----------------------------------------------------------------------
# fillna, many columns


setup = common_setup + """
values = np.random.randn(1000, 1000)
values[::2] = np.nan
df = DataFrame(values)
"""

frame_fillna_many_columns_pad = Benchmark("df.fillna(method='pad')",
                                          setup,
                                          start_date=datetime(2011, 3, 1))

#----------------------------------------------------------------------
# blog "pandas escaped the zoo"

setup = common_setup + """
n = 50000
indices = tm.makeStringIndex(n)

def sample(values, k):
    from random import shuffle
    sampler = np.arange(len(values))
    shuffle(sampler)
    return values.take(sampler[:k])

subsample_size = 40000

x = Series(np.random.randn(50000), indices)
y = Series(np.random.randn(subsample_size),
           index=sample(indices, subsample_size))
"""

series_align_irregular_string = Benchmark("x + y", setup,
                                          start_date=datetime(2010, 6, 1))
