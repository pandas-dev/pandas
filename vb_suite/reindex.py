from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

#----------------------------------------------------------------------
# DataFrame reindex columns

setup = common_setup + """
df = DataFrame(index=range(10000), data=np.random.rand(10000,30),
               columns=range(30))
"""
statement = "df.reindex(columns=df.columns[1:5])"

reindex_frame_columns = Benchmark(statement, setup,
                                  name='dataframe_reindex_columns')

#----------------------------------------------------------------------

setup = common_setup + """
rng = DateRange('1/1/1970', periods=10000, offset=datetools.Minute())
df = DataFrame(np.random.rand(10000, 10), index=rng,
               columns=range(10))
df['foo'] = 'bar'
rng2 = Index(rng[::2])
"""
statement = "df.reindex(rng2)"
reindex_frame_daterange = Benchmark(statement, setup,
                                    name='dataframe_reindex_daterange')

#----------------------------------------------------------------------
# multiindex reindexing

setup = common_setup + """
N = 1000
K = 20

level1 = np.array([tm.rands(10) for _ in xrange(N)], dtype='O').repeat(K)
level2 = np.tile(np.array([tm.rands(10) for _ in xrange(K)], dtype='O'),
                 N)
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

setup = common_setup + """
rng = DateRange('1/1/2000', periods=10000, offset=datetools.Minute())

ts = Series(np.random.randn(len(rng)), index=rng)
ts2 = ts[::2]
ts3 = ts2.reindex(ts.index)

def pad():
    try:
        ts2.reindex(ts.index, method='pad')
    except:
        ts2.reindex(ts.index, fillMethod='pad')
def backfill():
    try:
        ts2.reindex(ts.index, method='backfill')
    except:
        ts2.reindex(ts.index, fillMethod='backfill')
"""

statement = "pad()"
reindex_daterange_pad = Benchmark(statement, setup,
                                  name="reindex_daterange_pad")

statement = "backfill()"
reindex_daterange_backfill = Benchmark(statement, setup,
                                       name="reindex_daterange_backfill")

reindex_fillna_pad = Benchmark("ts3.fillna(method='pad')", setup,
                               name="reindex_fillna_pad",
                               start_date=datetime(2011, 3, 1))

reindex_fillna_backfill = Benchmark("ts3.fillna(method='backfill')", setup,
                                    name="reindex_fillna_backfill",
                                    start_date=datetime(2011, 3, 1))

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
import pandas._tseries as lib
N = 10000
K = 10

key1 = np.array([rands(10) for _ in xrange(N)], dtype='O').repeat(K)
key2 = np.array([rands(10) for _ in xrange(N)], dtype='O').repeat(K)

df = DataFrame({'key1' : key1, 'key2' : key2,
                'value' : np.random.randn(N * K)})
"""
statement = "df.sort_index(by=['key1', 'key2'])"
frame_sort_index_by_columns = Benchmark(statement, setup,
                                        name='frame_sort_index_by_columns',
                                        start_date=datetime(2011, 11, 1))

# drop_duplicates

statement = "df.drop_duplicates(['key1', 'key2'])"
frame_drop_duplicates = Benchmark(statement, setup,
                                  name='frame_drop_duplicates',
                                  start_date=datetime(2011, 11, 15))

statement = "df.drop_duplicates(['key1', 'key2'], inplace=True)"
frame_drop_dup_inplace = Benchmark(statement, setup,
                                  name='frame_drop_dup_inplace',
                                  start_date=datetime(2012, 5, 16))

lib_fast_zip = Benchmark('lib.fast_zip(df.values.T)', setup,
                         name='lib_fast_zip',
                         start_date=datetime(2012, 1, 1))

setup = setup + """
df.ix[:10000, :] = np.nan
"""
statement2 = "df.drop_duplicates(['key1', 'key2'])"
frame_drop_duplicates_na = Benchmark(statement2, setup,
                                     name='frame_drop_duplicates_na',
                                     start_date=datetime(2012, 5, 15))

lib_fast_zip_fillna = Benchmark('lib.fast_zip_fillna(df.values.T)', setup,
                                name='lib_fast_zip_fillna',
                                start_date=datetime(2012, 5, 15))

statement2 = "df.drop_duplicates(['key1', 'key2'], inplace=True)"
frame_drop_dup_na_inplace = Benchmark(statement2, setup,
                                  name='frame_drop_dup_na_inplace',
                                  start_date=datetime(2012, 5, 16))

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
