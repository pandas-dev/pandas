from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

setup = common_setup + """
level1 = np.array([rands(10) for _ in xrange(10)], dtype='O')
level2 = np.array([rands(10) for _ in xrange(1000)], dtype='O')
label1 = np.arange(10).repeat(1000)
label2 = np.tile(np.arange(1000), 10)

key1 = np.tile(level1.take(label1), 10)
key2 = np.tile(level2.take(label2), 10)

shuf = np.arange(100000)
random.shuffle(shuf)
try:
    index2 = MultiIndex(levels=[level1, level2], labels=[label1, label2])
    index3 = MultiIndex(levels=[np.arange(10), np.arange(100), np.arange(100)],
                        labels=[np.arange(10).repeat(10000),
                                np.tile(np.arange(100).repeat(100), 10),
                                np.tile(np.tile(np.arange(100), 100), 10)])
    df_multi = DataFrame(np.random.randn(len(index2), 4), index=index2,
                         columns=['A', 'B', 'C', 'D'])
except:  # pre-MultiIndex
    pass

try:
    DataFrame = DataMatrix
except:
    pass

df = DataFrame({'data1' : np.random.randn(100000),
                'data2' : np.random.randn(100000),
                'key1' : key1,
                'key2' : key2})


df_key1 = DataFrame(np.random.randn(len(level1), 4), index=level1,
                    columns=['A', 'B', 'C', 'D'])
df_key2 = DataFrame(np.random.randn(len(level2), 4), index=level2,
                    columns=['A', 'B', 'C', 'D'])

df_shuf = df.reindex(df.index[shuf])
"""

#----------------------------------------------------------------------
# DataFrame joins on key

join_dataframe_index_single_key_small = \
    Benchmark("df.join(df_key1, on='key1')", setup,
              name='join_dataframe_index_single_key_small')

join_dataframe_index_single_key_bigger = \
    Benchmark("df.join(df_key2, on='key2')", setup,
              name='join_dataframe_index_single_key_bigger')

join_dataframe_index_single_key_bigger_sort = \
    Benchmark("df_shuf.join(df_key2, on='key2', sort=True)", setup,
              name='join_dataframe_index_single_key_bigger',
              start_date=datetime(2012, 2, 5))

join_dataframe_index_multi = \
    Benchmark("df.join(df_multi, on=['key1', 'key2'])", setup,
              name='join_dataframe_index_multi',
              start_date=datetime(2011, 10, 20))

#----------------------------------------------------------------------
# Joins on integer keys

join_dataframe_integer_key = Benchmark("merge(df, df2, on='key')", setup,
                                       start_date=datetime(2011, 10, 20))

#----------------------------------------------------------------------
# DataFrame joins on index



#----------------------------------------------------------------------
# Merges

#----------------------------------------------------------------------
# Appending DataFrames

setup = common_setup + """
df1 = DataFrame(np.random.randn(10000, 4), columns=['A', 'B', 'C', 'D'])
df2 = df1.copy()
df2.index = np.arange(10000, 20000)
mdf1 = df1.copy()
mdf1['obj1'] = 'bar'
mdf1['obj2'] = 'bar'
mdf1['int1'] = 5
try:
    mdf1.consolidate(inplace=True)
except:
    pass
mdf2 = mdf1.copy()
mdf2.index = df2.index
"""

stmt = "df1.append(df2)"
append_frame_single_homogenous = \
    Benchmark(stmt, setup, name='append_frame_single_homogenous',
              ncalls=500, repeat=1)

stmt = "mdf1.append(mdf2)"
append_frame_single_mixed = Benchmark(stmt, setup,
                                      name='append_frame_single_mixed',
                                      ncalls=500, repeat=1)

#----------------------------------------------------------------------
# data alignment

setup = common_setup + """n = 1000000
# indices = Index([rands(10) for _ in xrange(n)])
def sample(values, k):
    sampler = np.random.permutation(len(values))
    return values.take(sampler[:k])
sz = 500000
rng = np.arange(0, 10000000000000, 10000000)
stamps = np.datetime64(datetime.now()).view('i8') + rng
idx1 = np.sort(sample(stamps, sz))
idx2 = np.sort(sample(stamps, sz))
ts1 = Series(np.random.randn(sz), idx1)
ts2 = Series(np.random.randn(sz), idx2)
"""
stmt = "ts1 + ts2"
series_align_int64_index = \
    Benchmark(stmt, setup,
              name="series_align_int64_index",
              start_date=datetime(2010, 6, 1), logy=True)

stmt = "ts1.align(ts2, join='left')"
series_align_left_monotonic = \
    Benchmark(stmt, setup,
              name="series_align_left_monotonic",
              start_date=datetime(2011, 12, 1), logy=True)

#----------------------------------------------------------------------
# Concat Series axis=1

setup = common_setup + """
n = 1000
indices = Index([rands(10) for _ in xrange(1000)])
s = Series(n, index=indices)
pieces = [s[i:-i] for i in range(1, 10)]
pieces = pieces * 50
"""

concat_series_axis1 = Benchmark('concat(pieces, axis=1)', setup,
                                start_date=datetime(2012, 2, 27))

#----------------------------------------------------------------------
# Ordered merge

setup = common_setup + """
groups = np.array([rands(10) for _ in xrange(10)], dtype='O')

left = DataFrame({'group': groups.repeat(5000),
                  'key' : np.tile(np.arange(0, 10000, 2), 10),
                  'lvalue': np.random.randn(50000)})

right = DataFrame({'key' : np.arange(10000),
                   'rvalue' : np.random.randn(10000)})

"""

stmt = "ordered_merge(left, right, on='key', left_by='group')"
