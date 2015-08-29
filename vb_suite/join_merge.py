from vbench.benchmark import Benchmark
from datetime import datetime

common_setup = """from .pandas_vb_common import *
"""

setup = common_setup + """
level1 = tm.makeStringIndex(10).values
level2 = tm.makeStringIndex(1000).values
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

df = pd.DataFrame({'data1' : np.random.randn(100000),
                'data2' : np.random.randn(100000),
                'key1' : key1,
                'key2' : key2})


df_key1 = pd.DataFrame(np.random.randn(len(level1), 4), index=level1,
                    columns=['A', 'B', 'C', 'D'])
df_key2 = pd.DataFrame(np.random.randn(len(level2), 4), index=level2,
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
              name='join_dataframe_index_single_key_bigger_sort',
              start_date=datetime(2012, 2, 5))

join_dataframe_index_multi = \
    Benchmark("df.join(df_multi, on=['key1', 'key2'])", setup,
              name='join_dataframe_index_multi',
              start_date=datetime(2011, 10, 20))

#----------------------------------------------------------------------
# Joins on integer keys
setup = common_setup + """
df = pd.DataFrame({'key1': np.tile(np.arange(500).repeat(10), 2),
                'key2': np.tile(np.arange(250).repeat(10), 4),
                'value': np.random.randn(10000)})
df2 = pd.DataFrame({'key1': np.arange(500), 'value2': randn(500)})
df3 = df[:5000]
"""


join_dataframe_integer_key = Benchmark("merge(df, df2, on='key1')", setup,
                                       start_date=datetime(2011, 10, 20))
join_dataframe_integer_2key = Benchmark("merge(df, df3)", setup,
                                        start_date=datetime(2011, 10, 20))

#----------------------------------------------------------------------
# DataFrame joins on index


#----------------------------------------------------------------------
# Merges
setup = common_setup + """
N = 10000

indices = tm.makeStringIndex(N).values
indices2 = tm.makeStringIndex(N).values
key = np.tile(indices[:8000], 10)
key2 = np.tile(indices2[:8000], 10)

left = pd.DataFrame({'key' : key, 'key2':key2,
                  'value' : np.random.randn(80000)})
right = pd.DataFrame({'key': indices[2000:], 'key2':indices2[2000:],
                   'value2' : np.random.randn(8000)})
"""

merge_2intkey_nosort = Benchmark('merge(left, right, sort=False)', setup,
                                 start_date=datetime(2011, 10, 20))

merge_2intkey_sort = Benchmark('merge(left, right, sort=True)', setup,
                               start_date=datetime(2011, 10, 20))

#----------------------------------------------------------------------
# Appending DataFrames

setup = common_setup + """
df1 = pd.DataFrame(np.random.randn(10000, 4), columns=['A', 'B', 'C', 'D'])
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
# indices = tm.makeStringIndex(n)
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
indices = tm.makeStringIndex(1000)
s = Series(n, index=indices)
pieces = [s[i:-i] for i in range(1, 10)]
pieces = pieces * 50
"""

concat_series_axis1 = Benchmark('concat(pieces, axis=1)', setup,
                                start_date=datetime(2012, 2, 27))

setup = common_setup + """
df = pd.DataFrame(randn(5, 4))
"""

concat_small_frames = Benchmark('concat([df] * 1000)', setup,
                                start_date=datetime(2012, 1, 1))


#----------------------------------------------------------------------
# Concat empty

setup = common_setup + """
df = pd.DataFrame(dict(A = range(10000)),index=date_range('20130101',periods=10000,freq='s'))
empty = pd.DataFrame()
"""

concat_empty_frames1 = Benchmark('concat([df,empty])', setup,
                                start_date=datetime(2012, 1, 1))
concat_empty_frames2 = Benchmark('concat([empty,df])', setup,
                                start_date=datetime(2012, 1, 1))


#----------------------------------------------------------------------
# Ordered merge

setup = common_setup + """
groups = tm.makeStringIndex(10).values

left = pd.DataFrame({'group': groups.repeat(5000),
                  'key' : np.tile(np.arange(0, 10000, 2), 10),
                  'lvalue': np.random.randn(50000)})

right = pd.DataFrame({'key' : np.arange(10000),
                   'rvalue' : np.random.randn(10000)})

"""

stmt = "ordered_merge(left, right, on='key', left_by='group')"

#----------------------------------------------------------------------
# outer join of non-unique
# GH 6329

setup = common_setup + """
date_index = date_range('01-Jan-2013', '23-Jan-2013', freq='T')
daily_dates = date_index.to_period('D').to_timestamp('S','S')
fracofday = date_index.view(np.ndarray) - daily_dates.view(np.ndarray)
fracofday = fracofday.astype('timedelta64[ns]').astype(np.float64)/864e11
fracofday = TimeSeries(fracofday, daily_dates)
index = date_range(date_index.min().to_period('A').to_timestamp('D','S'),
                      date_index.max().to_period('A').to_timestamp('D','E'),
                      freq='D')
temp = TimeSeries(1.0, index)
"""

join_non_unique_equal = Benchmark('fracofday * temp[fracofday.index]', setup,
                                   start_date=datetime(2013, 1, 1))


setup = common_setup + '''
np.random.seed(2718281)
n = 50000

left = pd.DataFrame(np.random.randint(1, n/500, (n, 2)),
        columns=['jim', 'joe'])

right = pd.DataFrame(np.random.randint(1, n/500, (n, 2)),
        columns=['jolie', 'jolia']).set_index('jolie')
'''

left_outer_join_index = Benchmark("left.join(right, on='jim')", setup,
                                  name='left_outer_join_index')


setup = common_setup + """
low, high, n = -1 << 10, 1 << 10, 1 << 20
left = pd.DataFrame(np.random.randint(low, high, (n, 7)),
                    columns=list('ABCDEFG'))
left['left'] = left.sum(axis=1)

i = np.random.permutation(len(left))
right = left.iloc[i].copy()
right.columns = right.columns[:-1].tolist() + ['right']
right.index = np.arange(len(right))
right['right'] *= -1
"""

i8merge = Benchmark("merge(left, right, how='outer')", setup,
                    name='i8merge')
