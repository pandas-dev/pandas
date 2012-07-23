from vbench.api import Benchmark
from datetime import datetime

common_setup = """from pandas_vb_common import *
"""

setup = common_setup + """
N = 100000
ngroups = 100

def get_test_data(ngroups=100, n=N):
    unique_groups = range(ngroups)
    arr = np.asarray(np.tile(unique_groups, n / ngroups), dtype=object)

    if len(arr) < n:
        arr = np.asarray(list(arr) + unique_groups[:n - len(arr)],
                         dtype=object)

    random.shuffle(arr)
    return arr

# aggregate multiple columns
df = DataFrame({'key1' : get_test_data(ngroups=ngroups),
                'key2' : get_test_data(ngroups=ngroups),
                'data1' : np.random.randn(N),
                'data2' : np.random.randn(N)})
def f():
    df.groupby(['key1', 'key2']).agg(lambda x: x.values.sum())

simple_series = Series(np.random.randn(N))
key1 = df['key1']
"""

stmt1 = "df.groupby(['key1', 'key2'])['data1'].agg(lambda x: x.values.sum())"
groupby_multi_python = Benchmark(stmt1, setup,
                                 start_date=datetime(2011, 7, 1))

stmt3 = "df.groupby(['key1', 'key2']).sum()"
groupby_multi_cython = Benchmark(stmt3, setup,
                                 start_date=datetime(2011, 7, 1))

stmt = "df.groupby(['key1', 'key2'])['data1'].agg(np.std)"
groupby_multi_series_op = Benchmark(stmt, setup,
                                    start_date=datetime(2011, 8, 1))

groupby_series_simple_cython = \
    Benchmark('simple_series.groupby(key1).sum()', setup,
              start_date=datetime(2011, 3, 1))

#----------------------------------------------------------------------
# 2d grouping, aggregate many columns

setup = common_setup + """
labels = np.random.randint(0, 100, size=1000)
df = DataFrame(randn(1000, 1000))
"""

groupby_frame_cython_many_columns = Benchmark('df.groupby(labels).sum()', setup,
                                              start_date=datetime(2011, 8, 1),
                                              logy=True)

#----------------------------------------------------------------------
# single key, long, integer key

setup = common_setup + """
data = np.random.randn(100000, 1)
labels = np.random.randint(0, 1000, size=100000)
df = DataFrame(data)
"""

groupby_frame_singlekey_integer = \
    Benchmark('df.groupby(labels).sum()', setup,
              start_date=datetime(2011, 8, 1), logy=True)

#----------------------------------------------------------------------
# group with different functions per column

setup = common_setup + """
fac1 = np.array(['A', 'B', 'C'], dtype='O')
fac2 = np.array(['one', 'two'], dtype='O')

df = DataFrame({'key1': fac1.take(np.random.randint(0, 3, size=100000)),
                'key2': fac2.take(np.random.randint(0, 2, size=100000)),
                'value1' : np.random.randn(100000),
                'value2' : np.random.randn(100000),
                'value3' : np.random.randn(100000)})
"""

groupby_multi_different_functions = \
    Benchmark("""df.groupby(['key1', 'key2']).agg({'value1' : 'mean',
                                                   'value2' : 'var',
                                                   'value3' : 'sum'})""",
              setup, start_date=datetime(2011, 9, 1))

groupby_multi_different_numpy_functions = \
    Benchmark("""df.groupby(['key1', 'key2']).agg({'value1' : np.mean,
                                                   'value2' : np.var,
                                                   'value3' : np.sum})""",
              setup, start_date=datetime(2011, 9, 1))

#----------------------------------------------------------------------
# size() speed

setup = common_setup + """
df = DataFrame({'key1': np.random.randint(0, 500, size=100000),
                'key2': np.random.randint(0, 100, size=100000),
                'value1' : np.random.randn(100000),
                'value2' : np.random.randn(100000),
                'value3' : np.random.randn(100000)})
"""

groupby_multi_size = Benchmark("df.groupby(['key1', 'key2']).size()",
                               setup, start_date=datetime(2011, 10, 1))

#----------------------------------------------------------------------
# Series.value_counts

setup = common_setup + """
s = Series(np.random.randint(0, 1000, size=100000))
"""

series_value_counts_int64 = Benchmark('s.value_counts()', setup,
                                      start_date=datetime(2011, 10, 21))

#----------------------------------------------------------------------
# pivot_table

setup = common_setup + """
fac1 = np.array(['A', 'B', 'C'], dtype='O')
fac2 = np.array(['one', 'two'], dtype='O')

ind1 = np.random.randint(0, 3, size=100000)
ind2 = np.random.randint(0, 2, size=100000)

df = DataFrame({'key1': fac1.take(ind1),
                'key2': fac2.take(ind2),
                'key3': fac2.take(ind2),
                'value1' : np.random.randn(100000),
                'value2' : np.random.randn(100000),
                'value3' : np.random.randn(100000)})
"""

stmt = "df.pivot_table(rows='key1', cols=['key2', 'key3'])"
groupby_pivot_table = Benchmark(stmt, setup, start_date=datetime(2011, 12, 15))


#----------------------------------------------------------------------
# dict return values

setup = common_setup + """
labels = np.arange(1000).repeat(10)
data = Series(randn(len(labels)))
f = lambda x: {'first': x.values[0], 'last': x.values[-1]}
"""

groupby_apply_dict_return = Benchmark('data.groupby(labels).apply(f)',
                                      setup, start_date=datetime(2011, 12, 15))

#----------------------------------------------------------------------
# First / last functions

setup = common_setup + """
labels = np.arange(10000).repeat(10)
data = Series(randn(len(labels)))
data[::3] = np.nan
data[1::3] = np.nan
labels = labels.take(np.random.permutation(len(labels)))
"""

groupby_first = Benchmark('data.groupby(labels).first()', setup,
                          start_date=datetime(2012, 5, 1))

groupby_last = Benchmark('data.groupby(labels).last()', setup,
                          start_date=datetime(2012, 5, 1))


#----------------------------------------------------------------------
# groupby_indices replacement, chop up Series

setup = common_setup + """
try:
    rng = date_range('1/1/2000', '12/31/2005', freq='H')
    year, month, day = rng.year, rng.month, rng.day
except:
    rng = date_range('1/1/2000', '12/31/2000', offset=datetools.Hour())
    year = rng.map(lambda x: x.year)
    month = rng.map(lambda x: x.month)
    day = rng.map(lambda x: x.day)

ts = Series(np.random.randn(len(rng)), index=rng)
"""

groupby_indices = Benchmark('len(ts.groupby([year, month, day]))',
                            setup, start_date=datetime(2012, 1, 1))

#----------------------------------------------------------------------
# median

#----------------------------------------------------------------------
# single key, long, integer key

setup = common_setup + """
data = np.random.randn(100000, 2)
labels = np.random.randint(0, 1000, size=100000)
df = DataFrame(data)
"""

groupby_frame_median = \
    Benchmark('df.groupby(labels).median()', setup,
              start_date=datetime(2011, 8, 1), logy=True)


setup = common_setup + """
data = np.random.randn(1000000, 2)
labels = np.random.randint(0, 1000, size=1000000)
df = DataFrame(data)
"""

groupby_simple_compress_timing = \
    Benchmark('df.groupby(labels).mean()', setup,
              start_date=datetime(2011, 8, 1))

