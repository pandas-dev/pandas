from vbench.api import Benchmark
from datetime import datetime

common_setup = """from .pandas_vb_common import *
"""

setup = common_setup + """
N = 100000
ngroups = 100

def get_test_data(ngroups=100, n=100000):
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


stmt4 = "df.groupby('key1').rank(pct=True)"
groupby_series_simple_cython = Benchmark(stmt4, setup,
                                    start_date=datetime(2014, 1, 16))

#----------------------------------------------------------------------
# 2d grouping, aggregate many columns

setup = common_setup + """
labels = np.random.randint(0, 100, size=1000)
df = DataFrame(randn(1000, 1000))
"""

groupby_frame_cython_many_columns = Benchmark(
    'df.groupby(labels).sum()', setup,
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
n = 100000
offsets = np.random.randint(n, size=n).astype('timedelta64[ns]')
dates = np.datetime64('now') + offsets
df = DataFrame({'key1': np.random.randint(0, 500, size=n),
                'key2': np.random.randint(0, 100, size=n),
                'value1' : np.random.randn(n),
                'value2' : np.random.randn(n),
                'value3' : np.random.randn(n),
                'dates' : dates})
"""

groupby_multi_size = Benchmark("df.groupby(['key1', 'key2']).size()",
                               setup, start_date=datetime(2011, 10, 1))

groupby_dt_size = Benchmark("df.groupby(['dates']).size()",
                            setup, start_date=datetime(2011, 10, 1))

groupby_dt_timegrouper_size = Benchmark("df.groupby(TimeGrouper(key='dates', freq='M')).size()",
                                        setup, start_date=datetime(2011, 10, 1))

#----------------------------------------------------------------------
# count() speed

setup = common_setup + """
n = 10000
offsets = np.random.randint(n, size=n).astype('timedelta64[ns]')

dates = np.datetime64('now') + offsets
dates[np.random.rand(n) > 0.5] = np.datetime64('nat')

offsets[np.random.rand(n) > 0.5] = np.timedelta64('nat')

value2 = np.random.randn(n)
value2[np.random.rand(n) > 0.5] = np.nan

obj = np.random.choice(list('ab'), size=n).astype(object)
obj[np.random.randn(n) > 0.5] = np.nan

df = DataFrame({'key1': np.random.randint(0, 500, size=n),
                'key2': np.random.randint(0, 100, size=n),
                'dates': dates,
                'value2' : value2,
                'value3' : np.random.randn(n),
                'ints': np.random.randint(0, 1000, size=n),
                'obj': obj,
                'offsets': offsets})
"""

groupby_multi_count = Benchmark("df.groupby(['key1', 'key2']).count()",
                                setup, name='groupby_multi_count',
                                start_date=datetime(2014, 5, 5))

setup = common_setup + """
n = 10000

df = DataFrame({'key1': randint(0, 500, size=n),
                'key2': randint(0, 100, size=n),
                'ints': randint(0, 1000, size=n),
                'ints2': randint(0, 1000, size=n)})
"""

groupby_int_count = Benchmark("df.groupby(['key1', 'key2']).count()",
                              setup, name='groupby_int_count',
                              start_date=datetime(2014, 5, 6))
#----------------------------------------------------------------------
# Series.value_counts

setup = common_setup + """
s = Series(np.random.randint(0, 1000, size=100000))
"""

series_value_counts_int64 = Benchmark('s.value_counts()', setup,
                                      start_date=datetime(2011, 10, 21))

# value_counts on lots of strings

setup = common_setup + """
K = 1000
N = 100000
uniques = tm.makeStringIndex(K).values
s = Series(np.tile(uniques, N // K))
"""

series_value_counts_strings = Benchmark('s.value_counts()', setup,
                                        start_date=datetime(2011, 10, 21))

#value_counts on float dtype

setup = common_setup + """
s = Series(np.random.randint(0, 1000, size=100000)).astype(float)
"""

series_value_counts_float64 = Benchmark('s.value_counts()', setup,
                                      start_date=datetime(2015, 8, 17))

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

stmt = "df.pivot_table(index='key1', columns=['key2', 'key3'])"
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
data2 = Series(randn(len(labels)),dtype='float32')
data2[::3] = np.nan
data2[1::3] = np.nan
labels = labels.take(np.random.permutation(len(labels)))
"""

groupby_first_float64 = Benchmark('data.groupby(labels).first()', setup,
                                  start_date=datetime(2012, 5, 1))

groupby_first_float32 = Benchmark('data2.groupby(labels).first()', setup,
                                  start_date=datetime(2013, 1, 1))

groupby_last_float64 = Benchmark('data.groupby(labels).last()', setup,
                                 start_date=datetime(2012, 5, 1))

groupby_last_float32 = Benchmark('data2.groupby(labels).last()', setup,
                                 start_date=datetime(2013, 1, 1))

groupby_nth_float64_none = Benchmark('data.groupby(labels).nth(0)', setup,
                                     start_date=datetime(2012, 5, 1))
groupby_nth_float32_none = Benchmark('data2.groupby(labels).nth(0)', setup,
                                     start_date=datetime(2013, 1, 1))
groupby_nth_float64_any = Benchmark('data.groupby(labels).nth(0,dropna="all")', setup,
                                    start_date=datetime(2012, 5, 1))
groupby_nth_float32_any = Benchmark('data2.groupby(labels).nth(0,dropna="all")', setup,
                                    start_date=datetime(2013, 1, 1))

# with datetimes (GH7555)
setup = common_setup + """
df = DataFrame({'a' : date_range('1/1/2011',periods=100000,freq='s'),'b' : range(100000)})
"""

groupby_first_datetimes = Benchmark('df.groupby("b").first()', setup,
                                    start_date=datetime(2013, 5, 1))
groupby_last_datetimes = Benchmark('df.groupby("b").last()', setup,
                                   start_date=datetime(2013, 5, 1))
groupby_nth_datetimes_none = Benchmark('df.groupby("b").nth(0)', setup,
                                       start_date=datetime(2013, 5, 1))
groupby_nth_datetimes_any = Benchmark('df.groupby("b").nth(0,dropna="all")', setup,
                                      start_date=datetime(2013, 5, 1))

# with object
setup = common_setup + """
df = DataFrame({'a' : ['foo']*100000,'b' : range(100000)})
"""

groupby_first_object = Benchmark('df.groupby("b").first()', setup,
                                 start_date=datetime(2013, 5, 1))
groupby_last_object = Benchmark('df.groupby("b").last()', setup,
                                 start_date=datetime(2013, 5, 1))
groupby_nth_object_none = Benchmark('df.groupby("b").nth(0)', setup,
                                    start_date=datetime(2013, 5, 1))
groupby_nth_object_any = Benchmark('df.groupby("b").nth(0,dropna="any")', setup,
                                   start_date=datetime(2013, 5, 1))

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


#----------------------------------------------------------------------
# DataFrame Apply overhead

setup = common_setup + """
N = 10000
labels = np.random.randint(0, 2000, size=N)
labels2 = np.random.randint(0, 3, size=N)
df = DataFrame({'key': labels,
'key2': labels2,
'value1': randn(N),
'value2': ['foo', 'bar', 'baz', 'qux'] * (N / 4)})
def f(g):
    return 1
"""

groupby_frame_apply_overhead = Benchmark("df.groupby('key').apply(f)", setup,
                                         start_date=datetime(2011, 10, 1))

groupby_frame_apply = Benchmark("df.groupby(['key', 'key2']).apply(f)", setup,
                                start_date=datetime(2011, 10, 1))


#----------------------------------------------------------------------
# DataFrame nth

setup = common_setup + """
df = DataFrame(np.random.randint(1, 100, (10000, 2)))
"""

# Not really a fair test as behaviour has changed!
groupby_frame_nth_none = Benchmark("df.groupby(0).nth(0)", setup,
                                   start_date=datetime(2014, 3, 1))

groupby_series_nth_none = Benchmark("df[1].groupby(df[0]).nth(0)", setup,
                                    start_date=datetime(2014, 3, 1))
groupby_frame_nth_any= Benchmark("df.groupby(0).nth(0,dropna='any')", setup,
                                 start_date=datetime(2014, 3, 1))

groupby_series_nth_any = Benchmark("df[1].groupby(df[0]).nth(0,dropna='any')", setup,
                                   start_date=datetime(2014, 3, 1))


#----------------------------------------------------------------------
# Sum booleans #2692

setup = common_setup + """
N = 500
df = DataFrame({'ii':range(N),'bb':[True for x in range(N)]})
"""

groupby_sum_booleans = Benchmark("df.groupby('ii').sum()", setup)


#----------------------------------------------------------------------
# multi-indexed group sum #9049

setup = common_setup + """
N = 50
df = DataFrame({'A': range(N) * 2, 'B': range(N*2), 'C': 1}).set_index(["A", "B"])
"""

groupby_sum_multiindex = Benchmark("df.groupby(level=[0, 1]).sum()", setup)


#----------------------------------------------------------------------
# Transform testing

setup = common_setup + """
n_dates = 400
n_securities = 250
n_columns = 3
share_na = 0.1

dates = date_range('1997-12-31', periods=n_dates, freq='B')
dates = Index(map(lambda x: x.year * 10000 + x.month * 100 + x.day, dates))

secid_min = int('10000000', 16)
secid_max = int('F0000000', 16)
step = (secid_max - secid_min) // (n_securities - 1)
security_ids = map(lambda x: hex(x)[2:10].upper(), range(secid_min, secid_max + 1, step))

data_index = MultiIndex(levels=[dates.values, security_ids],
    labels=[[i for i in range(n_dates) for _ in xrange(n_securities)], range(n_securities) * n_dates],
    names=['date', 'security_id'])
n_data = len(data_index)

columns = Index(['factor{}'.format(i) for i in range(1, n_columns + 1)])

data = DataFrame(np.random.randn(n_data, n_columns), index=data_index, columns=columns)

step = int(n_data * share_na)
for column_index in range(n_columns):
    index = column_index
    while index < n_data:
        data.set_value(data_index[index], columns[column_index], np.nan)
        index += step

f_fillna = lambda x: x.fillna(method='pad')
"""

groupby_transform = Benchmark("data.groupby(level='security_id').transform(f_fillna)", setup)
groupby_transform_ufunc = Benchmark("data.groupby(level='date').transform(np.max)", setup)

setup = common_setup + """
np.random.seed(0)

N = 120000
N_TRANSITIONS = 1400

# generate groups
transition_points = np.random.permutation(np.arange(N))[:N_TRANSITIONS]
transition_points.sort()
transitions = np.zeros((N,), dtype=np.bool)
transitions[transition_points] = True
g = transitions.cumsum()

df = DataFrame({ 'signal' : np.random.rand(N)})
"""
groupby_transform_series = Benchmark("df['signal'].groupby(g).transform(np.mean)", setup)

setup = common_setup + """
np.random.seed(0)

df=DataFrame( { 'id' : np.arange( 100000 ) / 3,
                'val': np.random.randn( 100000) } )
"""

groupby_transform_series2 = Benchmark("df.groupby('id')['val'].transform(np.mean)", setup)

setup = common_setup + '''
np.random.seed(2718281)
n = 20000
df = DataFrame(np.random.randint(1, n, (n, 3)),
        columns=['jim', 'joe', 'jolie'])
'''

stmt = "df.groupby(['jim', 'joe'])['jolie'].transform('max')";
groupby_transform_multi_key1 = Benchmark(stmt, setup)
groupby_transform_multi_key2 = Benchmark(stmt, setup + "df['jim'] = df['joe']")

setup = common_setup + '''
np.random.seed(2718281)
n = 200000
df = DataFrame(np.random.randint(1, n / 10, (n, 3)),
        columns=['jim', 'joe', 'jolie'])
'''
groupby_transform_multi_key3 = Benchmark(stmt, setup)
groupby_transform_multi_key4 = Benchmark(stmt, setup + "df['jim'] = df['joe']")

setup = common_setup + '''
np.random.seed(27182)
n = 100000
df = DataFrame(np.random.randint(1, n / 100, (n, 3)),
        columns=['jim', 'joe', 'jolie'])
'''

groupby_agg_builtins1 = Benchmark("df.groupby('jim').agg([sum, min, max])", setup)
groupby_agg_builtins2 = Benchmark("df.groupby(['jim', 'joe']).agg([sum, min, max])", setup)


setup = common_setup + '''
arr = np.random.randint(- 1 << 12, 1 << 12, (1 << 17, 5))
i = np.random.choice(len(arr), len(arr) * 5)
arr = np.vstack((arr, arr[i]))  # add sume duplicate rows

i = np.random.permutation(len(arr))
arr = arr[i]  # shuffle rows

df = DataFrame(arr, columns=list('abcde'))
df['jim'], df['joe'] = np.random.randn(2, len(df)) * 10
'''

groupby_int64_overflow = Benchmark("df.groupby(list('abcde')).max()", setup,
                                   name='groupby_int64_overflow')


setup = common_setup + '''
from itertools import product
from string import ascii_letters, digits

n = 5 * 7 * 11 * (1 << 9)
alpha = list(map(''.join, product(ascii_letters + digits, repeat=4)))
f = lambda k: np.repeat(np.random.choice(alpha, n // k), k)

df = DataFrame({'a': f(11), 'b': f(7), 'c': f(5), 'd': f(1)})
df['joe'] = (np.random.randn(len(df)) * 10).round(3)

i = np.random.permutation(len(df))
df = df.iloc[i].reset_index(drop=True).copy()
'''

groupby_multi_index = Benchmark("df.groupby(list('abcd')).max()", setup,
                                name='groupby_multi_index')

#----------------------------------------------------------------------
# groupby with a variable value for ngroups


ngroups_list = [100, 10000]
no_arg_func_list = [
    'all',
    'any',
    'count',
    'cumcount',
    'cummax',
    'cummin',
    'cumprod',
    'cumsum',
    'describe',
    'diff',
    'first',
    'head',
    'last',
    'mad',
    'max',
    'mean',
    'median',
    'min',
    'nunique',
    'pct_change',
    'prod',
    'rank',
    'sem',
    'size',
    'skew',
    'std',
    'sum',
    'tail',
    'unique',
    'var',
    'value_counts',
]


_stmt_template = "df.groupby('value')['timestamp'].%s"
_setup_template = common_setup + """
np.random.seed(1234)
ngroups = %s
size = ngroups * 2
rng = np.arange(ngroups)
df = DataFrame(dict(
    timestamp=rng.take(np.random.randint(0, ngroups, size=size)),
    value=np.random.randint(0, size, size=size)
))
"""
START_DATE = datetime(2011, 7, 1)


def make_large_ngroups_bmark(ngroups, func_name, func_args=''):
    bmark_name = 'groupby_ngroups_%s_%s' % (ngroups, func_name)
    stmt = _stmt_template % ('%s(%s)' % (func_name, func_args))
    setup = _setup_template % ngroups
    bmark = Benchmark(stmt, setup, start_date=START_DATE)
    # MUST set name
    bmark.name = bmark_name
    return bmark


def inject_bmark_into_globals(bmark):
    if not bmark.name:
        raise AssertionError('benchmark must have a name')
    globals()[bmark.name] = bmark


for ngroups in ngroups_list:
    for func_name in no_arg_func_list:
        bmark = make_large_ngroups_bmark(ngroups, func_name)
        inject_bmark_into_globals(bmark)

# avoid bmark to be collected as Benchmark object
del bmark
