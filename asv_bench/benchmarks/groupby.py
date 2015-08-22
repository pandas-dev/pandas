from .pandas_vb_common import *
from string import ascii_letters, digits
from itertools import product


class groupby_agg_builtins(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(27182)
        self.n = 100000
        self.df = DataFrame(np.random.randint(1, (self.n / 100), (self.n, 3)), columns=['jim', 'joe', 'jolie'])

    def time_groupby_agg_builtins1(self):
        self.df.groupby('jim').agg([sum, min, max])

    def time_groupby_agg_builtins2(self):
        self.df.groupby(['jim', 'joe']).agg([sum, min, max])

#----------------------------------------------------------------------
# dict return values

class groupby_apply_dict_return(object):
    goal_time = 0.2

    def setup(self):
        self.labels = np.arange(1000).repeat(10)
        self.data = Series(randn(len(self.labels)))
        self.f = (lambda x: {'first': x.values[0], 'last': x.values[(-1)], })

    def time_groupby_apply_dict_return(self):
        self.data.groupby(self.labels).apply(self.f)


#----------------------------------------------------------------------
# First / last functions

class groupby_first_last(object):
    goal_time = 0.2

    def setup(self):
        self.labels = np.arange(10000).repeat(10)
        self.data = Series(randn(len(self.labels)))
        self.data[::3] = np.nan
        self.data[1::3] = np.nan
        self.data2 = Series(randn(len(self.labels)), dtype='float32')
        self.data2[::3] = np.nan
        self.data2[1::3] = np.nan
        self.labels = self.labels.take(np.random.permutation(len(self.labels)))

    def time_groupby_first_float32(self):
        self.data2.groupby(self.labels).first()

    def time_groupby_first_float64(self):
        self.data.groupby(self.labels).first()

    def time_groupby_last_float32(self):
        self.data2.groupby(self.labels).last()

    def time_groupby_last_float64(self):
        self.data.groupby(self.labels).last()

    def time_groupby_nth_float32_any(self):
        self.data2.groupby(self.labels).nth(0, dropna='all')

    def time_groupby_nth_float32_none(self):
        self.data2.groupby(self.labels).nth(0)

    def time_groupby_nth_float64_any(self):
        self.data.groupby(self.labels).nth(0, dropna='all')

    def time_groupby_nth_float64_none(self):
        self.data.groupby(self.labels).nth(0)

# with datetimes (GH7555)

class groupby_first_last_datetimes(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'a': date_range('1/1/2011', periods=100000, freq='s'), 'b': range(100000), })

    def time_groupby_first_datetimes(self):
        self.df.groupby('b').first()

    def time_groupby_last_datetimes(self):
        self.df.groupby('b').last()

    def time_groupby_nth_datetimes_any(self):
        self.df.groupby('b').nth(0, dropna='all')

    def time_groupby_nth_datetimes_none(self):
        self.df.groupby('b').nth(0)


class groupby_first_last_object(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'a': (['foo'] * 100000), 'b': range(100000)})

    def time_groupby_first_object(self):
        self.df.groupby('b').first()

    def time_groupby_last_object(self):
        self.df.groupby('b').last()

    def time_groupby_nth_object_any(self):
        self.df.groupby('b').nth(0, dropna='any')

    def time_groupby_nth_object_none(self):
        self.df.groupby('b').nth(0)


#----------------------------------------------------------------------
# DataFrame Apply overhead

class groupby_frame_apply(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.labels = np.random.randint(0, 2000, size=self.N)
        self.labels2 = np.random.randint(0, 3, size=self.N)
        self.df = DataFrame({'key': self.labels, 'key2': self.labels2, 'value1': randn(self.N), 'value2': (['foo', 'bar', 'baz', 'qux'] * (self.N / 4)), })

    def f(self, g):
        return 1

    def time_groupby_frame_apply(self):
        self.df.groupby(['key', 'key2']).apply(self.f)

    def time_groupby_frame_apply_overhead(self):
        self.df.groupby('key').apply(self.f)


#----------------------------------------------------------------------
# 2d grouping, aggregate many columns

class groupby_frame_cython_many_columns(object):
    goal_time = 0.2

    def setup(self):
        self.labels = np.random.randint(0, 100, size=1000)
        self.df = DataFrame(randn(1000, 1000))

    def time_sum(self):
        self.df.groupby(self.labels).sum()


#----------------------------------------------------------------------
# single key, long, integer key

class groupby_frame_singlekey_integer(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(100000, 1)
        self.labels = np.random.randint(0, 1000, size=100000)
        self.df = DataFrame(self.data)

    def time_sum(self):
        self.df.groupby(self.labels).sum()


#----------------------------------------------------------------------
# median

class groupby_frame(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(100000, 2)
        self.labels = np.random.randint(0, 1000, size=100000)
        self.df = DataFrame(self.data)

    def time_groupby_frame_median(self):
        self.df.groupby(self.labels).median()

    def time_groupby_simple_compress_timing(self):
        self.df.groupby(self.labels).mean()


#----------------------------------------------------------------------
# DataFrame nth

class groupby_nth(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randint(1, 100, (10000, 2)))

    def time_groupby_frame_nth_any(self):
        self.df.groupby(0).nth(0, dropna='any')

    def time_groupby_frame_nth_none(self):
        self.df.groupby(0).nth(0)

    def time_groupby_series_nth_any(self):
        self.df[1].groupby(self.df[0]).nth(0, dropna='any')

    def time_groupby_series_nth_none(self):
        self.df[1].groupby(self.df[0]).nth(0)


#----------------------------------------------------------------------
# groupby_indices replacement, chop up Series

class groupby_indices(object):
    goal_time = 0.2

    def setup(self):
        try:
            self.rng = date_range('1/1/2000', '12/31/2005', freq='H')
            (self.year, self.month, self.day) = (self.rng.year, self.rng.month, self.rng.day)
        except:
            self.rng = date_range('1/1/2000', '12/31/2000', offset=datetools.Hour())
            self.year = self.rng.map((lambda x: x.year))
            self.month = self.rng.map((lambda x: x.month))
            self.day = self.rng.map((lambda x: x.day))
        self.ts = Series(np.random.randn(len(self.rng)), index=self.rng)

    def time_groupby_indices(self):
        len(self.ts.groupby([self.year, self.month, self.day]))


class groupby_int64_overflow(object):
    goal_time = 0.2

    def setup(self):
        self.arr = np.random.randint(((-1) << 12), (1 << 12), ((1 << 17), 5))
        self.i = np.random.choice(len(self.arr), (len(self.arr) * 5))
        self.arr = np.vstack((self.arr, self.arr[self.i]))
        self.i = np.random.permutation(len(self.arr))
        self.arr = self.arr[self.i]
        self.df = DataFrame(self.arr, columns=list('abcde'))
        (self.df['jim'], self.df['joe']) = (np.random.randn(2, len(self.df)) * 10)

    def time_groupby_int64_overflow(self):
        self.df.groupby(list('abcde')).max()


#----------------------------------------------------------------------
# count() speed

class groupby_multi_count(object):
    goal_time = 0.2

    def setup(self):
        self.n = 10000
        self.offsets = np.random.randint(self.n, size=self.n).astype('timedelta64[ns]')
        self.dates = (np.datetime64('now') + self.offsets)
        self.dates[(np.random.rand(self.n) > 0.5)] = np.datetime64('nat')
        self.offsets[(np.random.rand(self.n) > 0.5)] = np.timedelta64('nat')
        self.value2 = np.random.randn(self.n)
        self.value2[(np.random.rand(self.n) > 0.5)] = np.nan
        self.obj = tm.choice(list('ab'), size=self.n).astype(object)
        self.obj[(np.random.randn(self.n) > 0.5)] = np.nan
        self.df = DataFrame({'key1': np.random.randint(0, 500, size=self.n),
                             'key2': np.random.randint(0, 100, size=self.n),
                             'dates': self.dates,
                             'value2': self.value2,
                             'value3': np.random.randn(self.n),
                             'ints': np.random.randint(0, 1000, size=self.n),
                             'obj': self.obj,
                             'offsets': self.offsets, })

    def time_groupby_multi_count(self):
        self.df.groupby(['key1', 'key2']).count()


class groupby_int_count(object):
    goal_time = 0.2

    def setup(self):
        self.n = 10000
        self.df = DataFrame({'key1': randint(0, 500, size=self.n),
                             'key2': randint(0, 100, size=self.n),
                             'ints': randint(0, 1000, size=self.n),
                             'ints2': randint(0, 1000, size=self.n), })

    def time_groupby_int_count(self):
        self.df.groupby(['key1', 'key2']).count()


#----------------------------------------------------------------------
# group with different functions per column

class groupby_agg_multi(object):
    goal_time = 0.2

    def setup(self):
        self.fac1 = np.array(['A', 'B', 'C'], dtype='O')
        self.fac2 = np.array(['one', 'two'], dtype='O')
        self.df = DataFrame({'key1': self.fac1.take(np.random.randint(0, 3, size=100000)), 'key2': self.fac2.take(np.random.randint(0, 2, size=100000)), 'value1': np.random.randn(100000), 'value2': np.random.randn(100000), 'value3': np.random.randn(100000), })

    def time_groupby_multi_different_functions(self):
        self.df.groupby(['key1', 'key2']).agg({'value1': 'mean', 'value2': 'var', 'value3': 'sum'})

    def time_groupby_multi_different_numpy_functions(self):
        self.df.groupby(['key1', 'key2']).agg({'value1': np.mean, 'value2': np.var, 'value3': np.sum})


class groupby_multi_index(object):
    goal_time = 0.2

    def setup(self):
        self.n = (((5 * 7) * 11) * (1 << 9))
        self.alpha = list(map(''.join, product((ascii_letters + digits), repeat=4)))
        self.f = (lambda k: np.repeat(np.random.choice(self.alpha, (self.n // k)), k))
        self.df = DataFrame({'a': self.f(11), 'b': self.f(7), 'c': self.f(5), 'd': self.f(1), })
        self.df['joe'] = (np.random.randn(len(self.df)) * 10).round(3)
        self.i = np.random.permutation(len(self.df))
        self.df = self.df.iloc[self.i].reset_index(drop=True).copy()

    def time_groupby_multi_index(self):
        self.df.groupby(list('abcd')).max()


class groupby_multi(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.ngroups = 100
        self.df = DataFrame({'key1': self.get_test_data(ngroups=self.ngroups), 'key2': self.get_test_data(ngroups=self.ngroups), 'data1': np.random.randn(self.N), 'data2': np.random.randn(self.N), })
        self.simple_series = Series(np.random.randn(self.N))
        self.key1 = self.df['key1']

    def get_test_data(self, ngroups=100, n=100000):
        self.unique_groups = range(self.ngroups)
        self.arr = np.asarray(np.tile(self.unique_groups, (n / self.ngroups)), dtype=object)
        if (len(self.arr) < n):
            self.arr = np.asarray((list(self.arr) + self.unique_groups[:(n - len(self.arr))]), dtype=object)
        random.shuffle(self.arr)
        return self.arr

    def f(self):
        self.df.groupby(['key1', 'key2']).agg((lambda x: x.values.sum()))

    def time_groupby_multi_cython(self):
        self.df.groupby(['key1', 'key2']).sum()

    def time_groupby_multi_python(self):
        self.df.groupby(['key1', 'key2'])['data1'].agg((lambda x: x.values.sum()))

    def time_groupby_multi_series_op(self):
        self.df.groupby(['key1', 'key2'])['data1'].agg(np.std)

    def time_groupby_series_simple_cython(self):
        self.simple_series.groupby(self.key1).sum()

    def time_groupby_series_simple_rank(self):
        self.df.groupby('key1').rank(pct=True)


#----------------------------------------------------------------------
# size() speed

class groupby_size(object):
    goal_time = 0.2

    def setup(self):
        self.n = 100000
        self.offsets = np.random.randint(self.n, size=self.n).astype('timedelta64[ns]')
        self.dates = (np.datetime64('now') + self.offsets)
        self.df = DataFrame({'key1': np.random.randint(0, 500, size=self.n), 'key2': np.random.randint(0, 100, size=self.n), 'value1': np.random.randn(self.n), 'value2': np.random.randn(self.n), 'value3': np.random.randn(self.n), 'dates': self.dates, })

    def time_groupby_multi_size(self):
        self.df.groupby(['key1', 'key2']).size()

    def time_groupby_dt_size(self):
        self.df.groupby(['dates']).size()

    def time_groupby_dt_timegrouper_size(self):
        self.df.groupby(TimeGrouper(key='dates', freq='M')).size()


#----------------------------------------------------------------------
# groupby with a variable value for ngroups

class groupby_ngroups_10000(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_all(self):
        self.df.groupby('value')['timestamp'].all()

    def time_any(self):
        self.df.groupby('value')['timestamp'].any()

    def time_count(self):
        self.df.groupby('value')['timestamp'].count()

    def time_cumcount(self):
        self.df.groupby('value')['timestamp'].cumcount()

    def time_cummax(self):
        self.df.groupby('value')['timestamp'].cummax()

    def time_cummin(self):
        self.df.groupby('value')['timestamp'].cummin()

    def time_cumprod(self):
        self.df.groupby('value')['timestamp'].cumprod()

    def time_cumsum(self):
        self.df.groupby('value')['timestamp'].cumsum()

    def time_describe(self):
        self.df.groupby('value')['timestamp'].describe()

    def time_diff(self):
        self.df.groupby('value')['timestamp'].diff()

    def time_first(self):
        self.df.groupby('value')['timestamp'].first()

    def time_head(self):
        self.df.groupby('value')['timestamp'].head()

    def time_last(self):
        self.df.groupby('value')['timestamp'].last()

    def time_mad(self):
        self.df.groupby('value')['timestamp'].mad()

    def time_max(self):
        self.df.groupby('value')['timestamp'].max()

    def time_mean(self):
        self.df.groupby('value')['timestamp'].mean()

    def time_median(self):
        self.df.groupby('value')['timestamp'].median()

    def time_min(self):
        self.df.groupby('value')['timestamp'].min()

    def time_nunique(self):
        self.df.groupby('value')['timestamp'].nunique()

    def time_pct_change(self):
        self.df.groupby('value')['timestamp'].pct_change()

    def time_prod(self):
        self.df.groupby('value')['timestamp'].prod()

    def time_rank(self):
        self.df.groupby('value')['timestamp'].rank()

    def time_sem(self):
        self.df.groupby('value')['timestamp'].sem()

    def time_size(self):
        self.df.groupby('value')['timestamp'].size()

    def time_skew(self):
        self.df.groupby('value')['timestamp'].skew()

    def time_std(self):
        self.df.groupby('value')['timestamp'].std()

    def time_sum(self):
        self.df.groupby('value')['timestamp'].sum()

    def time_tail(self):
        self.df.groupby('value')['timestamp'].tail()

    def time_unique(self):
        self.df.groupby('value')['timestamp'].unique()

    def time_value_counts(self):
        self.df.groupby('value')['timestamp'].value_counts()

    def time_var(self):
        self.df.groupby('value')['timestamp'].var()


class groupby_ngroups_100(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_all(self):
        self.df.groupby('value')['timestamp'].all()

    def time_any(self):
        self.df.groupby('value')['timestamp'].any()

    def time_count(self):
        self.df.groupby('value')['timestamp'].count()

    def time_cumcount(self):
        self.df.groupby('value')['timestamp'].cumcount()

    def time_cummax(self):
        self.df.groupby('value')['timestamp'].cummax()

    def time_cummin(self):
        self.df.groupby('value')['timestamp'].cummin()

    def time_cumprod(self):
        self.df.groupby('value')['timestamp'].cumprod()

    def time_cumsum(self):
        self.df.groupby('value')['timestamp'].cumsum()

    def time_describe(self):
        self.df.groupby('value')['timestamp'].describe()

    def time_diff(self):
        self.df.groupby('value')['timestamp'].diff()

    def time_first(self):
        self.df.groupby('value')['timestamp'].first()

    def time_head(self):
        self.df.groupby('value')['timestamp'].head()

    def time_last(self):
        self.df.groupby('value')['timestamp'].last()

    def time_mad(self):
        self.df.groupby('value')['timestamp'].mad()

    def time_max(self):
        self.df.groupby('value')['timestamp'].max()

    def time_mean(self):
        self.df.groupby('value')['timestamp'].mean()

    def time_median(self):
        self.df.groupby('value')['timestamp'].median()

    def time_min(self):
        self.df.groupby('value')['timestamp'].min()

    def time_nunique(self):
        self.df.groupby('value')['timestamp'].nunique()

    def time_pct_change(self):
        self.df.groupby('value')['timestamp'].pct_change()

    def time_prod(self):
        self.df.groupby('value')['timestamp'].prod()

    def time_rank(self):
        self.df.groupby('value')['timestamp'].rank()

    def time_sem(self):
        self.df.groupby('value')['timestamp'].sem()

    def time_size(self):
        self.df.groupby('value')['timestamp'].size()

    def time_skew(self):
        self.df.groupby('value')['timestamp'].skew()

    def time_std(self):
        self.df.groupby('value')['timestamp'].std()

    def time_sum(self):
        self.df.groupby('value')['timestamp'].sum()

    def time_tail(self):
        self.df.groupby('value')['timestamp'].tail()

    def time_unique(self):
        self.df.groupby('value')['timestamp'].unique()

    def time_value_counts(self):
        self.df.groupby('value')['timestamp'].value_counts()

    def time_var(self):
        self.df.groupby('value')['timestamp'].var()


#----------------------------------------------------------------------
# Series.value_counts

class series_value_counts(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.randint(0, 1000, size=100000))
        self.s2 = self.s.astype(float)

        self.K = 1000
        self.N = 100000
        self.uniques = tm.makeStringIndex(self.K).values
        self.s3 = Series(np.tile(self.uniques, (self.N // self.K)))

    def time_value_counts_int64(self):
        self.s.value_counts()

    def time_value_counts_float64(self):
        self.s2.value_counts()

    def time_value_counts_strings(self):
        self.s.value_counts()


#----------------------------------------------------------------------
# pivot_table

class groupby_pivot_table(object):
    goal_time = 0.2

    def setup(self):
        self.fac1 = np.array(['A', 'B', 'C'], dtype='O')
        self.fac2 = np.array(['one', 'two'], dtype='O')
        self.ind1 = np.random.randint(0, 3, size=100000)
        self.ind2 = np.random.randint(0, 2, size=100000)
        self.df = DataFrame({'key1': self.fac1.take(self.ind1), 'key2': self.fac2.take(self.ind2), 'key3': self.fac2.take(self.ind2), 'value1': np.random.randn(100000), 'value2': np.random.randn(100000), 'value3': np.random.randn(100000), })

    def time_groupby_pivot_table(self):
        self.df.pivot_table(index='key1', columns=['key2', 'key3'])


#----------------------------------------------------------------------
# Sum booleans #2692

class groupby_sum_booleans(object):
    goal_time = 0.2

    def setup(self):
        self.N = 500
        self.df = DataFrame({'ii': range(self.N), 'bb': [True for x in range(self.N)], })

    def time_groupby_sum_booleans(self):
        self.df.groupby('ii').sum()


#----------------------------------------------------------------------
# multi-indexed group sum #9049

class groupby_sum_multiindex(object):
    goal_time = 0.2

    def setup(self):
        self.N = 50
        self.df = DataFrame({'A': (range(self.N) * 2), 'B': range((self.N * 2)), 'C': 1, }).set_index(['A', 'B'])

    def time_groupby_sum_multiindex(self):
        self.df.groupby(level=[0, 1]).sum()


#-------------------------------------------------------------------------------
# Transform testing

class groupby_transform(object):
    goal_time = 0.2

    def setup(self):
        self.n_dates = 400
        self.n_securities = 250
        self.n_columns = 3
        self.share_na = 0.1
        self.dates = date_range('1997-12-31', periods=self.n_dates, freq='B')
        self.dates = Index(map((lambda x: (((x.year * 10000) + (x.month * 100)) + x.day)), self.dates))
        self.secid_min = int('10000000', 16)
        self.secid_max = int('F0000000', 16)
        self.step = ((self.secid_max - self.secid_min) // (self.n_securities - 1))
        self.security_ids = map((lambda x: hex(x)[2:10].upper()), range(self.secid_min, (self.secid_max + 1), self.step))
        self.data_index = MultiIndex(levels=[self.dates.values, self.security_ids],
                                     labels=[[i for i in range(self.n_dates) for _ in range(self.n_securities)], (range(self.n_securities) * self.n_dates)],
                                     names=['date', 'security_id'])
        self.n_data = len(self.data_index)
        self.columns = Index(['factor{}'.format(i) for i in range(1, (self.n_columns + 1))])
        self.data = DataFrame(np.random.randn(self.n_data, self.n_columns), index=self.data_index, columns=self.columns)
        self.step = int((self.n_data * self.share_na))
        for column_index in range(self.n_columns):
            self.index = column_index
            while (self.index < self.n_data):
                self.data.set_value(self.data_index[self.index], self.columns[column_index], np.nan)
                self.index += self.step
        self.f_fillna = (lambda x: x.fillna(method='pad'))

    def time_groupby_transform(self):
        self.data.groupby(level='security_id').transform(self.f_fillna)

    def time_groupby_transform_ufunc(self):
        self.data.groupby(level='date').transform(np.max)


class groupby_transform_multi_key(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(2718281)
        self.n = 20000
        self.df = DataFrame(np.random.randint(1, self.n, (self.n, 3)), columns=['jim', 'joe', 'jolie'])

    def time_groupby_transform_multi_key1(self):
        self.df.groupby(['jim', 'joe'])['jolie'].transform('max')


class groupby_transform_multi_key2(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(2718281)
        self.n = 20000
        self.df = DataFrame(np.random.randint(1, self.n, (self.n, 3)), columns=['jim', 'joe', 'jolie'])
        self.df['jim'] = self.df['joe']

    def time_groupby_transform_multi_key2(self):
        self.df.groupby(['jim', 'joe'])['jolie'].transform('max')


class groupby_transform_multi_key3(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(2718281)
        self.n = 200000
        self.df = DataFrame(np.random.randint(1, (self.n / 10), (self.n, 3)), columns=['jim', 'joe', 'jolie'])

    def time_groupby_transform_multi_key3(self):
        self.df.groupby(['jim', 'joe'])['jolie'].transform('max')


class groupby_transform_multi_key4(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(2718281)
        self.n = 200000
        self.df = DataFrame(np.random.randint(1, (self.n / 10), (self.n, 3)), columns=['jim', 'joe', 'jolie'])
        self.df['jim'] = self.df['joe']

    def time_groupby_transform_multi_key4(self):
        self.df.groupby(['jim', 'joe'])['jolie'].transform('max')


class groupby_transform_series(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(0)
        self.N = 120000
        self.N_TRANSITIONS = 1400
        self.transition_points = np.random.permutation(np.arange(self.N))[:self.N_TRANSITIONS]
        self.transition_points.sort()
        self.transitions = np.zeros((self.N,), dtype=np.bool)
        self.transitions[self.transition_points] = True
        self.g = self.transitions.cumsum()
        self.df = DataFrame({'signal': np.random.rand(self.N), })

    def time_groupby_transform_series(self):
        self.df['signal'].groupby(self.g).transform(np.mean)


class groupby_transform_series2(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(0)
        self.df = DataFrame({'id': (np.arange(100000) / 3), 'val': np.random.randn(100000), })

    def time_groupby_transform_series2(self):
        self.df.groupby('id')['val'].transform(np.mean)

class groupby_transform_cythonized(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(0)
        self.df = DataFrame({'id': (np.arange(100000) / 3), 'val': np.random.randn(100000), })

    def time_groupby_transform_cumprod(self):
        self.df.groupby('id').cumprod()

    def time_groupby_transform_cumsum(self):
        self.df.groupby('id').cumsum()

    def time_groupby_transform_shift(self):
        self.df.groupby('id').shift()
