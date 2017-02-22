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
# groups

class Groups(object):
    goal_time = 0.1

    size = 2 ** 22
    data = {
        'int64_small': Series(np.random.randint(0, 100, size=size)),
        'int64_large' : Series(np.random.randint(0, 10000, size=size)),
        'object_small': Series(tm.makeStringIndex(100).take(np.random.randint(0, 100, size=size))),
        'object_large': Series(tm.makeStringIndex(10000).take(np.random.randint(0, 10000, size=size)))
    }

    param_names = ['df']
    params = ['int64_small', 'int64_large', 'object_small', 'object_large']

    def setup(self, df):
        self.df = self.data[df]

    def time_groupby_groups(self, df):
        self.df.groupby(self.df).groups


#----------------------------------------------------------------------
# First / last functions

class FirstLast(object):
    goal_time = 0.2

    param_names = ['dtype']
    params = ['float32', 'float64', 'datetime', 'object']

    # with datetimes (GH7555)

    def setup(self, dtype):

        if dtype == 'datetime':
            self.df = DataFrame(
                {'values': date_range('1/1/2011', periods=100000, freq='s'),
                 'key': range(100000),})
        elif dtype == 'object':
            self.df = DataFrame(
                {'values': (['foo'] * 100000),
                 'key': range(100000)})
        else:
            labels = np.arange(10000).repeat(10)
            data = Series(randn(len(labels)), dtype=dtype)
            data[::3] = np.nan
            data[1::3] = np.nan
            labels = labels.take(np.random.permutation(len(labels)))
            self.df = DataFrame({'values': data, 'key': labels})

    def time_groupby_first(self, dtype):
        self.df.groupby('key').first()

    def time_groupby_last(self, dtype):
        self.df.groupby('key').last()

    def time_groupby_nth_any(self, dtype):
        self.df.groupby('key').nth(0, dropna='all')

    def time_groupby_nth_none(self, dtype):
        self.df.groupby('key').nth(0)


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
        self.obj = np.random.choice(list('ab'), size=self.n).astype(object)
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
# nunique() speed

class groupby_nunique(object):

    def setup(self):
        self.n = 10000
        self.df = DataFrame({'key1': randint(0, 500, size=self.n),
                             'key2': randint(0, 100, size=self.n),
                             'ints': randint(0, 1000, size=self.n),
                             'ints2': randint(0, 1000, size=self.n), })

    def time_groupby_nunique(self):
        self.df.groupby(['key1', 'key2']).nunique()


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

class GroupBySuite(object):
    goal_time = 0.2

    param_names = ['dtype', 'ngroups']
    params = [['int', 'float'], [100, 10000]]

    def setup(self, dtype, ngroups):
        np.random.seed(1234)
        size = ngroups * 2
        rng = np.arange(ngroups)
        values = rng.take(np.random.randint(0, ngroups, size=size))
        if dtype == 'int':
            key = np.random.randint(0, size, size=size)
        else:
            key = np.concatenate([np.random.random(ngroups) * 0.1,
                                  np.random.random(ngroups) * 10.0])

        self.df = DataFrame({'values': values,
                             'key': key})

    def time_all(self, dtype, ngroups):
        self.df.groupby('key')['values'].all()

    def time_any(self, dtype, ngroups):
        self.df.groupby('key')['values'].any()

    def time_count(self, dtype, ngroups):
        self.df.groupby('key')['values'].count()

    def time_cumcount(self, dtype, ngroups):
        self.df.groupby('key')['values'].cumcount()

    def time_cummax(self, dtype, ngroups):
        self.df.groupby('key')['values'].cummax()

    def time_cummin(self, dtype, ngroups):
        self.df.groupby('key')['values'].cummin()

    def time_cumprod(self, dtype, ngroups):
        self.df.groupby('key')['values'].cumprod()

    def time_cumsum(self, dtype, ngroups):
        self.df.groupby('key')['values'].cumsum()

    def time_describe(self, dtype, ngroups):
        self.df.groupby('key')['values'].describe()

    def time_diff(self, dtype, ngroups):
        self.df.groupby('key')['values'].diff()

    def time_first(self, dtype, ngroups):
        self.df.groupby('key')['values'].first()

    def time_head(self, dtype, ngroups):
        self.df.groupby('key')['values'].head()

    def time_last(self, dtype, ngroups):
        self.df.groupby('key')['values'].last()

    def time_mad(self, dtype, ngroups):
        self.df.groupby('key')['values'].mad()

    def time_max(self, dtype, ngroups):
        self.df.groupby('key')['values'].max()

    def time_mean(self, dtype, ngroups):
        self.df.groupby('key')['values'].mean()

    def time_median(self, dtype, ngroups):
        self.df.groupby('key')['values'].median()

    def time_min(self, dtype, ngroups):
        self.df.groupby('key')['values'].min()

    def time_nunique(self, dtype, ngroups):
        self.df.groupby('key')['values'].nunique()

    def time_pct_change(self, dtype, ngroups):
        self.df.groupby('key')['values'].pct_change()

    def time_prod(self, dtype, ngroups):
        self.df.groupby('key')['values'].prod()

    def time_rank(self, dtype, ngroups):
        self.df.groupby('key')['values'].rank()

    def time_sem(self, dtype, ngroups):
        self.df.groupby('key')['values'].sem()

    def time_size(self, dtype, ngroups):
        self.df.groupby('key')['values'].size()

    def time_skew(self, dtype, ngroups):
        self.df.groupby('key')['values'].skew()

    def time_std(self, dtype, ngroups):
        self.df.groupby('key')['values'].std()

    def time_sum(self, dtype, ngroups):
        self.df.groupby('key')['values'].sum()

    def time_tail(self, dtype, ngroups):
        self.df.groupby('key')['values'].tail()

    def time_unique(self, dtype, ngroups):
        self.df.groupby('key')['values'].unique()

    def time_value_counts(self, dtype, ngroups):
        self.df.groupby('key')['values'].value_counts()

    def time_var(self, dtype, ngroups):
        self.df.groupby('key')['values'].var()


class groupby_float32(object):
    # GH 13335
    goal_time = 0.2

    def setup(self):
        tmp1 = (np.random.random(10000) * 0.1).astype(np.float32)
        tmp2 = (np.random.random(10000) * 10.0).astype(np.float32)
        tmp = np.concatenate((tmp1, tmp2))
        arr = np.repeat(tmp, 10)
        self.df = DataFrame(dict(a=arr, b=arr))

    def time_groupby_sum(self):
        self.df.groupby(['a'])['b'].sum()


class groupby_categorical(object):
    goal_time = 0.2

    def setup(self):
        N = 100000
        arr = np.random.random(N)

        self.df = DataFrame(dict(
            a=Categorical(np.random.randint(10000, size=N)),
            b=arr))
        self.df_ordered = DataFrame(dict(
            a=Categorical(np.random.randint(10000, size=N), ordered=True),
            b=arr))
        self.df_extra_cat = DataFrame(dict(
            a=Categorical(np.random.randint(100, size=N),
                          categories=np.arange(10000)),
            b=arr))

    def time_groupby_sort(self):
        self.df.groupby('a')['b'].count()

    def time_groupby_nosort(self):
        self.df.groupby('a', sort=False)['b'].count()

    def time_groupby_ordered_sort(self):
        self.df_ordered.groupby('a')['b'].count()

    def time_groupby_ordered_nosort(self):
        self.df_ordered.groupby('a', sort=False)['b'].count()

    def time_groupby_extra_cat_sort(self):
        self.df_extra_cat.groupby('a')['b'].count()

    def time_groupby_extra_cat_nosort(self):
        self.df_extra_cat.groupby('a', sort=False)['b'].count()


class groupby_period(object):
    # GH 14338
    goal_time = 0.2

    def make_grouper(self, N):
        return pd.period_range('1900-01-01', freq='D', periods=N)

    def setup(self):
        N = 10000
        self.grouper = self.make_grouper(N)
        self.df = pd.DataFrame(np.random.randn(N, 2))

    def time_groupby_sum(self):
        self.df.groupby(self.grouper).sum()


class groupby_datetime(groupby_period):
    def make_grouper(self, N):
        return pd.date_range('1900-01-01', freq='D', periods=N)


class groupby_datetimetz(groupby_period):
    def make_grouper(self, N):
        return pd.date_range('1900-01-01', freq='D', periods=N,
                             tz='US/Central')

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
        self.df = DataFrame({'A': (list(range(self.N)) * 2), 'B': list(range((self.N * 2))), 'C': 1, }).set_index(['A', 'B'])

    def time_groupby_sum_multiindex(self):
        self.df.groupby(level=[0, 1]).sum()


#-------------------------------------------------------------------------------
# Transform testing

class Transform(object):
    goal_time = 0.2

    def setup(self):
        n1 = 400
        n2 = 250

        index = MultiIndex(
            levels=[np.arange(n1), pd.util.testing.makeStringIndex(n2)],
            labels=[[i for i in range(n1) for _ in range(n2)],
                    (list(range(n2)) * n1)],
            names=['lev1', 'lev2'])

        data = DataFrame(np.random.randn(n1 * n2, 3),
                         index=index, columns=['col1', 'col20', 'col3'])
        step = int((n1 * n2 * 0.1))
        for col in range(len(data.columns)):
            idx = col
            while (idx < len(data)):
                data.set_value(data.index[idx], data.columns[col], np.nan)
                idx += step
        self.df = data
        self.f_fillna = (lambda x: x.fillna(method='pad'))

        np.random.seed(2718281)
        n = 20000
        self.df1 = DataFrame(np.random.randint(1, n, (n, 3)),
                             columns=['jim', 'joe', 'jolie'])
        self.df2 = self.df1.copy()
        self.df2['jim'] = self.df2['joe']

        self.df3 = DataFrame(np.random.randint(1, (n / 10), (n, 3)),
                             columns=['jim', 'joe', 'jolie'])
        self.df4 = self.df3.copy()
        self.df4['jim'] = self.df4['joe']

    def time_transform_func(self):
        self.df.groupby(level='lev2').transform(self.f_fillna)

    def time_transform_ufunc(self):
        self.df.groupby(level='lev1').transform(np.max)

    def time_transform_multi_key1(self):
        self.df1.groupby(['jim', 'joe'])['jolie'].transform('max')

    def time_transform_multi_key2(self):
        self.df2.groupby(['jim', 'joe'])['jolie'].transform('max')

    def time_transform_multi_key3(self):
        self.df3.groupby(['jim', 'joe'])['jolie'].transform('max')

    def time_transform_multi_key4(self):
        self.df4.groupby(['jim', 'joe'])['jolie'].transform('max')




np.random.seed(0)
N = 120000
N_TRANSITIONS = 1400
transition_points = np.random.permutation(np.arange(N))[:N_TRANSITIONS]
transition_points.sort()
transitions = np.zeros((N,), dtype=np.bool)
transitions[transition_points] = True
g = transitions.cumsum()
df = DataFrame({'signal': np.random.rand(N), })





class groupby_transform_series(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(0)
        N = 120000
        transition_points = np.sort(np.random.choice(np.arange(N), 1400))
        transitions = np.zeros((N,), dtype=np.bool)
        transitions[transition_points] = True
        self.g = transitions.cumsum()
        self.df = DataFrame({'signal': np.random.rand(N)})

    def time_groupby_transform_series(self):
        self.df['signal'].groupby(self.g).transform(np.mean)


class groupby_transform_series2(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(0)
        self.df = DataFrame({'key': (np.arange(100000) // 3),
                             'val': np.random.randn(100000)})

        self.df_nans = pd.DataFrame({'key': np.repeat(np.arange(1000), 10),
                                     'B': np.nan,
                                     'C': np.nan})
        self.df_nans.ix[4::10, 'B':'C'] = 5

    def time_transform_series2(self):
        self.df.groupby('key')['val'].transform(np.mean)

    def time_cumprod(self):
        self.df.groupby('key').cumprod()

    def time_cumsum(self):
        self.df.groupby('key').cumsum()

    def time_shift(self):
        self.df.groupby('key').shift()

    def time_transform_dataframe(self):
        # GH 12737
        self.df_nans.groupby('key').transform('first')
