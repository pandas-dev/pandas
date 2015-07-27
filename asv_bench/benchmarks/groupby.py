from pandas_vb_common import *
from itertools import product
from string import ascii_letters, digits


class groupby_agg_builtins1(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(27182)
        self.n = 100000
        self.df = DataFrame(np.random.randint(1, (self.n / 100), (self.n, 3)), columns=['jim', 'joe', 'jolie'])

    def time_groupby_agg_builtins1(self):
        self.df.groupby('jim').agg([sum, min, max])


class groupby_agg_builtins2(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(27182)
        self.n = 100000
        self.df = DataFrame(np.random.randint(1, (self.n / 100), (self.n, 3)), columns=['jim', 'joe', 'jolie'])

    def time_groupby_agg_builtins2(self):
        self.df.groupby(['jim', 'joe']).agg([sum, min, max])


class groupby_apply_dict_return(object):
    goal_time = 0.2

    def setup(self):
        self.labels = np.arange(1000).repeat(10)
        self.data = Series(randn(len(self.labels)))
        self.f = (lambda x: {'first': x.values[0], 'last': x.values[(-1)], })

    def time_groupby_apply_dict_return(self):
        self.data.groupby(self.labels).apply(self.f)


class groupby_dt_size(object):
    goal_time = 0.2

    def setup(self):
        self.n = 100000
        self.offsets = np.random.randint(self.n, size=self.n).astype('timedelta64[ns]')
        self.dates = (np.datetime64('now') + self.offsets)
        self.df = DataFrame({'key1': np.random.randint(0, 500, size=self.n), 'key2': np.random.randint(0, 100, size=self.n), 'value1': np.random.randn(self.n), 'value2': np.random.randn(self.n), 'value3': np.random.randn(self.n), 'dates': self.dates, })

    def time_groupby_dt_size(self):
        self.df.groupby(['dates']).size()


class groupby_dt_timegrouper_size(object):
    goal_time = 0.2

    def setup(self):
        self.n = 100000
        self.offsets = np.random.randint(self.n, size=self.n).astype('timedelta64[ns]')
        self.dates = (np.datetime64('now') + self.offsets)
        self.df = DataFrame({'key1': np.random.randint(0, 500, size=self.n), 'key2': np.random.randint(0, 100, size=self.n), 'value1': np.random.randn(self.n), 'value2': np.random.randn(self.n), 'value3': np.random.randn(self.n), 'dates': self.dates, })

    def time_groupby_dt_timegrouper_size(self):
        self.df.groupby(TimeGrouper(key='dates', freq='M')).size()


class groupby_first_datetimes(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'a': date_range('1/1/2011', periods=100000, freq='s'), 'b': range(100000), })

    def time_groupby_first_datetimes(self):
        self.df.groupby('b').first()


class groupby_first_float32(object):
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


class groupby_first_float64(object):
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

    def time_groupby_first_float64(self):
        self.data.groupby(self.labels).first()


class groupby_first_object(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'a': (['foo'] * 100000), 'b': range(100000), })

    def time_groupby_first_object(self):
        self.df.groupby('b').first()


class groupby_frame_apply(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.labels = np.random.randint(0, 2000, size=self.N)
        self.labels2 = np.random.randint(0, 3, size=self.N)
        self.df = DataFrame({'key': self.labels, 'key2': self.labels2, 'value1': randn(self.N), 'value2': (['foo', 'bar', 'baz', 'qux'] * (self.N / 4)), })

        def f(g):
            return 1

    def time_groupby_frame_apply(self):
        self.df.groupby(['key', 'key2']).apply(f)


class groupby_frame_apply_overhead(object):
    goal_time = 0.2

    def setup(self):
        self.N = 10000
        self.labels = np.random.randint(0, 2000, size=self.N)
        self.labels2 = np.random.randint(0, 3, size=self.N)
        self.df = DataFrame({'key': self.labels, 'key2': self.labels2, 'value1': randn(self.N), 'value2': (['foo', 'bar', 'baz', 'qux'] * (self.N / 4)), })

        def f(g):
            return 1

    def time_groupby_frame_apply_overhead(self):
        self.df.groupby('key').apply(f)


class groupby_frame_cython_many_columns(object):
    goal_time = 0.2

    def setup(self):
        self.labels = np.random.randint(0, 100, size=1000)
        self.df = DataFrame(randn(1000, 1000))

    def time_groupby_frame_cython_many_columns(self):
        self.df.groupby(self.labels).sum()


class groupby_frame_median(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(100000, 2)
        self.labels = np.random.randint(0, 1000, size=100000)
        self.df = DataFrame(self.data)

    def time_groupby_frame_median(self):
        self.df.groupby(self.labels).median()


class groupby_frame_nth_any(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randint(1, 100, (10000, 2)))

    def time_groupby_frame_nth_any(self):
        self.df.groupby(0).nth(0, dropna='any')


class groupby_frame_nth_none(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randint(1, 100, (10000, 2)))

    def time_groupby_frame_nth_none(self):
        self.df.groupby(0).nth(0)


class groupby_frame_singlekey_integer(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(100000, 1)
        self.labels = np.random.randint(0, 1000, size=100000)
        self.df = DataFrame(self.data)

    def time_groupby_frame_singlekey_integer(self):
        self.df.groupby(self.labels).sum()


class groupby_indices(object):
    goal_time = 0.2

    def setup(self):
        try:
            self.rng = date_range('1/1/2000', '12/31/2005', freq='H')
            (year, month, day) = (self.rng.year, self.rng.month, self.rng.day)
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


class groupby_int_count(object):
    goal_time = 0.2

    def setup(self):
        self.n = 10000
        self.df = DataFrame({'key1': randint(0, 500, size=self.n), 'key2': randint(0, 100, size=self.n), 'ints': randint(0, 1000, size=self.n), 'ints2': randint(0, 1000, size=self.n), })

    def time_groupby_int_count(self):
        self.df.groupby(['key1', 'key2']).count()


class groupby_last_datetimes(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'a': date_range('1/1/2011', periods=100000, freq='s'), 'b': range(100000), })

    def time_groupby_last_datetimes(self):
        self.df.groupby('b').last()


class groupby_last_float32(object):
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

    def time_groupby_last_float32(self):
        self.data2.groupby(self.labels).last()


class groupby_last_float64(object):
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

    def time_groupby_last_float64(self):
        self.data.groupby(self.labels).last()


class groupby_last_object(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'a': (['foo'] * 100000), 'b': range(100000), })

    def time_groupby_last_object(self):
        self.df.groupby('b').last()


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
        self.df = DataFrame({'key1': np.random.randint(0, 500, size=self.n), 'key2': np.random.randint(0, 100, size=self.n), 'dates': self.dates, 'value2': self.value2, 'value3': np.random.randn(self.n), 'ints': np.random.randint(0, 1000, size=self.n), 'obj': self.obj, 'offsets': self.offsets, })

    def time_groupby_multi_count(self):
        self.df.groupby(['key1', 'key2']).count()


class groupby_multi_cython(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.ngroups = 100

        def get_test_data(ngroups=100, n=self.N):
            self.unique_groups = range(self.ngroups)
            self.arr = np.asarray(np.tile(self.unique_groups, (n / self.ngroups)), dtype=object)
            if (len(self.arr) < n):
                self.arr = np.asarray((list(self.arr) + self.unique_groups[:(n - len(self.arr))]), dtype=object)
            random.shuffle(self.arr)
            return self.arr
        self.df = DataFrame({'key1': get_test_data(ngroups=self.ngroups), 'key2': get_test_data(ngroups=self.ngroups), 'data1': np.random.randn(self.N), 'data2': np.random.randn(self.N), })

        def f():
            self.df.groupby(['key1', 'key2']).agg((lambda x: x.values.sum()))
        self.simple_series = Series(np.random.randn(self.N))
        self.key1 = self.df['key1']

    def time_groupby_multi_cython(self):
        self.df.groupby(['key1', 'key2']).sum()


class groupby_multi_different_functions(object):
    goal_time = 0.2

    def setup(self):
        self.fac1 = np.array(['A', 'B', 'C'], dtype='O')
        self.fac2 = np.array(['one', 'two'], dtype='O')
        self.df = DataFrame({'key1': self.fac1.take(np.random.randint(0, 3, size=100000)), 'key2': self.fac2.take(np.random.randint(0, 2, size=100000)), 'value1': np.random.randn(100000), 'value2': np.random.randn(100000), 'value3': np.random.randn(100000), })

    def time_groupby_multi_different_functions(self):
        self.df.groupby(['key1', 'key2']).agg({'value1': 'mean', 'value2': 'var', 'value3': 'sum', })


class groupby_multi_different_numpy_functions(object):
    goal_time = 0.2

    def setup(self):
        self.fac1 = np.array(['A', 'B', 'C'], dtype='O')
        self.fac2 = np.array(['one', 'two'], dtype='O')
        self.df = DataFrame({'key1': self.fac1.take(np.random.randint(0, 3, size=100000)), 'key2': self.fac2.take(np.random.randint(0, 2, size=100000)), 'value1': np.random.randn(100000), 'value2': np.random.randn(100000), 'value3': np.random.randn(100000), })

    def time_groupby_multi_different_numpy_functions(self):
        self.df.groupby(['key1', 'key2']).agg({'value1': np.mean, 'value2': np.var, 'value3': np.sum, })


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


class groupby_multi_python(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.ngroups = 100

        def get_test_data(ngroups=100, n=self.N):
            self.unique_groups = range(self.ngroups)
            self.arr = np.asarray(np.tile(self.unique_groups, (n / self.ngroups)), dtype=object)
            if (len(self.arr) < n):
                self.arr = np.asarray((list(self.arr) + self.unique_groups[:(n - len(self.arr))]), dtype=object)
            random.shuffle(self.arr)
            return self.arr
        self.df = DataFrame({'key1': get_test_data(ngroups=self.ngroups), 'key2': get_test_data(ngroups=self.ngroups), 'data1': np.random.randn(self.N), 'data2': np.random.randn(self.N), })

        def f():
            self.df.groupby(['key1', 'key2']).agg((lambda x: x.values.sum()))
        self.simple_series = Series(np.random.randn(self.N))
        self.key1 = self.df['key1']

    def time_groupby_multi_python(self):
        self.df.groupby(['key1', 'key2'])['data1'].agg((lambda x: x.values.sum()))


class groupby_multi_series_op(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.ngroups = 100

        def get_test_data(ngroups=100, n=self.N):
            self.unique_groups = range(self.ngroups)
            self.arr = np.asarray(np.tile(self.unique_groups, (n / self.ngroups)), dtype=object)
            if (len(self.arr) < n):
                self.arr = np.asarray((list(self.arr) + self.unique_groups[:(n - len(self.arr))]), dtype=object)
            random.shuffle(self.arr)
            return self.arr
        self.df = DataFrame({'key1': get_test_data(ngroups=self.ngroups), 'key2': get_test_data(ngroups=self.ngroups), 'data1': np.random.randn(self.N), 'data2': np.random.randn(self.N), })

        def f():
            self.df.groupby(['key1', 'key2']).agg((lambda x: x.values.sum()))
        self.simple_series = Series(np.random.randn(self.N))
        self.key1 = self.df['key1']

    def time_groupby_multi_series_op(self):
        self.df.groupby(['key1', 'key2'])['data1'].agg(np.std)


class groupby_multi_size(object):
    goal_time = 0.2

    def setup(self):
        self.n = 100000
        self.offsets = np.random.randint(self.n, size=self.n).astype('timedelta64[ns]')
        self.dates = (np.datetime64('now') + self.offsets)
        self.df = DataFrame({'key1': np.random.randint(0, 500, size=self.n), 'key2': np.random.randint(0, 100, size=self.n), 'value1': np.random.randn(self.n), 'value2': np.random.randn(self.n), 'value3': np.random.randn(self.n), 'dates': self.dates, })

    def time_groupby_multi_size(self):
        self.df.groupby(['key1', 'key2']).size()


class groupby_ngroups_10000_all(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_all(self):
        self.df.groupby('value')['timestamp'].all()


class groupby_ngroups_10000_any(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_any(self):
        self.df.groupby('value')['timestamp'].any()


class groupby_ngroups_10000_count(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_count(self):
        self.df.groupby('value')['timestamp'].count()


class groupby_ngroups_10000_cumcount(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_cumcount(self):
        self.df.groupby('value')['timestamp'].cumcount()


class groupby_ngroups_10000_cummax(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_cummax(self):
        self.df.groupby('value')['timestamp'].cummax()


class groupby_ngroups_10000_cummin(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_cummin(self):
        self.df.groupby('value')['timestamp'].cummin()


class groupby_ngroups_10000_cumprod(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_cumprod(self):
        self.df.groupby('value')['timestamp'].cumprod()


class groupby_ngroups_10000_cumsum(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_cumsum(self):
        self.df.groupby('value')['timestamp'].cumsum()


class groupby_ngroups_10000_describe(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_describe(self):
        self.df.groupby('value')['timestamp'].describe()


class groupby_ngroups_10000_diff(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_diff(self):
        self.df.groupby('value')['timestamp'].diff()


class groupby_ngroups_10000_first(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_first(self):
        self.df.groupby('value')['timestamp'].first()


class groupby_ngroups_10000_head(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_head(self):
        self.df.groupby('value')['timestamp'].head()


class groupby_ngroups_10000_last(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_last(self):
        self.df.groupby('value')['timestamp'].last()


class groupby_ngroups_10000_mad(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_mad(self):
        self.df.groupby('value')['timestamp'].mad()


class groupby_ngroups_10000_max(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_max(self):
        self.df.groupby('value')['timestamp'].max()


class groupby_ngroups_10000_mean(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_mean(self):
        self.df.groupby('value')['timestamp'].mean()


class groupby_ngroups_10000_median(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_median(self):
        self.df.groupby('value')['timestamp'].median()


class groupby_ngroups_10000_min(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_min(self):
        self.df.groupby('value')['timestamp'].min()


class groupby_ngroups_10000_nunique(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_nunique(self):
        self.df.groupby('value')['timestamp'].nunique()


class groupby_ngroups_10000_pct_change(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_pct_change(self):
        self.df.groupby('value')['timestamp'].pct_change()


class groupby_ngroups_10000_prod(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_prod(self):
        self.df.groupby('value')['timestamp'].prod()


class groupby_ngroups_10000_rank(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_rank(self):
        self.df.groupby('value')['timestamp'].rank()


class groupby_ngroups_10000_sem(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_sem(self):
        self.df.groupby('value')['timestamp'].sem()


class groupby_ngroups_10000_size(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_size(self):
        self.df.groupby('value')['timestamp'].size()


class groupby_ngroups_10000_skew(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_skew(self):
        self.df.groupby('value')['timestamp'].skew()


class groupby_ngroups_10000_std(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_std(self):
        self.df.groupby('value')['timestamp'].std()


class groupby_ngroups_10000_sum(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_sum(self):
        self.df.groupby('value')['timestamp'].sum()


class groupby_ngroups_10000_tail(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_tail(self):
        self.df.groupby('value')['timestamp'].tail()


class groupby_ngroups_10000_unique(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_unique(self):
        self.df.groupby('value')['timestamp'].unique()


class groupby_ngroups_10000_value_counts(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_value_counts(self):
        self.df.groupby('value')['timestamp'].value_counts()


class groupby_ngroups_10000_var(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 10000
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_10000_var(self):
        self.df.groupby('value')['timestamp'].var()


class groupby_ngroups_100_all(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_all(self):
        self.df.groupby('value')['timestamp'].all()


class groupby_ngroups_100_any(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_any(self):
        self.df.groupby('value')['timestamp'].any()


class groupby_ngroups_100_count(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_count(self):
        self.df.groupby('value')['timestamp'].count()


class groupby_ngroups_100_cumcount(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_cumcount(self):
        self.df.groupby('value')['timestamp'].cumcount()


class groupby_ngroups_100_cummax(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_cummax(self):
        self.df.groupby('value')['timestamp'].cummax()


class groupby_ngroups_100_cummin(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_cummin(self):
        self.df.groupby('value')['timestamp'].cummin()


class groupby_ngroups_100_cumprod(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_cumprod(self):
        self.df.groupby('value')['timestamp'].cumprod()


class groupby_ngroups_100_cumsum(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_cumsum(self):
        self.df.groupby('value')['timestamp'].cumsum()


class groupby_ngroups_100_describe(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_describe(self):
        self.df.groupby('value')['timestamp'].describe()


class groupby_ngroups_100_diff(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_diff(self):
        self.df.groupby('value')['timestamp'].diff()


class groupby_ngroups_100_first(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_first(self):
        self.df.groupby('value')['timestamp'].first()


class groupby_ngroups_100_head(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_head(self):
        self.df.groupby('value')['timestamp'].head()


class groupby_ngroups_100_last(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_last(self):
        self.df.groupby('value')['timestamp'].last()


class groupby_ngroups_100_mad(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_mad(self):
        self.df.groupby('value')['timestamp'].mad()


class groupby_ngroups_100_max(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_max(self):
        self.df.groupby('value')['timestamp'].max()


class groupby_ngroups_100_mean(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_mean(self):
        self.df.groupby('value')['timestamp'].mean()


class groupby_ngroups_100_median(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_median(self):
        self.df.groupby('value')['timestamp'].median()


class groupby_ngroups_100_min(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_min(self):
        self.df.groupby('value')['timestamp'].min()


class groupby_ngroups_100_nunique(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_nunique(self):
        self.df.groupby('value')['timestamp'].nunique()


class groupby_ngroups_100_pct_change(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_pct_change(self):
        self.df.groupby('value')['timestamp'].pct_change()


class groupby_ngroups_100_prod(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_prod(self):
        self.df.groupby('value')['timestamp'].prod()


class groupby_ngroups_100_rank(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_rank(self):
        self.df.groupby('value')['timestamp'].rank()


class groupby_ngroups_100_sem(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_sem(self):
        self.df.groupby('value')['timestamp'].sem()


class groupby_ngroups_100_size(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_size(self):
        self.df.groupby('value')['timestamp'].size()


class groupby_ngroups_100_skew(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_skew(self):
        self.df.groupby('value')['timestamp'].skew()


class groupby_ngroups_100_std(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_std(self):
        self.df.groupby('value')['timestamp'].std()


class groupby_ngroups_100_sum(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_sum(self):
        self.df.groupby('value')['timestamp'].sum()


class groupby_ngroups_100_tail(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_tail(self):
        self.df.groupby('value')['timestamp'].tail()


class groupby_ngroups_100_unique(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_unique(self):
        self.df.groupby('value')['timestamp'].unique()


class groupby_ngroups_100_value_counts(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_value_counts(self):
        self.df.groupby('value')['timestamp'].value_counts()


class groupby_ngroups_100_var(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1234)
        self.ngroups = 100
        self.size = (self.ngroups * 2)
        self.rng = np.arange(self.ngroups)
        self.df = DataFrame(dict(timestamp=self.rng.take(np.random.randint(0, self.ngroups, size=self.size)), value=np.random.randint(0, self.size, size=self.size)))

    def time_groupby_ngroups_100_var(self):
        self.df.groupby('value')['timestamp'].var()


class groupby_nth_datetimes_any(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'a': date_range('1/1/2011', periods=100000, freq='s'), 'b': range(100000), })

    def time_groupby_nth_datetimes_any(self):
        self.df.groupby('b').nth(0, dropna='all')


class groupby_nth_datetimes_none(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'a': date_range('1/1/2011', periods=100000, freq='s'), 'b': range(100000), })

    def time_groupby_nth_datetimes_none(self):
        self.df.groupby('b').nth(0)


class groupby_nth_float32_any(object):
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

    def time_groupby_nth_float32_any(self):
        self.data2.groupby(self.labels).nth(0, dropna='all')


class groupby_nth_float32_none(object):
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

    def time_groupby_nth_float32_none(self):
        self.data2.groupby(self.labels).nth(0)


class groupby_nth_float64_any(object):
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

    def time_groupby_nth_float64_any(self):
        self.data.groupby(self.labels).nth(0, dropna='all')


class groupby_nth_float64_none(object):
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

    def time_groupby_nth_float64_none(self):
        self.data.groupby(self.labels).nth(0)


class groupby_nth_object_any(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'a': (['foo'] * 100000), 'b': range(100000), })

    def time_groupby_nth_object_any(self):
        self.df.groupby('b').nth(0, dropna='any')


class groupby_nth_object_none(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame({'a': (['foo'] * 100000), 'b': range(100000), })

    def time_groupby_nth_object_none(self):
        self.df.groupby('b').nth(0)


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


class groupby_series_nth_any(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randint(1, 100, (10000, 2)))

    def time_groupby_series_nth_any(self):
        self.df[1].groupby(self.df[0]).nth(0, dropna='any')


class groupby_series_nth_none(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randint(1, 100, (10000, 2)))

    def time_groupby_series_nth_none(self):
        self.df[1].groupby(self.df[0]).nth(0)


class groupby_series_simple_cython(object):
    goal_time = 0.2

    def setup(self):
        self.N = 100000
        self.ngroups = 100

        def get_test_data(ngroups=100, n=self.N):
            self.unique_groups = range(self.ngroups)
            self.arr = np.asarray(np.tile(self.unique_groups, (n / self.ngroups)), dtype=object)
            if (len(self.arr) < n):
                self.arr = np.asarray((list(self.arr) + self.unique_groups[:(n - len(self.arr))]), dtype=object)
            random.shuffle(self.arr)
            return self.arr
        self.df = DataFrame({'key1': get_test_data(ngroups=self.ngroups), 'key2': get_test_data(ngroups=self.ngroups), 'data1': np.random.randn(self.N), 'data2': np.random.randn(self.N), })

        def f():
            self.df.groupby(['key1', 'key2']).agg((lambda x: x.values.sum()))
        self.simple_series = Series(np.random.randn(self.N))
        self.key1 = self.df['key1']

    def time_groupby_series_simple_cython(self):
        self.df.groupby('key1').rank(pct=True)


class groupby_simple_compress_timing(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(1000000, 2)
        self.labels = np.random.randint(0, 1000, size=1000000)
        self.df = DataFrame(self.data)

    def time_groupby_simple_compress_timing(self):
        self.df.groupby(self.labels).mean()


class groupby_sum_booleans(object):
    goal_time = 0.2

    def setup(self):
        self.N = 500
        self.df = DataFrame({'ii': range(self.N), 'bb': [True for x in range(self.N)], })

    def time_groupby_sum_booleans(self):
        self.df.groupby('ii').sum()


class groupby_sum_multiindex(object):
    goal_time = 0.2

    def setup(self):
        self.N = 50
        self.df = DataFrame({'A': (range(self.N) * 2), 'B': range((self.N * 2)), 'C': 1, }).set_index(['A', 'B'])

    def time_groupby_sum_multiindex(self):
        self.df.groupby(level=[0, 1]).sum()


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
        self.data_index = MultiIndex(levels=[self.dates.values, self.security_ids], labels=[[i for i in xrange(self.n_dates) for _ in xrange(self.n_securities)], (range(self.n_securities) * self.n_dates)], names=['date', 'security_id'])
        self.n_data = len(self.data_index)
        self.columns = Index(['factor{}'.format(i) for i in xrange(1, (self.n_columns + 1))])
        self.data = DataFrame(np.random.randn(self.n_data, self.n_columns), index=self.data_index, columns=self.columns)
        self.step = int((self.n_data * self.share_na))
        for column_index in xrange(self.n_columns):
            self.index = column_index
            while (self.index < self.n_data):
                self.data.set_value(self.data_index[self.index], self.columns[column_index], np.nan)
                self.index += self.step
        self.f_fillna = (lambda x: x.fillna(method='pad'))

    def time_groupby_transform(self):
        self.data.groupby(level='security_id').transform(self.f_fillna)


class groupby_transform_multi_key1(object):
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


class groupby_transform_ufunc(object):
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
        self.data_index = MultiIndex(levels=[self.dates.values, self.security_ids], labels=[[i for i in xrange(self.n_dates) for _ in xrange(self.n_securities)], (range(self.n_securities) * self.n_dates)], names=['date', 'security_id'])
        self.n_data = len(self.data_index)
        self.columns = Index(['factor{}'.format(i) for i in xrange(1, (self.n_columns + 1))])
        self.data = DataFrame(np.random.randn(self.n_data, self.n_columns), index=self.data_index, columns=self.columns)
        self.step = int((self.n_data * self.share_na))
        for column_index in xrange(self.n_columns):
            self.index = column_index
            while (self.index < self.n_data):
                self.data.set_value(self.data_index[self.index], self.columns[column_index], np.nan)
                self.index += self.step
        self.f_fillna = (lambda x: x.fillna(method='pad'))

    def time_groupby_transform_ufunc(self):
        self.data.groupby(level='date').transform(np.max)


class series_value_counts_int64(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(np.random.randint(0, 1000, size=100000))

    def time_series_value_counts_int64(self):
        self.s.value_counts()


class series_value_counts_strings(object):
    goal_time = 0.2

    def setup(self):
        self.K = 1000
        self.N = 100000
        self.uniques = tm.makeStringIndex(self.K).values
        self.s = Series(np.tile(self.uniques, (self.N // self.K)))

    def time_series_value_counts_strings(self):
        self.s.value_counts()