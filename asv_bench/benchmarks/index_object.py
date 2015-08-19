from pandas_vb_common import *


class datetime_index_intersection(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=10000, freq='T')
        self.rng2 = self.rng[:(-1)]

    def time_datetime_index_intersection(self):
        self.rng.intersection(self.rng2)


class datetime_index_repr(object):
    goal_time = 0.2

    def setup(self):
        self.dr = pd.date_range('20000101', freq='D', periods=100000)

    def time_datetime_index_repr(self):
        self.dr._is_dates_only


class datetime_index_union(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=10000, freq='T')
        self.rng2 = self.rng[:(-1)]

    def time_datetime_index_union(self):
        self.rng.union(self.rng2)


class index_datetime_intersection(object):
    goal_time = 0.2

    def setup(self):
        self.rng = DatetimeIndex(start='1/1/2000', periods=10000, freq=datetools.Minute())
        if (self.rng.dtype == object):
            self.rng = self.rng.view(Index)
        else:
            self.rng = self.rng.asobject
        self.rng2 = self.rng[:(-1)]

    def time_index_datetime_intersection(self):
        self.rng.intersection(self.rng2)


class index_datetime_union(object):
    goal_time = 0.2

    def setup(self):
        self.rng = DatetimeIndex(start='1/1/2000', periods=10000, freq=datetools.Minute())
        if (self.rng.dtype == object):
            self.rng = self.rng.view(Index)
        else:
            self.rng = self.rng.asobject
        self.rng2 = self.rng[:(-1)]

    def time_index_datetime_union(self):
        self.rng.union(self.rng2)


class index_float64_boolean_indexer(object):
    goal_time = 0.2

    def setup(self):
        self.idx = tm.makeFloatIndex(1000000)
        self.mask = ((np.arange(self.idx.size) % 3) == 0)
        self.series_mask = Series(self.mask)

    def time_index_float64_boolean_indexer(self):
        self.idx[self.mask]


class index_float64_boolean_series_indexer(object):
    goal_time = 0.2

    def setup(self):
        self.idx = tm.makeFloatIndex(1000000)
        self.mask = ((np.arange(self.idx.size) % 3) == 0)
        self.series_mask = Series(self.mask)

    def time_index_float64_boolean_series_indexer(self):
        self.idx[self.series_mask]


class index_float64_construct(object):
    goal_time = 0.2

    def setup(self):
        self.baseidx = np.arange(1000000.0)

    def time_index_float64_construct(self):
        Index(self.baseidx)


class index_float64_div(object):
    goal_time = 0.2

    def setup(self):
        self.idx = tm.makeFloatIndex(1000000)
        self.mask = ((np.arange(self.idx.size) % 3) == 0)
        self.series_mask = Series(self.mask)

    def time_index_float64_div(self):
        (self.idx / 2)


class index_float64_get(object):
    goal_time = 0.2

    def setup(self):
        self.idx = tm.makeFloatIndex(1000000)
        self.mask = ((np.arange(self.idx.size) % 3) == 0)
        self.series_mask = Series(self.mask)

    def time_index_float64_get(self):
        self.idx[1]


class index_float64_mul(object):
    goal_time = 0.2

    def setup(self):
        self.idx = tm.makeFloatIndex(1000000)
        self.mask = ((np.arange(self.idx.size) % 3) == 0)
        self.series_mask = Series(self.mask)

    def time_index_float64_mul(self):
        (self.idx * 2)


class index_float64_slice_indexer_basic(object):
    goal_time = 0.2

    def setup(self):
        self.idx = tm.makeFloatIndex(1000000)
        self.mask = ((np.arange(self.idx.size) % 3) == 0)
        self.series_mask = Series(self.mask)

    def time_index_float64_slice_indexer_basic(self):
        self.idx[:(-1)]


class index_float64_slice_indexer_even(object):
    goal_time = 0.2

    def setup(self):
        self.idx = tm.makeFloatIndex(1000000)
        self.mask = ((np.arange(self.idx.size) % 3) == 0)
        self.series_mask = Series(self.mask)

    def time_index_float64_slice_indexer_even(self):
        self.idx[::2]


class index_int64_intersection(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.options = np.arange(self.N)
        self.left = Index(self.options.take(np.random.permutation(self.N)[:(self.N // 2)]))
        self.right = Index(self.options.take(np.random.permutation(self.N)[:(self.N // 2)]))

    def time_index_int64_intersection(self):
        self.left.intersection(self.right)


class index_int64_union(object):
    goal_time = 0.2

    def setup(self):
        self.N = 1000000
        self.options = np.arange(self.N)
        self.left = Index(self.options.take(np.random.permutation(self.N)[:(self.N // 2)]))
        self.right = Index(self.options.take(np.random.permutation(self.N)[:(self.N // 2)]))

    def time_index_int64_union(self):
        self.left.union(self.right)


class index_str_boolean_indexer(object):
    goal_time = 0.2

    def setup(self):
        self.idx = tm.makeStringIndex(1000000)
        self.mask = ((np.arange(1000000) % 3) == 0)
        self.series_mask = Series(self.mask)

    def time_index_str_boolean_indexer(self):
        self.idx[self.mask]


class index_str_boolean_series_indexer(object):
    goal_time = 0.2

    def setup(self):
        self.idx = tm.makeStringIndex(1000000)
        self.mask = ((np.arange(1000000) % 3) == 0)
        self.series_mask = Series(self.mask)

    def time_index_str_boolean_series_indexer(self):
        self.idx[self.series_mask]


class index_str_slice_indexer_basic(object):
    goal_time = 0.2

    def setup(self):
        self.idx = tm.makeStringIndex(1000000)
        self.mask = ((np.arange(1000000) % 3) == 0)
        self.series_mask = Series(self.mask)

    def time_index_str_slice_indexer_basic(self):
        self.idx[:(-1)]


class index_str_slice_indexer_even(object):
    goal_time = 0.2

    def setup(self):
        self.idx = tm.makeStringIndex(1000000)
        self.mask = ((np.arange(1000000) % 3) == 0)
        self.series_mask = Series(self.mask)

    def time_index_str_slice_indexer_even(self):
        self.idx[::2]


class multiindex_duplicated(object):
    goal_time = 0.2

    def setup(self):
        (n, k) = (200, 5000)
        self.levels = [np.arange(n), tm.makeStringIndex(n).values, (1000 + np.arange(n))]
        self.labels = [np.random.choice(n, (k * n)) for lev in self.levels]
        self.mi = MultiIndex(levels=self.levels, labels=self.labels)

    def time_multiindex_duplicated(self):
        self.mi.duplicated()


class multiindex_from_product(object):
    goal_time = 0.2

    def setup(self):
        self.iterables = [tm.makeStringIndex(10000), xrange(20)]

    def time_multiindex_from_product(self):
        MultiIndex.from_product(self.iterables)


class multiindex_sortlevel_int64(object):
    goal_time = 0.2

    def setup(self):
        self.n = ((((3 * 5) * 7) * 11) * (1 << 10))
        (low, high) = (((-1) << 12), (1 << 12))
        self.f = (lambda k: np.repeat(np.random.randint(low, high, (self.n // k)), k))
        self.i = np.random.permutation(self.n)
        self.mi = MultiIndex.from_arrays([self.f(11), self.f(7), self.f(5), self.f(3), self.f(1)])[self.i]

    def time_multiindex_sortlevel_int64(self):
        self.mi.sortlevel()


class multiindex_with_datetime_level_full(object):
    goal_time = 0.2

    def setup(self):
        self.level1 = range(1000)
        self.level2 = date_range(start='1/1/2012', periods=100)
        self.mi = MultiIndex.from_product([self.level1, self.level2])

    def time_multiindex_with_datetime_level_full(self):
        self.mi.copy().values


class multiindex_with_datetime_level_sliced(object):
    goal_time = 0.2

    def setup(self):
        self.level1 = range(1000)
        self.level2 = date_range(start='1/1/2012', periods=100)
        self.mi = MultiIndex.from_product([self.level1, self.level2])

    def time_multiindex_with_datetime_level_sliced(self):
        self.mi[:10].values