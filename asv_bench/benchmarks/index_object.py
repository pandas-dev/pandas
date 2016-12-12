from .pandas_vb_common import *


class SetOperations(object):
    goal_time = 0.2

    def setup(self):
        self.rng = date_range('1/1/2000', periods=10000, freq='T')
        self.rng2 = self.rng[:(-1)]

        # object index with datetime values
        if (self.rng.dtype == object):
            self.idx_rng = self.rng.view(Index)
        else:
            self.idx_rng = self.rng.asobject
        self.idx_rng2 = self.idx_rng[:(-1)]

        # other datetime
        N = 100000
        A = N - 20000
        B = N + 20000
        self.dtidx1 = DatetimeIndex(range(N))
        self.dtidx2 = DatetimeIndex(range(A, B))
        self.dtidx3 = DatetimeIndex(range(N, B))

        # integer
        self.N = 1000000
        self.options = np.arange(self.N)
        self.left = Index(
            self.options.take(np.random.permutation(self.N)[:(self.N // 2)]))
        self.right = Index(
            self.options.take(np.random.permutation(self.N)[:(self.N // 2)]))

        # strings
        N = 10000
        strs = tm.rands_array(10, N)
        self.leftstr = Index(strs[:N * 2 // 3])
        self.rightstr = Index(strs[N // 3:])

    def time_datetime_intersection(self):
        self.rng.intersection(self.rng2)

    def time_datetime_union(self):
        self.rng.union(self.rng2)

    def time_datetime_difference(self):
        self.dtidx1.difference(self.dtidx2)

    def time_datetime_difference_disjoint(self):
        self.dtidx1.difference(self.dtidx3)

    def time_datetime_symmetric_difference(self):
        self.dtidx1.symmetric_difference(self.dtidx2)

    def time_index_datetime_intersection(self):
        self.idx_rng.intersection(self.idx_rng2)

    def time_index_datetime_union(self):
        self.idx_rng.union(self.idx_rng2)

    def time_int64_intersection(self):
        self.left.intersection(self.right)

    def time_int64_union(self):
        self.left.union(self.right)

    def time_int64_difference(self):
        self.left.difference(self.right)

    def time_int64_symmetric_difference(self):
        self.left.symmetric_difference(self.right)

    def time_str_difference(self):
        self.leftstr.difference(self.rightstr)

    def time_str_symmetric_difference(self):
        self.leftstr.symmetric_difference(self.rightstr)


class Datetime(object):
    goal_time = 0.2

    def setup(self):
        self.dr = pd.date_range('20000101', freq='D', periods=10000)

    def time_is_dates_only(self):
        self.dr._is_dates_only


class Float64(object):
    goal_time = 0.2

    def setup(self):
        self.idx = tm.makeFloatIndex(1000000)
        self.mask = ((np.arange(self.idx.size) % 3) == 0)
        self.series_mask = Series(self.mask)

        self.baseidx = np.arange(1000000.0)

    def time_boolean_indexer(self):
        self.idx[self.mask]

    def time_boolean_series_indexer(self):
        self.idx[self.series_mask]

    def time_construct(self):
        Index(self.baseidx)

    def time_div(self):
        (self.idx / 2)

    def time_get(self):
        self.idx[1]

    def time_mul(self):
        (self.idx * 2)

    def time_slice_indexer_basic(self):
        self.idx[:(-1)]

    def time_slice_indexer_even(self):
        self.idx[::2]


class StringIndex(object):
    goal_time = 0.2

    def setup(self):
        self.idx = tm.makeStringIndex(1000000)
        self.mask = ((np.arange(1000000) % 3) == 0)
        self.series_mask = Series(self.mask)

    def time_boolean_indexer(self):
        self.idx[self.mask]

    def time_boolean_series_indexer(self):
        self.idx[self.series_mask]

    def time_slice_indexer_basic(self):
        self.idx[:(-1)]

    def time_slice_indexer_even(self):
        self.idx[::2]


class Multi1(object):
    goal_time = 0.2

    def setup(self):
        (n, k) = (200, 5000)
        self.levels = [np.arange(n), tm.makeStringIndex(n).values, (1000 + np.arange(n))]
        self.labels = [np.random.choice(n, (k * n)) for lev in self.levels]
        self.mi = MultiIndex(levels=self.levels, labels=self.labels)

        self.iterables = [tm.makeStringIndex(10000), range(20)]

    def time_duplicated(self):
        self.mi.duplicated()

    def time_from_product(self):
        MultiIndex.from_product(self.iterables)


class Multi2(object):
    goal_time = 0.2

    def setup(self):
        self.n = ((((3 * 5) * 7) * 11) * (1 << 10))
        (low, high) = (((-1) << 12), (1 << 12))
        self.f = (lambda k: np.repeat(np.random.randint(low, high, (self.n // k)), k))
        self.i = np.random.permutation(self.n)
        self.mi = MultiIndex.from_arrays([self.f(11), self.f(7), self.f(5), self.f(3), self.f(1)])[self.i]

        self.a = np.repeat(np.arange(100), 1000)
        self.b = np.tile(np.arange(1000), 100)
        self.midx2 = MultiIndex.from_arrays([self.a, self.b])
        self.midx2 = self.midx2.take(np.random.permutation(np.arange(100000)))

    def time_sortlevel_int64(self):
        self.mi.sortlevel()

    def time_sortlevel_zero(self):
        self.midx2.sortlevel(0)

    def time_sortlevel_one(self):
        self.midx2.sortlevel(1)


class Multi3(object):
    goal_time = 0.2

    def setup(self):
        self.level1 = range(1000)
        self.level2 = date_range(start='1/1/2012', periods=100)
        self.mi = MultiIndex.from_product([self.level1, self.level2])

    def time_datetime_level_values_full(self):
        self.mi.copy().values

    def time_datetime_level_values_sliced(self):
        self.mi[:10].values
