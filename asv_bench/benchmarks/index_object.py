import numpy as np
import pandas.util.testing as tm
from pandas import date_range, DatetimeIndex, Index, MultiIndex, RangeIndex

from .pandas_vb_common import setup  # noqa


class SetOperations(object):

    goal_time = 0.2

    def setup(self):
        self.dates_left = date_range('1/1/2000', periods=10000, freq='T')
        self.dates_right = self.dates_left[:(-1)]

        fmt = '%Y-%m-%d %H:%M:%S'
        self.date_str_left = Index(self.dates_left.strftime(fmt))
        self.date_str_right = self.date_str_left[:-1]

        # other datetime
        N = 100000
        A = N - 20000
        B = N + 20000
        self.datetime_left = DatetimeIndex(range(N))
        self.datetime_right = DatetimeIndex(range(A, B))
        self.datetime_right2 = DatetimeIndex(range(N, B))

        options = np.arange(N)
        self.int_left = Index(options.take(np.random.permutation(N)[:N // 2]))
        self.int_right = Index(options.take(np.random.permutation(N)[:N // 2]))

        strs = tm.rands_array(10, N / 10)
        self.str_left = Index(strs[:N / 10 * 2 // 3])
        self.str_right = Index(strs[N / 10 // 3:])

    def time_datetime_intersection(self):
        self.dates_left.intersection(self.dates_right)

    def time_datetime_union(self):
        self.dates_left.union(self.dates_right)

    def time_datetime_difference(self):
        self.datetime_left.difference(self.datetime_right)

    def time_datetime_difference_disjoint(self):
        self.datetime_left.difference(self.datetime_right2)

    def time_datetime_symmetric_difference(self):
        self.datetime_left.symmetric_difference(self.datetime_right)

    def time_index_datetime_intersection(self):
        self.date_str_left.intersection(self.date_str_right)

    def time_index_datetime_union(self):
        self.date_str_left.union(self.date_str_right)

    def time_int64_intersection(self):
        self.int_left.intersection(self.int_right)

    def time_int64_union(self):
        self.int_left.union(self.int_right)

    def time_int64_difference(self):
        self.int_left.difference(self.int_right)

    def time_int64_symmetric_difference(self):
        self.int_left.symmetric_difference(self.int_right)

    def time_str_difference(self):
        self.str_left.difference(self.str_right)

    def time_str_symmetric_difference(self):
        self.str_left.symmetric_difference(self.str_right)


class Datetime(object):

    goal_time = 0.2

    def setup(self):
        self.dr = date_range('20000101', freq='D', periods=10000)

    def time_is_dates_only(self):
        self.dr._is_dates_only


class Ops(object):

    sample_time = 0.2
    params = ['float', 'int']
    param_names = ['dtype']

    def setup(self, dtype):
        N = 10**6
        indexes = {'int': 'makeIntIndex', 'float': 'makeFloatIndex'}
        self.index = getattr(tm, indexes[dtype])(N)

    def time_add(self, dtype):
        self.index + 2

    def time_subtract(self, dtype):
        self.index - 2

    def time_multiply(self, dtype):
        self.index * 2

    def time_divide(self, dtype):
        self.index / 2

    def time_modulo(self, dtype):
        self.index % 2


class Duplicated(object):

    goal_time = 0.2

    def setup(self):
        n, k = 200, 5000
        levels = [np.arange(n),
                  tm.makeStringIndex(n).values,
                  1000 + np.arange(n)]
        labels = [np.random.choice(n, (k * n)) for lev in levels]
        self.mi = MultiIndex(levels=levels, labels=labels)

    def time_duplicated(self):
        self.mi.duplicated()


class Sortlevel(object):

    goal_time = 0.2

    def setup(self):
        n = 1182720
        low, high = -4096, 4096
        arrs = [np.repeat(np.random.randint(low, high, (n // k)), k)
                for k in [11, 7, 5, 3, 1]]
        self.mi_int = MultiIndex.from_arrays(arrs)[np.random.permutation(n)]

        a = np.repeat(np.arange(100), 1000)
        b = np.tile(np.arange(1000), 100)
        self.mi = MultiIndex.from_arrays([a, b])
        self.mi = self.mi.take(np.random.permutation(np.arange(100000)))

    def time_sortlevel_int64(self):
        self.mi_int.sortlevel()

    def time_sortlevel_zero(self):
        self.mi.sortlevel(0)

    def time_sortlevel_one(self):
        self.mi.sortlevel(1)


class MultiIndexValues(object):

    goal_time = 0.2

    def setup_cache(self):

        level1 = range(1000)
        level2 = date_range(start='1/1/2012', periods=100)
        mi = MultiIndex.from_product([level1, level2])
        return mi

    def time_datetime_level_values_copy(self, mi):
        mi.copy().values

    def time_datetime_level_values_sliced(self, mi):
        mi[:10].values


class Range(object):

    goal_time = 0.2

    def setup(self):
        self.idx_inc = RangeIndex(start=0, stop=10**7, step=3)
        self.idx_dec = RangeIndex(start=10**7, stop=-1, step=-3)

    def time_max(self):
        self.idx_inc.max()

    def time_max_trivial(self):
        self.idx_dec.max()

    def time_min(self):
        self.idx_dec.min()

    def time_min_trivial(self):
        self.idx_inc.min()


class IndexAppend(object):

    goal_time = 0.2

    def setup(self):

        N = 10000
        self.range_idx = RangeIndex(0, 100)
        self.int_idx = self.range_idx.astype(int)
        self.obj_idx = self.int_idx.astype(str)
        self.range_idxs = []
        self.int_idxs = []
        self.object_idxs = []
        for i in range(1, N):
            r_idx = RangeIndex(i * 100, (i + 1) * 100)
            self.range_idxs.append(r_idx)
            i_idx = r_idx.astype(int)
            self.int_idxs.append(i_idx)
            o_idx = i_idx.astype(str)
            self.object_idxs.append(o_idx)

    def time_append_range_list(self):
        self.range_idx.append(self.range_idxs)

    def time_append_int_list(self):
        self.int_idx.append(self.int_idxs)

    def time_append_obj_list(self):
        self.obj_idx.append(self.object_idxs)
