import string

import numpy as np
import pandas.util.testing as tm
from pandas import (Series, date_range, DatetimeIndex, Index, MultiIndex,
                    RangeIndex, Float64Index)

from .pandas_vb_common import setup  # noqa


class SetOperations(object):

    goal_time = 0.2
    params = (['datetime', 'date_string', 'int', 'strings'],
              ['intersection', 'union', 'symmetric_difference'])
    param_names = ['dtype', 'method']

    def setup(self, dtype, method):
        N = 10**5
        dates_left = date_range('1/1/2000', periods=N, freq='T')
        fmt = '%Y-%m-%d %H:%M:%S'
        date_str_left = Index(dates_left.strftime(fmt))
        int_left = Index(np.arange(N))
        str_left = tm.makeStringIndex(N)
        data = {'datetime': {'left': dates_left, 'right': dates_left[:-1]},
                'date_string': {'left': date_str_left,
                                'right': date_str_left[:-1]},
                'int': {'left': int_left, 'right': int_left[:-1]},
                'strings': {'left': str_left, 'right': str_left[:-1]}}
        self.left = data[dtype]['left']
        self.right = data[dtype]['right']

    def time_operation(self, dtype, method):
        getattr(self.left, method)(self.right)


class SetDisjoint(object):

    goal_time = 0.2

    def setup(self):
        N = 10**5
        B = N + 20000
        self.datetime_left = DatetimeIndex(range(N))
        self.datetime_right = DatetimeIndex(range(N, B))

    def time_datetime_difference_disjoint(self):
        self.datetime_left.difference(self.datetime_right)


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


class Indexing(object):

    goal_time = 0.2
    params = ['String', 'Float', 'Int']
    param_names = ['dtype']

    def setup(self, dtype):
        N = 10**6
        self.idx = getattr(tm, 'make{}Index'.format(dtype))(N)
        self.array_mask = (np.arange(N) % 3) == 0
        self.series_mask = Series(self.array_mask)

    def time_boolean_array(self, dtype):
        self.idx[self.array_mask]

    def time_boolean_series(self, dtype):
        self.idx[self.series_mask]

    def time_get(self, dtype):
        self.idx[1]

    def time_slice(self, dtype):
        self.idx[:-1]

    def time_slice_step(self, dtype):
        self.idx[::2]


class Float64IndexMethod(object):
    # GH 13166
    goal_time = 0.2

    def setup(self):
        N = 100000
        a = np.arange(N)
        self.ind = Float64Index(a * 4.8000000418824129e-08)

    def time_get_loc(self):
        self.ind.get_loc(0)


class MultiIndexGet(object):

    goal_time = 0.2

    def setup(self):
        self.mi_large = MultiIndex.from_product(
            [np.arange(1000), np.arange(20), list(string.ascii_letters)],
            names=['one', 'two', 'three'])
        self.mi_med = MultiIndex.from_product(
            [np.arange(1000), np.arange(10), list('A')],
            names=['one', 'two', 'three'])
        self.mi_small = MultiIndex.from_product(
            [np.arange(100), list('A'), list('A')],
            names=['one', 'two', 'three'])

    def time_multiindex_large_get_loc(self):
        self.mi_large.get_loc((999, 19, 'Z'))

    def time_multiindex_large_get_loc_warm(self):
        for _ in range(1000):
            self.mi_large.get_loc((999, 19, 'Z'))

    def time_multiindex_med_get_loc(self):
        self.mi_med.get_loc((999, 9, 'A'))

    def time_multiindex_med_get_loc_warm(self):
        for _ in range(1000):
            self.mi_med.get_loc((999, 9, 'A'))

    def time_multiindex_string_get_loc(self):
        self.mi_small.get_loc((99, 'A', 'A'))

    def time_multiindex_small_get_loc_warm(self):
        for _ in range(1000):
            self.mi_small.get_loc((99, 'A', 'A'))


class MultiIndexDuplicates(object):

    goal_time = 0.2

    def setup(self):
        size = 65536
        arrays = [np.random.randint(0, 8192, size),
                  np.random.randint(0, 1024, size)]
        mask = np.random.rand(size) < 0.1
        self.mi_unused_levels = MultiIndex.from_arrays(arrays)
        self.mi_unused_levels = self.mi_unused_levels[mask]

    def time_remove_unused_levels(self):
        self.mi_unused_levels.remove_unused_levels()


class MultiIndexInteger(object):

    goal_time = 0.2

    def setup(self):
        self.mi_int = MultiIndex.from_product([np.arange(1000),
                                               np.arange(1000)],
                                              names=['one', 'two'])
        self.obj_index = np.array([(0, 10), (0, 11), (0, 12),
                                   (0, 13), (0, 14), (0, 15),
                                   (0, 16), (0, 17), (0, 18),
                                   (0, 19)], dtype=object)

    def time_get_indexer(self):
        self.mi_int.get_indexer(self.obj_index)

    def time_is_monotonic(self):
        self.mi_int.is_monotonic
