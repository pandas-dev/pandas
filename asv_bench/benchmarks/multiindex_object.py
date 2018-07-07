import string

import numpy as np
import pandas.util.testing as tm
from pandas import date_range, MultiIndex

from .pandas_vb_common import setup  # noqa


class GetLoc(object):

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

    def time_large_get_loc(self):
        self.mi_large.get_loc((999, 19, 'Z'))

    def time_large_get_loc_warm(self):
        for _ in range(1000):
            self.mi_large.get_loc((999, 19, 'Z'))

    def time_med_get_loc(self):
        self.mi_med.get_loc((999, 9, 'A'))

    def time_med_get_loc_warm(self):
        for _ in range(1000):
            self.mi_med.get_loc((999, 9, 'A'))

    def time_string_get_loc(self):
        self.mi_small.get_loc((99, 'A', 'A'))

    def time_small_get_loc_warm(self):
        for _ in range(1000):
            self.mi_small.get_loc((99, 'A', 'A'))


class Duplicates(object):

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


class Integer(object):

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


class Values(object):

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
