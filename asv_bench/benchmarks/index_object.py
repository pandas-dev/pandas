import numpy as np
import pandas.util.testing as tm
from pandas import (Series, date_range, DatetimeIndex, Index, RangeIndex,
                    Float64Index)

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
    params = (['first', 'last', False], [True, False])
    param_names = ['keep', 'return_inverse']

    def setup(self, keep, return_inverse):
        if keep is False and return_inverse:
            raise NotImplementedError

        n, k = 200, 1000
        base = tm.makeStringIndex(n)
        self.idx = Index(base[np.random.choice(n, k * n)])

    def time_duplicated(self, keep, return_inverse):
        self.idx.duplicated(keep=keep, return_inverse=return_inverse)


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
        self.sorted = self.idx.sort_values()
        half = N // 2
        self.non_unique = self.idx[:half].append(self.idx[:half])
        self.non_unique_sorted = self.sorted[:half].append(self.sorted[:half])
        self.key = self.sorted[N // 4]

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

    def time_get_loc(self, dtype):
        self.idx.get_loc(self.key)

    def time_get_loc_sorted(self, dtype):
        self.sorted.get_loc(self.key)

    def time_get_loc_non_unique(self, dtype):
        self.non_unique.get_loc(self.key)

    def time_get_loc_non_unique_sorted(self, dtype):
        self.non_unique_sorted.get_loc(self.key)


class Float64IndexMethod(object):
    # GH 13166
    goal_time = 0.2

    def setup(self):
        N = 100000
        a = np.arange(N)
        self.ind = Float64Index(a * 4.8000000418824129e-08)

    def time_get_loc(self):
        self.ind.get_loc(0)
