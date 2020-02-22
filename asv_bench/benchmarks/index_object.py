import gc

import numpy as np

from pandas import (
    DatetimeIndex,
    Float64Index,
    Index,
    IntervalIndex,
    MultiIndex,
    RangeIndex,
    Series,
    date_range,
)

from .pandas_vb_common import tm


class SetOperations:

    params = (
        ["datetime", "date_string", "int", "strings"],
        ["intersection", "union", "symmetric_difference"],
    )
    param_names = ["dtype", "method"]

    def setup(self, dtype, method):
        N = 10 ** 5
        dates_left = date_range("1/1/2000", periods=N, freq="T")
        fmt = "%Y-%m-%d %H:%M:%S"
        date_str_left = Index(dates_left.strftime(fmt))
        int_left = Index(np.arange(N))
        str_left = tm.makeStringIndex(N)
        data = {
            "datetime": {"left": dates_left, "right": dates_left[:-1]},
            "date_string": {"left": date_str_left, "right": date_str_left[:-1]},
            "int": {"left": int_left, "right": int_left[:-1]},
            "strings": {"left": str_left, "right": str_left[:-1]},
        }
        self.left = data[dtype]["left"]
        self.right = data[dtype]["right"]

    def time_operation(self, dtype, method):
        getattr(self.left, method)(self.right)


class SetDisjoint:
    def setup(self):
        N = 10 ** 5
        B = N + 20000
        self.datetime_left = DatetimeIndex(range(N))
        self.datetime_right = DatetimeIndex(range(N, B))

    def time_datetime_difference_disjoint(self):
        self.datetime_left.difference(self.datetime_right)


class Datetime:
    def setup(self):
        self.dr = date_range("20000101", freq="D", periods=10000)

    def time_is_dates_only(self):
        self.dr._is_dates_only


class Range:
    def setup(self):
        self.idx_inc = RangeIndex(start=0, stop=10 ** 7, step=3)
        self.idx_dec = RangeIndex(start=10 ** 7, stop=-1, step=-3)

    def time_max(self):
        self.idx_inc.max()

    def time_max_trivial(self):
        self.idx_dec.max()

    def time_min(self):
        self.idx_dec.min()

    def time_min_trivial(self):
        self.idx_inc.min()

    def time_get_loc_inc(self):
        self.idx_inc.get_loc(900000)

    def time_get_loc_dec(self):
        self.idx_dec.get_loc(100000)


class IndexEquals:
    def setup(self):
        idx_large_fast = RangeIndex(100000)
        idx_small_slow = date_range(start="1/1/2012", periods=1)
        self.mi_large_slow = MultiIndex.from_product([idx_large_fast, idx_small_slow])

        self.idx_non_object = RangeIndex(1)

    def time_non_object_equals_multiindex(self):
        self.idx_non_object.equals(self.mi_large_slow)


class IndexAppend:
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


class Indexing:

    params = ["String", "Float", "Int"]
    param_names = ["dtype"]

    def setup(self, dtype):
        N = 10 ** 6
        self.idx = getattr(tm, f"make{dtype}Index")(N)
        self.array_mask = (np.arange(N) % 3) == 0
        self.series_mask = Series(self.array_mask)
        self.sorted = self.idx.sort_values()
        half = N // 2
        self.non_unique = self.idx[:half].append(self.idx[:half])
        self.non_unique_sorted = (
            self.sorted[:half].append(self.sorted[:half]).sort_values()
        )
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


class Float64IndexMethod:
    # GH 13166
    def setup(self):
        N = 100000
        a = np.arange(N)
        self.ind = Float64Index(a * 4.8000000418824129e-08)

    def time_get_loc(self):
        self.ind.get_loc(0)


class IntervalIndexMethod:
    # GH 24813
    params = [10 ** 3, 10 ** 5]

    def setup(self, N):
        left = np.append(np.arange(N), np.array(0))
        right = np.append(np.arange(1, N + 1), np.array(1))
        self.intv = IntervalIndex.from_arrays(left, right)
        self.intv._engine

        self.intv2 = IntervalIndex.from_arrays(left + 1, right + 1)
        self.intv2._engine

        self.left = IntervalIndex.from_breaks(np.arange(N))
        self.right = IntervalIndex.from_breaks(np.arange(N - 3, 2 * N - 3))

    def time_monotonic_inc(self, N):
        self.intv.is_monotonic_increasing

    def time_is_unique(self, N):
        self.intv.is_unique

    def time_intersection(self, N):
        self.left.intersection(self.right)

    def time_intersection_one_duplicate(self, N):
        self.intv.intersection(self.right)

    def time_intersection_both_duplicate(self, N):
        self.intv.intersection(self.intv2)


class GC:
    params = [1, 2, 5]

    def create_use_drop(self):
        idx = Index(list(range(1000 * 1000)))
        idx._engine

    def peakmem_gc_instances(self, N):
        try:
            gc.disable()

            for _ in range(N):
                self.create_use_drop()
        finally:
            gc.enable()


from .pandas_vb_common import setup  # noqa: F401 isort:skip
