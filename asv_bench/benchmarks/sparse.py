import numpy as np
import scipy.sparse

import pandas as pd
from pandas import MultiIndex, Series, SparseArray, date_range


def make_array(size, dense_proportion, fill_value, dtype):
    dense_size = int(size * dense_proportion)
    arr = np.full(size, fill_value, dtype)
    indexer = np.random.choice(np.arange(size), dense_size, replace=False)
    arr[indexer] = np.random.choice(np.arange(100, dtype=dtype), dense_size)
    return arr


class SparseSeriesToFrame:
    def setup(self):
        K = 50
        N = 50001
        rng = date_range("1/1/2000", periods=N, freq="T")
        self.series = {}
        for i in range(1, K):
            data = np.random.randn(N)[:-i]
            idx = rng[:-i]
            data[100:] = np.nan
            self.series[i] = pd.Series(pd.SparseArray(data), index=idx)

    def time_series_to_frame(self):
        pd.DataFrame(self.series)


class SparseArrayConstructor:

    params = ([0.1, 0.01], [0, np.nan], [np.int64, np.float64, np.object])
    param_names = ["dense_proportion", "fill_value", "dtype"]

    def setup(self, dense_proportion, fill_value, dtype):
        N = 10 ** 6
        self.array = make_array(N, dense_proportion, fill_value, dtype)

    def time_sparse_array(self, dense_proportion, fill_value, dtype):
        SparseArray(self.array, fill_value=fill_value, dtype=dtype)


class SparseDataFrameConstructor:
    def setup(self):
        N = 1000
        self.arr = np.arange(N)
        self.sparse = scipy.sparse.rand(N, N, 0.005)

    def time_from_scipy(self):
        pd.DataFrame.sparse.from_spmatrix(self.sparse)


class FromCoo:
    def setup(self):
        self.matrix = scipy.sparse.coo_matrix(
            ([3.0, 1.0, 2.0], ([1, 0, 0], [0, 2, 3])), shape=(100, 100)
        )

    def time_sparse_series_from_coo(self):
        pd.Series.sparse.from_coo(self.matrix)


class ToCoo:
    def setup(self):
        s = Series([np.nan] * 10000)
        s[0] = 3.0
        s[100] = -1.0
        s[999] = 12.1
        s.index = MultiIndex.from_product([range(10)] * 4)
        self.ss = s.astype("Sparse")

    def time_sparse_series_to_coo(self):
        self.ss.sparse.to_coo(row_levels=[0, 1], column_levels=[2, 3], sort_labels=True)


class Arithmetic:

    params = ([0.1, 0.01], [0, np.nan])
    param_names = ["dense_proportion", "fill_value"]

    def setup(self, dense_proportion, fill_value):
        N = 10 ** 6
        arr1 = make_array(N, dense_proportion, fill_value, np.int64)
        self.array1 = SparseArray(arr1, fill_value=fill_value)
        arr2 = make_array(N, dense_proportion, fill_value, np.int64)
        self.array2 = SparseArray(arr2, fill_value=fill_value)

    def time_make_union(self, dense_proportion, fill_value):
        self.array1.sp_index.make_union(self.array2.sp_index)

    def time_intersect(self, dense_proportion, fill_value):
        self.array1.sp_index.intersect(self.array2.sp_index)

    def time_add(self, dense_proportion, fill_value):
        self.array1 + self.array2

    def time_divide(self, dense_proportion, fill_value):
        self.array1 / self.array2


class ArithmeticBlock:

    params = [np.nan, 0]
    param_names = ["fill_value"]

    def setup(self, fill_value):
        N = 10 ** 6
        self.arr1 = self.make_block_array(
            length=N, num_blocks=1000, block_size=10, fill_value=fill_value
        )
        self.arr2 = self.make_block_array(
            length=N, num_blocks=1000, block_size=10, fill_value=fill_value
        )

    def make_block_array(self, length, num_blocks, block_size, fill_value):
        arr = np.full(length, fill_value)
        indicies = np.random.choice(
            np.arange(0, length, block_size), num_blocks, replace=False
        )
        for ind in indicies:
            arr[ind : ind + block_size] = np.random.randint(0, 100, block_size)
        return SparseArray(arr, fill_value=fill_value)

    def time_make_union(self, fill_value):
        self.arr1.sp_index.make_union(self.arr2.sp_index)

    def time_intersect(self, fill_value):
        self.arr2.sp_index.intersect(self.arr2.sp_index)

    def time_addition(self, fill_value):
        self.arr1 + self.arr2

    def time_division(self, fill_value):
        self.arr1 / self.arr2


from .pandas_vb_common import setup  # noqa: F401 isort:skip
