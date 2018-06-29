import itertools

import numpy as np
import scipy.sparse
from pandas import (SparseSeries, SparseDataFrame, SparseArray, Series,
                    date_range, MultiIndex)

from .pandas_vb_common import setup  # noqa


def make_array(size, dense_proportion, fill_value, dtype):
    dense_size = int(size * dense_proportion)
    arr = np.full(size, fill_value, dtype)
    indexer = np.random.choice(np.arange(size), dense_size, replace=False)
    arr[indexer] = np.random.choice(np.arange(100, dtype=dtype), dense_size)
    return arr


class SparseSeriesToFrame(object):

    goal_time = 0.2

    def setup(self):
        K = 50
        N = 50001
        rng = date_range('1/1/2000', periods=N, freq='T')
        self.series = {}
        for i in range(1, K):
            data = np.random.randn(N)[:-i]
            idx = rng[:-i]
            data[100:] = np.nan
            self.series[i] = SparseSeries(data, index=idx)

    def time_series_to_frame(self):
        SparseDataFrame(self.series)


class SparseArrayConstructor(object):

    goal_time = 0.2
    params = ([0.1, 0.01], [0, np.nan],
              [np.int64, np.float64, np.object])
    param_names = ['dense_proportion', 'fill_value', 'dtype']

    def setup(self, dense_proportion, fill_value, dtype):
        N = 10**6
        self.array = make_array(N, dense_proportion, fill_value, dtype)

    def time_sparse_array(self, dense_proportion, fill_value, dtype):
        SparseArray(self.array, fill_value=fill_value, dtype=dtype)


class SparseDataFrameConstructor(object):

    goal_time = 0.2

    def setup(self):
        N = 1000
        self.arr = np.arange(N)
        self.sparse = scipy.sparse.rand(N, N, 0.005)
        self.dict = dict(zip(range(N), itertools.repeat([0])))

    def time_constructor(self):
        SparseDataFrame(columns=self.arr, index=self.arr)

    def time_from_scipy(self):
        SparseDataFrame(self.sparse)

    def time_from_dict(self):
        SparseDataFrame(self.dict)


class FromCoo(object):

    goal_time = 0.2

    def setup(self):
        self.matrix = scipy.sparse.coo_matrix(([3.0, 1.0, 2.0],
                                               ([1, 0, 0], [0, 2, 3])),
                                              shape=(100, 100))

    def time_sparse_series_from_coo(self):
        SparseSeries.from_coo(self.matrix)


class ToCoo(object):

    goal_time = 0.2

    def setup(self):
        s = Series([np.nan] * 10000)
        s[0] = 3.0
        s[100] = -1.0
        s[999] = 12.1
        s.index = MultiIndex.from_product([range(10)] * 4)
        self.ss = s.to_sparse()

    def time_sparse_series_to_coo(self):
        self.ss.to_coo(row_levels=[0, 1],
                       column_levels=[2, 3],
                       sort_labels=True)


class Arithmetic(object):

    goal_time = 0.2
    params = ([0.1, 0.01], [0, np.nan])
    param_names = ['dense_proportion', 'fill_value']

    def setup(self, dense_proportion, fill_value):
        N = 10**6
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


class ArithmeticBlock(object):

    goal_time = 0.2
    params = [np.nan, 0]
    param_names = ['fill_value']

    def setup(self, fill_value):
        N = 10**6
        self.arr1 = self.make_block_array(length=N, num_blocks=1000,
                                          block_size=10, fill_value=fill_value)
        self.arr2 = self.make_block_array(length=N, num_blocks=1000,
                                          block_size=10, fill_value=fill_value)

    def make_block_array(self, length, num_blocks, block_size, fill_value):
        arr = np.full(length, fill_value)
        indicies = np.random.choice(np.arange(0, length, block_size),
                                    num_blocks,
                                    replace=False)
        for ind in indicies:
            arr[ind:ind + block_size] = np.random.randint(0, 100, block_size)
        return SparseArray(arr, fill_value=fill_value)

    def time_make_union(self, fill_value):
        self.arr1.sp_index.make_union(self.arr2.sp_index)

    def time_intersect(self, fill_value):
        self.arr2.sp_index.intersect(self.arr2.sp_index)

    def time_addition(self, fill_value):
        self.arr1 + self.arr2

    def time_division(self, fill_value):
        self.arr1 / self.arr2
