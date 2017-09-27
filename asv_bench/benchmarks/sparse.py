import itertools

from .pandas_vb_common import *
import scipy.sparse
from pandas import SparseSeries, SparseDataFrame, SparseArray


class sparse_series_to_frame(object):
    goal_time = 0.2

    def setup(self):
        self.K = 50
        self.N = 50000
        self.rng = np.asarray(date_range('1/1/2000', periods=self.N, freq='T'))
        self.series = {}
        for i in range(1, (self.K + 1)):
            self.data = np.random.randn(self.N)[:(- i)]
            self.this_rng = self.rng[:(- i)]
            self.data[100:] = np.nan
            self.series[i] = SparseSeries(self.data, index=self.this_rng)

    def time_sparse_series_to_frame(self):
        SparseDataFrame(self.series)


class sparse_array_constructor(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1)
        self.int64_10percent = self.make_numeric_array(length=1000000, dense_size=100000, fill_value=0, dtype=np.int64)
        self.int64_1percent = self.make_numeric_array(length=1000000, dense_size=10000, fill_value=0, dtype=np.int64)

        self.float64_10percent = self.make_numeric_array(length=1000000, dense_size=100000, fill_value=np.nan, dtype=np.float64)
        self.float64_1percent = self.make_numeric_array(length=1000000, dense_size=10000, fill_value=np.nan, dtype=np.float64)

        self.object_nan_fill_value_10percent = self.make_object_array(length=1000000, dense_size=100000, fill_value=np.nan)
        self.object_nan_fill_value_1percent = self.make_object_array(length=1000000, dense_size=10000, fill_value=np.nan)

        self.object_non_nan_fill_value_10percent = self.make_object_array(length=1000000, dense_size=100000, fill_value=0)
        self.object_non_nan_fill_value_1percent = self.make_object_array(length=1000000, dense_size=10000, fill_value=0)

    def make_numeric_array(self, length, dense_size, fill_value, dtype):
        arr = np.array([fill_value] * length, dtype=dtype)
        indexer = np.unique(np.random.randint(0, length, dense_size))
        arr[indexer] = np.random.randint(0, 100, len(indexer))
        return (arr, fill_value, dtype)

    def make_object_array(self, length, dense_size, fill_value):
        elems = np.array(['a', 0.0, False, 1, 2], dtype=np.object)
        arr = np.array([fill_value] * length, dtype=np.object)
        indexer = np.unique(np.random.randint(0, length, dense_size))
        arr[indexer] = np.random.choice(elems, len(indexer))
        return (arr, fill_value, np.object)

    def time_sparse_array_constructor_int64_10percent(self):
        arr, fill_value, dtype = self.int64_10percent
        SparseArray(arr, fill_value=fill_value, dtype=dtype)

    def time_sparse_array_constructor_int64_1percent(self):
        arr, fill_value, dtype = self.int64_1percent
        SparseArray(arr, fill_value=fill_value, dtype=dtype)

    def time_sparse_array_constructor_float64_10percent(self):
        arr, fill_value, dtype = self.float64_10percent
        SparseArray(arr, fill_value=fill_value, dtype=dtype)

    def time_sparse_array_constructor_float64_1percent(self):
        arr, fill_value, dtype = self.float64_1percent
        SparseArray(arr, fill_value=fill_value, dtype=dtype)

    def time_sparse_array_constructor_object_nan_fill_value_10percent(self):
        arr, fill_value, dtype = self.object_nan_fill_value_10percent
        SparseArray(arr, fill_value=fill_value, dtype=dtype)

    def time_sparse_array_constructor_object_nan_fill_value_1percent(self):
        arr, fill_value, dtype = self.object_nan_fill_value_1percent
        SparseArray(arr, fill_value=fill_value, dtype=dtype)

    def time_sparse_array_constructor_object_non_nan_fill_value_10percent(self):
        arr, fill_value, dtype = self.object_non_nan_fill_value_10percent
        SparseArray(arr, fill_value=fill_value, dtype=dtype)

    def time_sparse_array_constructor_object_non_nan_fill_value_1percent(self):
        arr, fill_value, dtype = self.object_non_nan_fill_value_1percent
        SparseArray(arr, fill_value=fill_value, dtype=dtype)


class sparse_frame_constructor(object):
    goal_time = 0.2

    def time_sparse_frame_constructor(self):
        SparseDataFrame(columns=np.arange(100), index=np.arange(1000))

    def time_sparse_from_scipy(self):
        SparseDataFrame(scipy.sparse.rand(1000, 1000, 0.005))

    def time_sparse_from_dict(self):
        SparseDataFrame(dict(zip(range(1000), itertools.repeat([0]))))


class sparse_series_from_coo(object):
    goal_time = 0.2

    def setup(self):
        self.A = scipy.sparse.coo_matrix(([3.0, 1.0, 2.0], ([1, 0, 0], [0, 2, 3])), shape=(100, 100))

    def time_sparse_series_from_coo(self):
        self.ss = SparseSeries.from_coo(self.A)


class sparse_series_to_coo(object):
    goal_time = 0.2

    def setup(self):
        self.s = pd.Series(([np.nan] * 10000))
        self.s[0] = 3.0
        self.s[100] = (-1.0)
        self.s[999] = 12.1
        self.s.index = pd.MultiIndex.from_product((range(10), range(10), range(10), range(10)))
        self.ss = self.s.to_sparse()

    def time_sparse_series_to_coo(self):
        self.ss.to_coo(row_levels=[0, 1], column_levels=[2, 3], sort_labels=True)


class sparse_arithmetic_int(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1)
        self.a_10percent = self.make_sparse_array(length=1000000, dense_size=100000, fill_value=np.nan)
        self.b_10percent = self.make_sparse_array(length=1000000, dense_size=100000, fill_value=np.nan)

        self.a_10percent_zero = self.make_sparse_array(length=1000000, dense_size=100000, fill_value=0)
        self.b_10percent_zero = self.make_sparse_array(length=1000000, dense_size=100000, fill_value=0)

        self.a_1percent = self.make_sparse_array(length=1000000, dense_size=10000, fill_value=np.nan)
        self.b_1percent = self.make_sparse_array(length=1000000, dense_size=10000, fill_value=np.nan)

    def make_sparse_array(self, length, dense_size, fill_value):
        arr = np.array([fill_value] * length, dtype=np.float64)
        indexer = np.unique(np.random.randint(0, length, dense_size))
        arr[indexer] = np.random.randint(0, 100, len(indexer))
        return pd.SparseArray(arr, fill_value=fill_value)

    def time_sparse_make_union(self):
        self.a_10percent.sp_index.make_union(self.b_10percent.sp_index)

    def time_sparse_intersect(self):
        self.a_10percent.sp_index.intersect(self.b_10percent.sp_index)

    def time_sparse_addition_10percent(self):
        self.a_10percent + self.b_10percent

    def time_sparse_addition_10percent_zero(self):
        self.a_10percent_zero + self.b_10percent_zero

    def time_sparse_addition_1percent(self):
        self.a_1percent + self.b_1percent

    def time_sparse_division_10percent(self):
        self.a_10percent / self.b_10percent

    def time_sparse_division_10percent_zero(self):
        self.a_10percent_zero / self.b_10percent_zero

    def time_sparse_division_1percent(self):
        self.a_1percent / self.b_1percent



class sparse_arithmetic_block(object):
    goal_time = 0.2

    def setup(self):
        np.random.seed(1)
        self.a = self.make_sparse_array(length=1000000, num_blocks=1000,
                                        block_size=10, fill_value=np.nan)
        self.b = self.make_sparse_array(length=1000000, num_blocks=1000,
                                        block_size=10, fill_value=np.nan)

        self.a_zero = self.make_sparse_array(length=1000000, num_blocks=1000,
                                             block_size=10, fill_value=0)
        self.b_zero = self.make_sparse_array(length=1000000, num_blocks=1000,
                                             block_size=10, fill_value=np.nan)

    def make_sparse_array(self, length, num_blocks, block_size, fill_value):
        a = np.array([fill_value] * length)
        for block in range(num_blocks):
            i = np.random.randint(0, length)
            a[i:i + block_size] = np.random.randint(0, 100, len(a[i:i + block_size]))
        return pd.SparseArray(a, fill_value=fill_value)

    def time_sparse_make_union(self):
        self.a.sp_index.make_union(self.b.sp_index)

    def time_sparse_intersect(self):
        self.a.sp_index.intersect(self.b.sp_index)

    def time_sparse_addition(self):
        self.a + self.b

    def time_sparse_addition_zero(self):
        self.a_zero + self.b_zero

    def time_sparse_division(self):
        self.a / self.b

    def time_sparse_division_zero(self):
        self.a_zero / self.b_zero
