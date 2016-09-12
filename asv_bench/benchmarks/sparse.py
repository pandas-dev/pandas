from .pandas_vb_common import *
import pandas.sparse.series
import scipy.sparse
from pandas.core.sparse import SparseSeries, SparseDataFrame
from pandas.core.sparse import SparseDataFrame


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


class sparse_frame_constructor(object):
    goal_time = 0.2

    def time_sparse_frame_constructor(self):
        SparseDataFrame(columns=np.arange(100), index=np.arange(1000))


class sparse_series_from_coo(object):
    goal_time = 0.2

    def setup(self):
        self.A = scipy.sparse.coo_matrix(([3.0, 1.0, 2.0], ([1, 0, 0], [0, 2, 3])), shape=(100, 100))

    def time_sparse_series_from_coo(self):
        self.ss = pandas.sparse.series.SparseSeries.from_coo(self.A)


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
