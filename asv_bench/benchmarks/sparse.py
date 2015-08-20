from pandas_vb_common import *
import scipy.sparse
import pandas.sparse.series
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