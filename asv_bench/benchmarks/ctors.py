import numpy as np
from pandas import DataFrame, Series, Index, DatetimeIndex, Timestamp

from .pandas_vb_common import setup # noqa


class Constructors(object):

    goal_time = 0.2

    def setup(self):
        N = 10**2
        self.arr = np.random.randn(N, N)
        self.arr_str = np.array(['foo', 'bar', 'baz'], dtype=object)

        self.data = np.random.randn(N)
        self.index = Index(np.arange(N))

        self.s = Series([Timestamp('20110101'), Timestamp('20120101'),
                         Timestamp('20130101')] * N * 10)

    def time_frame_from_ndarray(self):
        DataFrame(self.arr)

    def time_series_from_ndarray(self):
        Series(self.data, index=self.index)

    def time_index_from_array_string(self):
        Index(self.arr_str)

    def time_dtindex_from_series(self):
        DatetimeIndex(self.s)

    def time_dtindex_from_index_with_series(self):
        Index(self.s)
