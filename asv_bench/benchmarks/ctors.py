from .pandas_vb_common import *


class Constructors(object):
    goal_time = 0.2

    def setup(self):
        self.arr = np.random.randn(100, 100)
        self.arr_str = np.array(['foo', 'bar', 'baz'], dtype=object)

        self.data = np.random.randn(100)
        self.index = Index(np.arange(100))

        self.s = Series(([Timestamp('20110101'), Timestamp('20120101'),
                          Timestamp('20130101')] * 1000))

    def time_frame_from_ndarray(self):
        DataFrame(self.arr)

    def time_series_from_ndarray(self):
        pd.Series(self.data, index=self.index)

    def time_index_from_array_string(self):
        Index(self.arr_str)

    def time_dtindex_from_series(self):
        DatetimeIndex(self.s)

    def time_dtindex_from_series2(self):
        Index(self.s)
