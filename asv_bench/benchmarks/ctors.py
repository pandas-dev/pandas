from pandas_vb_common import *


class frame_constructor_ndarray(object):
    goal_time = 0.2

    def setup(self):
        self.arr = np.random.randn(100, 100)

    def time_frame_constructor_ndarray(self):
        DataFrame(self.arr)


class ctor_index_array_string(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.array(['foo', 'bar', 'baz'], dtype=object)

    def time_ctor_index_array_string(self):
        Index(self.data)


class series_constructor_ndarray(object):
    goal_time = 0.2

    def setup(self):
        self.data = np.random.randn(100)
        self.index = Index(np.arange(100))

    def time_series_constructor_ndarray(self):
        Series(self.data, index=self.index)


class dtindex_from_series_ctor(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(([Timestamp('20110101'), Timestamp('20120101'), Timestamp('20130101')] * 1000))

    def time_dtindex_from_series_ctor(self):
        DatetimeIndex(self.s)


class index_from_series_ctor(object):
    goal_time = 0.2

    def setup(self):
        self.s = Series(([Timestamp('20110101'), Timestamp('20120101'), Timestamp('20130101')] * 1000))

    def time_index_from_series_ctor(self):
        Index(self.s)