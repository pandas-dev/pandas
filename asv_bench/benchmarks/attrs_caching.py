from pandas_vb_common import *


class getattr_dataframe_index(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(10, 6))
        self.cur_index = self.df.index

    def time_getattr_dataframe_index(self):
        self.foo = self.df.index


class setattr_dataframe_index(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(10, 6))
        self.cur_index = self.df.index

    def time_setattr_dataframe_index(self):
        self.df.index = self.cur_index