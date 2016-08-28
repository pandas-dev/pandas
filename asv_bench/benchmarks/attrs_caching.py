from .pandas_vb_common import *


class DataFrameAttributes(object):
    goal_time = 0.2

    def setup(self):
        self.df = DataFrame(np.random.randn(10, 6))
        self.cur_index = self.df.index

    def time_get_index(self):
        self.foo = self.df.index

    def time_set_index(self):
        self.df.index = self.cur_index
