from pandas_vb_common import *


class concat_categorical(object):
    goal_time = 0.2

    def setup(self):
        self.s = pd.Series((list('aabbcd') * 1000000)).astype('category')

    def time_concat_categorical(self):
        concat([self.s, self.s])