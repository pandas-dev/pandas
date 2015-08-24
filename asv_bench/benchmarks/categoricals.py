from pandas_vb_common import *


class concat_categorical(object):
    goal_time = 0.2

    def setup(self):
        self.s = pd.Series((list('aabbcd') * 1000000)).astype('category')

    def time_concat_categorical(self):
        concat([self.s, self.s])


class categorical_value_counts(object):
    goal_time = 1

    def setup(self):
        n = 500000
        np.random.seed(2718281)
        arr = ['s%04d' % i for i in np.random.randint(0, n // 10, size=n)]
        self.ts = Series(arr).astype('category')

    def time_value_counts(self):
        self.ts.value_counts(dropna=False)

    def time_value_counts_dropna(self):
        self.ts.value_counts(dropna=True)
