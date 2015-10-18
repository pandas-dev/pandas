from .pandas_vb_common import *
import string


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


class categorical_constructor(object):
    goal_time = 0.2

    def setup(self):
        n = 5
        N = 1e6
        self.categories = list(string.ascii_letters[:n])
        self.cat_idx = Index(self.categories)
        self.values = np.tile(self.categories, N)
        self.codes = np.tile(range(n), N)

    def time_regular_constructor(self):
        Categorical(self.values, self.categories)

    def time_fastpath(self):
        Categorical(self.codes, self.cat_idx, fastpath=True)


class categorical_rendering(object):
    goal_time = 3e-3

    def setup(self):
        n = 1000
        items = [str(i) for i in range(n)]
        s = pd.Series(items, dtype='category')
        df = pd.DataFrame({'C': s, 'data': np.random.randn(n)})
        self.data = df[df.C == '20']

    def time_rendering(self):
        str(self.data.C)
