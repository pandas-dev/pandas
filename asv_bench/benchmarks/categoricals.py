from .pandas_vb_common import *
try:
    from pandas.api.types import union_categoricals
except ImportError:
    try:
        from pandas.types.concat import union_categoricals
    except ImportError:
        pass


class Categoricals(object):
    goal_time = 0.2

    def setup(self):
        N = 100000
        self.s = pd.Series((list('aabbcd') * N)).astype('category')

        self.a = pd.Categorical((list('aabbcd') * N))
        self.b = pd.Categorical((list('bbcdjk') * N))

        self.categories = list('abcde')
        self.cat_idx = Index(self.categories)
        self.values = np.tile(self.categories, N)
        self.codes = np.tile(range(len(self.categories)), N)

        self.datetimes = pd.Series(pd.date_range(
            '1995-01-01 00:00:00', periods=10000, freq='s'))

    def time_concat(self):
        concat([self.s, self.s])

    def time_union(self):
        union_categoricals([self.a, self.b])

    def time_constructor_regular(self):
        Categorical(self.values, self.categories)

    def time_constructor_fastpath(self):
        Categorical(self.codes, self.cat_idx, fastpath=True)

    def time_constructor_datetimes(self):
        Categorical(self.datetimes)

    def time_constructor_datetimes_with_nat(self):
        t = self.datetimes
        t.iloc[-1] = pd.NaT
        Categorical(t)


class Categoricals2(object):
    goal_time = 0.2

    def setup(self):
        n = 500000
        np.random.seed(2718281)
        arr = ['s%04d' % i for i in np.random.randint(0, n // 10, size=n)]
        self.ts = Series(arr).astype('category')

        self.sel = self.ts.loc[[0]]

    def time_value_counts(self):
        self.ts.value_counts(dropna=False)

    def time_value_counts_dropna(self):
        self.ts.value_counts(dropna=True)

    def time_rendering(self):
        str(self.sel)


class Categoricals3(object):
    goal_time = 0.2

    def setup(self):
        N = 100000
        ncats = 100

        self.s1 = Series(np.array(tm.makeCategoricalIndex(N, ncats)))
        self.s1_cat = self.s1.astype('category')
        self.s1_cat_ordered = self.s1.astype('category', ordered=True)

        self.s2 = Series(np.random.randint(0, ncats, size=N))
        self.s2_cat = self.s2.astype('category')
        self.s2_cat_ordered = self.s2.astype('category', ordered=True)

    def time_rank_string(self):
        self.s1.rank()

    def time_rank_string_cat(self):
        self.s1_cat.rank()

    def time_rank_string_cat_ordered(self):
        self.s1_cat_ordered.rank()

    def time_rank_int(self):
        self.s2.rank()

    def time_rank_int_cat(self):
        self.s2_cat.rank()

    def time_rank_int_cat_ordered(self):
        self.s2_cat_ordered.rank()
