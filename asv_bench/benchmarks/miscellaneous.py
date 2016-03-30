from .pandas_vb_common import *
from pandas.util.decorators import cache_readonly


class match_strings(object):
    goal_time = 0.2

    def setup(self):
        self.uniques = tm.makeStringIndex(1000).values
        self.all = self.uniques.repeat(10)

    def time_match_strings(self):
        match(self.all, self.uniques)


class misc_cache_readonly(object):
    goal_time = 0.2

    def setup(self):


        class Foo:

            @cache_readonly
            def prop(self):
                return 5
        self.obj = Foo()

    def time_misc_cache_readonly(self):
        self.obj.prop


class to_numeric(object):
    goal_time = 0.2

    def setup(self):
        self.n = 10000
        self.float = Series(np.random.randn(self.n * 100))
        self.numstr = self.float.astype('str')
        self.str = Series(tm.makeStringIndex(self.n))

    def time_from_float(self):
        pd.to_numeric(self.float)

    def time_from_numeric_str(self):
        pd.to_numeric(self.numstr)

    def time_from_str_ignore(self):
        pd.to_numeric(self.str, errors='ignore')

    def time_from_str_coerce(self):
        pd.to_numeric(self.str, errors='coerce')
