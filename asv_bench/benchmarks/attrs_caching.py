import numpy as np

import pandas as pd
from pandas import DataFrame

try:
    from pandas.util import cache_readonly
except ImportError:
    from pandas.util.decorators import cache_readonly

try:
    from pandas.core.construction import extract_array
except ImportError:
    extract_array = None


class DataFrameAttributes:
    def setup(self):
        self.df = DataFrame(np.random.randn(10, 6))
        self.cur_index = self.df.index

    def time_get_index(self):
        self.foo = self.df.index

    def time_set_index(self):
        self.df.index = self.cur_index


class SeriesArrayAttribute:

    params = [["numeric", "object", "category", "datetime64", "datetime64tz"]]
    param_names = ["dtype"]

    def setup(self, dtype):
        if dtype == "numeric":
            self.series = pd.Series([1, 2, 3])
        elif dtype == "object":
            self.series = pd.Series(["a", "b", "c"], dtype=object)
        elif dtype == "category":
            self.series = pd.Series(["a", "b", "c"], dtype="category")
        elif dtype == "datetime64":
            self.series = pd.Series(pd.date_range("2013", periods=3))
        elif dtype == "datetime64tz":
            self.series = pd.Series(pd.date_range("2013", periods=3, tz="UTC"))

    def time_array(self, dtype):
        self.series.array

    def time_extract_array(self, dtype):
        extract_array(self.series)

    def time_extract_array_numpy(self, dtype):
        extract_array(self.series, extract_numpy=True)


class CacheReadonly:
    def setup(self):
        class Foo:
            @cache_readonly
            def prop(self):
                return 5

        self.obj = Foo()

    def time_cache_readonly(self):
        self.obj.prop


from .pandas_vb_common import setup  # noqa: F401 isort:skip
