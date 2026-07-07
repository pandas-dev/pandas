import numpy as np

import pandas as pd
from pandas import DataFrame

try:
    from pandas.core.construction import extract_array
except ImportError:
    extract_array = None


class DataFrameAttributes:
    def setup(self):
        self.df = DataFrame(np.random.randn(10, 6))
        self.cur_index = self.df.index

    def time_get_index(self):
        self.df.index

    def time_set_index(self):
        self.df.index = self.cur_index


class Dir:
    # GH#18587 measure uncached tab-completion cost over a large index.
    number = 1
    repeat = (3, 100, 20)

    params = ["unique", "non_unique"]
    param_names = ["kind"]

    def setup(self, kind):
        N = 10**5
        if kind == "unique":
            labels = [f"col_{i}" for i in range(N)]
        else:
            labels = [f"col_{i % 50}" for i in range(N)]
        self.df = DataFrame(np.zeros((1, N)), columns=pd.Index(labels))
        self.ser = pd.Series(np.zeros(N), index=pd.Index(labels))

    def time_dir_frame(self, kind):
        self.df.columns._cache = {}
        dir(self.df)

    def time_dir_series(self, kind):
        self.ser.index._cache = {}
        dir(self.ser)


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


from .pandas_vb_common import setup  # noqa: F401 isort:skip
