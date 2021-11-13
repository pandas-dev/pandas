import numpy as np

import pandas as pd


class UniqueForLargePyObjectInts:
    def setup(self):
        lst = [x << 32 for x in range(5000)]
        self.arr = np.array(lst, dtype=np.object_)

    def time_unique(self):
        pd.unique(self.arr)


class Float64GroupIndex:
    # GH28303
    def setup(self):
        self.df = pd.date_range(
            start="1/1/2018", end="1/2/2018", periods=10 ** 6
        ).to_frame()
        self.group_index = np.round(self.df.index.astype(int) / 10 ** 9)

    def time_groupby(self):
        self.df.groupby(self.group_index).last()


class UniqueAndFactorizeArange:
    params = range(4, 16)
    param_names = ["exponent"]

    def setup(self, exponent):
        a = np.arange(10 ** 4, dtype="float64")
        self.a2 = (a + 10 ** exponent).repeat(100)

    def time_factorize(self, exponent):
        pd.factorize(self.a2)

    def time_unique(self, exponent):
        pd.unique(self.a2)


class NumericSeriesIndexing:

    params = [
        (pd.Int64Index, pd.UInt64Index, pd.Float64Index),
        (10 ** 4, 10 ** 5, 5 * 10 ** 5, 10 ** 6, 5 * 10 ** 6),
    ]
    param_names = ["index_dtype", "N"]

    def setup(self, index, N):
        vals = np.array(list(range(55)) + [54] + list(range(55, N - 1)))
        indices = index(vals)
        self.data = pd.Series(np.arange(N), index=indices)

    def time_loc_slice(self, index, N):
        # trigger building of mapping
        self.data.loc[:800]


class NumericSeriesIndexingShuffled:

    params = [
        (pd.Int64Index, pd.UInt64Index, pd.Float64Index),
        (10 ** 4, 10 ** 5, 5 * 10 ** 5, 10 ** 6, 5 * 10 ** 6),
    ]
    param_names = ["index_dtype", "N"]

    def setup(self, index, N):
        vals = np.array(list(range(55)) + [54] + list(range(55, N - 1)))
        np.random.shuffle(vals)
        indices = index(vals)
        self.data = pd.Series(np.arange(N), index=indices)

    def time_loc_slice(self, index, N):
        # trigger building of mapping
        self.data.loc[:800]
