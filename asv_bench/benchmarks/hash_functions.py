import numpy as np

import pandas as pd


class IsinAlmostFullWithRandomInt:
    params = [
        [np.float64, np.int64, np.uint64, np.object],
        range(10, 21),
    ]
    param_names = ["dtype", "exponent"]

    def setup(self, dtype, exponent):
        M = 3 * 2 ** (exponent - 2)
        # 0.77-the maximal share of occupied buckets
        np.random.seed(42)
        self.s = pd.Series(np.random.randint(0, M, M)).astype(dtype)
        self.values = np.random.randint(0, M, M).astype(dtype)
        self.values_outside = self.values + M

    def time_isin(self, dtype, exponent):
        self.s.isin(self.values)

    def time_isin_outside(self, dtype, exponent):
        self.s.isin(self.values_outside)


class IsinWithRandomFloat:
    params = [
        [np.float64, np.object],
        [
            1_300,
            2_000,
            7_000,
            8_000,
            70_000,
            80_000,
            750_000,
            900_000,
        ],
    ]
    param_names = ["dtype", "M"]

    def setup(self, dtype, M):
        np.random.seed(42)
        self.values = np.random.rand(M)
        self.s = pd.Series(self.values).astype(dtype)
        np.random.shuffle(self.values)
        self.values_outside = self.values + 0.1

    def time_isin(self, dtype, M):
        self.s.isin(self.values)

    def time_isin_outside(self, dtype, M):
        self.s.isin(self.values_outside)


class IsinWithArangeSorted:
    params = [
        [np.float64, np.int64, np.uint64, np.object],
        [
            1_000,
            2_000,
            8_000,
            100_000,
            1_000_000,
        ],
    ]
    param_names = ["dtype", "M"]

    def setup(self, dtype, M):
        self.s = pd.Series(np.arange(M)).astype(dtype)
        self.values = np.arange(M).astype(dtype)

    def time_isin(self, dtype, M):
        self.s.isin(self.values)


class IsinWithArange:
    params = [
        [np.float64, np.int64, np.uint64, np.object],
        [
            1_000,
            2_000,
            8_000,
        ],
        [-2, 0, 2],
    ]
    param_names = ["dtype", "M", "offset_factor"]

    def setup(self, dtype, M, offset_factor):
        offset = int(M * offset_factor)
        np.random.seed(42)
        tmp = pd.Series(np.random.randint(offset, M + offset, 10 ** 6))
        self.s = tmp.astype(dtype)
        self.values = np.arange(M).astype(dtype)

    def time_isin(self, dtype, M, offset_factor):
        self.s.isin(self.values)


class Float64GroupIndex:
    # GH28303
    def setup(self):
        self.df = pd.date_range(
            start="1/1/2018", end="1/2/2018", periods=1e6
        ).to_frame()
        self.group_index = np.round(self.df.index.astype(int) / 1e9)

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
        np.random.seed(42)
        np.random.shuffle(vals)
        indices = index(vals)
        self.data = pd.Series(np.arange(N), index=indices)

    def time_loc_slice(self, index, N):
        # trigger building of mapping
        self.data.loc[:800]
