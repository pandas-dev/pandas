import numpy as np

import pandas as pd

ops = ["mean", "sum", "median", "std", "skew", "kurt", "prod", "sem", "var"]


class FrameOps:
    params = [ops, ["float", "int", "Int64"], [0, 1, None]]
    param_names = ["op", "dtype", "axis"]

    def setup(self, op, dtype, axis):
        values = np.random.randn(100000, 4)
        if dtype == "Int64":
            values = values.astype(int)
        df = pd.DataFrame(values).astype(dtype)
        self.df_func = getattr(df, op)

    def time_op(self, op, dtype, axis):
        self.df_func(axis=axis)


class FrameMixedDtypesOps:
    params = [ops, [0, 1, None]]
    param_names = ["op", "axis"]

    def setup(self, op, axis):
        if op in ("sum", "skew", "kurt", "prod", "sem", "var") or (
            (op, axis)
            in (
                ("mean", 1),
                ("mean", None),
                ("median", 1),
                ("median", None),
                ("std", 1),
                ("std", None),
            )
        ):
            # Skipping cases where datetime aggregations are not implemented
            raise NotImplementedError

        N = 1_000_000
        df = pd.DataFrame(
            {
                "f": np.random.normal(0.0, 1.0, N),
                "i": np.random.randint(0, N, N),
                "ts": pd.date_range(start="1/1/2000", periods=N, freq="h"),
            }
        )

        self.df_func = getattr(df, op)

    def time_op(self, op, axis):
        self.df_func(axis=axis)


class FrameAxis1BlockFusion:
    params = [
        ["sum", "prod", "min", "max"],
        [2048, 32_768, 100_000],
        [2, 8],
        ["float64", "mixed_float", "int64"],
    ]
    param_names = ["op", "nrows", "nblocks", "dtype_mix"]

    def setup(self, op, nrows, nblocks, dtype_mix):
        rng = np.random.default_rng(2)
        frame = pd.DataFrame(index=range(nrows))
        for i in range(nblocks):
            if dtype_mix == "mixed_float":
                dtype = "float32" if i % 2 == 0 else "float64"
            else:
                dtype = dtype_mix

            if dtype == "int64":
                values = ((np.arange(nrows) + i) % 3).astype(dtype)
            else:
                values = rng.uniform(0.9, 1.1, nrows).astype(dtype)
                values[i::1000] = np.nan
            frame[f"column-{i}"] = values

        assert len(frame._mgr.blocks) == nblocks
        self.frame_op = getattr(frame, op)

    def time_axis1(self, op, nrows, nblocks, dtype_mix):
        self.frame_op(axis=1)


class FrameAxis1BlockFusionPeakMemory:
    params = [["sum", "min"]]
    param_names = ["op"]

    def setup(self, op):
        rng = np.random.default_rng(2)
        frame = pd.DataFrame(index=range(1_000_000))
        for i in range(8):
            values = rng.uniform(0.9, 1.1, len(frame))
            values[i::1000] = np.nan
            frame[f"column-{i}"] = values

        assert len(frame._mgr.blocks) == 8
        self.frame_op = getattr(frame, op)

    def peakmem_axis1(self, op):
        self.frame_op(axis=1)


class FrameMultiIndexOps:
    params = [ops]
    param_names = ["op"]

    def setup(self, op):
        levels = [np.arange(10), np.arange(100), np.arange(100)]
        codes = [
            np.arange(10).repeat(10000),
            np.tile(np.arange(100).repeat(100), 10),
            np.tile(np.tile(np.arange(100), 100), 10),
        ]
        index = pd.MultiIndex(levels=levels, codes=codes)
        df = pd.DataFrame(np.random.randn(len(index), 4), index=index)
        self.df_func = getattr(df, op)

    def time_op(self, op):
        self.df_func()


class SeriesOps:
    params = [ops, ["float", "int"]]
    param_names = ["op", "dtype"]

    def setup(self, op, dtype):
        s = pd.Series(np.random.randn(100000)).astype(dtype)
        self.s_func = getattr(s, op)

    def time_op(self, op, dtype):
        self.s_func()


class SeriesMultiIndexOps:
    params = [ops]
    param_names = ["op"]

    def setup(self, op):
        levels = [np.arange(10), np.arange(100), np.arange(100)]
        codes = [
            np.arange(10).repeat(10000),
            np.tile(np.arange(100).repeat(100), 10),
            np.tile(np.tile(np.arange(100), 100), 10),
        ]
        index = pd.MultiIndex(levels=levels, codes=codes)
        s = pd.Series(np.random.randn(len(index)), index=index)
        self.s_func = getattr(s, op)

    def time_op(self, op):
        self.s_func()


class Rank:
    params = [["DataFrame", "Series"], [True, False]]
    param_names = ["constructor", "pct"]

    def setup(self, constructor, pct):
        values = np.random.randn(10**5)
        self.data = getattr(pd, constructor)(values)

    def time_rank(self, constructor, pct):
        self.data.rank(pct=pct)

    def time_average_old(self, constructor, pct):
        self.data.rank(pct=pct) / len(self.data)


class Correlation:
    params = [["spearman", "kendall", "pearson"]]
    param_names = ["method"]

    def setup(self, method):
        self.df = pd.DataFrame(np.random.randn(500, 15))
        self.df2 = pd.DataFrame(np.random.randn(500, 15))
        self.df_wide = pd.DataFrame(np.random.randn(500, 100))
        self.df_wide_nans = self.df_wide.where(np.random.random((500, 100)) < 0.9)
        self.s = pd.Series(np.random.randn(500))
        self.s2 = pd.Series(np.random.randn(500))

    def time_corr(self, method):
        self.df.corr(method=method)

    def time_corr_wide(self, method):
        self.df_wide.corr(method=method)

    def time_corr_wide_nans(self, method):
        self.df_wide_nans.corr(method=method)

    def peakmem_corr_wide(self, method):
        self.df_wide.corr(method=method)

    def time_corr_series(self, method):
        self.s.corr(self.s2, method=method)

    def time_corrwith_cols(self, method):
        self.df.corrwith(self.df2, method=method)

    def time_corrwith_rows(self, method):
        self.df.corrwith(self.df2, axis=1, method=method)


class PearsonCorrelation:
    params = [
        ["clean", "constant", "sorted_nans", "unsorted_nans"],
        [1_000, 1_000_000],
    ]
    param_names = ["scenario", "size"]

    def setup(self, scenario, size):
        rng = np.random.default_rng(2)
        left = rng.standard_normal(size)
        right = rng.standard_normal(size)

        if scenario == "constant":
            left[:] = 1.0
        elif scenario == "sorted_nans":
            left[: size // 10] = np.nan
            right[-size // 10 :] = np.nan
        elif scenario == "unsorted_nans":
            left[size // 2] = np.nan
            right[size // 3] = np.nan

        self.s = pd.Series(left)
        self.s2 = pd.Series(right)

    def time_corr_series(self, scenario, size):
        self.s.corr(self.s2, method="pearson")


class Covariance:
    params = []
    param_names = []

    def setup(self):
        self.s = pd.Series(np.random.randn(100000))
        self.s2 = pd.Series(np.random.randn(100000))

    def time_cov_series(self):
        self.s.cov(self.s2)


from .pandas_vb_common import setup  # noqa: F401 isort:skip
