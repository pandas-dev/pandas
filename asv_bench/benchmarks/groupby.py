from functools import partial
from itertools import product
from string import ascii_letters

import numpy as np

from pandas import (
    Categorical,
    DataFrame,
    MultiIndex,
    Series,
    Timestamp,
    date_range,
    period_range,
)

from .pandas_vb_common import tm

method_blocklist = {
    "object": {
        "median",
        "prod",
        "sem",
        "cumsum",
        "sum",
        "cummin",
        "mean",
        "max",
        "skew",
        "cumprod",
        "cummax",
        "pct_change",
        "min",
        "var",
        "mad",
        "describe",
        "std",
        "quantile",
    },
    "datetime": {
        "median",
        "prod",
        "sem",
        "cumsum",
        "sum",
        "mean",
        "skew",
        "cumprod",
        "cummax",
        "pct_change",
        "var",
        "mad",
        "describe",
        "std",
    },
}


class ApplyDictReturn:
    def setup(self):
        self.labels = np.arange(1000).repeat(10)
        self.data = Series(np.random.randn(len(self.labels)))

    def time_groupby_apply_dict_return(self):
        self.data.groupby(self.labels).apply(
            lambda x: {"first": x.values[0], "last": x.values[-1]}
        )


class Apply:

    param_names = ["factor"]
    params = [4, 5]

    def setup(self, factor):
        N = 10 ** factor
        # two cases:
        # - small groups: small data (N**4) + many labels (2000) -> average group
        #   size of 5 (-> larger overhead of slicing method)
        # - larger groups: larger data (N**5) + fewer labels (20) -> average group
        #   size of 5000
        labels = np.random.randint(0, 2000 if factor == 4 else 20, size=N)
        labels2 = np.random.randint(0, 3, size=N)
        df = DataFrame(
            {
                "key": labels,
                "key2": labels2,
                "value1": np.random.randn(N),
                "value2": ["foo", "bar", "baz", "qux"] * (N // 4),
            }
        )
        self.df = df

    def time_scalar_function_multi_col(self, factor):
        self.df.groupby(["key", "key2"]).apply(lambda x: 1)

    def time_scalar_function_single_col(self, factor):
        self.df.groupby("key").apply(lambda x: 1)

    @staticmethod
    def df_copy_function(g):
        # ensure that the group name is available (see GH #15062)
        g.name
        return g.copy()

    def time_copy_function_multi_col(self, factor):
        self.df.groupby(["key", "key2"]).apply(self.df_copy_function)

    def time_copy_overhead_single_col(self, factor):
        self.df.groupby("key").apply(self.df_copy_function)


class Groups:

    param_names = ["key"]
    params = ["int64_small", "int64_large", "object_small", "object_large"]

    def setup_cache(self):
        size = 10 ** 6
        data = {
            "int64_small": Series(np.random.randint(0, 100, size=size)),
            "int64_large": Series(np.random.randint(0, 10000, size=size)),
            "object_small": Series(
                tm.makeStringIndex(100).take(np.random.randint(0, 100, size=size))
            ),
            "object_large": Series(
                tm.makeStringIndex(10000).take(np.random.randint(0, 10000, size=size))
            ),
        }
        return data

    def setup(self, data, key):
        self.ser = data[key]

    def time_series_groups(self, data, key):
        self.ser.groupby(self.ser).groups

    def time_series_indices(self, data, key):
        self.ser.groupby(self.ser).indices


class GroupManyLabels:

    params = [1, 1000]
    param_names = ["ncols"]

    def setup(self, ncols):
        N = 1000
        data = np.random.randn(N, ncols)
        self.labels = np.random.randint(0, 100, size=N)
        self.df = DataFrame(data)

    def time_sum(self, ncols):
        self.df.groupby(self.labels).sum()


class Nth:

    param_names = ["dtype"]
    params = ["float32", "float64", "datetime", "object"]

    def setup(self, dtype):
        N = 10 ** 5
        # with datetimes (GH7555)
        if dtype == "datetime":
            values = date_range("1/1/2011", periods=N, freq="s")
        elif dtype == "object":
            values = ["foo"] * N
        else:
            values = np.arange(N).astype(dtype)

        key = np.arange(N)
        self.df = DataFrame({"key": key, "values": values})
        self.df.iloc[1, 1] = np.nan  # insert missing data

    def time_frame_nth_any(self, dtype):
        self.df.groupby("key").nth(0, dropna="any")

    def time_groupby_nth_all(self, dtype):
        self.df.groupby("key").nth(0, dropna="all")

    def time_frame_nth(self, dtype):
        self.df.groupby("key").nth(0)

    def time_series_nth_any(self, dtype):
        self.df["values"].groupby(self.df["key"]).nth(0, dropna="any")

    def time_series_nth_all(self, dtype):
        self.df["values"].groupby(self.df["key"]).nth(0, dropna="all")

    def time_series_nth(self, dtype):
        self.df["values"].groupby(self.df["key"]).nth(0)


class DateAttributes:
    def setup(self):
        rng = date_range("1/1/2000", "12/31/2005", freq="H")
        self.year, self.month, self.day = rng.year, rng.month, rng.day
        self.ts = Series(np.random.randn(len(rng)), index=rng)

    def time_len_groupby_object(self):
        len(self.ts.groupby([self.year, self.month, self.day]))


class Int64:
    def setup(self):
        arr = np.random.randint(-1 << 12, 1 << 12, (1 << 17, 5))
        i = np.random.choice(len(arr), len(arr) * 5)
        arr = np.vstack((arr, arr[i]))
        i = np.random.permutation(len(arr))
        arr = arr[i]
        self.cols = list("abcde")
        self.df = DataFrame(arr, columns=self.cols)
        self.df["jim"], self.df["joe"] = np.random.randn(2, len(self.df)) * 10

    def time_overflow(self):
        self.df.groupby(self.cols).max()


class CountMultiDtype:
    def setup_cache(self):
        n = 10000
        offsets = np.random.randint(n, size=n).astype("timedelta64[ns]")
        dates = np.datetime64("now") + offsets
        dates[np.random.rand(n) > 0.5] = np.datetime64("nat")
        offsets[np.random.rand(n) > 0.5] = np.timedelta64("nat")
        value2 = np.random.randn(n)
        value2[np.random.rand(n) > 0.5] = np.nan
        obj = np.random.choice(list("ab"), size=n).astype(object)
        obj[np.random.randn(n) > 0.5] = np.nan
        df = DataFrame(
            {
                "key1": np.random.randint(0, 500, size=n),
                "key2": np.random.randint(0, 100, size=n),
                "dates": dates,
                "value2": value2,
                "value3": np.random.randn(n),
                "ints": np.random.randint(0, 1000, size=n),
                "obj": obj,
                "offsets": offsets,
            }
        )
        return df

    def time_multi_count(self, df):
        df.groupby(["key1", "key2"]).count()


class CountMultiInt:
    def setup_cache(self):
        n = 10000
        df = DataFrame(
            {
                "key1": np.random.randint(0, 500, size=n),
                "key2": np.random.randint(0, 100, size=n),
                "ints": np.random.randint(0, 1000, size=n),
                "ints2": np.random.randint(0, 1000, size=n),
            }
        )
        return df

    def time_multi_int_count(self, df):
        df.groupby(["key1", "key2"]).count()

    def time_multi_int_nunique(self, df):
        df.groupby(["key1", "key2"]).nunique()


class AggFunctions:
    def setup_cache(self):
        N = 10 ** 5
        fac1 = np.array(["A", "B", "C"], dtype="O")
        fac2 = np.array(["one", "two"], dtype="O")
        df = DataFrame(
            {
                "key1": fac1.take(np.random.randint(0, 3, size=N)),
                "key2": fac2.take(np.random.randint(0, 2, size=N)),
                "value1": np.random.randn(N),
                "value2": np.random.randn(N),
                "value3": np.random.randn(N),
            }
        )
        return df

    def time_different_str_functions(self, df):
        df.groupby(["key1", "key2"]).agg(
            {"value1": "mean", "value2": "var", "value3": "sum"}
        )

    def time_different_numpy_functions(self, df):
        df.groupby(["key1", "key2"]).agg(
            {"value1": np.mean, "value2": np.var, "value3": np.sum}
        )

    def time_different_python_functions_multicol(self, df):
        df.groupby(["key1", "key2"]).agg([sum, min, max])

    def time_different_python_functions_singlecol(self, df):
        df.groupby("key1").agg([sum, min, max])


class GroupStrings:
    def setup(self):
        n = 2 * 10 ** 5
        alpha = list(map("".join, product(ascii_letters, repeat=4)))
        data = np.random.choice(alpha, (n // 5, 4), replace=False)
        data = np.repeat(data, 5, axis=0)
        self.df = DataFrame(data, columns=list("abcd"))
        self.df["joe"] = (np.random.randn(len(self.df)) * 10).round(3)
        self.df = self.df.sample(frac=1).reset_index(drop=True)

    def time_multi_columns(self):
        self.df.groupby(list("abcd")).max()


class MultiColumn:
    def setup_cache(self):
        N = 10 ** 5
        key1 = np.tile(np.arange(100, dtype=object), 1000)
        key2 = key1.copy()
        np.random.shuffle(key1)
        np.random.shuffle(key2)
        df = DataFrame(
            {
                "key1": key1,
                "key2": key2,
                "data1": np.random.randn(N),
                "data2": np.random.randn(N),
            }
        )
        return df

    def time_lambda_sum(self, df):
        df.groupby(["key1", "key2"]).agg(lambda x: x.values.sum())

    def time_cython_sum(self, df):
        df.groupby(["key1", "key2"]).sum()

    def time_col_select_lambda_sum(self, df):
        df.groupby(["key1", "key2"])["data1"].agg(lambda x: x.values.sum())

    def time_col_select_numpy_sum(self, df):
        df.groupby(["key1", "key2"])["data1"].agg(np.sum)


class Size:
    def setup(self):
        n = 10 ** 5
        offsets = np.random.randint(n, size=n).astype("timedelta64[ns]")
        dates = np.datetime64("now") + offsets
        self.df = DataFrame(
            {
                "key1": np.random.randint(0, 500, size=n),
                "key2": np.random.randint(0, 100, size=n),
                "value1": np.random.randn(n),
                "value2": np.random.randn(n),
                "value3": np.random.randn(n),
                "dates": dates,
            }
        )
        self.draws = Series(np.random.randn(n))
        labels = Series(["foo", "bar", "baz", "qux"] * (n // 4))
        self.cats = labels.astype("category")

    def time_multi_size(self):
        self.df.groupby(["key1", "key2"]).size()

    def time_category_size(self):
        self.draws.groupby(self.cats).size()


class FillNA:
    def setup(self):
        N = 100
        self.df = DataFrame(
            {"group": [1] * N + [2] * N, "value": [np.nan, 1.0] * N}
        ).set_index("group")

    def time_df_ffill(self):
        self.df.groupby("group").fillna(method="ffill")

    def time_df_bfill(self):
        self.df.groupby("group").fillna(method="bfill")

    def time_srs_ffill(self):
        self.df.groupby("group")["value"].fillna(method="ffill")

    def time_srs_bfill(self):
        self.df.groupby("group")["value"].fillna(method="bfill")


class GroupByMethods:

    param_names = ["dtype", "method", "application"]
    params = [
        ["int", "float", "object", "datetime", "uint"],
        [
            "all",
            "any",
            "bfill",
            "count",
            "cumcount",
            "cummax",
            "cummin",
            "cumprod",
            "cumsum",
            "describe",
            "ffill",
            "first",
            "head",
            "last",
            "mad",
            "max",
            "min",
            "median",
            "mean",
            "nunique",
            "pct_change",
            "prod",
            "quantile",
            "rank",
            "sem",
            "shift",
            "size",
            "skew",
            "std",
            "sum",
            "tail",
            "unique",
            "value_counts",
            "var",
        ],
        ["direct", "transformation"],
    ]

    def setup(self, dtype, method, application):
        if method in method_blocklist.get(dtype, {}):
            raise NotImplementedError  # skip benchmark
        ngroups = 1000
        size = ngroups * 2
        rng = np.arange(ngroups)
        values = rng.take(np.random.randint(0, ngroups, size=size))
        if dtype == "int":
            key = np.random.randint(0, size, size=size)
        elif dtype == "uint":
            key = np.random.randint(0, size, size=size, dtype=dtype)
        elif dtype == "float":
            key = np.concatenate(
                [np.random.random(ngroups) * 0.1, np.random.random(ngroups) * 10.0]
            )
        elif dtype == "object":
            key = ["foo"] * size
        elif dtype == "datetime":
            key = date_range("1/1/2011", periods=size, freq="s")

        df = DataFrame({"values": values, "key": key})

        if application == "transform":
            if method == "describe":
                raise NotImplementedError

            self.as_group_method = lambda: df.groupby("key")["values"].transform(method)
            self.as_field_method = lambda: df.groupby("values")["key"].transform(method)
        else:
            self.as_group_method = getattr(df.groupby("key")["values"], method)
            self.as_field_method = getattr(df.groupby("values")["key"], method)

    def time_dtype_as_group(self, dtype, method, application):
        self.as_group_method()

    def time_dtype_as_field(self, dtype, method, application):
        self.as_field_method()


class GroupByCythonAgg:
    """
    Benchmarks specifically targetting our cython aggregation algorithms
    (using a big enough dataframe with simple key, so a large part of the
    time is actually spent in the grouped aggregation).
    """

    param_names = ["dtype", "method"]
    params = [
        ["float64"],
        [
            "sum",
            "prod",
            "min",
            "max",
            "mean",
            "median",
            "var",
            "first",
            "last",
            "any",
            "all",
        ],
    ]

    def setup(self, dtype, method):
        N = 1_000_000
        df = DataFrame(np.random.randn(N, 10), columns=list("abcdefghij"))
        df["key"] = np.random.randint(0, 100, size=N)
        self.df = df

    def time_frame_agg(self, dtype, method):
        self.df.groupby("key").agg(method)


class Cumulative:
    param_names = ["dtype", "method"]
    params = [
        ["float64", "int64", "Float64", "Int64"],
        ["cummin", "cummax", "cumsum"],
    ]

    def setup(self, dtype, method):
        N = 500_000
        vals = np.random.randint(-10, 10, (N, 5))
        null_vals = vals.astype(float, copy=True)
        null_vals[::2, :] = np.nan
        null_vals[::3, :] = np.nan
        df = DataFrame(vals, columns=list("abcde"), dtype=dtype)
        null_df = DataFrame(null_vals, columns=list("abcde"), dtype=dtype)
        keys = np.random.randint(0, 100, size=N)
        df["key"] = keys
        null_df["key"] = keys
        self.df = df
        self.null_df = null_df

    def time_frame_transform(self, dtype, method):
        self.df.groupby("key").transform(method)

    def time_frame_transform_many_nulls(self, dtype, method):
        self.null_df.groupby("key").transform(method)


class RankWithTies:
    # GH 21237
    param_names = ["dtype", "tie_method"]
    params = [
        ["float64", "float32", "int64", "datetime64"],
        ["first", "average", "dense", "min", "max"],
    ]

    def setup(self, dtype, tie_method):
        N = 10 ** 4
        if dtype == "datetime64":
            data = np.array([Timestamp("2011/01/01")] * N, dtype=dtype)
        else:
            data = np.array([1] * N, dtype=dtype)
        self.df = DataFrame({"values": data, "key": ["foo"] * N})

    def time_rank_ties(self, dtype, tie_method):
        self.df.groupby("key").rank(method=tie_method)


class Float32:
    # GH 13335
    def setup(self):
        tmp1 = (np.random.random(10000) * 0.1).astype(np.float32)
        tmp2 = (np.random.random(10000) * 10.0).astype(np.float32)
        tmp = np.concatenate((tmp1, tmp2))
        arr = np.repeat(tmp, 10)
        self.df = DataFrame({"a": arr, "b": arr})

    def time_sum(self):
        self.df.groupby(["a"])["b"].sum()


class String:
    # GH#41596
    param_names = ["dtype", "method"]
    params = [
        ["str", "string[python]"],
        [
            "sum",
            "prod",
            "min",
            "max",
            "mean",
            "median",
            "var",
            "first",
            "last",
            "any",
            "all",
        ],
    ]

    def setup(self, dtype, method):
        cols = list("abcdefghjkl")
        self.df = DataFrame(
            np.random.randint(0, 100, size=(1_000_000, len(cols))),
            columns=cols,
            dtype=dtype,
        )

    def time_str_func(self, dtype, method):
        self.df.groupby("a")[self.df.columns[1:]].agg(method)


class Categories:
    def setup(self):
        N = 10 ** 5
        arr = np.random.random(N)
        data = {"a": Categorical(np.random.randint(10000, size=N)), "b": arr}
        self.df = DataFrame(data)
        data = {
            "a": Categorical(np.random.randint(10000, size=N), ordered=True),
            "b": arr,
        }
        self.df_ordered = DataFrame(data)
        data = {
            "a": Categorical(
                np.random.randint(100, size=N), categories=np.arange(10000)
            ),
            "b": arr,
        }
        self.df_extra_cat = DataFrame(data)

    def time_groupby_sort(self):
        self.df.groupby("a")["b"].count()

    def time_groupby_nosort(self):
        self.df.groupby("a", sort=False)["b"].count()

    def time_groupby_ordered_sort(self):
        self.df_ordered.groupby("a")["b"].count()

    def time_groupby_ordered_nosort(self):
        self.df_ordered.groupby("a", sort=False)["b"].count()

    def time_groupby_extra_cat_sort(self):
        self.df_extra_cat.groupby("a")["b"].count()

    def time_groupby_extra_cat_nosort(self):
        self.df_extra_cat.groupby("a", sort=False)["b"].count()


class Datelike:
    # GH 14338
    params = ["period_range", "date_range", "date_range_tz"]
    param_names = ["grouper"]

    def setup(self, grouper):
        N = 10 ** 4
        rng_map = {
            "period_range": period_range,
            "date_range": date_range,
            "date_range_tz": partial(date_range, tz="US/Central"),
        }
        self.grouper = rng_map[grouper]("1900-01-01", freq="D", periods=N)
        self.df = DataFrame(np.random.randn(10 ** 4, 2))

    def time_sum(self, grouper):
        self.df.groupby(self.grouper).sum()


class SumBools:
    # GH 2692
    def setup(self):
        N = 500
        self.df = DataFrame({"ii": range(N), "bb": [True] * N})

    def time_groupby_sum_booleans(self):
        self.df.groupby("ii").sum()


class SumMultiLevel:
    # GH 9049
    timeout = 120.0

    def setup(self):
        N = 50
        self.df = DataFrame(
            {"A": list(range(N)) * 2, "B": range(N * 2), "C": 1}
        ).set_index(["A", "B"])

    def time_groupby_sum_multiindex(self):
        self.df.groupby(level=[0, 1]).sum()


class Transform:
    def setup(self):
        n1 = 400
        n2 = 250
        index = MultiIndex(
            levels=[np.arange(n1), tm.makeStringIndex(n2)],
            codes=[np.repeat(range(n1), n2).tolist(), list(range(n2)) * n1],
            names=["lev1", "lev2"],
        )
        arr = np.random.randn(n1 * n2, 3)
        arr[::10000, 0] = np.nan
        arr[1::10000, 1] = np.nan
        arr[2::10000, 2] = np.nan
        data = DataFrame(arr, index=index, columns=["col1", "col20", "col3"])
        self.df = data

        n = 20000
        self.df1 = DataFrame(
            np.random.randint(1, n, (n, 3)), columns=["jim", "joe", "jolie"]
        )
        self.df2 = self.df1.copy()
        self.df2["jim"] = self.df2["joe"]

        self.df3 = DataFrame(
            np.random.randint(1, (n / 10), (n, 3)), columns=["jim", "joe", "jolie"]
        )
        self.df4 = self.df3.copy()
        self.df4["jim"] = self.df4["joe"]

    def time_transform_lambda_max(self):
        self.df.groupby(level="lev1").transform(lambda x: max(x))

    def time_transform_ufunc_max(self):
        self.df.groupby(level="lev1").transform(np.max)

    def time_transform_multi_key1(self):
        self.df1.groupby(["jim", "joe"])["jolie"].transform("max")

    def time_transform_multi_key2(self):
        self.df2.groupby(["jim", "joe"])["jolie"].transform("max")

    def time_transform_multi_key3(self):
        self.df3.groupby(["jim", "joe"])["jolie"].transform("max")

    def time_transform_multi_key4(self):
        self.df4.groupby(["jim", "joe"])["jolie"].transform("max")


class TransformBools:
    def setup(self):
        N = 120000
        transition_points = np.sort(np.random.choice(np.arange(N), 1400))
        transitions = np.zeros(N, dtype=np.bool_)
        transitions[transition_points] = True
        self.g = transitions.cumsum()
        self.df = DataFrame({"signal": np.random.rand(N)})

    def time_transform_mean(self):
        self.df["signal"].groupby(self.g).transform(np.mean)


class TransformNaN:
    # GH 12737
    def setup(self):
        self.df_nans = DataFrame(
            {"key": np.repeat(np.arange(1000), 10), "B": np.nan, "C": np.nan}
        )
        self.df_nans.loc[4::10, "B":"C"] = 5

    def time_first(self):
        self.df_nans.groupby("key").transform("first")


class TransformEngine:

    param_names = ["parallel"]
    params = [[True, False]]

    def setup(self, parallel):
        N = 10 ** 3
        data = DataFrame(
            {0: [str(i) for i in range(100)] * N, 1: list(range(100)) * N},
            columns=[0, 1],
        )
        self.parallel = parallel
        self.grouper = data.groupby(0)

    def time_series_numba(self, parallel):
        def function(values, index):
            return values * 5

        self.grouper[1].transform(
            function, engine="numba", engine_kwargs={"parallel": self.parallel}
        )

    def time_series_cython(self, parallel):
        def function(values):
            return values * 5

        self.grouper[1].transform(function, engine="cython")

    def time_dataframe_numba(self, parallel):
        def function(values, index):
            return values * 5

        self.grouper.transform(
            function, engine="numba", engine_kwargs={"parallel": self.parallel}
        )

    def time_dataframe_cython(self, parallel):
        def function(values):
            return values * 5

        self.grouper.transform(function, engine="cython")


class AggEngine:

    param_names = ["parallel"]
    params = [[True, False]]

    def setup(self, parallel):
        N = 10 ** 3
        data = DataFrame(
            {0: [str(i) for i in range(100)] * N, 1: list(range(100)) * N},
            columns=[0, 1],
        )
        self.parallel = parallel
        self.grouper = data.groupby(0)

    def time_series_numba(self, parallel):
        def function(values, index):
            total = 0
            for i, value in enumerate(values):
                if i % 2:
                    total += value + 5
                else:
                    total += value * 2
            return total

        self.grouper[1].agg(
            function, engine="numba", engine_kwargs={"parallel": self.parallel}
        )

    def time_series_cython(self, parallel):
        def function(values):
            total = 0
            for i, value in enumerate(values):
                if i % 2:
                    total += value + 5
                else:
                    total += value * 2
            return total

        self.grouper[1].agg(function, engine="cython")

    def time_dataframe_numba(self, parallel):
        def function(values, index):
            total = 0
            for i, value in enumerate(values):
                if i % 2:
                    total += value + 5
                else:
                    total += value * 2
            return total

        self.grouper.agg(
            function, engine="numba", engine_kwargs={"parallel": self.parallel}
        )

    def time_dataframe_cython(self, parallel):
        def function(values):
            total = 0
            for i, value in enumerate(values):
                if i % 2:
                    total += value + 5
                else:
                    total += value * 2
            return total

        self.grouper.agg(function, engine="cython")


from .pandas_vb_common import setup  # noqa: F401 isort:skip
