import string
import warnings

import numpy as np

from pandas import (
    DataFrame,
    Index,
    MultiIndex,
    NaT,
    Series,
    date_range,
    isnull,
    period_range,
    timedelta_range,
)


class AsType:
    params = [
        [
            # from_dtype == to_dtype
            ("Float64", "Float64"),
            ("float64[pyarrow]", "float64[pyarrow]"),
            # from non-EA to EA
            ("float64", "Float64"),
            ("float64", "float64[pyarrow]"),
            # from EA to non-EA
            ("Float64", "float64"),
            ("float64[pyarrow]", "float64"),
            # from EA to EA
            ("Int64", "Float64"),
            ("int64[pyarrow]", "float64[pyarrow]"),
        ],
        [False, True],
    ]
    param_names = ["from_to_dtypes", "copy"]

    def setup(self, from_to_dtypes, copy):
        from_dtype = from_to_dtypes[0]
        if from_dtype in ("float64", "Float64", "float64[pyarrow]"):
            data = np.random.randn(100, 100)
        elif from_dtype in ("int64", "Int64", "int64[pyarrow]"):
            data = np.random.randint(0, 1000, (100, 100))
        else:
            raise NotImplementedError
        self.df = DataFrame(data, dtype=from_dtype)

    def time_astype(self, from_to_dtypes, copy):
        self.df.astype(from_to_dtypes[1], copy=copy)


class Clip:
    params = [
        ["float64", "Float64", "float64[pyarrow]"],
    ]
    param_names = ["dtype"]

    def setup(self, dtype):
        data = np.random.randn(100_000, 10)
        df = DataFrame(data, dtype=dtype)
        self.df = df

    def time_clip(self, dtype):
        self.df.clip(-1.0, 1.0)


class GetNumericData:
    def setup(self):
        self.df = DataFrame(np.random.randn(10000, 25))
        self.df["foo"] = "bar"
        self.df["bar"] = "baz"
        self.df = self.df._consolidate()

    def time_frame_get_numeric_data(self):
        self.df._get_numeric_data()


class Reindex:
    def setup(self):
        N = 10**3
        self.df = DataFrame(np.random.randn(N * 10, N))
        self.idx = np.arange(4 * N, 7 * N)
        self.idx_cols = np.random.randint(0, N, N)
        self.df2 = DataFrame(
            {
                c: {
                    0: np.random.randint(0, 2, N).astype(np.bool_),
                    1: np.random.randint(0, N, N).astype(np.int16),
                    2: np.random.randint(0, N, N).astype(np.int32),
                    3: np.random.randint(0, N, N).astype(np.int64),
                }[np.random.randint(0, 4)]
                for c in range(N)
            }
        )

    def time_reindex_axis0(self):
        self.df.reindex(self.idx)

    def time_reindex_axis1(self):
        self.df.reindex(columns=self.idx_cols)

    def time_reindex_axis1_missing(self):
        self.df.reindex(columns=self.idx)

    def time_reindex_both_axes(self):
        self.df.reindex(index=self.idx, columns=self.idx_cols)

    def time_reindex_upcast(self):
        self.df2.reindex(np.random.permutation(range(1200)))


class Rename:
    def setup(self):
        N = 10**3
        self.df = DataFrame(np.random.randn(N * 10, N))
        self.idx = np.arange(4 * N, 7 * N)
        self.dict_idx = {k: k for k in self.idx}
        self.df2 = DataFrame(
            {
                c: {
                    0: np.random.randint(0, 2, N).astype(np.bool_),
                    1: np.random.randint(0, N, N).astype(np.int16),
                    2: np.random.randint(0, N, N).astype(np.int32),
                    3: np.random.randint(0, N, N).astype(np.int64),
                }[np.random.randint(0, 4)]
                for c in range(N)
            }
        )

    def time_rename_single(self):
        self.df.rename({0: 0})

    def time_rename_axis0(self):
        self.df.rename(self.dict_idx)

    def time_rename_axis1(self):
        self.df.rename(columns=self.dict_idx)

    def time_rename_both_axes(self):
        self.df.rename(index=self.dict_idx, columns=self.dict_idx)

    def time_dict_rename_both_axes(self):
        self.df.rename(index=self.dict_idx, columns=self.dict_idx)


class Iteration:
    # mem_itertuples_* benchmarks are slow
    timeout = 120

    def setup(self):
        N = 1000
        self.df = DataFrame(np.random.randn(N * 10, N))
        self.df2 = DataFrame(np.random.randn(N * 50, 10))
        self.df3 = DataFrame(
            np.random.randn(N, 5 * N), columns=["C" + str(c) for c in range(N * 5)]
        )
        self.df4 = DataFrame(np.random.randn(N * 1000, 10))

    def time_items(self):
        # (monitor no-copying behaviour)
        for name, col in self.df.items():
            pass

    def time_iteritems_indexing(self):
        for col in self.df3:
            self.df3[col]

    def time_itertuples_start(self):
        self.df4.itertuples()

    def time_itertuples_read_first(self):
        next(self.df4.itertuples())

    def time_itertuples(self):
        for row in self.df4.itertuples():
            pass

    def time_itertuples_to_list(self):
        list(self.df4.itertuples())

    def mem_itertuples_start(self):
        return self.df4.itertuples()

    def peakmem_itertuples_start(self):
        self.df4.itertuples()

    def mem_itertuples_read_first(self):
        return next(self.df4.itertuples())

    def peakmem_itertuples(self):
        for row in self.df4.itertuples():
            pass

    def mem_itertuples_to_list(self):
        return list(self.df4.itertuples())

    def peakmem_itertuples_to_list(self):
        list(self.df4.itertuples())

    def time_itertuples_raw_start(self):
        self.df4.itertuples(index=False, name=None)

    def time_itertuples_raw_read_first(self):
        next(self.df4.itertuples(index=False, name=None))

    def time_itertuples_raw_tuples(self):
        for row in self.df4.itertuples(index=False, name=None):
            pass

    def time_itertuples_raw_tuples_to_list(self):
        list(self.df4.itertuples(index=False, name=None))

    def mem_itertuples_raw_start(self):
        return self.df4.itertuples(index=False, name=None)

    def peakmem_itertuples_raw_start(self):
        self.df4.itertuples(index=False, name=None)

    def peakmem_itertuples_raw_read_first(self):
        next(self.df4.itertuples(index=False, name=None))

    def peakmem_itertuples_raw(self):
        for row in self.df4.itertuples(index=False, name=None):
            pass

    def mem_itertuples_raw_to_list(self):
        return list(self.df4.itertuples(index=False, name=None))

    def peakmem_itertuples_raw_to_list(self):
        list(self.df4.itertuples(index=False, name=None))

    def time_iterrows(self):
        for row in self.df.iterrows():
            pass


class ToString:
    def setup(self):
        self.df = DataFrame(np.random.randn(100, 10))

    def time_to_string_floats(self):
        self.df.to_string()


class ToHTML:
    def setup(self):
        nrows = 500
        self.df2 = DataFrame(np.random.randn(nrows, 10))
        self.df2[0] = period_range("2000", periods=nrows)
        self.df2[1] = range(nrows)

    def time_to_html_mixed(self):
        self.df2.to_html()


class ToDict:
    params = [["dict", "list", "series", "split", "records", "index"]]
    param_names = ["orient"]

    def setup(self, orient):
        data = np.random.randint(0, 1000, size=(10000, 4))
        self.int_df = DataFrame(data)
        self.datetimelike_df = self.int_df.astype("timedelta64[ns]")

    def time_to_dict_ints(self, orient):
        self.int_df.to_dict(orient=orient)

    def time_to_dict_datetimelike(self, orient):
        self.datetimelike_df.to_dict(orient=orient)


class ToNumpy:
    def setup(self):
        N = 10000
        M = 10
        self.df_tall = DataFrame(np.random.randn(N, M))
        self.df_wide = DataFrame(np.random.randn(M, N))
        self.df_mixed_tall = self.df_tall.copy()
        self.df_mixed_tall["foo"] = "bar"
        self.df_mixed_tall[0] = period_range("2000", periods=N)
        self.df_mixed_tall[1] = range(N)
        self.df_mixed_wide = self.df_wide.copy()
        self.df_mixed_wide["foo"] = "bar"
        self.df_mixed_wide[0] = period_range("2000", periods=M)
        self.df_mixed_wide[1] = range(M)

    def time_to_numpy_tall(self):
        self.df_tall.to_numpy()

    def time_to_numpy_wide(self):
        self.df_wide.to_numpy()

    def time_to_numpy_mixed_tall(self):
        self.df_mixed_tall.to_numpy()

    def time_to_numpy_mixed_wide(self):
        self.df_mixed_wide.to_numpy()

    def time_values_tall(self):
        self.df_tall.values

    def time_values_wide(self):
        self.df_wide.values

    def time_values_mixed_tall(self):
        self.df_mixed_tall.values

    def time_values_mixed_wide(self):
        self.df_mixed_wide.values


class ToRecords:
    def setup(self):
        N = 100_000
        data = np.random.randn(N, 2)
        mi = MultiIndex.from_arrays(
            [
                np.arange(N),
                date_range("1970-01-01", periods=N, freq="ms"),
            ]
        )
        self.df = DataFrame(data)
        self.df_mi = DataFrame(data, index=mi)

    def time_to_records(self):
        self.df.to_records(index=True)

    def time_to_records_multiindex(self):
        self.df_mi.to_records(index=True)


class Repr:
    def setup(self):
        nrows = 10000
        data = np.random.randn(nrows, 10)
        arrays = np.tile(np.random.randn(3, nrows // 100), 100)
        idx = MultiIndex.from_arrays(arrays)
        self.df3 = DataFrame(data, index=idx)
        self.df4 = DataFrame(data, index=np.random.randn(nrows))
        self.df_tall = DataFrame(np.random.randn(nrows, 10))
        self.df_wide = DataFrame(np.random.randn(10, nrows))

    def time_html_repr_trunc_mi(self):
        self.df3._repr_html_()

    def time_html_repr_trunc_si(self):
        self.df4._repr_html_()

    def time_repr_tall(self):
        repr(self.df_tall)

    def time_frame_repr_wide(self):
        repr(self.df_wide)


class MaskBool:
    def setup(self):
        data = np.random.randn(1000, 500)
        df = DataFrame(data)
        df = df.where(df > 0)
        self.bools = df > 0
        self.mask = isnull(df)

    def time_frame_mask_bools(self):
        self.bools.mask(self.mask)

    def time_frame_mask_floats(self):
        self.bools.astype(float).mask(self.mask)


class Isnull:
    def setup(self):
        N = 10**3
        self.df_no_null = DataFrame(np.random.randn(N, N))

        sample = np.array([np.nan, 1.0])
        data = np.random.choice(sample, (N, N))
        self.df = DataFrame(data)

        sample = np.array(list(string.ascii_letters + string.whitespace))
        data = np.random.choice(sample, (N, N))
        self.df_strings = DataFrame(data)

        sample = np.array(
            [
                NaT,
                np.nan,
                None,
                np.datetime64("NaT"),
                np.timedelta64("NaT"),
                0,
                1,
                2.0,
                "",
                "abcd",
            ]
        )
        data = np.random.choice(sample, (N, N))
        self.df_obj = DataFrame(data)

    def time_isnull_floats_no_null(self):
        isnull(self.df_no_null)

    def time_isnull(self):
        isnull(self.df)

    def time_isnull_strngs(self):
        isnull(self.df_strings)

    def time_isnull_obj(self):
        isnull(self.df_obj)


class Fillna:
    params = (
        [True, False],
        [
            "float64",
            "float32",
            "object",
            "Int64",
            "Float64",
            "datetime64[ns]",
            "datetime64[ns, tz]",
            "timedelta64[ns]",
        ],
    )
    param_names = ["inplace", "dtype"]

    def setup(self, inplace, dtype):
        N, M = 10000, 100
        if dtype in ("datetime64[ns]", "datetime64[ns, tz]", "timedelta64[ns]"):
            data = {
                "datetime64[ns]": date_range("2011-01-01", freq="h", periods=N),
                "datetime64[ns, tz]": date_range(
                    "2011-01-01", freq="h", periods=N, tz="Asia/Tokyo"
                ),
                "timedelta64[ns]": timedelta_range(start="1 day", periods=N, freq="1D"),
            }
            self.df = DataFrame({f"col_{i}": data[dtype] for i in range(M)})
            self.df[::2] = None
        else:
            values = np.random.randn(N, M)
            values[::2] = np.nan
            if dtype == "Int64":
                values = values.round()
            self.df = DataFrame(values, dtype=dtype)
        self.fill_values = self.df.iloc[self.df.first_valid_index()].to_dict()

    def time_fillna(self, inplace, dtype):
        self.df.fillna(value=self.fill_values, inplace=inplace)

    def time_ffill(self, inplace, dtype):
        self.df.ffill(inplace=inplace)

    def time_bfill(self, inplace, dtype):
        self.df.bfill(inplace=inplace)


class Dropna:
    params = (["all", "any"], [0, 1])
    param_names = ["how", "axis"]

    def setup(self, how, axis):
        self.df = DataFrame(np.random.randn(10000, 1000))
        self.df.iloc[50:1000, 20:50] = np.nan
        self.df.iloc[2000:3000] = np.nan
        self.df.iloc[:, 60:70] = np.nan
        self.df_mixed = self.df.copy()
        self.df_mixed["foo"] = "bar"

    def time_dropna(self, how, axis):
        self.df.dropna(how=how, axis=axis)

    def time_dropna_axis_mixed_dtypes(self, how, axis):
        self.df_mixed.dropna(how=how, axis=axis)


class Isna:
    params = ["float64", "Float64", "float64[pyarrow]"]
    param_names = ["dtype"]

    def setup(self, dtype):
        data = np.random.randn(10000, 1000)
        # all-na columns
        data[:, 600:800] = np.nan
        # partial-na columns
        data[800:1000, 4000:5000] = np.nan
        self.df = DataFrame(data, dtype=dtype)

    def time_isna(self, dtype):
        self.df.isna()


class Count:
    params = [0, 1]
    param_names = ["axis"]

    def setup(self, axis):
        self.df = DataFrame(np.random.randn(10000, 1000))
        self.df.iloc[50:1000, 20:50] = np.nan
        self.df.iloc[2000:3000] = np.nan
        self.df.iloc[:, 60:70] = np.nan
        self.df_mixed = self.df.copy()
        self.df_mixed["foo"] = "bar"

    def time_count(self, axis):
        self.df.count(axis=axis)

    def time_count_mixed_dtypes(self, axis):
        self.df_mixed.count(axis=axis)


class Apply:
    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 100))

        self.s = Series(np.arange(1028.0))
        self.df2 = DataFrame(dict.fromkeys(range(1028), self.s))
        self.df3 = DataFrame(np.random.randn(1000, 3), columns=list("ABC"))

    def time_apply_user_func(self):
        self.df2.apply(lambda x: np.corrcoef(x, self.s)[(0, 1)])

    def time_apply_axis_1(self):
        self.df.apply(lambda x: x + 1, axis=1)

    def time_apply_lambda_mean(self):
        self.df.apply(lambda x: x.mean())

    def time_apply_str_mean(self):
        self.df.apply("mean")

    def time_apply_pass_thru(self):
        self.df.apply(lambda x: x)

    def time_apply_ref_by_name(self):
        self.df3.apply(lambda x: x["A"] + x["B"], axis=1)


class Dtypes:
    def setup(self):
        self.df = DataFrame(np.random.randn(1000, 1000))

    def time_frame_dtypes(self):
        self.df.dtypes


class Equals:
    def setup(self):
        N = 10**3
        self.float_df = DataFrame(np.random.randn(N, N))
        self.float_df_nan = self.float_df.copy()
        self.float_df_nan.iloc[-1, -1] = np.nan

        self.object_df = DataFrame("foo", index=range(N), columns=range(N))
        self.object_df_nan = self.object_df.copy()
        self.object_df_nan.iloc[-1, -1] = np.nan

        self.nonunique_cols = self.object_df.copy()
        self.nonunique_cols.columns = ["A"] * len(self.nonunique_cols.columns)
        self.nonunique_cols_nan = self.nonunique_cols.copy()
        self.nonunique_cols_nan.iloc[-1, -1] = np.nan

    def time_frame_float_equal(self):
        self.float_df.equals(self.float_df)

    def time_frame_float_unequal(self):
        self.float_df.equals(self.float_df_nan)

    def time_frame_nonunique_equal(self):
        self.nonunique_cols.equals(self.nonunique_cols)

    def time_frame_nonunique_unequal(self):
        self.nonunique_cols.equals(self.nonunique_cols_nan)

    def time_frame_object_equal(self):
        self.object_df.equals(self.object_df)

    def time_frame_object_unequal(self):
        self.object_df.equals(self.object_df_nan)


class Interpolate:
    def setup(self):
        N = 10000
        # this is the worst case, where every column has NaNs.
        arr = np.random.randn(N, 100)
        arr[::2] = np.nan

        self.df = DataFrame(arr)

        self.df2 = DataFrame(
            {
                "A": np.arange(0, N),
                "B": np.random.randint(0, 100, N),
                "C": np.random.randn(N),
                "D": np.random.randn(N),
            }
        )
        self.df2.loc[1::5, "A"] = np.nan
        self.df2.loc[1::5, "C"] = np.nan

    def time_interpolate(self):
        self.df.interpolate()

    def time_interpolate_some_good(self):
        self.df2.interpolate()


class Shift:
    # frame shift speedup issue-5609
    params = [0, 1]
    param_names = ["axis"]

    def setup(self, axis):
        self.df = DataFrame(np.random.rand(10000, 500))

    def time_shift(self, axis):
        self.df.shift(1, axis=axis)


class Nunique:
    def setup(self):
        self.df = DataFrame(np.random.randn(10000, 1000))

    def time_frame_nunique(self):
        self.df.nunique()


class SeriesNuniqueWithNan:
    def setup(self):
        values = 100 * [np.nan] + list(range(100))
        self.ser = Series(np.tile(values, 10000), dtype=float)

    def time_series_nunique_nan(self):
        self.ser.nunique()


class Duplicated:
    def setup(self):
        n = 1 << 20
        t = date_range("2015-01-01", freq="s", periods=(n // 64))
        xs = np.random.randn(n // 64).round(2)
        self.df = DataFrame(
            {
                "a": np.random.randint(-1 << 8, 1 << 8, n),
                "b": np.random.choice(t, n),
                "c": np.random.choice(xs, n),
            }
        )
        self.df2 = DataFrame(np.random.randn(1000, 100).astype(str)).T

    def time_frame_duplicated(self):
        self.df.duplicated()

    def time_frame_duplicated_wide(self):
        self.df2.duplicated()

    def time_frame_duplicated_subset(self):
        self.df.duplicated(subset=["a"])


class XS:
    params = [0, 1]
    param_names = ["axis"]

    def setup(self, axis):
        self.N = 10**4
        self.df = DataFrame(np.random.randn(self.N, self.N))

    def time_frame_xs(self, axis):
        self.df.xs(self.N / 2, axis=axis)


class SortValues:
    params = [True, False]
    param_names = ["ascending"]

    def setup(self, ascending):
        self.df = DataFrame(np.random.randn(1000000, 2), columns=list("AB"))

    def time_frame_sort_values(self, ascending):
        self.df.sort_values(by="A", ascending=ascending)


class SortMultiKey:
    params = [True, False]
    param_names = ["monotonic"]

    def setup(self, monotonic):
        N = 10000
        K = 10
        df = DataFrame(
            {
                "key1": Index([f"i-{i}" for i in range(N)], dtype=object).values.repeat(
                    K
                ),
                "key2": Index([f"i-{i}" for i in range(N)], dtype=object).values.repeat(
                    K
                ),
                "value": np.random.randn(N * K),
            }
        )
        if monotonic:
            df = df.sort_values(["key1", "key2"])
        self.df_by_columns = df
        self.df_by_index = df.set_index(["key1", "key2"])

    def time_sort_values(self, monotonic):
        self.df_by_columns.sort_values(by=["key1", "key2"])

    def time_sort_index(self, monotonic):
        self.df_by_index.sort_index()


class Quantile:
    params = [0, 1]
    param_names = ["axis"]

    def setup(self, axis):
        self.df = DataFrame(np.random.randn(1000, 3), columns=list("ABC"))

    def time_frame_quantile(self, axis):
        self.df.quantile([0.1, 0.5], axis=axis)


class Rank:
    param_names = ["dtype"]
    params = [
        ["int", "uint", "float", "object"],
    ]

    def setup(self, dtype):
        self.df = DataFrame(
            np.random.randn(10000, 10).astype(dtype), columns=range(10), dtype=dtype
        )

    def time_rank(self, dtype):
        self.df.rank()


class GetDtypeCounts:
    # 2807
    def setup(self):
        self.df = DataFrame(np.random.randn(10, 10000))

    def time_frame_get_dtype_counts(self):
        with warnings.catch_warnings(record=True):
            self.df.dtypes.value_counts()

    def time_info(self):
        self.df.info()


class NSort:
    params = ["first", "last", "all"]
    param_names = ["keep"]

    def setup(self, keep):
        self.df = DataFrame(np.random.randn(100000, 3), columns=list("ABC"))

    def time_nlargest_one_column(self, keep):
        self.df.nlargest(100, "A", keep=keep)

    def time_nlargest_two_columns(self, keep):
        self.df.nlargest(100, ["A", "B"], keep=keep)

    def time_nsmallest_one_column(self, keep):
        self.df.nsmallest(100, "A", keep=keep)

    def time_nsmallest_two_columns(self, keep):
        self.df.nsmallest(100, ["A", "B"], keep=keep)


class Describe:
    def setup(self):
        self.df = DataFrame(
            {
                "a": np.random.randint(0, 100, 10**6),
                "b": np.random.randint(0, 100, 10**6),
                "c": np.random.randint(0, 100, 10**6),
            }
        )

    def time_series_describe(self):
        self.df["a"].describe()

    def time_dataframe_describe(self):
        self.df.describe()


class MemoryUsage:
    def setup(self):
        self.df = DataFrame(np.random.randn(100000, 2), columns=list("AB"))
        self.df2 = self.df.copy()
        self.df2["A"] = self.df2["A"].astype("object")

    def time_memory_usage(self):
        self.df.memory_usage(deep=True)

    def time_memory_usage_object_dtype(self):
        self.df2.memory_usage(deep=True)


class Round:
    def setup(self):
        self.df = DataFrame(np.random.randn(10000, 10))
        self.df_t = self.df.transpose(copy=True)

    def time_round(self):
        self.df.round()

    def time_round_transposed(self):
        self.df_t.round()

    def peakmem_round(self):
        self.df.round()

    def peakmem_round_transposed(self):
        self.df_t.round()


class Where:
    params = (
        [True, False],
        ["float64", "Float64", "float64[pyarrow]"],
    )
    param_names = ["dtype"]

    def setup(self, inplace, dtype):
        self.df = DataFrame(np.random.randn(100_000, 10), dtype=dtype)
        self.mask = self.df < 0

    def time_where(self, inplace, dtype):
        self.df.where(self.mask, other=0.0, inplace=inplace)


class FindValidIndex:
    param_names = ["dtype"]
    params = [
        ["float", "Float64", "float64[pyarrow]"],
    ]

    def setup(self, dtype):
        df = DataFrame(
            np.random.randn(100000, 2),
            columns=list("AB"),
            dtype=dtype,
        )
        df.iloc[:100, 0] = None
        df.iloc[:200, 1] = None
        df.iloc[-100:, 0] = None
        df.iloc[-200:, 1] = None
        self.df = df

    def time_first_valid_index(self, dtype):
        self.df.first_valid_index()

    def time_last_valid_index(self, dtype):
        self.df.last_valid_index()


class Update:
    def setup(self):
        rng = np.random.default_rng()
        self.df = DataFrame(rng.uniform(size=(1_000_000, 10)))

        idx = rng.choice(range(1_000_000), size=1_000_000, replace=False)
        self.df_random = DataFrame(self.df, index=idx)

        idx = rng.choice(range(1_000_000), size=100_000, replace=False)
        cols = rng.choice(range(10), size=2, replace=False)
        self.df_sample = DataFrame(
            rng.uniform(size=(100_000, 2)), index=idx, columns=cols
        )

    def time_to_update_big_frame_small_arg(self):
        self.df.update(self.df_sample)

    def time_to_update_random_indices(self):
        self.df_random.update(self.df_sample)

    def time_to_update_small_frame_big_arg(self):
        self.df_sample.update(self.df)


from .pandas_vb_common import setup  # noqa: F401 isort:skip
