"""
These benchmarks are for Series and DataFrame indexing methods.  For the
lower-level methods directly on Index and subclasses, see index_object.py,
indexing_engine.py, and index_cached.py
"""

from datetime import datetime
import warnings

import numpy as np

from pandas import (
    NA,
    CategoricalIndex,
    DataFrame,
    Index,
    IntervalIndex,
    MultiIndex,
    Series,
    concat,
    date_range,
    option_context,
    period_range,
)


class NumericSeriesIndexing:
    params = [
        (np.int64, np.uint64, np.float64),
        ("unique_monotonic_inc", "nonunique_monotonic_inc"),
    ]
    param_names = ["dtype", "index_structure"]

    def setup(self, dtype, index_structure):
        N = 10**6
        indices = {
            "unique_monotonic_inc": Index(range(N), dtype=dtype),
            "nonunique_monotonic_inc": Index(
                list(range(55)) + [54] + list(range(55, N - 1)), dtype=dtype
            ),
        }
        self.data = Series(np.random.rand(N), index=indices[index_structure])
        self.array = np.arange(10000)
        self.array_list = self.array.tolist()

    def time_getitem_scalar(self, index, index_structure):
        self.data[800000]

    def time_getitem_slice(self, index, index_structure):
        self.data[:800000]

    def time_getitem_list_like(self, index, index_structure):
        self.data[[800000]]

    def time_getitem_array(self, index, index_structure):
        self.data[self.array]

    def time_getitem_lists(self, index, index_structure):
        self.data[self.array_list]

    def time_iloc_array(self, index, index_structure):
        self.data.iloc[self.array]

    def time_iloc_list_like(self, index, index_structure):
        self.data.iloc[[800000]]

    def time_iloc_scalar(self, index, index_structure):
        self.data.iloc[800000]

    def time_iloc_slice(self, index, index_structure):
        self.data.iloc[:800000]

    def time_loc_array(self, index, index_structure):
        self.data.loc[self.array]

    def time_loc_list_like(self, index, index_structure):
        self.data.loc[[800000]]

    def time_loc_scalar(self, index, index_structure):
        self.data.loc[800000]

    def time_loc_slice(self, index, index_structure):
        self.data.loc[:800000]


class NumericMaskedIndexing:
    monotonic_list = list(range(10**6))
    non_monotonic_list = list(range(50)) + [54, 53, 52, 51] + list(range(55, 10**6 - 1))

    params = [
        ("Int64", "UInt64", "Float64"),
        (True, False),
    ]
    param_names = ["dtype", "monotonic"]

    def setup(self, dtype, monotonic):
        indices = {
            True: Index(self.monotonic_list, dtype=dtype),
            False: Index(self.non_monotonic_list, dtype=dtype).append(
                Index([NA], dtype=dtype)
            ),
        }
        self.data = indices[monotonic]
        self.indexer = np.arange(300, 1_000)
        self.data_dups = self.data.append(self.data)

    def time_get_indexer(self, dtype, monotonic):
        self.data.get_indexer(self.indexer)

    def time_get_indexer_dups(self, dtype, monotonic):
        self.data.get_indexer_for(self.indexer)


class NonNumericSeriesIndexing:
    params = [
        ("string", "datetime", "period"),
        ("unique_monotonic_inc", "nonunique_monotonic_inc", "non_monotonic"),
    ]
    param_names = ["index_dtype", "index_structure"]

    def setup(self, index, index_structure):
        N = 10**6
        if index == "string":
            index = Index([f"i-{i}" for i in range(N)], dtype=object)
        elif index == "datetime":
            index = date_range("1900", periods=N, freq="s")
        elif index == "period":
            index = period_range("1900", periods=N, freq="s")
        index = index.sort_values()
        assert index.is_unique and index.is_monotonic_increasing
        if index_structure == "nonunique_monotonic_inc":
            index = index.insert(item=index[2], loc=2)[:-1]
        elif index_structure == "non_monotonic":
            index = index[::2].append(index[1::2])
            assert len(index) == N
        self.s = Series(np.random.rand(N), index=index)
        self.lbl = index[80000]
        # warm up index mapping
        self.s[self.lbl]

    def time_getitem_label_slice(self, index, index_structure):
        self.s[: self.lbl]

    def time_getitem_pos_slice(self, index, index_structure):
        self.s[:80000]

    def time_getitem_scalar(self, index, index_structure):
        self.s[self.lbl]

    def time_getitem_list_like(self, index, index_structure):
        self.s[[self.lbl]]


class DataFrameStringIndexing:
    def setup(self):
        index = Index([f"i-{i}" for i in range(1000)], dtype=object)
        columns = Index([f"i-{i}" for i in range(30)], dtype=object)
        with warnings.catch_warnings(record=True):
            self.df = DataFrame(np.random.randn(1000, 30), index=index, columns=columns)
        self.idx_scalar = index[100]
        self.col_scalar = columns[10]
        self.bool_indexer = self.df[self.col_scalar] > 0
        self.bool_obj_indexer = self.bool_indexer.astype(object)
        self.boolean_indexer = (self.df[self.col_scalar] > 0).astype("boolean")

    def time_loc(self):
        self.df.loc[self.idx_scalar, self.col_scalar]

    def time_at(self):
        self.df.at[self.idx_scalar, self.col_scalar]

    def time_at_setitem(self):
        self.df.at[self.idx_scalar, self.col_scalar] = 0.0

    def time_getitem_scalar(self):
        self.df[self.col_scalar][self.idx_scalar]

    def time_boolean_rows(self):
        self.df[self.bool_indexer]

    def time_boolean_rows_object(self):
        self.df[self.bool_obj_indexer]

    def time_boolean_rows_boolean(self):
        self.df[self.boolean_indexer]


class DataFrameNumericIndexing:
    params = [
        (np.int64, np.uint64, np.float64),
        ("unique_monotonic_inc", "nonunique_monotonic_inc"),
    ]
    param_names = ["dtype", "index_structure"]

    def setup(self, dtype, index_structure):
        N = 10**5
        indices = {
            "unique_monotonic_inc": Index(range(N), dtype=dtype),
            "nonunique_monotonic_inc": Index(
                list(range(55)) + [54] + list(range(55, N - 1)), dtype=dtype
            ),
        }
        self.idx_dupe = np.array(range(30)) * 99
        self.df = DataFrame(np.random.randn(N, 5), index=indices[index_structure])
        self.df_dup = concat([self.df, 2 * self.df, 3 * self.df])
        self.bool_indexer = [True] * (N // 2) + [False] * (N - N // 2)

    def time_iloc_dups(self, index, index_structure):
        self.df_dup.iloc[self.idx_dupe]

    def time_loc_dups(self, index, index_structure):
        self.df_dup.loc[self.idx_dupe]

    def time_iloc(self, index, index_structure):
        self.df.iloc[:100, 0]

    def time_loc(self, index, index_structure):
        self.df.loc[:100, 0]

    def time_bool_indexer(self, index, index_structure):
        self.df[self.bool_indexer]


class Take:
    params = ["int", "datetime"]
    param_names = ["index"]

    def setup(self, index):
        N = 100000
        indexes = {
            "int": Index(np.arange(N), dtype=np.int64),
            "datetime": date_range("2011-01-01", freq="s", periods=N),
        }
        index = indexes[index]
        self.s = Series(np.random.rand(N), index=index)
        self.indexer = np.random.randint(0, N, size=N)

    def time_take(self, index):
        self.s.take(self.indexer)


class MultiIndexing:
    params = [True, False]
    param_names = ["unique_levels"]

    def setup(self, unique_levels):
        self.nlevels = 2
        if unique_levels:
            mi = MultiIndex.from_arrays([range(1000000)] * self.nlevels)
        else:
            mi = MultiIndex.from_product([range(1000)] * self.nlevels)
        self.df = DataFrame(np.random.randn(len(mi)), index=mi)

        self.tgt_slice = slice(200, 800)
        self.tgt_null_slice = slice(None)
        self.tgt_list = list(range(0, 1000, 10))
        self.tgt_scalar = 500

        bool_indexer = np.zeros(len(mi), dtype=np.bool_)
        bool_indexer[slice(0, len(mi), 100)] = True
        self.tgt_bool_indexer = bool_indexer

    def time_loc_partial_key_slice(self, unique_levels):
        self.df.loc[self.tgt_slice, :]

    def time_loc_partial_key_null_slice(self, unique_levels):
        self.df.loc[self.tgt_null_slice, :]

    def time_loc_partial_key_list(self, unique_levels):
        self.df.loc[self.tgt_list, :]

    def time_loc_partial_key_scalar(self, unique_levels):
        self.df.loc[self.tgt_scalar, :]

    def time_loc_partial_key_bool_indexer(self, unique_levels):
        self.df.loc[self.tgt_bool_indexer, :]

    def time_loc_all_slices(self, unique_levels):
        target = tuple([self.tgt_slice] * self.nlevels)
        self.df.loc[target, :]

    def time_loc_all_null_slices(self, unique_levels):
        target = tuple([self.tgt_null_slice] * self.nlevels)
        self.df.loc[target, :]

    def time_loc_all_lists(self, unique_levels):
        target = tuple([self.tgt_list] * self.nlevels)
        self.df.loc[target, :]

    def time_loc_all_scalars(self, unique_levels):
        target = tuple([self.tgt_scalar] * self.nlevels)
        self.df.loc[target, :]

    def time_loc_all_bool_indexers(self, unique_levels):
        target = tuple([self.tgt_bool_indexer] * self.nlevels)
        self.df.loc[target, :]

    def time_loc_slice_plus_null_slice(self, unique_levels):
        target = (self.tgt_slice, self.tgt_null_slice)
        self.df.loc[target, :]

    def time_loc_null_slice_plus_slice(self, unique_levels):
        target = (self.tgt_null_slice, self.tgt_slice)
        self.df.loc[target, :]

    def time_loc_multiindex(self, unique_levels):
        target = self.df.index[::10]
        self.df.loc[target]

    def time_xs_level_0(self, unique_levels):
        target = self.tgt_scalar
        self.df.xs(target, level=0)

    def time_xs_level_1(self, unique_levels):
        target = self.tgt_scalar
        self.df.xs(target, level=1)

    def time_xs_full_key(self, unique_levels):
        target = tuple([self.tgt_scalar] * self.nlevels)
        self.df.xs(target)


class IntervalIndexing:
    def setup_cache(self):
        idx = IntervalIndex.from_breaks(np.arange(1000001))
        monotonic = Series(np.arange(1000000), index=idx)
        return monotonic

    def time_getitem_scalar(self, monotonic):
        monotonic[80000]

    def time_loc_scalar(self, monotonic):
        monotonic.loc[80000]

    def time_getitem_list(self, monotonic):
        monotonic[80000:]

    def time_loc_list(self, monotonic):
        monotonic.loc[80000:]


class DatetimeIndexIndexing:
    def setup(self):
        dti = date_range("2016-01-01", periods=10000, tz="US/Pacific")
        dti2 = dti.tz_convert("UTC")
        self.dti = dti
        self.dti2 = dti2

    def time_get_indexer_mismatched_tz(self):
        # reached via e.g.
        #  ser = Series(range(len(dti)), index=dti)
        #  ser[dti2]
        self.dti.get_indexer(self.dti2)


class SortedAndUnsortedDatetimeIndexLoc:
    def setup(self):
        dti = date_range("2016-01-01", periods=10000, tz="US/Pacific")
        index = np.array(dti)

        unsorted_index = index.copy()
        unsorted_index[10] = unsorted_index[20]

        self.df_unsorted = DataFrame(index=unsorted_index, data={"a": 1})
        self.df_sort = DataFrame(index=index, data={"a": 1})

    def time_loc_unsorted(self):
        self.df_unsorted.loc["2016-6-11"]

    def time_loc_sorted(self):
        self.df_sort.loc["2016-6-11"]


class CategoricalIndexIndexing:
    params = ["monotonic_incr", "monotonic_decr", "non_monotonic"]
    param_names = ["index"]

    def setup(self, index):
        N = 10**5
        values = list("a" * N + "b" * N + "c" * N)
        indices = {
            "monotonic_incr": CategoricalIndex(values),
            "monotonic_decr": CategoricalIndex(reversed(values)),
            "non_monotonic": CategoricalIndex(list("abc" * N)),
        }
        self.data = indices[index]
        self.data_unique = CategoricalIndex([str(i) for i in range(N * 3)])

        self.int_scalar = 10000
        self.int_list = list(range(10000))

        self.cat_scalar = "b"
        self.cat_list = ["1", "3"]

    def time_getitem_scalar(self, index):
        self.data[self.int_scalar]

    def time_getitem_slice(self, index):
        self.data[: self.int_scalar]

    def time_getitem_list_like(self, index):
        self.data[[self.int_scalar]]

    def time_getitem_list(self, index):
        self.data[self.int_list]

    def time_getitem_bool_array(self, index):
        self.data[self.data == self.cat_scalar]

    def time_get_loc_scalar(self, index):
        self.data.get_loc(self.cat_scalar)

    def time_get_indexer_list(self, index):
        self.data_unique.get_indexer(self.cat_list)


class MethodLookup:
    def setup_cache(self):
        s = Series()
        return s

    def time_lookup_iloc(self, s):
        s.iloc

    def time_lookup_loc(self, s):
        s.loc


class GetItemSingleColumn:
    def setup(self):
        self.df_string_col = DataFrame(np.random.randn(3000, 1), columns=["A"])
        self.df_int_col = DataFrame(np.random.randn(3000, 1))

    def time_frame_getitem_single_column_label(self):
        self.df_string_col["A"]

    def time_frame_getitem_single_column_int(self):
        self.df_int_col[0]


class IndexSingleRow:
    params = [True, False]
    param_names = ["unique_cols"]

    def setup(self, unique_cols):
        arr = np.arange(10**7).reshape(-1, 10)
        df = DataFrame(arr)
        dtypes = ["u1", "u2", "u4", "u8", "i1", "i2", "i4", "i8", "f8", "f4"]
        for i, d in enumerate(dtypes):
            df[i] = df[i].astype(d)

        if not unique_cols:
            # GH#33032 single-row lookups with non-unique columns were
            #  15x slower than with unique columns
            df.columns = ["A", "A"] + list(df.columns[2:])

        self.df = df

    def time_iloc_row(self, unique_cols):
        self.df.iloc[10000]

    def time_loc_row(self, unique_cols):
        self.df.loc[10000]


class AssignTimeseriesIndex:
    def setup(self):
        N = 100000
        idx = date_range("1/1/2000", periods=N, freq="h")
        self.df = DataFrame(np.random.randn(N, 1), columns=["A"], index=idx)

    def time_frame_assign_timeseries_index(self):
        self.df["date"] = self.df.index


class InsertColumns:
    def setup(self):
        self.N = 10**3
        self.df = DataFrame(index=range(self.N))
        self.df2 = DataFrame(np.random.randn(self.N, 2))

    def time_insert(self):
        for i in range(100):
            self.df.insert(0, i, np.random.randn(self.N), allow_duplicates=True)

    def time_insert_middle(self):
        # same as time_insert but inserting to a middle column rather than
        #  front or back (which have fast-paths)
        for i in range(100):
            self.df2.insert(
                1, "colname", np.random.randn(self.N), allow_duplicates=True
            )

    def time_assign_with_setitem(self):
        for i in range(100):
            self.df[i] = np.random.randn(self.N)

    def time_assign_list_like_with_setitem(self):
        self.df[list(range(100))] = np.random.randn(self.N, 100)

    def time_assign_list_of_columns_concat(self):
        df = DataFrame(np.random.randn(self.N, 100))
        concat([self.df, df], axis=1)


class Setitem:
    def setup(self):
        N = 500_000
        cols = 500
        self.df = DataFrame(np.random.rand(N, cols))

    def time_setitem(self):
        self.df[100] = 100

    def time_setitem_list(self):
        self.df[[100, 200, 300]] = 100


class SetitemObjectDtype:
    # GH#19299

    def setup(self):
        N = 1000
        cols = 500
        self.df = DataFrame(index=range(N), columns=range(cols), dtype=object)

    def time_setitem_object_dtype(self):
        self.df.loc[0, 1] = 1.0


class ChainIndexing:
    params = [None, "warn"]
    param_names = ["mode"]

    def setup(self, mode):
        self.N = 1000000
        self.df = DataFrame({"A": np.arange(self.N), "B": "foo"})

    def time_chained_indexing(self, mode):
        df = self.df
        N = self.N
        with warnings.catch_warnings(record=True):
            with option_context("mode.chained_assignment", mode):
                df2 = df[df.A > N // 2]
                df2["C"] = 1.0


class Block:
    params = [
        (True, "True"),
        (np.array(True), "np.array(True)"),
    ]

    def setup(self, true_value, mode):
        self.df = DataFrame(
            False,
            columns=np.arange(500).astype(str),
            index=date_range("2010-01-01", "2011-01-01"),
        )

        self.true_value = true_value

    def time_test(self, true_value, mode):
        start = datetime(2010, 5, 1)
        end = datetime(2010, 9, 1)
        self.df.loc[start:end, :] = true_value


from .pandas_vb_common import setup  # noqa: F401 isort:skip
