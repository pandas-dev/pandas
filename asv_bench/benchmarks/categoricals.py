import string
import sys
import warnings

import numpy as np

import pandas as pd

from .pandas_vb_common import tm

try:
    from pandas.api.types import union_categoricals
except ImportError:
    try:
        from pandas.types.concat import union_categoricals
    except ImportError:
        pass


class Constructor:
    def setup(self):
        N = 10 ** 5
        self.categories = list("abcde")
        self.cat_idx = pd.Index(self.categories)
        self.values = np.tile(self.categories, N)
        self.codes = np.tile(range(len(self.categories)), N)

        self.datetimes = pd.Series(
            pd.date_range("1995-01-01 00:00:00", periods=N / 10, freq="s")
        )
        self.datetimes_with_nat = self.datetimes.copy()
        self.datetimes_with_nat.iloc[-1] = pd.NaT

        self.values_some_nan = list(np.tile(self.categories + [np.nan], N))
        self.values_all_nan = [np.nan] * len(self.values)
        self.values_all_int8 = np.ones(N, "int8")
        self.categorical = pd.Categorical(self.values, self.categories)
        self.series = pd.Series(self.categorical)
        self.intervals = pd.interval_range(0, 1, periods=N // 10)

    def time_regular(self):
        pd.Categorical(self.values, self.categories)

    def time_fastpath(self):
        pd.Categorical(self.codes, self.cat_idx, fastpath=True)

    def time_datetimes(self):
        pd.Categorical(self.datetimes)

    def time_interval(self):
        pd.Categorical(self.datetimes, categories=self.datetimes)

    def time_datetimes_with_nat(self):
        pd.Categorical(self.datetimes_with_nat)

    def time_with_nan(self):
        pd.Categorical(self.values_some_nan)

    def time_all_nan(self):
        pd.Categorical(self.values_all_nan)

    def time_from_codes_all_int8(self):
        pd.Categorical.from_codes(self.values_all_int8, self.categories)

    def time_existing_categorical(self):
        pd.Categorical(self.categorical)

    def time_existing_series(self):
        pd.Categorical(self.series)


class AsType:
    def setup(self):
        N = 10 ** 5

        random_pick = np.random.default_rng().choice

        categories = {
            "str": list(string.ascii_letters),
            "int": np.random.randint(2 ** 16, size=154),
            "float": sys.maxsize * np.random.random((38,)),
            "timestamp": [
                pd.Timestamp(x, unit="s") for x in np.random.randint(2 ** 18, size=578)
            ],
        }

        self.df = pd.DataFrame(
            {col: random_pick(cats, N) for col, cats in categories.items()}
        )

        for col in ("int", "float", "timestamp"):
            self.df[col + "_as_str"] = self.df[col].astype(str)

        for col in self.df.columns:
            self.df[col] = self.df[col].astype("category")

    def astype_str(self):
        [self.df[col].astype("str") for col in "int float timestamp".split()]

    def astype_int(self):
        [self.df[col].astype("int") for col in "int_as_str timestamp".split()]

    def astype_float(self):
        [
            self.df[col].astype("float")
            for col in "float_as_str int int_as_str timestamp".split()
        ]

    def astype_datetime(self):
        self.df["float"].astype(pd.DatetimeTZDtype(tz="US/Pacific"))


class Concat:
    def setup(self):
        N = 10 ** 5
        self.s = pd.Series(list("aabbcd") * N).astype("category")

        self.a = pd.Categorical(list("aabbcd") * N)
        self.b = pd.Categorical(list("bbcdjk") * N)

        self.idx_a = pd.CategoricalIndex(range(N), range(N))
        self.idx_b = pd.CategoricalIndex(range(N + 1), range(N + 1))
        self.df_a = pd.DataFrame(range(N), columns=["a"], index=self.idx_a)
        self.df_b = pd.DataFrame(range(N + 1), columns=["a"], index=self.idx_b)

    def time_concat(self):
        pd.concat([self.s, self.s])

    def time_union(self):
        union_categoricals([self.a, self.b])

    def time_append_overlapping_index(self):
        self.idx_a.append(self.idx_a)

    def time_append_non_overlapping_index(self):
        self.idx_a.append(self.idx_b)

    def time_concat_overlapping_index(self):
        pd.concat([self.df_a, self.df_a])

    def time_concat_non_overlapping_index(self):
        pd.concat([self.df_a, self.df_b])


class ValueCounts:

    params = [True, False]
    param_names = ["dropna"]

    def setup(self, dropna):
        n = 5 * 10 ** 5
        arr = [f"s{i:04d}" for i in np.random.randint(0, n // 10, size=n)]
        self.ts = pd.Series(arr).astype("category")

    def time_value_counts(self, dropna):
        self.ts.value_counts(dropna=dropna)


class Repr:
    def setup(self):
        self.sel = pd.Series(["s1234"]).astype("category")

    def time_rendering(self):
        str(self.sel)


class SetCategories:
    def setup(self):
        n = 5 * 10 ** 5
        arr = [f"s{i:04d}" for i in np.random.randint(0, n // 10, size=n)]
        self.ts = pd.Series(arr).astype("category")

    def time_set_categories(self):
        self.ts.cat.set_categories(self.ts.cat.categories[::2])


class RemoveCategories:
    def setup(self):
        n = 5 * 10 ** 5
        arr = [f"s{i:04d}" for i in np.random.randint(0, n // 10, size=n)]
        self.ts = pd.Series(arr).astype("category")

    def time_remove_categories(self):
        self.ts.cat.remove_categories(self.ts.cat.categories[::2])


class Rank:
    def setup(self):
        N = 10 ** 5
        ncats = 100

        self.s_str = pd.Series(tm.makeCategoricalIndex(N, ncats)).astype(str)
        self.s_str_cat = pd.Series(self.s_str, dtype="category")
        with warnings.catch_warnings(record=True):
            str_cat_type = pd.CategoricalDtype(set(self.s_str), ordered=True)
            self.s_str_cat_ordered = self.s_str.astype(str_cat_type)

        self.s_int = pd.Series(np.random.randint(0, ncats, size=N))
        self.s_int_cat = pd.Series(self.s_int, dtype="category")
        with warnings.catch_warnings(record=True):
            int_cat_type = pd.CategoricalDtype(set(self.s_int), ordered=True)
            self.s_int_cat_ordered = self.s_int.astype(int_cat_type)

    def time_rank_string(self):
        self.s_str.rank()

    def time_rank_string_cat(self):
        self.s_str_cat.rank()

    def time_rank_string_cat_ordered(self):
        self.s_str_cat_ordered.rank()

    def time_rank_int(self):
        self.s_int.rank()

    def time_rank_int_cat(self):
        self.s_int_cat.rank()

    def time_rank_int_cat_ordered(self):
        self.s_int_cat_ordered.rank()


class IsMonotonic:
    def setup(self):
        N = 1000
        self.c = pd.CategoricalIndex(list("a" * N + "b" * N + "c" * N))
        self.s = pd.Series(self.c)

    def time_categorical_index_is_monotonic_increasing(self):
        self.c.is_monotonic_increasing

    def time_categorical_index_is_monotonic_decreasing(self):
        self.c.is_monotonic_decreasing

    def time_categorical_series_is_monotonic_increasing(self):
        self.s.is_monotonic_increasing

    def time_categorical_series_is_monotonic_decreasing(self):
        self.s.is_monotonic_decreasing


class Contains:
    def setup(self):
        N = 10 ** 5
        self.ci = tm.makeCategoricalIndex(N)
        self.c = self.ci.values
        self.key = self.ci.categories[0]

    def time_categorical_index_contains(self):
        self.key in self.ci

    def time_categorical_contains(self):
        self.key in self.c


class CategoricalSlicing:

    params = ["monotonic_incr", "monotonic_decr", "non_monotonic"]
    param_names = ["index"]

    def setup(self, index):
        N = 10 ** 6
        categories = ["a", "b", "c"]
        values = [0] * N + [1] * N + [2] * N
        if index == "monotonic_incr":
            self.data = pd.Categorical.from_codes(values, categories=categories)
        elif index == "monotonic_decr":
            self.data = pd.Categorical.from_codes(
                list(reversed(values)), categories=categories
            )
        elif index == "non_monotonic":
            self.data = pd.Categorical.from_codes([0, 1, 2] * N, categories=categories)
        else:
            raise ValueError(f"Invalid index param: {index}")

        self.scalar = 10000
        self.list = list(range(10000))
        self.cat_scalar = "b"

    def time_getitem_scalar(self, index):
        self.data[self.scalar]

    def time_getitem_slice(self, index):
        self.data[: self.scalar]

    def time_getitem_list_like(self, index):
        self.data[[self.scalar]]

    def time_getitem_list(self, index):
        self.data[self.list]

    def time_getitem_bool_array(self, index):
        self.data[self.data == self.cat_scalar]


class Indexing:
    def setup(self):
        N = 10 ** 5
        self.index = pd.CategoricalIndex(range(N), range(N))
        self.series = pd.Series(range(N), index=self.index).sort_index()
        self.category = self.index[500]

    def time_get_loc(self):
        self.index.get_loc(self.category)

    def time_shallow_copy(self):
        self.index._view()

    def time_align(self):
        pd.DataFrame({"a": self.series, "b": self.series[:500]})

    def time_intersection(self):
        self.index[:750].intersection(self.index[250:])

    def time_unique(self):
        self.index.unique()

    def time_reindex(self):
        self.index.reindex(self.index[:500])

    def time_reindex_missing(self):
        self.index.reindex(["a", "b", "c", "d"])

    def time_sort_values(self):
        self.index.sort_values(ascending=False)


class SearchSorted:
    def setup(self):
        N = 10 ** 5
        self.ci = tm.makeCategoricalIndex(N).sort_values()
        self.c = self.ci.values
        self.key = self.ci.categories[1]

    def time_categorical_index_contains(self):
        self.ci.searchsorted(self.key)

    def time_categorical_contains(self):
        self.c.searchsorted(self.key)


from .pandas_vb_common import setup  # noqa: F401 isort:skip
