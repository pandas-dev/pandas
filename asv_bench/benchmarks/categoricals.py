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

    def time_regular(self):
        pd.Categorical(self.values, self.categories)

    def time_fastpath(self):
        pd.Categorical(self.codes, self.cat_idx, fastpath=True)

    def time_datetimes(self):
        pd.Categorical(self.datetimes)

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


class Concat:
    def setup(self):
        N = 10 ** 5
        self.s = pd.Series(list("aabbcd") * N).astype("category")

        self.a = pd.Categorical(list("aabbcd") * N)
        self.b = pd.Categorical(list("bbcdjk") * N)

    def time_concat(self):
        pd.concat([self.s, self.s])

    def time_union(self):
        union_categoricals([self.a, self.b])


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


class Isin:

    params = ["object", "int64"]
    param_names = ["dtype"]

    def setup(self, dtype):
        np.random.seed(1234)
        n = 5 * 10 ** 5
        sample_size = 100
        arr = list(np.random.randint(0, n // 10, size=n))
        if dtype == "object":
            arr = [f"s{i:04d}" for i in arr]
        self.sample = np.random.choice(arr, sample_size)
        self.series = pd.Series(arr).astype("category")

    def time_isin_categorical(self, dtype):
        self.series.isin(self.sample)


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
        self.index._shallow_copy()

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
