from itertools import product
import string

import numpy as np

import pandas as pd
from pandas import DataFrame, MultiIndex, date_range, melt, wide_to_long


class Melt:
    def setup(self):
        self.df = DataFrame(np.random.randn(10000, 3), columns=["A", "B", "C"])
        self.df["id1"] = np.random.randint(0, 10, 10000)
        self.df["id2"] = np.random.randint(100, 1000, 10000)

    def time_melt_dataframe(self):
        melt(self.df, id_vars=["id1", "id2"])


class Pivot:
    def setup(self):
        N = 10000
        index = date_range("1/1/2000", periods=N, freq="h")
        data = {
            "value": np.random.randn(N * 50),
            "variable": np.arange(50).repeat(N),
            "date": np.tile(index.values, 50),
        }
        self.df = DataFrame(data)

    def time_reshape_pivot_time_series(self):
        self.df.pivot("date", "variable", "value")


class SimpleReshape:
    def setup(self):
        arrays = [np.arange(100).repeat(100), np.roll(np.tile(np.arange(100), 100), 25)]
        index = MultiIndex.from_arrays(arrays)
        self.df = DataFrame(np.random.randn(10000, 4), index=index)
        self.udf = self.df.unstack(1)

    def time_stack(self):
        self.udf.stack()

    def time_unstack(self):
        self.df.unstack(1)


class Unstack:

    params = ["int", "category"]

    def setup(self, dtype):
        m = 100
        n = 1000

        levels = np.arange(m)
        index = MultiIndex.from_product([levels] * 2)
        columns = np.arange(n)
        if dtype == "int":
            values = np.arange(m * m * n).reshape(m * m, n)
        else:
            # the category branch is ~20x slower than int. So we
            # cut down the size a bit. Now it's only ~3x slower.
            n = 50
            columns = columns[:n]
            indices = np.random.randint(0, 52, size=(m * m, n))
            values = np.take(list(string.ascii_letters), indices)
            values = [pd.Categorical(v) for v in values.T]

        self.df = DataFrame(values, index, columns)
        self.df2 = self.df.iloc[:-1]

    def time_full_product(self, dtype):
        self.df.unstack()

    def time_without_last_row(self, dtype):
        self.df2.unstack()


class SparseIndex:
    def setup(self):
        NUM_ROWS = 1000
        self.df = DataFrame(
            {
                "A": np.random.randint(50, size=NUM_ROWS),
                "B": np.random.randint(50, size=NUM_ROWS),
                "C": np.random.randint(-10, 10, size=NUM_ROWS),
                "D": np.random.randint(-10, 10, size=NUM_ROWS),
                "E": np.random.randint(10, size=NUM_ROWS),
                "F": np.random.randn(NUM_ROWS),
            }
        )
        self.df = self.df.set_index(["A", "B", "C", "D", "E"])

    def time_unstack(self):
        self.df.unstack()


class WideToLong:
    def setup(self):
        nyrs = 20
        nidvars = 20
        N = 5000
        self.letters = list("ABCD")
        yrvars = [l + str(num) for l, num in product(self.letters, range(1, nyrs + 1))]
        columns = [str(i) for i in range(nidvars)] + yrvars
        self.df = DataFrame(np.random.randn(N, nidvars + len(yrvars)), columns=columns)
        self.df["id"] = self.df.index

    def time_wide_to_long_big(self):
        wide_to_long(self.df, self.letters, i="id", j="year")


class PivotTable:
    def setup(self):
        N = 100000
        fac1 = np.array(["A", "B", "C"], dtype="O")
        fac2 = np.array(["one", "two"], dtype="O")
        ind1 = np.random.randint(0, 3, size=N)
        ind2 = np.random.randint(0, 2, size=N)
        self.df = DataFrame(
            {
                "key1": fac1.take(ind1),
                "key2": fac2.take(ind2),
                "key3": fac2.take(ind2),
                "value1": np.random.randn(N),
                "value2": np.random.randn(N),
                "value3": np.random.randn(N),
            }
        )
        self.df2 = DataFrame(
            {"col1": list("abcde"), "col2": list("fghij"), "col3": [1, 2, 3, 4, 5]}
        )
        self.df2.col1 = self.df2.col1.astype("category")
        self.df2.col2 = self.df2.col2.astype("category")

    def time_pivot_table(self):
        self.df.pivot_table(index="key1", columns=["key2", "key3"])

    def time_pivot_table_agg(self):
        self.df.pivot_table(
            index="key1", columns=["key2", "key3"], aggfunc=["sum", "mean"]
        )

    def time_pivot_table_margins(self):
        self.df.pivot_table(index="key1", columns=["key2", "key3"], margins=True)

    def time_pivot_table_categorical(self):
        self.df2.pivot_table(
            index="col1", values="col3", columns="col2", aggfunc=np.sum, fill_value=0
        )

    def time_pivot_table_categorical_observed(self):
        self.df2.pivot_table(
            index="col1",
            values="col3",
            columns="col2",
            aggfunc=np.sum,
            fill_value=0,
            observed=True,
        )

    def time_pivot_table_margins_only_column(self):
        self.df.pivot_table(columns=["key2", "key3"], margins=True)


class Crosstab:
    def setup(self):
        N = 100000
        fac1 = np.array(["A", "B", "C"], dtype="O")
        fac2 = np.array(["one", "two"], dtype="O")
        self.ind1 = np.random.randint(0, 3, size=N)
        self.ind2 = np.random.randint(0, 2, size=N)
        self.vec1 = fac1.take(self.ind1)
        self.vec2 = fac2.take(self.ind2)

    def time_crosstab(self):
        pd.crosstab(self.vec1, self.vec2)

    def time_crosstab_values(self):
        pd.crosstab(self.vec1, self.vec2, values=self.ind1, aggfunc="sum")

    def time_crosstab_normalize(self):
        pd.crosstab(self.vec1, self.vec2, normalize=True)

    def time_crosstab_normalize_margins(self):
        pd.crosstab(self.vec1, self.vec2, normalize=True, margins=True)


class GetDummies:
    def setup(self):
        categories = list(string.ascii_letters[:12])
        s = pd.Series(
            np.random.choice(categories, size=1000000),
            dtype=pd.api.types.CategoricalDtype(categories),
        )
        self.s = s

    def time_get_dummies_1d(self):
        pd.get_dummies(self.s, sparse=False)

    def time_get_dummies_1d_sparse(self):
        pd.get_dummies(self.s, sparse=True)


class Cut:
    params = [[4, 10, 1000]]
    param_names = ["bins"]

    def setup(self, bins):
        N = 10 ** 5
        self.int_series = pd.Series(np.arange(N).repeat(5))
        self.float_series = pd.Series(np.random.randn(N).repeat(5))
        self.timedelta_series = pd.Series(
            np.random.randint(N, size=N), dtype="timedelta64[ns]"
        )
        self.datetime_series = pd.Series(
            np.random.randint(N, size=N), dtype="datetime64[ns]"
        )
        self.interval_bins = pd.IntervalIndex.from_breaks(np.linspace(0, N, bins))

    def time_cut_int(self, bins):
        pd.cut(self.int_series, bins)

    def time_cut_float(self, bins):
        pd.cut(self.float_series, bins)

    def time_cut_timedelta(self, bins):
        pd.cut(self.timedelta_series, bins)

    def time_cut_datetime(self, bins):
        pd.cut(self.datetime_series, bins)

    def time_qcut_int(self, bins):
        pd.qcut(self.int_series, bins)

    def time_qcut_float(self, bins):
        pd.qcut(self.float_series, bins)

    def time_qcut_timedelta(self, bins):
        pd.qcut(self.timedelta_series, bins)

    def time_qcut_datetime(self, bins):
        pd.qcut(self.datetime_series, bins)

    def time_cut_interval(self, bins):
        # GH 27668
        pd.cut(self.int_series, self.interval_bins)

    def peakmem_cut_interval(self, bins):
        # GH 27668
        pd.cut(self.int_series, self.interval_bins)


class Explode:
    param_names = ["n_rows", "max_list_length"]
    params = [[100, 1000, 10000], [3, 5, 10]]

    def setup(self, n_rows, max_list_length):

        data = [np.arange(np.random.randint(max_list_length)) for _ in range(n_rows)]
        self.series = pd.Series(data)

    def time_explode(self, n_rows, max_list_length):
        self.series.explode()


from .pandas_vb_common import setup  # noqa: F401 isort:skip
