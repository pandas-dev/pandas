import numpy as np

from pandas import (
    DataFrame,
    HDFStore,
    Index,
    date_range,
    read_hdf,
)

from ..pandas_vb_common import BaseIO


class HDFStoreDataFrame(BaseIO):
    def setup(self):
        N = 25000
        index = Index([f"i-{i}" for i in range(N)], dtype=object)
        self.df = DataFrame(
            {"float1": np.random.randn(N), "float2": np.random.randn(N)}, index=index
        )
        self.df_mixed = DataFrame(
            {
                "float1": np.random.randn(N),
                "float2": np.random.randn(N),
                "string1": ["foo"] * N,
                "bool1": [True] * N,
                "int1": np.random.randint(0, N, size=N),
            },
            index=index,
        )
        self.df_wide = DataFrame(np.random.randn(N, 100))
        self.start_wide = self.df_wide.index[10000]
        self.stop_wide = self.df_wide.index[15000]
        self.df2 = DataFrame(
            {"float1": np.random.randn(N), "float2": np.random.randn(N)},
            index=date_range("1/1/2000", periods=N),
        )
        self.start = self.df2.index[10000]
        self.stop = self.df2.index[15000]
        self.df_wide2 = DataFrame(
            np.random.randn(N, 100), index=date_range("1/1/2000", periods=N)
        )
        self.df_dc = DataFrame(
            np.random.randn(N, 10), columns=[f"C{i:03d}" for i in range(10)]
        )

        self.fname = "__test__.h5"

        self.store = HDFStore(self.fname)
        self.store.put("fixed", self.df)
        self.store.put("fixed_mixed", self.df_mixed)
        self.store.append("table", self.df2)
        self.store.append("table_mixed", self.df_mixed)
        self.store.append("table_wide", self.df_wide)
        self.store.append("table_wide2", self.df_wide2)

    def teardown(self):
        self.store.close()
        self.remove(self.fname)

    def time_read_store(self):
        self.store.get("fixed")

    def time_read_store_mixed(self):
        self.store.get("fixed_mixed")

    def time_write_store(self):
        self.store.put("fixed_write", self.df)

    def time_write_store_mixed(self):
        self.store.put("fixed_mixed_write", self.df_mixed)

    def time_read_store_table_mixed(self):
        self.store.select("table_mixed")

    def time_write_store_table_mixed(self):
        self.store.append("table_mixed_write", self.df_mixed)

    def time_read_store_table(self):
        self.store.select("table")

    def time_write_store_table(self):
        self.store.append("table_write", self.df)

    def time_read_store_table_wide(self):
        self.store.select("table_wide")

    def time_write_store_table_wide(self):
        self.store.append("table_wide_write", self.df_wide)

    def time_write_store_table_dc(self):
        self.store.append("table_dc_write", self.df_dc, data_columns=True)

    def time_query_store_table_wide(self):
        self.store.select(
            "table_wide", where="index > self.start_wide and index < self.stop_wide"
        )

    def time_query_store_table(self):
        self.store.select("table", where="index > self.start and index < self.stop")

    def time_store_repr(self):
        repr(self.store)

    def time_store_str(self):
        str(self.store)

    def time_store_info(self):
        self.store.info()


class HDF(BaseIO):
    params = ["table", "fixed"]
    param_names = ["format"]

    def setup(self, format):
        self.fname = "__test__.h5"
        N = 100000
        C = 5
        self.df = DataFrame(
            np.random.randn(N, C),
            columns=[f"float{i}" for i in range(C)],
            index=date_range("20000101", periods=N, freq="h"),
        )
        self.df["object"] = Index([f"i-{i}" for i in range(N)], dtype=object)
        self.df.to_hdf(self.fname, key="df", format=format)

        # Numeric df
        self.df1 = self.df.copy()
        self.df1 = self.df1.reset_index()
        self.df1.to_hdf(self.fname, key="df1", format=format)

    def time_read_hdf(self, format):
        read_hdf(self.fname, "df")

    def peakmem_read_hdf(self, format):
        read_hdf(self.fname, "df")

    def time_write_hdf(self, format):
        self.df.to_hdf(self.fname, key="df", format=format)


from ..pandas_vb_common import setup  # noqa: F401 isort:skip
