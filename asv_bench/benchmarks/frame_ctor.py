import numpy as np

import pandas as pd
from pandas import (
    DataFrame,
    MultiIndex,
    Series,
    Timestamp,
    date_range,
)

from .pandas_vb_common import tm

try:
    from pandas.tseries.offsets import (
        Hour,
        Nano,
    )
except ImportError:
    # For compatibility with older versions
    from pandas.core.datetools import *  # noqa


class FromDicts:
    def setup(self):
        N, K = 5000, 50
        self.index = tm.makeStringIndex(N)
        self.columns = tm.makeStringIndex(K)
        frame = DataFrame(np.random.randn(N, K), index=self.index, columns=self.columns)
        self.data = frame.to_dict()
        self.dict_list = frame.to_dict(orient="records")
        self.data2 = {i: {j: float(j) for j in range(100)} for i in range(2000)}

    def time_list_of_dict(self):
        DataFrame(self.dict_list)

    def time_nested_dict(self):
        DataFrame(self.data)

    def time_nested_dict_index(self):
        DataFrame(self.data, index=self.index)

    def time_nested_dict_columns(self):
        DataFrame(self.data, columns=self.columns)

    def time_nested_dict_index_columns(self):
        DataFrame(self.data, index=self.index, columns=self.columns)

    def time_nested_dict_int64(self):
        # nested dict, integer indexes, regression described in #621
        DataFrame(self.data2)


class FromSeries:
    def setup(self):
        mi = MultiIndex.from_product([range(100), range(100)])
        self.s = Series(np.random.randn(10000), index=mi)

    def time_mi_series(self):
        DataFrame(self.s)


class FromDictwithTimestamp:

    params = [Nano(1), Hour(1)]
    param_names = ["offset"]

    def setup(self, offset):
        N = 10 ** 3
        idx = date_range(Timestamp("1/1/1900"), freq=offset, periods=N)
        df = DataFrame(np.random.randn(N, 10), index=idx)
        self.d = df.to_dict()

    def time_dict_with_timestamp_offsets(self, offset):
        DataFrame(self.d)


class FromRecords:

    params = [None, 1000]
    param_names = ["nrows"]

    # Generators get exhausted on use, so run setup before every call
    number = 1
    repeat = (3, 250, 10)

    def setup(self, nrows):
        N = 100000
        self.gen = ((x, (x * 20), (x * 100)) for x in range(N))

    def time_frame_from_records_generator(self, nrows):
        # issue-6700
        self.df = DataFrame.from_records(self.gen, nrows=nrows)


class FromNDArray:
    def setup(self):
        N = 100000
        self.data = np.random.randn(N)

    def time_frame_from_ndarray(self):
        self.df = DataFrame(self.data)


class FromLists:

    goal_time = 0.2

    def setup(self):
        N = 1000
        M = 100
        self.data = [list(range(M)) for i in range(N)]

    def time_frame_from_lists(self):
        self.df = DataFrame(self.data)


class FromRange:

    goal_time = 0.2

    def setup(self):
        N = 1_000_000
        self.data = range(N)

    def time_frame_from_range(self):
        self.df = DataFrame(self.data)


class FromArrays:

    goal_time = 0.2

    def setup(self):
        N_rows = 1000
        N_cols = 1000
        self.float_arrays = [np.random.randn(N_rows) for _ in range(N_cols)]
        self.sparse_arrays = [
            pd.arrays.SparseArray(np.random.randint(0, 2, N_rows), dtype="float64")
            for _ in range(N_cols)
        ]
        self.int_arrays = [
            pd.array(np.random.randint(1000, size=N_rows), dtype="Int64")
            for _ in range(N_cols)
        ]
        self.index = pd.Index(range(N_rows))
        self.columns = pd.Index(range(N_cols))

    def time_frame_from_arrays_float(self):
        self.df = DataFrame._from_arrays(
            self.float_arrays,
            index=self.index,
            columns=self.columns,
            verify_integrity=False,
        )

    def time_frame_from_arrays_int(self):
        self.df = DataFrame._from_arrays(
            self.int_arrays,
            index=self.index,
            columns=self.columns,
            verify_integrity=False,
        )

    def time_frame_from_arrays_sparse(self):
        self.df = DataFrame._from_arrays(
            self.sparse_arrays,
            index=self.index,
            columns=self.columns,
            verify_integrity=False,
        )


from .pandas_vb_common import setup  # noqa: F401 isort:skip
