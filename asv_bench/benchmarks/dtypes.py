import string

import numpy as np

import pandas as pd
from pandas import DataFrame
import pandas._testing as tm
from pandas.api.types import (
    is_extension_array_dtype,
    pandas_dtype,
)

from .pandas_vb_common import (
    datetime_dtypes,
    extension_dtypes,
    numeric_dtypes,
    string_dtypes,
)

_numpy_dtypes = [
    np.dtype(dtype) for dtype in (numeric_dtypes + datetime_dtypes + string_dtypes)
]
_dtypes = _numpy_dtypes + extension_dtypes


class Dtypes:
    params = _dtypes + list(map(lambda dt: dt.name, _dtypes))
    param_names = ["dtype"]

    def time_pandas_dtype(self, dtype):
        pandas_dtype(dtype)


class DtypesInvalid:
    param_names = ["dtype"]
    params = ["scalar-string", "scalar-int", "list-string", "array-string"]
    data_dict = {
        "scalar-string": "foo",
        "scalar-int": 1,
        "list-string": ["foo"] * 1000,
        "array-string": np.array(["foo"] * 1000),
    }

    def time_pandas_dtype_invalid(self, dtype):
        try:
            pandas_dtype(self.data_dict[dtype])
        except TypeError:
            pass


class SelectDtypes:

    params = [
        tm.ALL_INT_DTYPES
        + tm.ALL_EA_INT_DTYPES
        + tm.FLOAT_DTYPES
        + tm.COMPLEX_DTYPES
        + tm.DATETIME64_DTYPES
        + tm.TIMEDELTA64_DTYPES
        + tm.BOOL_DTYPES
    ]
    param_names = ["dtype"]

    def setup(self, dtype):
        N, K = 5000, 50
        self.index = tm.makeStringIndex(N)
        self.columns = tm.makeStringIndex(K)

        def create_df(data):
            return DataFrame(data, index=self.index, columns=self.columns)

        self.df_int = create_df(np.random.randint(low=100, size=(N, K)))
        self.df_float = create_df(np.random.randn(N, K))
        self.df_bool = create_df(np.random.choice([True, False], size=(N, K)))
        self.df_string = create_df(
            np.random.choice(list(string.ascii_letters), size=(N, K))
        )

    def time_select_dtype_int_include(self, dtype):
        self.df_int.select_dtypes(include=dtype)

    def time_select_dtype_int_exclude(self, dtype):
        self.df_int.select_dtypes(exclude=dtype)

    def time_select_dtype_float_include(self, dtype):
        self.df_float.select_dtypes(include=dtype)

    def time_select_dtype_float_exclude(self, dtype):
        self.df_float.select_dtypes(exclude=dtype)

    def time_select_dtype_bool_include(self, dtype):
        self.df_bool.select_dtypes(include=dtype)

    def time_select_dtype_bool_exclude(self, dtype):
        self.df_bool.select_dtypes(exclude=dtype)

    def time_select_dtype_string_include(self, dtype):
        self.df_string.select_dtypes(include=dtype)

    def time_select_dtype_string_exclude(self, dtype):
        self.df_string.select_dtypes(exclude=dtype)


class CheckDtypes:
    def setup(self):
        self.ext_dtype = pd.Int64Dtype()
        self.np_dtype = np.dtype("int64")

    def time_is_extension_array_dtype_true(self):
        is_extension_array_dtype(self.ext_dtype)

    def time_is_extension_array_dtype_false(self):
        is_extension_array_dtype(self.np_dtype)


from .pandas_vb_common import setup  # noqa: F401 isort:skip
