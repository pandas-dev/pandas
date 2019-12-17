import numpy as np

from pandas.api.types import pandas_dtype

from .pandas_vb_common import (
    datetime_dtypes,
    extension_dtypes,
    lib,
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


class InferDtypes:
    param_names = ["dtype"]
    data_dict = {
        "np-object": np.array([1] * 100000, dtype="O"),
        "py-object": [1] * 100000,
        "np-null": np.array([1] * 50000 + [np.nan] * 50000),
        "py-null": [1] * 50000 + [None] * 50000,
        "np-int": np.array([1] * 100000, dtype=int),
        "np-floating": np.array([1.0] * 100000, dtype=float),
        "empty": [],
        "bytes": [b"a"] * 100000,
    }
    params = list(data_dict.keys())

    def time_infer_skipna(self, dtype):
        lib.infer_dtype(self.data_dict[dtype], skipna=True)

    def time_infer(self, dtype):
        lib.infer_dtype(self.data_dict[dtype], skipna=False)


from .pandas_vb_common import setup  # noqa: F401 isort:skip
