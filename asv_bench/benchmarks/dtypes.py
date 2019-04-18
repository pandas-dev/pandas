from pandas.api.types import pandas_dtype

import numpy as np
from .pandas_vb_common import (
    numeric_dtypes, datetime_dtypes, string_dtypes, extension_dtypes)


_numpy_dtypes = [np.dtype(dtype)
                 for dtype in (numeric_dtypes +
                               datetime_dtypes +
                               string_dtypes)]
_dtypes = _numpy_dtypes + extension_dtypes


class Dtypes:
    params = (_dtypes +
              list(map(lambda dt: dt.name, _dtypes)))
    param_names = ['dtype']

    def time_pandas_dtype(self, dtype):
        pandas_dtype(dtype)


class DtypesInvalid:
    param_names = ['dtype']
    params = ['scalar-string', 'scalar-int', 'list-string', 'array-string']
    data_dict = {'scalar-string': 'foo',
                 'scalar-int': 1,
                 'list-string': ['foo'] * 1000,
                 'array-string': np.array(['foo'] * 1000)}

    def time_pandas_dtype_invalid(self, dtype):
        try:
            pandas_dtype(self.data_dict[dtype])
        except TypeError:
            pass


from .pandas_vb_common import setup  # noqa: F401
