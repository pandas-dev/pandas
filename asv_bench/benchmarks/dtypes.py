from pandas.api.types import pandas_dtype

import numpy as np
from .pandas_vb_common import (
    numeric_dtypes, datetime_dtypes, string_dtypes, extension_dtypes)


_numpy_dtypes = list(map(np.dtype, (numeric_dtypes +
                                    datetime_dtypes +
                                    string_dtypes)))
_dtypes = _numpy_dtypes + extension_dtypes


class Dtypes(object):
    params = (_dtypes +
              list(map(lambda dt: dt.name, _dtypes)))
    param_names = ['dtype']

    def time_pandas_dtype(self, dtype):
        pandas_dtype(dtype)


class DtypesInvalid(object):
    params = ['foo', 1, ['foo'] * 1000, np.array(['foo'] * 1000)]
    param_names = ['dtype']

    def time_pandas_dtype_invalid(self, dtype):
        try:
            pandas_dtype(dtype)
        except TypeError:
            pass

    param_names = ['dtype']
    params = ['scalar-string', 'scalar-int', 'list-string', 'array-string']
    data_dict = {'scalar-string': 'foo',
                 'scalar-int': 1,
                 'list-string': ['foo'] * 1000,
                 'array-string': np.array(['foo'] * 1000)}

    def setup(self, dtype):
        self.data = self.data_dict[dtype]


from .pandas_vb_common import setup  # noqa: F401
