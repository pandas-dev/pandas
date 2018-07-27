import numpy as np

from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.dtypes.dtypes import Registry
from pandas import compat


class SparseDtype(ExtensionDtype):

    def __init__(self, dtype=np.float64):
        self._dtype = np.dtype(dtype)

    @property
    def kind(self):
        return self.dtype.kind

    @property
    def dtype(self):
        return self._dtype

    @property
    def type(self):
        return self.dtype.type

    @property
    def subdtype(self):
        return self.type

    @property
    def name(self):
        return 'Sparse[{}]'.format(self.dtype.name)

    def __repr__(self):
        return self.name

    @classmethod
    def construct_array_type(cls):
        from .array import SparseArray
        return SparseArray

    @classmethod
    def construct_from_string(cls, string):
        if string.startswith("Sparse"):
            sub_type = cls._parse_subtype(string)
            try:
                return SparseDtype(sub_type)
            except Exception:
                raise TypeError
        else:
            raise TypeError

    @staticmethod
    def _parse_subtype(dtype):
        if dtype.startswith("Sparse["):
            sub_type = dtype[7:-1]
        elif dtype == "Sparse":
            sub_type = 'float64'
        else:
            raise ValueError
        return sub_type

    @classmethod
    def is_dtype(cls, dtype):
        dtype = getattr(dtype, 'dtype', dtype)
        if isinstance(dtype, compat.string_types) and dtype.startswith("Sparse"):
            dtype = np.dtype(cls._parse_subtype(dtype))
        elif isinstance(dtype, cls):
            return True
        return isinstance(dtype, np.dtype) or dtype == 'Sparse'


Registry.register(SparseDtype)
