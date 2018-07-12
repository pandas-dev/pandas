import numpy as np

from pandas.core.dtypes.base import ExtensionDtype


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
    def name(self):
        return 'sparse'

    @classmethod
    def construct_array_type(cls):
        from .array import SparseArray
        return SparseArray

    @classmethod
    def construct_from_string(cls, string):
        if string == 'sparse':
            string = 'float64'
        try:
            return SparseDtype(string)
        except:
            raise TypeError

