"""
This module provide a numpy-boolean boolean array
"""

import numpy as np

from pandas.api.extensions import take
from pandas.core.arrays.mask._base import MaskArray, MaskDtype


class NumpyMaskDtype(MaskDtype):

    na_value = np.nan

    @classmethod
    def construct_array_type(cls):
        return NumpyMaskArray


class NumpyMaskArray(MaskArray):
    """Generic class which can be used to represent missing data.
    """

    dtype = NumpyMaskDtype()

    @classmethod
    def from_scalars(cls, values):
        arr = np.asarray(values).astype(np.bool_, copy=False)
        return cls(arr, copy=False)

    def __init__(self, mask, copy=True):
        """
        Parameters
        ----------
        mask : numpy array
            Mask of missing values.
        """
        assert isinstance(mask, np.ndarray)
        assert mask.dtype == np.bool_

        if copy:
            mask = mask.copy()
        self._data = mask

    def __setitem__(self, key, value):
        self._data[key] = value

    def __array__(self, dtype=None):
        return self._data

    def __iter__(self):
        return iter(self._data)

    @property
    def nbytes(self):
        return self._data.nbytes

    def reshape(self, shape, **kwargs):
        return np.array(self, copy=False).reshape(shape, **kwargs)

    def astype(self, dtype, copy=True):
        # needed to fix this astype for the Series constructor.
        if isinstance(dtype, type(self.dtype)) and dtype == self.dtype:
            if copy:
                return self.copy()
            return self
        return super(NumpyMaskArray, self).astype(dtype, copy)

    def take(self, indices, allow_fill=False, fill_value=None, axis=None):
        # TODO: had to add axis here
        data = self._data

        if allow_fill and fill_value is None:
            fill_value = self.dtype.na_value

        result = take(data, indices, fill_value=fill_value,
                      allow_fill=allow_fill)
        return self._from_sequence(result, dtype=self.dtype)

    def _concat_same_type(cls, to_concat):
        concat = np.concatenate(to_concat)
        return cls.from_scalars(concat)
