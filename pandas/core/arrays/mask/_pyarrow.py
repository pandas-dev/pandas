"""Rudimentary Apache Arrow-backed ExtensionArray.

At the moment, just a boolean array / type is implemented.
Eventually, we'll want to parametrize the type and support
multiple dtypes. Not all methods are implemented yet, and the
current implementation is not efficient.
"""
from distutils.version import LooseVersion
import itertools

import numpy as np

from pandas.api.extensions import take
from pandas.core.arrays.mask._base import MaskArray, MaskDtype

# we require pyarrow >= 0.10.0

try:
    import pyarrow as pa
    if pa.__version__ < LooseVersion('0.10.0'):
        raise ImportError("pyarrow minimum for bool suppport is 0.10.0")
except ImportError:
    raise


class ArrowMaskDtype(MaskDtype):

    na_value = pa.NULL

    @classmethod
    def construct_array_type(cls):
        return ArrowMaskArray


class ArrowMaskArray(MaskArray):

    dtype = ArrowMaskDtype()

    @classmethod
    def from_scalars(cls, values):
        values = np.asarray(values).astype(np.bool_, copy=False)
        arr = pa.chunked_array([values])
        return cls(arr)

    def __init__(self, values, copy=False):

        # TODO: we need to rationalize the return types from
        # various ops, we oftentimes return boolean array arrays
        # but not chunked ones
        if not isinstance(values, pa.ChunkedArray):
            values = pa.chunked_array([values])
        assert values.type == pa.bool_()
        if copy:
            values = values.copy()

        self._data = values

    def __setitem__(self, key, value):
        # TODO: hack-a-minute
        data = np.array(self._data)
        data[key] = value
        self._data = pa.array(data)

    def astype(self, dtype, copy=True):
        # needed to fix this astype for the Series constructor.
        if isinstance(dtype, type(self.dtype)) and dtype == self.dtype:
            if copy:
                return self.copy()
            return self
        return super(ArrowMaskArray, self).astype(dtype, copy)

    @property
    def nbytes(self):
        return sum(x.size for chunk in self._data.chunks
                   for x in chunk.buffers()
                   if x is not None)

    def take(self, indices, allow_fill=False, fill_value=None, axis=None):
        # TODO: had to add axis here
        data = self._data.to_pandas()

        if allow_fill and fill_value is None:
            fill_value = self.dtype.na_value

        result = take(data, indices, fill_value=fill_value,
                      allow_fill=allow_fill)
        return self._from_sequence(result, dtype=self.dtype)

    def _concat_same_type(cls, to_concat):
        chunks = list(itertools.chain.from_iterable(x._data.chunks
                                                    for x in to_concat))
        arr = pa.chunked_array(chunks)
        return cls(arr)

    def __array__(self, dtype=None):
        return np.array(self._data, copy=False)
