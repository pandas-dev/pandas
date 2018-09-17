import copy
import itertools

import numpy as np
import pyarrow as pa
import pandas as pd
from pandas.api.extensions import (
    ExtensionDtype, ExtensionArray
)


# @register_extension_dtype
class ArrowBoolDtype(ExtensionDtype):

    type = np.bool_
    kind = 'b'
    name = 'arrow_bool'
    na_value = pa.NULL

    @classmethod
    def construct_from_string(cls, string):
        if string == cls.name:
            return cls()
        else:
            raise TypeError("Cannot construct a '{}' from "
                            "'{}'".format(cls, string))

    @classmethod
    def construct_array_type(cls):
        return ArrowBoolArray

    def _is_boolean(self):
        return True


class ArrowBoolArray(ExtensionArray):
    def __init__(self, values):
        if not isinstance(values, pa.ChunkedArray):
            raise ValueError

        assert values.type == pa.bool_()
        self._data = values
        self._dtype = ArrowBoolDtype()

    def __repr__(self):
        return "ArrowBoolArray({})".format(repr(self._data))

    @classmethod
    def from_scalars(cls, values):
        arr = pa.chunked_array([pa.array(np.asarray(values))])
        return cls(arr)

    @classmethod
    def from_array(cls, arr):
        assert isinstance(arr, pa.Array)
        return cls(pa.chunked_array([arr]))

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        return cls.from_scalars(scalars)

    def __getitem__(self, item):
        return self._data.to_pandas()[item]

    def __len__(self):
        return len(self._data)

    @property
    def dtype(self):
        return self._dtype

    @property
    def nbytes(self):
        return sum(x.size for chunk in self._data.chunks
                   for x in chunk.buffers()
                   if x is not None)

    def isna(self):
        return pd.isna(self._data.to_pandas())

    def take(self, indices, allow_fill=False, fill_value=None):
        from pandas.core.algorithms import take
        data = self._data.to_pandas()

        if allow_fill and fill_value is None:
            fill_value = self.dtype.na_value

        result = take(data, indices, fill_value=fill_value,
                      allow_fill=allow_fill)
        return self._from_sequence(result, dtype=self.dtype)

    def copy(self, deep=False):
        if deep:
            return copy.deepcopy(self._data)
        else:
            return copy.copy(self._data)

    def _concat_same_type(cls, to_concat):
        chunks = list(itertools.chain.from_iterable(x._data.chunks
                                                    for x in to_concat))
        arr = pa.chunked_array(chunks)
        return cls(arr)
