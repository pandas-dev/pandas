from __future__ import annotations

import numpy as np

from pandas._typing import Dtype

from pandas.core.arrays import ExtensionArray
from pandas.core.indexes.range import RangeIndex


class NoIndex(RangeIndex):
    def __new__(
        cls,
        len: int,
        dtype: Dtype | None = None,
        copy: bool = False,
    ):
        rng = range(0, len)
        return cls._simple_new(rng)

    def append(self, other):
        if not isinstance(other, list):
            other = [other]
        length = len(self)
        for _other in other:
            if not isinstance(_other, NoIndex):
                raise TypeError(f'Can only concatenate NoIndex to NoIndex - got {_other}')
            length += len(_other)
        return NoIndex(length)

    def __getitem__(self, key):
        _super = super().__getitem__(key)
        try:
            return NoIndex(len(_super))
        except TypeError:
            return _super

    def get_loc(self, key, method=None, tolerance=None):
        breakpoint()
        from pandas.core import common as com

        if not com.is_bool_indexer(key):
            raise TypeError("Cannot use label-based indexing on NoIndex!")

    @property
    def _constructor(self):  # type: ignore[override]
        """return the class to use for construction"""
        return NoIndex

    @classmethod
    def _simple_new(cls, values, name=None):
        result = object.__new__(cls)
        assert isinstance(values, (range, np.ndarray, ExtensionArray))
        values = range(len(values))

        result._range = values
        result._cache = {}
        result._reset_identity()
        return result

    def sort_values(
        self,
        return_indexer: bool = False,
        ascending: bool = True,
        na_position: str = "last",
        key=None,
    ):
        raise NotImplementedError()

    def delete(self, loc):
        raise NotImplementedError()

    def insert(self, loc):
        raise NotImplementedError()

    def reindex(self, *args, **kwargs):
        raise NotImplementedError(
            "Can't reindex a DataFrame without an index. "
            "First, give it an index."
        )