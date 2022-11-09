from __future__ import annotations

import numpy as np

from pandas._typing import Dtype

from pandas.core.arrays import ExtensionArray
from pandas.core.indexes.range import RangeIndex


class NoIndex(RangeIndex):
    _typ = "noindex"

    def __new__(
        cls,
        len: int,
        dtype: Dtype | None = None,
        copy: bool = False,
        name=None,
    ):
        rng = range(0, len)
        return cls._simple_new(rng, name=name)

    def __mul__(self, other):
        raise NotImplementedError()

    def __rmul__(self, other):
        raise NotImplementedError()

    def __add__(self, other):
        raise NotImplementedError()

    def __radd__(self, other):
        raise NotImplementedError()

    def __div__(self, other):
        raise NotImplementedError()

    def __rdiv__(self, other):
        raise NotImplementedError()

    def __sub__(self, other):
        raise NotImplementedError()

    def __rsub__(self, other):
        raise NotImplementedError()

    def __pow__(self, other):
        raise NotImplementedError()

    @property
    def name(self):
        return None

    @name.setter
    def name(self, new_name):
        if new_name is not None:
            raise TypeError("Can't set name of NoIndex!")

    def _set_names(self, values, *, level=None) -> None:
        raise TypeError("Can't set name of NoIndex!")

    def __repr__(self) -> str:
        return f"NoIndex(len={self.stop})"

    def append(self, other):
        if not isinstance(other, list):
            other = [other]
        length = len(self)
        for _other in other:
            if not isinstance(_other, NoIndex):
                raise TypeError(
                    f"Can only concatenate NoIndex to NoIndex - got {_other}"
                )
            length += len(_other)
        return NoIndex(length)

    def __getitem__(self, key):
        _super = super().__getitem__(key)
        try:
            return NoIndex(len(_super))
        except TypeError:
            return _super

    def get_loc(self, key, method=None, tolerance=None):
        from pandas.core import common as com

        if not com.is_bool_indexer(key):
            raise IndexError("Cannot use label-based indexing on NoIndex!")
        return super().get_loc(key, method, tolerance)

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
        result.name = name
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

    def insert(self, loc, item):
        raise NotImplementedError()

    def reindex(self, *args, **kwargs):
        raise NotImplementedError(
            "Can't reindex a DataFrame without an index. First, give it an index."
        )

    def join(
        self,
        other,
        *,
        how: str = "left",
        level=None,
        return_indexers: bool = False,
        sort: bool = False,
    ) -> NoIndex:
        if not isinstance(other, NoIndex):
            raise TypeError("Can't join NoIndex with Index")
        if not len(self) == len(other):
            raise TypeError("Can't join NoIndex of different lengths")
        return super().join(other, how=how, return_indexers=return_indexers, sort=sort)
