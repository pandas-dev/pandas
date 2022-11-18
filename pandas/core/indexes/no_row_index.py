from __future__ import annotations

import numpy as np

from pandas._typing import Dtype

from pandas.core import common as com
from pandas.core.arrays import ExtensionArray
from pandas.core.indexes.range import RangeIndex


class NoRowIndex(RangeIndex):
    _typ = "NoRowIndex"

    def __new__(
        cls,
        len: int,
        dtype: Dtype | None = None,
        copy: bool = False,
        name=None,
    ) -> NoRowIndex:
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
    def name(self) -> None:
        return None

    @name.setter
    def name(self, new_name):
        if new_name is not None:
            raise TypeError("Can't set name of NoRowIndex!")

    def _set_names(self, values, *, level=None):
        raise TypeError("Can't set name of NoRowIndex!")

    def __repr__(self) -> str:
        return f"NoRowIndex(len={self.stop})"

    def append(self, other) -> NoRowIndex:
        if not isinstance(other, list):
            other = [other]
        length = len(self)
        for _other in other:
            if not isinstance(_other, NoRowIndex):
                raise TypeError(
                    f"Can only concatenate NoRowIndex to NoRowIndex - got {_other}"
                )
            length += len(_other)
        return NoRowIndex(length)

    def __getitem__(self, key):
        _super = super().__getitem__(key)
        try:
            return NoRowIndex(len(_super))
        except TypeError:
            # wait, why do we need this?
            return _super

    def get_loc(self, key, method=None, tolerance=None):
        if not com.is_bool_indexer(key):
            raise IndexError("Cannot use label-based indexing on NoRowIndex!")
        return super().get_loc(key, method, tolerance)

    @property
    def _constructor(self) -> type[NoRowIndex]:  # type: ignore[override]
        """return the class to use for construction"""
        return NoRowIndex

    @classmethod
    def _simple_new(cls, values, name=None) -> NoRowIndex:
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
    ) -> NoRowIndex:
        if not isinstance(other, NoRowIndex):
            raise TypeError("Can't join NoRowIndex with Index")
        if not len(self) == len(other):
            raise TypeError("Can't join NoRowIndex of different lengths")
        return super().join(other, how=how, return_indexers=return_indexers, sort=sort)
