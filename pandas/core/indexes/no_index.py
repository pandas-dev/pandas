from __future__ import annotations
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
            assert isinstance(_other, NoIndex)
            length += len(_other)
        return NoIndex(length)

    def __getitem__(self, key):
        return NoIndex(len(super().__getitem__(key)))
