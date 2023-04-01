from typing import Generic

_IndexingMixinT = ...

class NDFrameIndexerBase(Generic[_IndexingMixinT]):
    name: str
    # in practice obj is either a DataFrame or a Series
    obj: _IndexingMixinT

    def __init__(self, name: str, obj: _IndexingMixinT) -> None: ...
    @property
    def ndim(self) -> int: ...
