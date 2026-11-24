from pandas.core.indexing import IndexingMixin

class NDFrameIndexerBase[IndexingMixinT: IndexingMixin]:
    name: str
    # in practice obj is either a DataFrame or a Series
    obj: IndexingMixinT

    def __init__(self, name: str, obj: IndexingMixinT) -> None: ...
    @property
    def ndim(self) -> int: ...
