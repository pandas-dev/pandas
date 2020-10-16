"""
Base class for the internal managers. Both BlockManager and ArrayManager
inherit from this class.
"""
from pandas.core.base import PandasObject
from pandas.core.indexes.api import ensure_index


class DataManager(PandasObject):

    # TODO share more methods/attributes

    def __len__(self) -> int:
        return len(self.items)

    @property
    def ndim(self) -> int:
        return len(self.axes)

    def reindex_axis(
        self,
        new_index,
        axis: int,
        method=None,
        limit=None,
        fill_value=None,
        copy: bool = True,
    ):
        """
        Conform block manager to new index.
        """
        new_index = ensure_index(new_index)
        new_index, indexer = self.axes[axis].reindex(
            new_index, method=method, limit=limit
        )

        return self.reindex_indexer(
            new_index, indexer, axis=axis, fill_value=fill_value, copy=copy
        )
