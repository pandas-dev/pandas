"""
Base class for the internal managers. Both BlockManager and ArrayManager
inherit from this class.
"""
from typing import (
    List,
    Optional,
    TypeVar,
)

from pandas._typing import (
    DtypeObj,
    Shape,
    final,
)
from pandas.errors import AbstractMethodError

from pandas.core.dtypes.cast import find_common_type

from pandas.core.base import PandasObject
from pandas.core.indexes.api import Index

T = TypeVar("T", bound="DataManager")


class DataManager(PandasObject):

    # TODO share more methods/attributes

    axes: List[Index]

    @property
    def items(self) -> Index:
        raise AbstractMethodError(self)

    def __len__(self) -> int:
        return len(self.items)

    @property
    def ndim(self) -> int:
        return len(self.axes)

    @property
    def shape(self) -> Shape:
        return tuple(len(ax) for ax in self.axes)

    def reindex_indexer(
        self: T,
        new_axis,
        indexer,
        axis: int,
        fill_value=None,
        allow_dups: bool = False,
        copy: bool = True,
        consolidate: bool = True,
        only_slice: bool = False,
    ) -> T:
        raise AbstractMethodError(self)

    @final
    def reindex_axis(
        self: T,
        new_index: Index,
        axis: int,
        fill_value=None,
        consolidate: bool = True,
        only_slice: bool = False,
    ) -> T:
        """
        Conform data manager to new index.
        """
        new_index, indexer = self.axes[axis].reindex(new_index)

        return self.reindex_indexer(
            new_index,
            indexer,
            axis=axis,
            fill_value=fill_value,
            copy=False,
            consolidate=consolidate,
            only_slice=only_slice,
        )

    def _equal_values(self: T, other: T) -> bool:
        """
        To be implemented by the subclasses. Only check the column values
        assuming shape and indexes have already been checked.
        """
        raise AbstractMethodError(self)

    def equals(self, other: object) -> bool:
        """
        Implementation for DataFrame.equals
        """
        if not isinstance(other, DataManager):
            return False

        self_axes, other_axes = self.axes, other.axes
        if len(self_axes) != len(other_axes):
            return False
        if not all(ax1.equals(ax2) for ax1, ax2 in zip(self_axes, other_axes)):
            return False

        return self._equal_values(other)

    def apply(
        self: T,
        f,
        align_keys: Optional[List[str]] = None,
        ignore_failures: bool = False,
        **kwargs,
    ) -> T:
        raise AbstractMethodError(self)

    def isna(self: T, func) -> T:
        return self.apply("apply", func=func)


class SingleDataManager(DataManager):
    ndim = 1

    @property
    def array(self):
        """
        Quick access to the backing array of the Block or SingleArrayManager.
        """
        return self.arrays[0]  # type: ignore[attr-defined]


def interleaved_dtype(dtypes: List[DtypeObj]) -> Optional[DtypeObj]:
    """
    Find the common dtype for `blocks`.

    Parameters
    ----------
    blocks : List[DtypeObj]

    Returns
    -------
    dtype : np.dtype, ExtensionDtype, or None
        None is returned when `blocks` is empty.
    """
    if not len(dtypes):
        return None

    return find_common_type(dtypes)
