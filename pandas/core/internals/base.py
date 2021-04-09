"""
Base class for the internal managers. Both BlockManager and ArrayManager
inherit from this class.
"""
from __future__ import annotations

from typing import (
    Any,
    TypeVar,
)

import numpy as np

from pandas._typing import (
    DtypeObj,
    Shape,
    final,
)
from pandas.errors import AbstractMethodError
from pandas.util._validators import validate_bool_kwarg

from pandas.core.dtypes.cast import (
    astype_array_safe,
    find_common_type,
)

from pandas.core.base import PandasObject
from pandas.core.construction import extract_array
from pandas.core.indexes.api import (
    Index,
    ensure_index,
)
from pandas.core.internals.blocks import to_native_types

T = TypeVar("T", bound="DataManager")


class DataManager(PandasObject):

    # TODO share more methods/attributes

    axes: list[Index]

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

    def reindex_axis(
        self,
        new_index,
        axis: int,
        method=None,
        limit=None,
        fill_value=None,
        copy: bool = True,
        consolidate: bool = True,
        only_slice: bool = False,
    ):
        """
        Conform data manager to new index.
        """
        new_index = ensure_index(new_index)
        new_index, indexer = self.axes[axis].reindex(
            new_index, method=method, limit=limit
        )

        return self.reindex_indexer(
            new_index,
            indexer,
            axis=axis,
            fill_value=fill_value,
            copy=copy,
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
        align_keys: list[str] | None = None,
        ignore_failures: bool = False,
        **kwargs,
    ) -> T:
        raise AbstractMethodError(self)

    def apply_with_block(
        self: T,
        f,
        align_keys: list[str] | None = None,
        ignore_failures: bool = False,
        **kwargs,
    ) -> T:
        raise AbstractMethodError(self)

    @final
    def isna(self: T, func) -> T:
        return self.apply("apply", func=func)

    @final
    def where(self: T, other, cond, align: bool, errors: str) -> T:
        if align:
            align_keys = ["other", "cond"]
        else:
            align_keys = ["cond"]
            other = extract_array(other, extract_numpy=True)

        return self.apply_with_block(
            "where",
            align_keys=align_keys,
            other=other,
            cond=cond,
            errors=errors,
        )

    @final
    def putmask(self: T, mask, new, align: bool = True) -> T:
        if align:
            align_keys = ["new", "mask"]
        else:
            align_keys = ["mask"]
            new = extract_array(new, extract_numpy=True)

        return self.apply_with_block(
            "putmask",
            align_keys=align_keys,
            mask=mask,
            new=new,
        )

    @final
    def fillna(self: T, value, limit, inplace: bool, downcast) -> T:
        return self.apply_with_block(
            "fillna", value=value, limit=limit, inplace=inplace, downcast=downcast
        )

    @final
    def downcast(self: T) -> T:
        return self.apply_with_block("downcast")

    @final
    def replace(self: T, to_replace, value, inplace: bool, regex: bool) -> T:
        assert np.ndim(value) == 0, value
        return self.apply_with_block(
            "replace", to_replace=to_replace, value=value, inplace=inplace, regex=regex
        )

    def replace_list(
        self: T,
        src_list: list[Any],
        dest_list: list[Any],
        inplace: bool = False,
        regex: bool = False,
    ) -> T:
        """ do a list replace """
        inplace = validate_bool_kwarg(inplace, "inplace")

        return self.apply_with_block(
            "_replace_list",
            src_list=src_list,
            dest_list=dest_list,
            inplace=inplace,
            regex=regex,
        )

    @final
    def astype(self: T, dtype, copy: bool = False, errors: str = "raise") -> T:
        return self.apply(astype_array_safe, dtype=dtype, copy=copy, errors=errors)

    @final
    def to_native_types(self: T, **kwargs) -> T:
        return self.apply(to_native_types, **kwargs)


class SingleDataManager(DataManager):
    ndim = 1

    @property
    def array(self):
        """
        Quick access to the backing array of the Block or SingleArrayManager.
        """
        return self.arrays[0]  # type: ignore[attr-defined]


def interleaved_dtype(dtypes: list[DtypeObj]) -> DtypeObj | None:
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
