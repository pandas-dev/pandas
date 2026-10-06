"""
Shared methods for Index subclasses backed by ExtensionArray.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from pandas.util._decorators import cache_readonly

from pandas.core.indexes.base import Index

if TYPE_CHECKING:
    import numpy as np

    from pandas._typing import (
        ArrayLike,
        npt,
    )

    from pandas.core.arrays import IntervalArray
    from pandas.core.arrays._mixins import NDArrayBackedExtensionArray


class ExtensionIndex(Index):
    """
    Index subclass for indexes backed by ExtensionArray.
    """

    # The base class already passes through to _data:
    #  size, __len__, dtype

    _data: IntervalArray | NDArrayBackedExtensionArray

    # ---------------------------------------------------------------------

    def _validate_fill_value(self, value):
        """
        Convert value to be insertable to underlying array.
        """
        return self._data._validate_setitem_value(value)

    @cache_readonly
    def _isnan(self) -> npt.NDArray[np.bool_]:
        # error: Incompatible return value type (got "ExtensionArray", expected
        # "ndarray")
        return self._data.isna()  # type: ignore[return-value]


class NDArrayBackedExtensionIndex(ExtensionIndex):
    """
    Index subclass for indexes backed by NDArrayBackedExtensionArray.
    """

    _data: NDArrayBackedExtensionArray

    def _get_engine_target(self) -> np.ndarray:
        return self._data._ndarray

    def _from_join_target(self, result: np.ndarray) -> ArrayLike:
        assert result.dtype == self._data._ndarray.dtype
        return self._data._from_backing_data(result)
