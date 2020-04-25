from typing import Any, Sequence, TypeVar

import numpy as np

from pandas.errors import AbstractMethodError

from pandas.core.algorithms import take
from pandas.core.arrays.base import ExtensionArray

_T = TypeVar("_T", bound="NDArrayBackedExtensionArray")


class NDArrayBackedExtensionArray(ExtensionArray):
    """
    ExtensionArray that is backed by a single NumPy ndarray.
    """

    _ndarray: np.ndarray

    def _from_backing_data(self: _T, arr: np.ndarray) -> _T:
        """
        Construct a new ExtensionArray `new_array` with `arr` as its _ndarray.

        This should round-trip:
            self == self._from_backing_data(self._ndarray)
        """
        raise AbstractMethodError(self)

    # ------------------------------------------------------------------------

    def take(
        self: _T,
        indices: Sequence[int],
        allow_fill: bool = False,
        fill_value: Any = None,
    ) -> _T:
        if allow_fill:
            fill_value = self._validate_fill_value(fill_value)

        new_data = take(
            self._ndarray, indices, allow_fill=allow_fill, fill_value=fill_value,
        )
        return self._from_backing_data(new_data)

    def _validate_fill_value(self, fill_value):
        """
        If a fill_value is passed to `take` convert it to a representation
        suitable for self._ndarray, raising ValueError if this is not possible.

        Parameters
        ----------
        fill_value : object

        Returns
        -------
        fill_value : native representation

        Raises
        ------
        ValueError
        """
        raise AbstractMethodError(self)
