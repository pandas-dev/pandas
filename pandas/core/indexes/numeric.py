from __future__ import annotations

import numpy as np

from pandas._typing import Dtype
from pandas.util._decorators import (
    cache_readonly,
    doc,
)

from pandas.core.indexes.base import Index


class NumericIndex(Index):
    """
    Immutable numeric sequence used for indexing and alignment.

    The basic object storing axis labels for all pandas objects.
    NumericIndex is a special case of `Index` with purely numpy int/uint/float labels.

    .. versionadded:: 1.4.0

    Parameters
    ----------
    data : array-like (1-dimensional)
    dtype : NumPy dtype (default: None)
    copy : bool
        Make a copy of input ndarray.
    name : object
        Name to be stored in the index.

    Attributes
    ----------
    None

    Methods
    -------
    None

    See Also
    --------
    Index : The base pandas Index type.

    Notes
    -----
    An NumericIndex instance can **only** contain numpy int64/32/16/8, uint64/32/16/8 or
    float64/32 dtype. In particular, ``NumericIndex`` *can not* hold numpy float16
    dtype or Pandas numeric dtypes (:class:`Int64Dtype`, :class:`Int32Dtype` etc.).
    """

    _typ = "numericindex"
    _default_dtype: np.dtype | None = None
    _can_hold_strings = False

    def __new__(
        cls, data=None, dtype: Dtype | None = None, copy: bool = False, name=None
    ) -> NumericIndex:
        # temporary scaffolding, will be removed soon.
        if isinstance(data, list) and len(data) == 0:
            data = np.array([], dtype=np.int64)
        elif isinstance(data, range):
            data = np.arange(data.start, data.stop, data.step, dtype=np.int64)
        return super().__new__(
            cls, data=data, dtype=dtype, copy=copy, name=name
        )  # type: ignore[return-value]

    # ----------------------------------------------------------------
    # Indexing Methods

    @cache_readonly
    @doc(Index._should_fallback_to_positional)
    def _should_fallback_to_positional(self) -> bool:
        return False
