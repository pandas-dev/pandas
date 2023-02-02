from __future__ import annotations

import numpy as np

from pandas._typing import Dtype
from pandas.util._decorators import (
    cache_readonly,
    doc,
)

from pandas.core.indexes.base import Index


class NumericIndex(Index):
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
