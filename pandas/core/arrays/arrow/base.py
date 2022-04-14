from __future__ import annotations

import numpy as np
import pyarrow as pa

from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.dtypes import BaseMaskedDtype


class BaseArrowDtype(BaseMaskedDtype):
    """
    Base class for dtypes for BaseArrowArray subclasses.
    """

    type: pa.DataType
    na_value = pa.NA

    @cache_readonly
    def numpy_dtype(self) -> np.dtype:
        """Return an instance of our numpy dtype"""
        return self.type.to_pandas_dtype()

    @classmethod
    def from_numpy_dtype(cls, dtype: np.dtype) -> BaseArrowDtype:
        """
        Construct the ArrowDtype corresponding to the given numpy dtype.
        """
        # TODO: Fill when the other ArrowDtyes are created
        raise NotImplementedError(dtype)
