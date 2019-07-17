from pathlib import Path
from typing import IO, TYPE_CHECKING, AnyStr, TypeVar, Union

import numpy as np

if TYPE_CHECKING:  # Use for any internal imports
    from pandas._libs import Period, Timedelta, Timestamp

    from pandas.core.arrays.base import ExtensionArray
    from pandas.core.dtypes.dtypes import ExtensionDtype
    from pandas.core.dtypes.generic import (
        ABCDataFrame,
        ABCExtensionArray,
        ABCIndexClass,
        ABCSeries,
        ABCSparseSeries,
    )
    from pandas.core.indexes.base import Index
    from pandas.core.frame import DataFrame
    from pandas.core.series import Series
    from pandas.core.sparse.series import SparseSeries


AnyArrayLike = TypeVar(
    "AnyArrayLike", "ExtensionArray", "Index", "Series", "SparseSeries", np.ndarray
)
ArrayLike = TypeVar("ArrayLike", "ExtensionArray", np.ndarray)
DatetimeLikeScalar = TypeVar("DatetimeLikeScalar", "Period", "Timestamp", "Timedelta")
Dtype = Union[str, np.dtype, "ExtensionDtype"]
FilePathOrBuffer = Union[str, Path, IO[AnyStr]]

FrameOrSeries = TypeVar("FrameOrSeries", "Series", "DataFrame")
Scalar = Union[str, int, float]
