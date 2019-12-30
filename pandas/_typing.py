from pathlib import Path
from typing import (
    IO,
    TYPE_CHECKING,
    AnyStr,
    Collection,
    Dict,
    List,
    Optional,
    TypeVar,
    Union,
)

import numpy as np

# To prevent import cycles place any internal imports in the branch below
# and use a string literal forward reference to it in subsequent types
# https://mypy.readthedocs.io/en/latest/common_issues.html#import-cycles
if TYPE_CHECKING:
    from pandas._libs import Period, Timedelta, Timestamp  # noqa: F401
    from pandas.core.arrays.base import ExtensionArray  # noqa: F401
    from pandas.core.dtypes.dtypes import ExtensionDtype  # noqa: F401
    from pandas.core.indexes.base import Index  # noqa: F401
    from pandas.core.series import Series  # noqa: F401
    from pandas.core.generic import NDFrame  # noqa: F401


AnyArrayLike = TypeVar("AnyArrayLike", "ExtensionArray", "Index", "Series", np.ndarray)
ArrayLike = TypeVar("ArrayLike", "ExtensionArray", np.ndarray)
DatetimeLikeScalar = TypeVar("DatetimeLikeScalar", "Period", "Timestamp", "Timedelta")
Dtype = Union[str, np.dtype, "ExtensionDtype"]
FilePathOrBuffer = Union[str, Path, IO[AnyStr]]

FrameOrSeries = TypeVar("FrameOrSeries", bound="NDFrame")
Scalar = Union[str, int, float, bool]
Axis = Union[str, int]
Ordered = Optional[bool]
JSONSerializable = Union[Scalar, List, Dict]

Axes = Collection

# to maintain type information across generic functions and parametrization
_T = TypeVar("_T")
