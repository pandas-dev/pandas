from pathlib import Path
from typing import (
    IO,
    TYPE_CHECKING,
    Any,
    AnyStr,
    Callable,
    Collection,
    Dict,
    Hashable,
    List,
    Mapping,
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
    from pandas.core.generic import NDFrame  # noqa: F401
    from pandas import Interval  # noqa: F401
    from pandas.core.series import Series  # noqa: F401
    from pandas.core.frame import DataFrame  # noqa: F401

# array-like

AnyArrayLike = TypeVar("AnyArrayLike", "ExtensionArray", "Index", "Series", np.ndarray)
ArrayLike = TypeVar("ArrayLike", "ExtensionArray", np.ndarray)

# scalars

PythonScalar = Union[str, int, float, bool]
DatetimeLikeScalar = TypeVar("DatetimeLikeScalar", "Period", "Timestamp", "Timedelta")
PandasScalar = Union["Period", "Timestamp", "Timedelta", "Interval"]
Scalar = Union[PythonScalar, PandasScalar]

# other

Dtype = Union[str, np.dtype, "ExtensionDtype"]
DtypeObj = Union[np.dtype, "ExtensionDtype"]
FilePathOrBuffer = Union[str, Path, IO[AnyStr]]

# FrameOrSeriesUnion  means either a DataFrame or a Series. E.g.
# `def func(a: FrameOrSeriesUnion) -> FrameOrSeriesUnion: ...` means that if a Series
# is passed in, either a Series or DataFrame is returned, and if a DataFrame is passed
# in, either a DataFrame or a Series is returned.
FrameOrSeriesUnion = Union["DataFrame", "Series"]

# FrameOrSeries is stricter and ensures that the same subclass of NDFrame always is
# used. E.g. `def func(a: FrameOrSeries) -> FrameOrSeries: ...` means that if a
# Series is passed into a function, a Series is always returned and if a DataFrame is
# passed in, a DataFrame is always returned.
FrameOrSeries = TypeVar("FrameOrSeries", bound="NDFrame")

Axis = Union[str, int]
Label = Optional[Hashable]
Level = Union[Label, int]
Ordered = Optional[bool]
JSONSerializable = Union[PythonScalar, List, Dict]
Axes = Collection

# For functions like rename that convert one label to another
Renamer = Union[Mapping[Label, Any], Callable[[Label], Label]]

# to maintain type information across generic functions and parametrization
T = TypeVar("T")
