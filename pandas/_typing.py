from pathlib import Path
from typing import IO, AnyStr, TypeVar, Union

import numpy as np

from pandas._libs import Timestamp
from pandas._libs.tslibs.period import Period
from pandas._libs.tslibs.timedeltas import Timedelta

from pandas.core.dtypes.dtypes import ExtensionDtype
from pandas.core.dtypes.generic import (
    ABCExtensionArray, ABCIndexClass, ABCSeries, ABCSparseSeries)

from pandas.core.arrays.datetimes import DatetimeArray
from pandas.core.arrays.period import PeriodArray
from pandas.core.arrays.timedeltas import TimedeltaArray

AnyArrayLike = TypeVar('AnyArrayLike',
                       ABCExtensionArray,
                       ABCIndexClass,
                       ABCSeries,
                       ABCSparseSeries,
                       np.ndarray)
DatetimeLikeArray = TypeVar('DatetimeLikeArray', DatetimeArray,
                            PeriodArray, TimedeltaArray)

ArrayLike = TypeVar('ArrayLike', ABCExtensionArray, np.ndarray)
DatetimeLikeScalar = TypeVar('DatetimeLikeScalar', Period, Timestamp,
                             Timedelta)
Dtype = Union[str, np.dtype, ExtensionDtype]
FilePathOrBuffer = Union[str, Path, IO[AnyStr]]
