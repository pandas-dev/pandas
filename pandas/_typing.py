from pathlib import Path
from typing import IO, AnyStr, Type, TypeVar, Union

import numpy as np

from pandas._libs import Timestamp
from pandas._libs.tslibs.period import Period
from pandas._libs.tslibs.timedeltas import Timedelta

from pandas.core.dtypes.dtypes import ExtensionDtype
from pandas.core.dtypes.generic import (
    ABCExtensionArray, ABCIndexClass, ABCSeries, ABCSparseSeries)

from pandas.core.arrays import DatetimeArray, PeriodArray, TimedeltaArray

AnyArrayLike = Union[ABCExtensionArray,
                     ABCIndexClass,
                     ABCSeries,
                     ABCSparseSeries,
                     np.ndarray]
ArrayLike = Union[ABCExtensionArray, np.ndarray]
DatetimeLikeArray = TypeVar('DatetimeLikeArray', DatetimeArray,
                            PeriodArray, TimedeltaArray)
DatetimeLikeScalar = Type[Union[Period, Timestamp, Timedelta]]
Dtype = Union[str, np.dtype, ExtensionDtype]
FilePathOrBuffer = Union[str, Path, IO[AnyStr]]
