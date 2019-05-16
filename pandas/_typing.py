from pathlib import Path
from typing import IO, AnyStr, Type, Union

import numpy as np

from pandas._libs import Timestamp
from pandas._libs.tslibs.period import Period
from pandas._libs.tslibs.timedeltas import Timedelta

from pandas.core.dtypes.dtypes import ExtensionDtype
from pandas.core.dtypes.generic import ABCExtensionArray

ArrayLike = Union[ABCExtensionArray, np.ndarray]
DatetimeLikeScalar = Type[Union[Period, Timestamp, Timedelta]]
Dtype = Union[str, np.dtype, ExtensionDtype]
FilePathOrBuffer = Union[str, Path, IO[AnyStr]]
