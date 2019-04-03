from pathlib import Path
from typing import IO, AnyStr, Type, Union

import numpy as np

from pandas.core.dtypes.dtypes import ExtensionDtype
from pandas.core.dtypes.generic import ABCExtensionArray

ArrayLike = Union[ABCExtensionArray, np.ndarray]
SparseDtype = Union[str, np.dtype, ExtensionDtype,
                    Type[float], Type[int], Type[object]]
FilePathOrBuffer = Union[str, Path, IO[AnyStr]]
