from typing import Generic, NamedTuple
from typing_extensions import TypeVar

import numpy as np

_FltT_co = TypeVar("_FltT_co", default=np.float64, covariant=True)

class ConfidenceInterval(NamedTuple, Generic[_FltT_co]):
    low: _FltT_co
    high: _FltT_co
