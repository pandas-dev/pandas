from datetime import (
    timedelta,
    tzinfo,
)
from typing import (
    Iterable,
    Optional,
    Union,
)

import numpy as np

def tz_convert_from_utc(
    vals: np.ndarray,  # const int64_t[:]
    tz: tzinfo,
) -> np.ndarray: ...  # np.ndarray[np.int64]

def tz_convert_from_utc_single(val: np.int64, tz: tzinfo) -> np.int64: ...

def tz_localize_to_utc(
    vals: np.ndarray,  # np.ndarray[np.int64]
    tz: Optional[tzinfo],
    ambiguous: Optional[Union[str, bool, Iterable[bool]]] = None,
    nonexistent: Optional[Union[str, timedelta, np.timedelta64]] = None,
) -> np.ndarray: ...  # np.ndarray[np.int64]
