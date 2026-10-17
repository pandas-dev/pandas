from collections.abc import Iterable
from datetime import (
    timedelta,
    tzinfo,
)

import numpy as np

from pandas._libs.tslibs.fields import RoundTo
from pandas._typing import npt

# tz_convert_from_utc_single exposed for testing
def tz_convert_from_utc_single(
    utc_val: np.int64, tz: tzinfo, creso: int = ...
) -> np.int64: ...
def tz_localize_to_utc(
    vals: npt.NDArray[np.int64],
    tz: tzinfo | None,
    ambiguous: str | bool | Iterable[bool] | None = ...,
    nonexistent: str | timedelta | np.timedelta64 | None = ...,
    creso: int = ...,  # NPY_DATETIMEUNIT
) -> npt.NDArray[np.int64]: ...
def tz_localize_rounded(
    rounded: npt.NDArray[np.int64],
    orig: npt.NDArray[np.int64],
    tz: tzinfo | None,
    mode: RoundTo,
    nonexistent: str | timedelta | None,
    creso: int,  # NPY_DATETIMEUNIT
) -> npt.NDArray[np.int64]: ...
