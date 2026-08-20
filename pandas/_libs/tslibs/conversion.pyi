from datetime import (
    datetime,
    tzinfo,
)

import numpy as np
import numpy.typing as npt

DT64NS_DTYPE: np.dtype
TD64NS_DTYPE: np.dtype

def localize_pydatetime(dt: datetime, tz: tzinfo | None) -> datetime: ...
def cast_from_unit_vectorized(
    values: npt.NDArray[np.float64], unit: str, out_unit: str = ...
) -> npt.NDArray[np.int64]: ...
def datetime_from_fields(
    years: npt.NDArray[np.int64],
    months: npt.NDArray[np.int64],
    days: npt.NDArray[np.int64],
    hours: npt.NDArray[np.int64],
    minutes: npt.NDArray[np.int64],
    seconds: npt.NDArray[np.int64],
) -> tuple[npt.NDArray[np.int64], int]: ...
