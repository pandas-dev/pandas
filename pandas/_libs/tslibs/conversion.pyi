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
