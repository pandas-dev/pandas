from datetime import tzinfo

import numpy as np

from pandas._typing import npt

def array_strptime(
    values: npt.NDArray[np.object_],
    fmt: str | None,
    exact: bool = ...,
    errors: str = ...,
    utc: bool = ...,
    creso: int = ...,  # NPY_DATETIMEUNIT
) -> tuple[npt.NDArray[np.datetime64], tzinfo | None]: ...
