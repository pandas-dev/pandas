from typing import Optional

import numpy as np

def array_strptime(
    values: np.ndarray,  # np.ndarray[object]
    fmt: Optional[str],
    exact: bool = True,
    errors: str = "raise"
) -> tuple[np.ndarray, np.ndarray]: ...
# first  ndarray is M8[ns], second is object ndarray of Optional[tzinfo]
