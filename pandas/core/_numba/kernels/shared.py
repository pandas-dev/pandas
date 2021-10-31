import numba
import numpy as np


# error: Untyped decorator makes function "is_monotonic_increasing" untyped
@numba.jit(  # type: ignore[misc]
    numba.boolean(numba.int64[:]), nopython=True, nogil=True, parallel=False
)
def is_monotonic_increasing(bounds: np.ndarray) -> bool:
    """Check if int64 values are monotonically increasing."""
    n = len(bounds)
    if n < 2:
        return True
    prev = bounds[0]
    for i in range(1, n):
        cur = bounds[i]
        if cur < prev:
            return False
        prev = cur
    return True
