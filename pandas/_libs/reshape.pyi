import numpy as np

from pandas._typing import npt

def unstack(
    values: np.ndarray,  # numeric_object_complex_t[:, :]
    mask: np.ndarray,  # const uint8_t[:]
    stride: int,
    length: int,
    width: int,
    new_values: np.ndarray,  # numeric_object_complex_t[:, :]
) -> None: ...
def explode(
    values: npt.NDArray[np.object_],
) -> tuple[npt.NDArray[np.object_], npt.NDArray[np.int64]]: ...
