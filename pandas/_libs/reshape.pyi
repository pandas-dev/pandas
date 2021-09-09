import numpy as np

import pandas._tying as npt

def unstack(
    values: np.ndarray,  # reshape_t[:, :]
    mask: np.ndarray,  # const uint8_t[:]
    stride: int,
    length: int,
    width: int,
    new_values: np.ndarray,  # reshape_t[:, :]
    new_mask: np.ndarray,  # uint8_t[:, :]
) -> None: ...
def explode(
    values: npt.NDArray[np.object_],
) -> tuple[npt.NDArray[np.object_], npt.NDArray[np.int64]]: ...
