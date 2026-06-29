import numpy as np

from pandas._typing import npt

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
def _get_dummies_dense_from_codes(
    dummy_mat: np.ndarray,
    codes: npt.NDArray[np.intp],
    n_rows: int,
) -> None: ...
