from typing import Final

import numpy as np
import optype.numpy as onp

__all__ = ["max_len_seq"]

_mls_taps: Final[dict[int, list[int]]] = ...

def max_len_seq(
    nbits: onp.ToJustInt,
    state: onp.ToInt1D | None = None,
    length: onp.ToJustInt | None = None,
    taps: onp.ToJustInt1D | None = None,
) -> tuple[onp.Array1D[np.int8], onp.Array1D[np.int8]]: ...
