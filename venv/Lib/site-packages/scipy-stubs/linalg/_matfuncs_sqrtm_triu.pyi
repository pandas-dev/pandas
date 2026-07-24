from collections.abc import Iterable
from typing import overload

import numpy as np
import optype.numpy as onp

@overload
def within_block_loop(
    R: onp.ArrayND[np.float64], T: onp.ArrayND[np.float64], start_stop_pairs: Iterable[tuple[int, int]], nblocks: int | np.intp
) -> None: ...
@overload
def within_block_loop(
    R: onp.ArrayND[np.complex128],
    T: onp.ArrayND[np.complex128],
    start_stop_pairs: Iterable[tuple[int, int]],
    nblocks: int | np.intp,
) -> None: ...
