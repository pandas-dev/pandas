from typing import Final, TypeVar

import numpy as np
import optype.numpy as onp

_InexactT = TypeVar("_InexactT", bound=np.float32 | np.float64 | np.complex64 | np.complex128)

###

__pythran__: Final[tuple[str, str]] = ...

def _funm_loops(
    F: onp.Array2D[_InexactT], T: onp.Array2D[_InexactT], n: int, minden: _InexactT
) -> tuple[onp.Array2D[_InexactT], _InexactT]: ...
