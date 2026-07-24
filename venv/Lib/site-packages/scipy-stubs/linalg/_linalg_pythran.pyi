from typing import Final

import numpy as np
import optype.numpy as onp

###

__pythran__: Final[tuple[str, str]] = ...

def _funm_loops[InexactT: np.float32 | np.float64 | np.complex64 | np.complex128](
    F: onp.Array2D[InexactT], T: onp.Array2D[InexactT], n: int, minden: InexactT
) -> tuple[onp.Array2D[InexactT], InexactT]: ...
