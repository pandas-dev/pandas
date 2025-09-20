from typing import TypeAlias

import numpy as np
import optype.numpy as onp

_InexactND: TypeAlias = onp.ArrayND[np.float32 | np.float64 | np.complex64 | np.complex128]

###

class error(Exception): ...  # undocumented

# `Am` must have a shape like `(5, n, n)`.
def pick_pade_structure(Am: _InexactND) -> tuple[int, int]: ...
def pade_UV_calc(Am: _InexactND, n: int) -> int: ...
