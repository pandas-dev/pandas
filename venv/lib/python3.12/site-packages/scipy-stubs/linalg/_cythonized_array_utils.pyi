from typing import TypeAlias

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["bandwidth", "ishermitian", "issymmetric"]

# see `scipy/linalg/_cythonized_array_utils.pxd`
_Numeric: TypeAlias = npc.integer | np.float32 | np.float64 | npc.floating80 | np.complex64 | np.complex128

def bandwidth(a: onp.ArrayND[_Numeric]) -> tuple[int, int]: ...
def issymmetric(a: onp.ArrayND[_Numeric], atol: float | None = None, rtol: float | None = None) -> bool: ...
def ishermitian(a: onp.ArrayND[_Numeric], atol: float | None = None, rtol: float | None = None) -> bool: ...
