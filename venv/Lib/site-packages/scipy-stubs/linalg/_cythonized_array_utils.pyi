import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["ishermitian", "issymmetric"]

###

# see `scipy/linalg/_cythonized_array_utils.pxd`
type _Numeric = npc.integer | np.float32 | np.float64 | npc.floating80 | np.complex64 | np.complex128

###

def issymmetric(a: onp.ArrayND[_Numeric], atol: float | None = None, rtol: float | None = None) -> bool: ...
def ishermitian(a: onp.ArrayND[_Numeric], atol: float | None = None, rtol: float | None = None) -> bool: ...
