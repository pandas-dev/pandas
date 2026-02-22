import numpy as np
import optype.numpy as onp

__all__: list[str] = []

class SqrtmError(np.linalg.LinAlgError): ...  # undocumented

def _sqrtm_triu(T: onp.ToComplex2D, blocksize: onp.ToJustInt = 64) -> onp.Array2D[np.float64 | np.complex128]: ...  # undocumented
