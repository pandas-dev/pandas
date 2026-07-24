from typing import overload

import numpy as np
import optype.numpy as onp

@overload
def lambertw(z: onp.ToComplex, k: onp.ToInt = 0, tol: float | np.float64 = 1e-8) -> np.complex128: ...
@overload
def lambertw(z: onp.ToComplex, k: onp.ToIntND, tol: float | np.float64 = 1e-8) -> onp.ArrayND[np.complex128]: ...
@overload
def lambertw(z: onp.ToComplexND, k: onp.ToInt | onp.ToIntND, tol: float | np.float64 = 1e-8) -> onp.ArrayND[np.complex128]: ...
