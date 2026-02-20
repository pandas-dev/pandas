from typing import Literal, TypeAlias, overload

import numpy as np
import optype.numpy as onp

__all__ = ["expm_cond", "expm_frechet"]

_Method: TypeAlias = Literal["SPS", "blockEnlarge"]
_ArrayF64: TypeAlias = onp.ArrayND[np.float64]
_ArrayC128: TypeAlias = onp.ArrayND[np.complex128]

###

@overload
def expm_frechet(
    A: onp.ToComplexND,
    E: onp.ToComplexND,
    method: _Method | None = None,
    compute_expm: onp.ToTrue = True,
    check_finite: bool = True,
) -> tuple[_ArrayF64, _ArrayF64] | tuple[_ArrayF64 | _ArrayC128, _ArrayC128]: ...
@overload
def expm_frechet(
    A: onp.ToComplexND, E: onp.ToComplexND, method: _Method | None, compute_expm: onp.ToFalse, check_finite: bool = True
) -> tuple[_ArrayF64, _ArrayF64] | tuple[_ArrayF64 | _ArrayC128, _ArrayC128]: ...
@overload
def expm_frechet(
    A: onp.ToComplexND, E: onp.ToComplexND, method: _Method | None = None, *, compute_expm: onp.ToFalse, check_finite: bool = True
) -> tuple[_ArrayF64, _ArrayF64] | tuple[_ArrayF64 | _ArrayC128, _ArrayC128]: ...

#
@overload
def expm_cond(A: onp.ToComplexStrict2D, check_finite: bool = True) -> np.float64: ...
@overload
def expm_cond(A: onp.ToComplexStrict3D, check_finite: bool = True) -> onp.Array1D[np.float64]: ...
@overload
def expm_cond(A: onp.ToComplexND, check_finite: bool = True) -> np.float64 | onp.ArrayND[np.float64]: ...
