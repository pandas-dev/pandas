from typing import TypeAlias, overload

import numpy as np
import optype as op
import optype.numpy as onp

__all__ = ["orthogonal_procrustes"]

_Float: TypeAlias = np.float32 | np.float64
_Complex: TypeAlias = np.complex64 | np.complex128

###

@overload
def orthogonal_procrustes(
    A: onp.ToFloat64, B: onp.ToIntND | onp.ToJustFloat64, check_finite: op.CanBool = True
) -> tuple[onp.ArrayND[np.float64], np.float64]: ...
@overload
def orthogonal_procrustes(
    A: onp.ToIntND | onp.ToJustFloat64, B: onp.ToFloat64, check_finite: op.CanBool = True
) -> tuple[onp.ArrayND[np.float64], np.float64]: ...
@overload
def orthogonal_procrustes(
    A: onp.ToFloatND, B: onp.ToFloatND, check_finite: op.CanBool = True
) -> tuple[onp.ArrayND[_Float], _Float]: ...
@overload
def orthogonal_procrustes(
    A: onp.ToComplex128_ND, B: onp.ToJustComplex128_ND, check_finite: onp.ToBool = True
) -> tuple[onp.ArrayND[np.complex128], np.float64]: ...
@overload
def orthogonal_procrustes(
    A: onp.ToJustComplex128_ND, B: onp.ToComplex128_ND, check_finite: onp.ToBool = True
) -> tuple[onp.ArrayND[np.complex128], np.float64]: ...
@overload
def orthogonal_procrustes(
    A: onp.ToJustComplexND, B: onp.ToComplexND, check_finite: onp.ToBool = True
) -> tuple[onp.ArrayND[_Complex], _Float]: ...
@overload
def orthogonal_procrustes(
    A: onp.ToComplexND, B: onp.ToJustComplexND, check_finite: onp.ToBool = True
) -> tuple[onp.ArrayND[_Complex], _Float]: ...
@overload
def orthogonal_procrustes(
    A: onp.ToComplexND, B: onp.ToComplexND, check_finite: onp.ToBool = True
) -> tuple[onp.ArrayND[_Float | _Complex], _Float]: ...
