from typing import overload

import numpy as np
import optype.numpy as onp

__all__ = ["orthogonal_procrustes"]

###

type _Float = np.float32 | np.float64
type _Complex = np.complex64 | np.complex128

###

@overload
def orthogonal_procrustes(
    A: onp.ToFloat64, B: onp.ToIntND | onp.ToJustFloat64, check_finite: bool = True
) -> tuple[onp.ArrayND[np.float64], np.float64]: ...
@overload
def orthogonal_procrustes(
    A: onp.ToIntND | onp.ToJustFloat64, B: onp.ToFloat64, check_finite: bool = True
) -> tuple[onp.ArrayND[np.float64], np.float64]: ...
@overload
def orthogonal_procrustes(
    A: onp.ToFloatND, B: onp.ToFloatND, check_finite: bool = True
) -> tuple[onp.ArrayND[_Float], _Float]: ...
@overload
def orthogonal_procrustes(
    A: onp.ToComplex128_ND, B: onp.ToJustComplex128_ND, check_finite: bool = True
) -> tuple[onp.ArrayND[np.complex128], np.float64]: ...
@overload
def orthogonal_procrustes(
    A: onp.ToJustComplex128_ND, B: onp.ToComplex128_ND, check_finite: bool = True
) -> tuple[onp.ArrayND[np.complex128], np.float64]: ...
@overload
def orthogonal_procrustes(
    A: onp.ToJustComplexND, B: onp.ToComplexND, check_finite: bool = True
) -> tuple[onp.ArrayND[_Complex], _Float]: ...
@overload
def orthogonal_procrustes(
    A: onp.ToComplexND, B: onp.ToJustComplexND, check_finite: bool = True
) -> tuple[onp.ArrayND[_Complex], _Float]: ...
@overload
def orthogonal_procrustes(
    A: onp.ToComplexND, B: onp.ToComplexND, check_finite: bool = True
) -> tuple[onp.ArrayND[_Float | _Complex], _Float]: ...
