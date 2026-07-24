from typing import overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["ldl"]

type _ISize1D = onp.Array1D[np.intp]
type _ISizeND = onp.ArrayND[np.intp]
type _Float2D = onp.Array2D[npc.floating]
type _FloatND = onp.ArrayND[npc.floating]
type _Complex2D = onp.Array2D[npc.complexfloating]
type _ComplexND = onp.ArrayND[npc.complexfloating]
type _InexactND = onp.ArrayND[npc.inexact]

###

@overload
def ldl(
    A: onp.ToFloatStrict2D, lower: bool = True, hermitian: bool = True, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[_Float2D, _Float2D, _ISize1D]: ...
@overload
def ldl(
    A: onp.ToFloatND, lower: bool = True, hermitian: bool = True, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[_FloatND, _FloatND, _ISizeND]: ...
@overload
def ldl(
    A: onp.ToJustComplexStrict2D, lower: bool = True, hermitian: bool = True, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[_Complex2D, _Complex2D, _ISize1D]: ...
@overload
def ldl(
    A: onp.ToJustComplexND, lower: bool = True, hermitian: bool = True, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[_ComplexND, _ComplexND, _ISizeND]: ...
@overload
def ldl(
    A: onp.ToComplexND, lower: bool = True, hermitian: bool = True, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[_InexactND, _InexactND, _ISizeND]: ...
