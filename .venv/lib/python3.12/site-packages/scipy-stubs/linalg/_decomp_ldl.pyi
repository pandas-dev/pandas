from typing import TypeAlias, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["ldl"]

_ISize1D: TypeAlias = onp.Array1D[np.intp]
_ISizeND: TypeAlias = onp.ArrayND[np.intp]
_Float2D: TypeAlias = onp.Array2D[npc.floating]
_FloatND: TypeAlias = onp.ArrayND[npc.floating]
_Complex2D: TypeAlias = onp.Array2D[npc.complexfloating]
_ComplexND: TypeAlias = onp.ArrayND[npc.complexfloating]
_InexactND: TypeAlias = onp.ArrayND[npc.inexact]

###

@overload
def ldl(
    A: onp.ToFloatStrict2D,
    lower: onp.ToBool = True,
    hermitian: onp.ToBool = True,
    overwrite_a: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> tuple[_Float2D, _Float2D, _ISize1D]: ...
@overload
def ldl(
    A: onp.ToFloatND,
    lower: onp.ToBool = True,
    hermitian: onp.ToBool = True,
    overwrite_a: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> tuple[_FloatND, _FloatND, _ISizeND]: ...
@overload
def ldl(
    A: onp.ToJustComplexStrict2D,
    lower: onp.ToBool = True,
    hermitian: onp.ToBool = True,
    overwrite_a: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> tuple[_Complex2D, _Complex2D, _ISize1D]: ...
@overload
def ldl(
    A: onp.ToJustComplexND,
    lower: onp.ToBool = True,
    hermitian: onp.ToBool = True,
    overwrite_a: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> tuple[_ComplexND, _ComplexND, _ISizeND]: ...
@overload
def ldl(
    A: onp.ToComplexND,
    lower: onp.ToBool = True,
    hermitian: onp.ToBool = True,
    overwrite_a: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> tuple[_InexactND, _InexactND, _ISizeND]: ...
