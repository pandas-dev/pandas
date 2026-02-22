from typing import Literal, TypeAlias, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["lu", "lu_factor", "lu_solve"]

_ISizeND: TypeAlias = onp.ArrayND[np.intp]

_Float1D: TypeAlias = onp.Array1D[npc.floating]
_Float2D: TypeAlias = onp.Array2D[npc.floating]
_FloatND: TypeAlias = onp.ArrayND[npc.floating]

_Complex1D: TypeAlias = onp.Array1D[npc.complexfloating]
_Complex2D: TypeAlias = onp.Array2D[npc.complexfloating]
_ComplexND: TypeAlias = onp.ArrayND[npc.complexfloating]

_InexactND: TypeAlias = onp.ArrayND[npc.inexact]

_Trans: TypeAlias = Literal[0, 1, 2]

@overload  # float[:, :] -> (float[:, :], float[:])
def lu_factor(
    a: onp.ToFloatStrict2D, overwrite_a: onp.ToBool = False, check_finite: onp.ToBool = True
) -> tuple[_Float2D, _Float1D]: ...
@overload  # float[:, :, ...] -> (float[:, :, ...], float[:, ...])
def lu_factor(
    a: onp.ToFloatND, overwrite_a: onp.ToBool = False, check_finite: onp.ToBool = True
) -> tuple[_FloatND, _FloatND]: ...
@overload  # complex[:, :] -> (complex[:, :], complex[:])
def lu_factor(
    a: onp.ToJustComplexStrict2D, overwrite_a: onp.ToBool = False, check_finite: onp.ToBool = True
) -> tuple[_Complex2D, _Complex1D]: ...
@overload  # complex[:, :, ...] -> (complex[:, :, ...], complex[:, ...])
def lu_factor(
    a: onp.ToJustComplexND, overwrite_a: onp.ToBool = False, check_finite: onp.ToBool = True
) -> tuple[_ComplexND, _ComplexND]: ...
@overload  # fallback
def lu_factor(
    a: onp.ToComplexND, overwrite_a: onp.ToBool = False, check_finite: onp.ToBool = True
) -> tuple[_InexactND, _InexactND]: ...

#
@overload  # (float[:, :], float[:]) -> float[:, :]
def lu_solve(
    lu_and_piv: tuple[onp.ToFloatStrict2D, onp.ToFloatStrict1D],
    b: onp.ToFloat1D,
    trans: _Trans = 0,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _Float2D: ...
@overload  # (float[:, :, ...], float[:, ...]) -> float[:, :, ...]
def lu_solve(
    lu_and_piv: tuple[onp.ToFloatND, onp.ToFloatND],
    b: onp.ToFloatND,
    trans: _Trans = 0,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _FloatND: ...
@overload  # (complex[:, :,], complex[:]) -> complex[:, :]
def lu_solve(
    lu_and_piv: tuple[onp.ToJustComplexStrict2D, onp.ToJustComplexStrict1D],
    b: onp.ToJustComplex1D,
    trans: _Trans = 0,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _Complex2D: ...
@overload  # (complex[:, :, ...], complex[:, ...]) -> complex[:, :, ...]
def lu_solve(
    lu_and_piv: tuple[onp.ToJustComplexND, onp.ToJustComplexND],
    b: onp.ToJustComplexND,
    trans: _Trans = 0,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _ComplexND: ...
@overload  # fallback
def lu_solve(
    lu_and_piv: tuple[onp.ToComplexND, onp.ToComplexND],
    b: onp.ToComplexND,
    trans: _Trans = 0,
    overwrite_b: onp.ToBool = False,
    check_finite: onp.ToBool = True,
) -> _InexactND: ...

#
@overload  # (float[:, :], permute_l=False, p_indices=False) -> (float[...], float[...], float[...])
def lu(
    a: onp.ToFloatND,
    permute_l: onp.ToFalse = False,
    overwrite_a: onp.ToBool = False,
    check_finite: onp.ToBool = True,
    p_indices: onp.ToFalse = False,
) -> tuple[_FloatND, _FloatND, _FloatND]: ...
@overload  # (float[:, :], permute_l=False, p_indices=True, /) -> (intp[...], float[...], float[...])
def lu(
    a: onp.ToFloatND, permute_l: onp.ToFalse, overwrite_a: onp.ToBool, check_finite: onp.ToBool, p_indices: onp.ToTrue
) -> tuple[_ISizeND, _FloatND, _FloatND]: ...
@overload  # (float[:, :], permute_l=False, *, p_indices=True) -> (intp[...], float[...], float[...])
def lu(
    a: onp.ToFloatND,
    permute_l: onp.ToFalse = False,
    overwrite_a: onp.ToBool = False,
    check_finite: onp.ToBool = True,
    *,
    p_indices: onp.ToTrue,
) -> tuple[_ISizeND, _FloatND, _FloatND]: ...
@overload  # (float[:, :], permute_l=True, p_indices=False) -> (intp[...], float[...], float[...])
def lu(
    a: onp.ToFloatND,
    permute_l: onp.ToTrue,
    overwrite_a: onp.ToBool = False,
    check_finite: onp.ToBool = True,
    p_indices: onp.ToBool = False,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # (complex[:, :], permute_l=False, p_indices=False) -> (complex[...], complex[...], complex[...])
def lu(
    a: onp.ToJustComplexND,
    permute_l: onp.ToFalse = False,
    overwrite_a: onp.ToBool = False,
    check_finite: onp.ToBool = True,
    p_indices: onp.ToFalse = False,
) -> tuple[_ComplexND, _ComplexND, _ComplexND]: ...
@overload  # (complex[:, :], permute_l=False, p_indices=True, /) -> (intp[...], complex[...], complex[...])
def lu(
    a: onp.ToJustComplexND, permute_l: onp.ToFalse, overwrite_a: onp.ToBool, check_finite: onp.ToBool, p_indices: onp.ToTrue
) -> tuple[_ISizeND, _ComplexND, _ComplexND]: ...
@overload  # (complex[:, :], permute_l=False, *, p_indices=True) -> (intp[...], complex[...], complex[...])
def lu(
    a: onp.ToJustComplexND,
    permute_l: onp.ToFalse = False,
    overwrite_a: onp.ToBool = False,
    check_finite: onp.ToBool = True,
    *,
    p_indices: onp.ToTrue,
) -> tuple[_ISizeND, _ComplexND, _ComplexND]: ...
@overload  # (complex[:, :], permute_l=True, p_indices=False) -> (intp[...], complex[...], complex[...])
def lu(
    a: onp.ToJustComplexND,
    permute_l: onp.ToTrue,
    overwrite_a: onp.ToBool = False,
    check_finite: onp.ToBool = True,
    p_indices: onp.ToBool = False,
) -> tuple[_ComplexND, _ComplexND]: ...
@overload  # fallback, permute_l=False, p_indices=False
def lu(
    a: onp.ToComplexND,
    permute_l: onp.ToFalse = False,
    overwrite_a: onp.ToBool = False,
    check_finite: onp.ToBool = True,
    p_indices: onp.ToFalse = False,
) -> tuple[_InexactND, _InexactND, _InexactND]: ...
@overload  # fallback, permute_l=False, p_indices=True
def lu(
    a: onp.ToComplexND, permute_l: onp.ToFalse, overwrite_a: onp.ToBool, check_finite: onp.ToBool, p_indices: onp.ToTrue
) -> tuple[_ISizeND, _InexactND, _InexactND]: ...
@overload  # fallback, permute_l=False, *, p_indices=True
def lu(
    a: onp.ToComplexND,
    permute_l: onp.ToFalse = False,
    overwrite_a: onp.ToBool = False,
    check_finite: onp.ToBool = True,
    *,
    p_indices: onp.ToTrue,
) -> tuple[_ISizeND, _InexactND, _InexactND]: ...
@overload  # fallback, permute_l=True, p_indices=False
def lu(
    a: onp.ToComplexND,
    permute_l: onp.ToTrue,
    overwrite_a: onp.ToBool = False,
    check_finite: onp.ToBool = True,
    p_indices: onp.ToBool = False,
) -> tuple[_InexactND, _InexactND]: ...
