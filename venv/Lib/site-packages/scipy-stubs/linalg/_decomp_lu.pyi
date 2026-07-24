from typing import Any, Literal, Never, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["lu", "lu_factor", "lu_solve"]

###

type _as_f32 = np.float32 | np.float16 | npc.integer16 | npc.integer8 | np.bool  # noqa: PYI042
type _as_f64 = npc.floating64 | npc.floating80 | npc.integer64 | npc.integer32  # noqa: PYI042
type _as_c128 = npc.complexfloating160 | npc.complexfloating128  # noqa: PYI042

type _PivND = onp.ToArrayND[int, npc.integer]
type _Trans = Literal[0, 1, 2]

# workaround for mypy & pyright's failure to conform to the overload typing specification
type _JustAnyShape = tuple[Never, Never, Never, Never]

###

# NOTE: The ignored mypy `overload-overlap` errors are false positives

@overload  # ?d f64
def lu_factor(  # type: ignore[overload-overlap]
    a: onp.ArrayND[_as_f64, _JustAnyShape], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.int32]]: ...
@overload  # ?d f32
def lu_factor(
    a: onp.ArrayND[_as_f32, _JustAnyShape], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.ArrayND[np.float32], onp.ArrayND[np.int32]]: ...
@overload  # ?d c128
def lu_factor(
    a: onp.ArrayND[_as_c128, _JustAnyShape], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.ArrayND[np.complex128], onp.ArrayND[np.int32]]: ...
@overload  # ?d c64
def lu_factor(
    a: onp.ArrayND[np.complex64, _JustAnyShape], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.ArrayND[np.complex64], onp.ArrayND[np.int32]]: ...
@overload  # 2d f64
def lu_factor(  # type: ignore[overload-overlap]
    a: onp.ToArrayStrict2D[float, _as_f64], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.Array2D[np.float64], onp.Array1D[np.int32]]: ...
@overload  # 2d f32
def lu_factor(
    a: onp.ToArrayStrict2D[np.float32, _as_f32], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.Array2D[np.float32], onp.Array1D[np.int32]]: ...
@overload  # 2d c128
def lu_factor(
    a: onp.ToArrayStrict2D[op.JustComplex, _as_c128], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.Array2D[np.complex128], onp.Array1D[np.int32]]: ...
@overload  # 2d c64
def lu_factor(
    a: onp.ToJustComplex64Strict2D, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.Array2D[np.complex64], onp.Array1D[np.int32]]: ...
@overload  # 3d f64
def lu_factor(  # type: ignore[overload-overlap]
    a: onp.ToArrayStrict3D[float, _as_f64], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.Array3D[np.float64], onp.Array2D[np.int32]]: ...
@overload  # 3d f32
def lu_factor(
    a: onp.ToArrayStrict3D[np.float32, _as_f32], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.Array3D[np.float32], onp.Array2D[np.int32]]: ...
@overload  # 3d c128
def lu_factor(
    a: onp.ToArrayStrict3D[op.JustComplex, _as_c128], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.Array3D[np.complex128], onp.Array2D[np.int32]]: ...
@overload  # 3d c64
def lu_factor(
    a: onp.ToJustComplex64Strict3D, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.Array3D[np.complex64], onp.Array2D[np.int32]]: ...
@overload  # nd f64
def lu_factor(  # type: ignore[overload-overlap]
    a: onp.ToArrayND[float, _as_f64], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.int32]]: ...
@overload  # nd f32
def lu_factor(
    a: onp.ToArrayND[np.float32, _as_f32], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.ArrayND[np.float32], onp.ArrayND[np.int32]]: ...
@overload  # nd c128
def lu_factor(
    a: onp.ToArrayND[op.JustComplex, _as_c128], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.ArrayND[np.complex128], onp.ArrayND[np.int32]]: ...
@overload  # nd c64
def lu_factor(
    a: onp.ToJustComplex64_ND, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.ArrayND[np.complex64], onp.ArrayND[np.int32]]: ...

# keep in sync with `cho_solve` in `_decomp_cholesky`
@overload  # nd +f64\+f32, ?d +f64
def lu_solve(  # type: ignore[overload-overlap]
    lu_and_piv: tuple[onp.ToArrayND[float, _as_f64], _PivND],
    b: onp.ToFloatND,
    trans: _Trans = 0,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # nd +f64, ?d +f64\+f32
def lu_solve(  # type: ignore[overload-overlap]
    lu_and_piv: tuple[onp.ToFloatND, _PivND],
    b: onp.ToArrayND[float, _as_f64],
    trans: _Trans = 0,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # nd +f32, ?d +f32
def lu_solve(
    lu_and_piv: tuple[onp.ToFloat32_ND, _PivND],
    b: onp.ToFloat32_ND,
    trans: _Trans = 0,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float32]: ...
@overload  # nd ~c128|c160, ?d +c128
def lu_solve(
    lu_and_piv: tuple[onp.ToArrayND[op.JustComplex, _as_c128], _PivND],
    b: onp.ToComplexND,
    trans: _Trans = 0,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload  # nd +c128, ?d ~c128|c160
def lu_solve(
    lu_and_piv: tuple[onp.ToComplexND, _PivND],
    b: onp.ToArrayND[op.JustComplex, _as_c128],
    trans: _Trans = 0,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload  # nd ~c64, ?d +c64
def lu_solve(
    lu_and_piv: tuple[onp.ToJustComplex64_ND, _PivND],
    b: onp.ToComplex64_ND,
    trans: _Trans = 0,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex64]: ...
@overload  # nd +c64, ?d ~c64
def lu_solve(
    lu_and_piv: tuple[onp.ToComplex64_ND, _PivND],
    b: onp.ToJustComplex64_ND,
    trans: _Trans = 0,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex64]: ...
@overload  # nd +cfloating, ?d ~cfloating (fallback)
def lu_solve(
    lu_and_piv: tuple[onp.ToComplexND, _PivND],
    b: onp.ToComplexND,
    trans: _Trans = 0,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[Any]: ...

#
@overload  # nd f64
def lu(  # type: ignore[overload-overlap]
    a: onp.ToArrayND[float, _as_f64],
    permute_l: Literal[False] = False,
    overwrite_a: bool = False,
    check_finite: bool = True,
    p_indices: Literal[False] = False,
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
@overload  # nd f64, p_indices=True
def lu(  # type: ignore[overload-overlap]
    a: onp.ToArrayND[float, _as_f64],
    permute_l: Literal[False] = False,
    overwrite_a: bool = False,
    check_finite: bool = True,
    *,
    p_indices: Literal[True],
) -> tuple[onp.ArrayND[np.int32], onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
@overload  # nd f64, permute_l=True
def lu(  # type: ignore[overload-overlap]
    a: onp.ToArrayND[float, _as_f64],
    permute_l: Literal[True],
    overwrite_a: bool = False,
    check_finite: bool = True,
    p_indices: bool = False,
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
@overload  # nd f32
def lu(
    a: onp.ToArrayND[np.float32, _as_f32],
    permute_l: Literal[False] = False,
    overwrite_a: bool = False,
    check_finite: bool = True,
    p_indices: Literal[False] = False,
) -> tuple[onp.ArrayND[np.float32], onp.ArrayND[np.float32], onp.ArrayND[np.float32]]: ...
@overload  # nd f32, p_indices=True
def lu(
    a: onp.ToArrayND[np.float32, _as_f32],
    permute_l: Literal[False] = False,
    overwrite_a: bool = False,
    check_finite: bool = True,
    *,
    p_indices: Literal[True],
) -> tuple[onp.ArrayND[np.int32], onp.ArrayND[np.float32], onp.ArrayND[np.float32]]: ...
@overload  # nd f32, permute_l=True
def lu(
    a: onp.ToArrayND[np.float32, _as_f32],
    permute_l: Literal[True],
    overwrite_a: bool = False,
    check_finite: bool = True,
    p_indices: bool = False,
) -> tuple[onp.ArrayND[np.float32], onp.ArrayND[np.float32]]: ...
@overload  # nd c128
def lu(
    a: onp.ToArrayND[op.JustComplex, _as_c128],
    permute_l: Literal[False] = False,
    overwrite_a: bool = False,
    check_finite: bool = True,
    p_indices: Literal[False] = False,
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.complex128], onp.ArrayND[np.complex128]]: ...
@overload  # nd c128, p_indices=True
def lu(
    a: onp.ToArrayND[op.JustComplex, _as_c128],
    permute_l: Literal[False] = False,
    overwrite_a: bool = False,
    check_finite: bool = True,
    *,
    p_indices: Literal[True],
) -> tuple[onp.ArrayND[np.int32], onp.ArrayND[np.complex128], onp.ArrayND[np.complex128]]: ...
@overload  # nd c128, permute_l=True
def lu(
    a: onp.ToArrayND[op.JustComplex, _as_c128],
    permute_l: Literal[True],
    overwrite_a: bool = False,
    check_finite: bool = True,
    p_indices: bool = False,
) -> tuple[onp.ArrayND[np.complex128], onp.ArrayND[np.complex128]]: ...
@overload  # nd c64
def lu(
    a: onp.ToJustComplex64_ND,
    permute_l: Literal[False] = False,
    overwrite_a: bool = False,
    check_finite: bool = True,
    p_indices: Literal[False] = False,
) -> tuple[onp.ArrayND[np.float32], onp.ArrayND[np.complex64], onp.ArrayND[np.complex64]]: ...
@overload  # nd c64, p_indices=True
def lu(
    a: onp.ToJustComplex64_ND,
    permute_l: Literal[False] = False,
    overwrite_a: bool = False,
    check_finite: bool = True,
    *,
    p_indices: Literal[True],
) -> tuple[onp.ArrayND[np.int32], onp.ArrayND[np.complex64], onp.ArrayND[np.complex64]]: ...
@overload  # nd c64, permute_l=True
def lu(
    a: onp.ToJustComplex64_ND,
    permute_l: Literal[True],
    overwrite_a: bool = False,
    check_finite: bool = True,
    p_indices: bool = False,
) -> tuple[onp.ArrayND[np.complex64], onp.ArrayND[np.complex64]]: ...
