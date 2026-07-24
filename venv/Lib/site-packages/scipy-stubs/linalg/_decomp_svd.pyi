from typing import Literal, SupportsIndex, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["diagsvd", "null_space", "orth", "subspace_angles", "svd", "svdvals"]

###

type _SVD_ND[ScalarT1: np.generic, ScalarT2: np.generic] = tuple[
    onp.ArrayND[ScalarT1],
    onp.ArrayND[ScalarT2],
    onp.ArrayND[ScalarT1],
]  # fmt: skip

type _LapackDriver = Literal["gesdd", "gesvd"]

###

@overload  # nd float64
def svd(
    a: onp.ToArrayND[float, npc.floating80 | npc.floating64 | npc.integer | np.bool],
    full_matrices: bool = True,
    compute_uv: Literal[True] = True,
    overwrite_a: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver = "gesdd",
) -> _SVD_ND[np.float64, np.float64]: ...
@overload  # nd float32
def svd(
    a: onp.CanArrayND[np.float32 | np.float16],
    full_matrices: bool = True,
    compute_uv: Literal[True] = True,
    overwrite_a: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver = "gesdd",
) -> _SVD_ND[np.float32, np.float32]: ...
@overload  # nd complex128
def svd(
    a: onp.ToArrayND[op.JustComplex, npc.complexfloating128 | npc.complexfloating160],
    full_matrices: bool = True,
    compute_uv: Literal[True] = True,
    overwrite_a: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver = "gesdd",
) -> _SVD_ND[np.complex128, np.float64]: ...
@overload  # nd complex64
def svd(
    a: onp.CanArrayND[np.complex64],
    full_matrices: bool = True,
    compute_uv: Literal[True] = True,
    overwrite_a: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver = "gesdd",
) -> _SVD_ND[np.complex64, np.float32]: ...
@overload  # nd float64 | complex128, compute_uv=False (keyword)
def svd(
    a: onp.ToArrayND[complex, npc.inexact80 | npc.inexact64 | npc.integer | np.bool],
    full_matrices: bool = True,
    *,
    compute_uv: Literal[False],
    overwrite_a: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver = "gesdd",
) -> onp.ArrayND[np.float64]: ...
@overload  # nd float32 | complex64, compute_uv=False (keyword)
def svd(
    a: onp.CanArrayND[npc.inexact32 | np.float16],
    full_matrices: bool = True,
    *,
    compute_uv: Literal[False],
    overwrite_a: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver = "gesdd",
) -> onp.ArrayND[np.float32]: ...

#
@overload
def svdvals(  # type: ignore[overload-overlap]
    a: onp.ToArrayND[complex, npc.number64 | npc.inexact80 | npc.integer32], overwrite_a: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.float64]: ...
@overload
def svdvals(
    a: onp.ToArrayND[np.float32, npc.inexact32 | npc.number16 | npc.integer8 | np.bool],
    overwrite_a: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float32]: ...
@overload
def svdvals(a: onp.ToComplexND, overwrite_a: bool = False, check_finite: bool = True) -> onp.ArrayND[np.float64 | np.float32]: ...

#
@overload
def diagsvd[RealT: np.bool | npc.integer | npc.floating](
    s: onp.ToArrayND[RealT, RealT], M: SupportsIndex, N: SupportsIndex
) -> onp.ArrayND[RealT]: ...
@overload
def diagsvd(s: onp.SequenceND[bool], M: SupportsIndex, N: SupportsIndex) -> onp.ArrayND[np.bool]: ...
@overload
def diagsvd(s: onp.SequenceND[op.JustInt], M: SupportsIndex, N: SupportsIndex) -> onp.ArrayND[np.int_]: ...
@overload
def diagsvd(s: onp.SequenceND[op.JustFloat], M: SupportsIndex, N: SupportsIndex) -> onp.ArrayND[np.float64]: ...

#
@overload
def orth(
    A: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool], rcond: float | None = None
) -> onp.ArrayND[np.float64]: ...
@overload
def orth(A: onp.ToJustComplex128_ND, rcond: float | None = None) -> onp.ArrayND[np.complex128]: ...
@overload
def orth[InexactT: np.float32 | np.float64 | np.complex64 | np.complex128](
    A: onp.ToArrayND[InexactT, InexactT], rcond: float | None = None
) -> onp.ArrayND[InexactT]: ...

#
@overload
def null_space(
    A: onp.ToArrayND[float, npc.floating64 | npc.integer | np.bool],
    rcond: float | None = None,
    *,
    overwrite_a: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver = "gesdd",
) -> onp.ArrayND[np.float64]: ...
@overload
def null_space(
    A: onp.ToJustComplex128_ND,
    rcond: float | None = None,
    *,
    overwrite_a: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver = "gesdd",
) -> onp.ArrayND[np.complex128]: ...
@overload
def null_space[InexactT: np.float32 | np.float64 | np.complex64 | np.complex128](
    A: onp.ToArrayND[InexactT, InexactT],
    rcond: float | None = None,
    *,
    overwrite_a: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver = "gesdd",
) -> onp.ArrayND[InexactT]: ...

#
@overload
def subspace_angles(  # type: ignore[overload-overlap]
    A: onp.ToArray2D[complex, npc.number64 | npc.inexact80 | npc.integer32], B: onp.ToComplex2D
) -> onp.Array1D[np.float64]: ...
@overload
def subspace_angles(  # type: ignore[overload-overlap]
    A: onp.ToComplex2D, B: onp.ToArray2D[complex, npc.number64 | npc.inexact80 | npc.integer32]
) -> onp.Array1D[np.float64]: ...
@overload
def subspace_angles(
    A: onp.ToArray2D[np.float32, npc.inexact32 | npc.number16 | npc.integer8 | np.bool],
    B: onp.ToArray2D[np.float32, npc.inexact32 | npc.number16 | npc.integer8 | np.bool],
) -> onp.Array1D[np.float32]: ...
@overload
def subspace_angles(A: onp.ToComplex2D, B: onp.ToComplex2D) -> onp.Array1D[np.float64 | np.float32]: ...
