from typing import Literal, TypeAlias, TypeVar, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["diagsvd", "null_space", "orth", "subspace_angles", "svd", "svdvals"]

_RealT = TypeVar("_RealT", bound=np.bool_ | npc.integer | npc.floating)
_InexactT = TypeVar("_InexactT", bound=_Float | _Complex)
_ScalarT = TypeVar("_ScalarT", bound=np.generic)
_ScalarT1 = TypeVar("_ScalarT1", bound=np.generic)

_SVD_ND: TypeAlias = tuple[onp.ArrayND[_ScalarT], onp.ArrayND[_ScalarT1], onp.ArrayND[_ScalarT]]

_Float: TypeAlias = np.float64 | np.float32
_Complex: TypeAlias = np.complex128 | np.complex64

_as_f32: TypeAlias = np.float32 | np.float16  # noqa: PYI042
_as_f64: TypeAlias = npc.floating80 | np.float64 | npc.integer | np.bool_  # noqa: PYI042
_as_c128: TypeAlias = np.complex128 | npc.complexfloating160  # noqa: PYI042

_ToSafeFloat64ND: TypeAlias = onp.ToArrayND[float, np.float64 | npc.integer | np.bool_]
_ToArrayND: TypeAlias = onp.CanArrayND[_ScalarT] | onp.SequenceND[_ScalarT]

_LapackDriver: TypeAlias = Literal["gesdd", "gesvd"]

###

@overload  # nd float64
def svd(
    a: onp.ToArrayND[float, _as_f64],
    full_matrices: bool = True,
    compute_uv: Literal[True] = True,
    overwrite_a: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver = "gesdd",
) -> _SVD_ND[np.float64, np.float64]: ...
@overload  # nd float32
def svd(
    a: onp.CanArrayND[_as_f32],
    full_matrices: bool = True,
    compute_uv: Literal[True] = True,
    overwrite_a: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver = "gesdd",
) -> _SVD_ND[np.float32, np.float32]: ...
@overload  # nd complex128
def svd(
    a: onp.ToArrayND[op.JustComplex, _as_c128],
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
    a: onp.ToArrayND[complex, _as_f64 | _as_c128],
    full_matrices: bool = True,
    *,
    compute_uv: Literal[False],
    overwrite_a: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver = "gesdd",
) -> onp.ArrayND[np.float64]: ...
@overload  # nd float32 | complex64, compute_uv=False (keyword)
def svd(
    a: onp.CanArrayND[_as_f32 | np.complex64],
    full_matrices: bool = True,
    *,
    compute_uv: Literal[False],
    overwrite_a: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver = "gesdd",
) -> onp.ArrayND[np.float32]: ...

#
def svdvals(a: onp.ToComplexND, overwrite_a: bool = False, check_finite: bool = True) -> onp.ArrayND[np.float64 | np.float32]: ...

#
@overload
def diagsvd(s: _ToArrayND[_RealT], M: op.CanIndex, N: op.CanIndex) -> onp.ArrayND[_RealT]: ...
@overload
def diagsvd(s: onp.SequenceND[bool], M: op.CanIndex, N: op.CanIndex) -> onp.ArrayND[np.bool_]: ...
@overload
def diagsvd(s: onp.SequenceND[op.JustInt], M: op.CanIndex, N: op.CanIndex) -> onp.ArrayND[np.int_]: ...
@overload
def diagsvd(s: onp.SequenceND[op.JustFloat], M: op.CanIndex, N: op.CanIndex) -> onp.ArrayND[np.float64]: ...

#
@overload
def orth(A: _ToSafeFloat64ND, rcond: float | None = None) -> onp.ArrayND[np.float64]: ...
@overload
def orth(A: onp.ToJustComplex128_ND, rcond: float | None = None) -> onp.ArrayND[np.complex128]: ...
@overload
def orth(A: _ToArrayND[_InexactT], rcond: float | None = None) -> onp.ArrayND[_InexactT]: ...

#
@overload
def null_space(
    A: _ToSafeFloat64ND,
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
def null_space(
    A: _ToArrayND[_InexactT],
    rcond: float | None = None,
    *,
    overwrite_a: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver = "gesdd",
) -> onp.ArrayND[_InexactT]: ...

#
def subspace_angles(A: onp.ToComplexND, B: onp.ToComplexND) -> onp.ArrayND[np.float64 | np.float32]: ...
