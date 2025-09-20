from collections.abc import Mapping
from typing import Any, Literal, Protocol, TypeAlias, overload, type_check_only
from typing_extensions import TypeVar, deprecated

import numpy as np
import numpy_typing_compat as nptc
import optype.numpy as onp
import optype.numpy.compat as npc

from ._superlu import SuperLU
from scipy.sparse._base import SparseEfficiencyWarning, _spbase
from scipy.sparse._bsr import _bsr_base
from scipy.sparse._lil import _lil_base
from scipy.sparse._typing import _Sparse2D

__all__ = [
    "MatrixRankWarning",
    "factorized",
    "is_sptriangular",
    "spbandwidth",
    "spilu",
    "splu",
    "spsolve",
    "spsolve_triangular",
    "use_solver",
]

_SparseT = TypeVar("_SparseT", bound=_spbase)
_NumberT_contra = TypeVar("_NumberT_contra", bound=npc.number, contravariant=True)
_InexactT_co = TypeVar("_InexactT_co", bound=np.float32 | np.float64 | np.complex64 | np.complex128, covariant=True)

_PermcSpec: TypeAlias = Literal["COLAMD", "NATURAL", "MMD_ATA", "MMD_AT_PLUS_A"]

_ToF32Mat: TypeAlias = _Sparse2D[np.float32] | nptc.CanArray[tuple[Any, ...], np.dtype[np.float32]]
_ToF64Mat: TypeAlias = _Sparse2D[np.float64 | npc.integer] | onp.ToInt2D | onp.ToJustFloat64_2D
_ToC64Mat: TypeAlias = _Sparse2D[np.complex64] | nptc.CanArray[tuple[Any, ...], np.dtype[np.complex64]]
_ToC128Mat: TypeAlias = _Sparse2D[np.complex128] | onp.ToJustComplex128_2D

_ToFloatMat: TypeAlias = _Sparse2D[npc.floating | npc.integer | np.bool_] | onp.ToFloat2D
_ToFloatMatStrict: TypeAlias = _Sparse2D[npc.floating | npc.integer | np.bool_] | onp.ToFloatStrict2D
_ToComplexMat: TypeAlias = _Sparse2D[npc.complexfloating] | onp.ToJustComplex2D
_ToInexactMat: TypeAlias = _Sparse2D[Any] | onp.ToComplex128_2D
_ToInexactMatStrict: TypeAlias = _Sparse2D[Any] | onp.ToComplex128Strict2D

_AsF32: TypeAlias = npc.integer8 | npc.number16 | np.int32 | np.float16 | np.float32
_AsF64: TypeAlias = npc.integer | np.float16 | np.float32 | np.float64
_AsC64: TypeAlias = npc.integer8 | npc.integer16 | np.int32 | np.float16 | npc.inexact32
_AsC128: TypeAlias = npc.integer | np.float16 | npc.inexact32 | npc.inexact64

@type_check_only
class _SuperLU_solve(Protocol[_NumberT_contra, _InexactT_co]):
    @overload
    def __call__(self, rhs: onp.Array1D[_NumberT_contra | np.bool_]) -> onp.Array1D[_InexactT_co]: ...
    @overload
    def __call__(self, rhs: onp.Array2D[_NumberT_contra | np.bool_]) -> onp.Array2D[_InexactT_co]: ...
    @overload
    def __call__(self, rhs: onp.ArrayND[_NumberT_contra | np.bool_]) -> onp.ArrayND[_InexactT_co]: ...

###

class MatrixRankWarning(UserWarning): ...

def use_solver(*, useUmfpack: bool = True, assumeSortedIndices: bool = False) -> None: ...

# NOTE: the mypy ignores work around a mypy bug in overload overlap checking
@overload
def factorized(A: _ToF64Mat) -> _SuperLU_solve[_AsF64, np.float64]: ...  # type: ignore[overload-overlap]
@overload
def factorized(A: _ToC128Mat) -> _SuperLU_solve[_AsC128, np.complex128]: ...  # type: ignore[overload-overlap]
@overload
def factorized(A: _ToF32Mat) -> _SuperLU_solve[_AsF32, np.float32]: ...
@overload
def factorized(A: _ToC64Mat) -> _SuperLU_solve[_AsC64, np.complex64]: ...

#
@overload  # 2d float, sparse 2d
def spsolve(A: _ToInexactMat, b: _SparseT, permc_spec: _PermcSpec | None = None, use_umfpack: bool = True) -> _SparseT: ...
@overload  # 2d float, 1d float
def spsolve(
    A: _ToFloatMat, b: onp.ToFloatStrict1D, permc_spec: _PermcSpec | None = None, use_umfpack: bool = True
) -> onp.Array1D[np.float64]: ...
@overload  # 2d float, 2d float
def spsolve(
    A: _ToFloatMat, b: onp.ToFloatStrict2D, permc_spec: _PermcSpec | None = None, use_umfpack: bool = True
) -> onp.Array2D[np.float64]: ...
@overload  # 2d float, 1d or 2d float
def spsolve(
    A: _ToFloatMat, b: onp.ToFloat2D | onp.ToFloat1D, permc_spec: _PermcSpec | None = None, use_umfpack: bool = True
) -> onp.ArrayND[np.float64]: ...
@overload  # 2d complex, 1d complex
def spsolve(
    A: _ToComplexMat, b: onp.ToComplexStrict1D, permc_spec: _PermcSpec | None = None, use_umfpack: bool = True
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d complex, 2d complex
def spsolve(
    A: _ToComplexMat, b: onp.ToComplexStrict2D, permc_spec: _PermcSpec | None = None, use_umfpack: bool = True
) -> onp.Array2D[np.complex128]: ...
@overload  # 2d complex, 1d or 2d complex
def spsolve(
    A: _ToComplexMat, b: onp.ToComplex2D | onp.ToComplex1D, permc_spec: _PermcSpec | None = None, use_umfpack: bool = True
) -> onp.ArrayND[np.complex128]: ...
@overload  # 2d inexact, 1d or 2d inexact
def spsolve(
    A: _ToInexactMat, b: onp.ToComplex2D | onp.ToComplex1D, permc_spec: _PermcSpec | None = None, use_umfpack: bool = True
) -> onp.ArrayND[np.float64 | np.complex128]: ...

#
@overload  # 2d float, 1d float
def spsolve_triangular(
    A: _ToFloatMat,
    b: onp.ToFloatStrict1D,
    lower: bool = True,
    overwrite_A: bool = False,
    overwrite_b: bool = False,
    unit_diagonal: bool = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d float, 2d float
def spsolve_triangular(
    A: _ToFloatMat,
    b: _ToFloatMatStrict,
    lower: bool = True,
    overwrite_A: bool = False,
    overwrite_b: bool = False,
    unit_diagonal: bool = False,
) -> onp.Array2D[np.float64]: ...
@overload  # 2d float, 1d or 2d float
def spsolve_triangular(
    A: _ToFloatMat,
    b: _ToFloatMat | onp.ToFloat1D,
    lower: bool = True,
    overwrite_A: bool = False,
    overwrite_b: bool = False,
    unit_diagonal: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload  # 2d complex, 1d complex
def spsolve_triangular(
    A: _ToComplexMat,
    b: onp.ToComplexStrict1D,
    lower: bool = True,
    overwrite_A: bool = False,
    overwrite_b: bool = False,
    unit_diagonal: bool = False,
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d complex, 2d complex
def spsolve_triangular(
    A: _ToComplexMat,
    b: _ToInexactMatStrict,
    lower: bool = True,
    overwrite_A: bool = False,
    overwrite_b: bool = False,
    unit_diagonal: bool = False,
) -> onp.Array2D[np.complex128]: ...
@overload  # 2d complex, 1d or 2d complex
def spsolve_triangular(
    A: _ToComplexMat,
    b: _ToInexactMat | onp.ToComplex1D,
    lower: bool = True,
    overwrite_A: bool = False,
    overwrite_b: bool = False,
    unit_diagonal: bool = False,
) -> onp.ArrayND[np.complex128]: ...
@overload  # 2d inexact, 1d or 2d inexact
def spsolve_triangular(
    A: _ToInexactMat,
    b: _ToInexactMat | onp.ToComplex1D,
    lower: bool = True,
    overwrite_A: bool = False,
    overwrite_b: bool = False,
    unit_diagonal: bool = False,
) -> onp.ArrayND[np.float64 | np.complex128]: ...

# NOTE: keep in sync with `spilu`
# NOTE: the mypy ignores work around a mypy bug in overload overlap checking
@overload  # float64
def splu(  # type: ignore[overload-overlap]
    A: _ToF64Mat,
    permc_spec: _PermcSpec | None = None,
    diag_pivot_thresh: onp.ToFloat | None = None,
    relax: int | None = None,
    panel_size: int | None = None,
    options: Mapping[str, object] | None = None,
) -> SuperLU[np.float64]: ...
@overload  # complex128
def splu(  # type: ignore[overload-overlap]
    A: _ToC128Mat,
    permc_spec: _PermcSpec | None = None,
    diag_pivot_thresh: onp.ToFloat | None = None,
    relax: int | None = None,
    panel_size: int | None = None,
    options: Mapping[str, object] | None = None,
) -> SuperLU[np.complex128]: ...
@overload  # float32
def splu(
    A: _ToF32Mat,
    permc_spec: _PermcSpec | None = None,
    diag_pivot_thresh: onp.ToFloat | None = None,
    relax: int | None = None,
    panel_size: int | None = None,
    options: Mapping[str, object] | None = None,
) -> SuperLU[np.float32]: ...
@overload  # complex64
def splu(
    A: _ToC64Mat,
    permc_spec: _PermcSpec | None = None,
    diag_pivot_thresh: onp.ToFloat | None = None,
    relax: int | None = None,
    panel_size: int | None = None,
    options: Mapping[str, object] | None = None,
) -> SuperLU[np.complex64]: ...
@overload  # floating
def splu(
    A: _ToFloatMat,
    permc_spec: _PermcSpec | None = None,
    diag_pivot_thresh: onp.ToFloat | None = None,
    relax: int | None = None,
    panel_size: int | None = None,
    options: Mapping[str, object] | None = None,
) -> SuperLU[np.float32 | np.float64]: ...
@overload  # complexfloating
def splu(
    A: _ToComplexMat,
    permc_spec: _PermcSpec | None = None,
    diag_pivot_thresh: onp.ToFloat | None = None,
    relax: int | None = None,
    panel_size: int | None = None,
    options: Mapping[str, object] | None = None,
) -> SuperLU[np.complex64 | np.complex128]: ...
@overload  # unknown
def splu(
    A: _ToInexactMat,
    permc_spec: _PermcSpec | None = None,
    diag_pivot_thresh: onp.ToFloat | None = None,
    relax: int | None = None,
    panel_size: int | None = None,
    options: Mapping[str, object] | None = None,
) -> SuperLU: ...

# NOTE: keep in sync with `splu`
@overload  # float64
def spilu(  # type: ignore[overload-overlap]
    A: _ToF64Mat,
    drop_tol: onp.ToFloat | None = None,
    fill_factor: onp.ToFloat | None = None,
    drop_rule: str | None = None,
    permc_spec: _PermcSpec | None = None,
    diag_pivot_thresh: onp.ToFloat | None = None,
    relax: int | None = None,
    panel_size: int | None = None,
    options: Mapping[str, object] | None = None,
) -> SuperLU[np.float64]: ...
@overload  # complex128
def spilu(  # type: ignore[overload-overlap]
    A: _ToC128Mat,
    drop_tol: onp.ToFloat | None = None,
    fill_factor: onp.ToFloat | None = None,
    drop_rule: str | None = None,
    permc_spec: _PermcSpec | None = None,
    diag_pivot_thresh: onp.ToFloat | None = None,
    relax: int | None = None,
    panel_size: int | None = None,
    options: Mapping[str, object] | None = None,
) -> SuperLU[np.complex128]: ...
@overload  # float32
def spilu(
    A: _ToF32Mat,
    drop_tol: onp.ToFloat | None = None,
    fill_factor: onp.ToFloat | None = None,
    drop_rule: str | None = None,
    permc_spec: _PermcSpec | None = None,
    diag_pivot_thresh: onp.ToFloat | None = None,
    relax: int | None = None,
    panel_size: int | None = None,
    options: Mapping[str, object] | None = None,
) -> SuperLU[np.float32]: ...
@overload  # complex64
def spilu(
    A: _ToC64Mat,
    drop_tol: onp.ToFloat | None = None,
    fill_factor: onp.ToFloat | None = None,
    drop_rule: str | None = None,
    permc_spec: _PermcSpec | None = None,
    diag_pivot_thresh: onp.ToFloat | None = None,
    relax: int | None = None,
    panel_size: int | None = None,
    options: Mapping[str, object] | None = None,
) -> SuperLU[np.complex64]: ...
@overload  # floating
def spilu(
    A: _ToFloatMat,
    drop_tol: onp.ToFloat | None = None,
    fill_factor: onp.ToFloat | None = None,
    drop_rule: str | None = None,
    permc_spec: _PermcSpec | None = None,
    diag_pivot_thresh: onp.ToFloat | None = None,
    relax: int | None = None,
    panel_size: int | None = None,
    options: Mapping[str, object] | None = None,
) -> SuperLU[np.float32 | np.float64]: ...
@overload  # complexfloating
def spilu(
    A: _ToComplexMat,
    drop_tol: onp.ToFloat | None = None,
    fill_factor: onp.ToFloat | None = None,
    drop_rule: str | None = None,
    permc_spec: _PermcSpec | None = None,
    diag_pivot_thresh: onp.ToFloat | None = None,
    relax: int | None = None,
    panel_size: int | None = None,
    options: Mapping[str, object] | None = None,
) -> SuperLU[np.complex64 | np.complex128]: ...
@overload
def spilu(
    A: _ToInexactMat,
    drop_tol: onp.ToFloat | None = None,
    fill_factor: onp.ToFloat | None = None,
    drop_rule: str | None = None,
    permc_spec: _PermcSpec | None = None,
    diag_pivot_thresh: onp.ToFloat | None = None,
    relax: int | None = None,
    panel_size: int | None = None,
    options: Mapping[str, object] | None = None,
) -> SuperLU: ...

#
@overload
@deprecated("is_sptriangular needs sparse and not BSR format. Converting to CSR.", category=SparseEfficiencyWarning)
def is_sptriangular(A: _bsr_base) -> tuple[bool, bool]: ...
@overload
def is_sptriangular(A: _spbase) -> tuple[bool, bool]: ...

#
@overload
@deprecated("spbandwidth needs sparse format not LIL and BSR. Converting to CSR.", category=SparseEfficiencyWarning)
def spbandwidth(A: _bsr_base | _lil_base) -> tuple[int, int]: ...
@overload
def spbandwidth(A: _spbase) -> tuple[int, int]: ...
