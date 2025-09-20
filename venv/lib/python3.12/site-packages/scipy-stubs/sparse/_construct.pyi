from _typeshed import Incomplete
from collections.abc import Callable, Iterable, Sequence as Seq
from typing import Any, Literal, Never, Protocol, TypeAlias, TypeVar, overload, type_check_only

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc
import optype.typing as opt

from ._base import _spbase, sparray
from ._bsr import bsr_array, bsr_matrix
from ._coo import coo_array, coo_matrix
from ._csc import csc_array, csc_matrix
from ._csr import csr_array, csr_matrix
from ._dia import dia_array, dia_matrix
from ._dok import dok_array, dok_matrix
from ._lil import lil_array, lil_matrix
from ._matrix import spmatrix
from ._typing import _CanStack, _CanStackAs, _Format, _ToShape2D

__all__ = [
    "block_array",
    "block_diag",
    "bmat",
    "diags",
    "diags_array",
    "eye",
    "eye_array",
    "hstack",
    "identity",
    "kron",
    "kronsum",
    "rand",
    "random",
    "random_array",
    "spdiags",
    "vstack",
]

_T = TypeVar("_T")
_SCT = TypeVar("_SCT", bound=_Numeric, default=_Numeric)
_SCT0 = TypeVar("_SCT0", bound=_Numeric)
_ShapeT = TypeVar("_ShapeT", bound=tuple[int, *tuple[int, ...]])

_Numeric: TypeAlias = npc.number | np.bool_

_ToDType: TypeAlias = type[complex | _Numeric] | np.dtype[_Numeric] | str

_SpBase2D: TypeAlias = _spbase[_SCT, tuple[int, int]]
_SpArray2D: TypeAlias = sparray[_SCT, tuple[int, int]]

_BSRArray: TypeAlias = bsr_array[_SCT]
_COOArray2D: TypeAlias = coo_array[_SCT, tuple[int, int]]
_CSCArray: TypeAlias = csc_array[_SCT]
_CSRArray2D: TypeAlias = csr_array[_SCT, tuple[int, int]]
_DIAArray: TypeAlias = dia_array[_SCT]
_DOKArray2D: TypeAlias = dok_array[_SCT, tuple[int, int]]
_LILArray: TypeAlias = lil_array[_SCT]

_ToArray1D2D: TypeAlias = onp.CanArray[tuple[int] | tuple[int, int], np.dtype[_SCT]] | Seq[_SCT | Seq[_SCT]]
_ToSpMatrix: TypeAlias = spmatrix[_SCT] | onp.ToArray2D[Never, _SCT]
_ToSparse2D: TypeAlias = _SpBase2D[_SCT] | onp.ToArray2D[Never, _SCT]

_FmtBSR: TypeAlias = Literal["bsr"]
_FmtCOO: TypeAlias = Literal["coo"]
_FmtCSC: TypeAlias = Literal["csc"]
_FmtCSR: TypeAlias = Literal["csr"]
_FmtDIA: TypeAlias = Literal["dia"]
_FmtDOK: TypeAlias = Literal["dok"]
_FmtLIL: TypeAlias = Literal["lil"]
_FmtNonCOO: TypeAlias = Literal["bsr", "csc", "csr", "dia", "dok", "lil"]

# TODO(julvandenbroeck): find a way to separate float and complex
_ComplexSeq1D2D: TypeAlias = Seq[Seq[complex] | complex]
_ToComplex1D2D: TypeAlias = onp.ToComplex1D | onp.ToComplex2D
_Offsets: TypeAlias = onp.ToInt | onp.ToInt1D

_DataRVS: TypeAlias = Callable[[int], onp.ArrayND[_Numeric]]

_ToBlocks: TypeAlias = Seq[Seq[_T | None]]
_ToBlocksSPArray: TypeAlias = _ToBlocks[_SpArray2D[_SCT]]
_ToBlocksSPMatrix: TypeAlias = _ToBlocks[spmatrix[_SCT]]
_ToBlocksCanStack: TypeAlias = _ToBlocks[_CanStack[_T]]
_ToBlocksCanStackAs: TypeAlias = _ToBlocks[_CanStackAs[_SCT0, _T]]
_ToBlocksArray: TypeAlias = onp.ArrayND[np.object_]
_ToBlocksUnkown: TypeAlias = _ToBlocksSPArray | _ToBlocksArray

_ToMatsDiag: TypeAlias = Iterable[_spbase[_SCT] | onp.ToArrayND[_SCT]]
_ToMatsDiagUnknown: TypeAlias = Iterable[_spbase | onp.ArrayND[_Numeric] | complex | Seq[onp.ToComplex] | Seq[onp.ToComplex1D]]

@type_check_only
class _DataSampler(Protocol):
    def __call__(self, /, *, size: int) -> onp.ArrayND[_Numeric]: ...

###
@overload  # diagonals: <complex>, format: "dia" | None, dtype: None
def diags_array(
    diagonals: _ComplexSeq1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDIA | None = None,
    dtype: None = None,
) -> _DIAArray[np.float64] | _DIAArray[np.complex128]: ...
@overload  # diagonals: <complex>, format: "bsr", dtype: None
def diags_array(
    diagonals: _ComplexSeq1D2D, /, *, offsets: _Offsets = 0, shape: _ToShape2D | None = None, format: _FmtBSR, dtype: None = None
) -> _BSRArray[np.float64] | _BSRArray[np.complex128]: ...
@overload  # diagonals: <complex>, format: "coo", dtype: None
def diags_array(
    diagonals: _ComplexSeq1D2D, /, *, offsets: _Offsets = 0, shape: _ToShape2D | None = None, format: _FmtCOO, dtype: None = None
) -> _COOArray2D[np.float64] | _COOArray2D[np.complex128]: ...
@overload  # diagonals: <complex>, format: "csc", dtype: None
def diags_array(
    diagonals: _ComplexSeq1D2D, /, *, offsets: _Offsets = 0, shape: _ToShape2D | None = None, format: _FmtCSC, dtype: None = None
) -> _CSCArray[np.float64] | _CSCArray[np.complex128]: ...
@overload  # diagonals: <complex>, format: "csr", dtype: None
def diags_array(
    diagonals: _ComplexSeq1D2D, /, *, offsets: _Offsets = 0, shape: _ToShape2D | None = None, format: _FmtCSR, dtype: None = None
) -> _CSRArray2D[np.float64] | _CSRArray2D[np.complex128]: ...
@overload  # diagonals: <complex>, format: "dok", dtype: None
def diags_array(
    diagonals: _ComplexSeq1D2D, /, *, offsets: _Offsets = 0, shape: _ToShape2D | None = None, format: _FmtDOK, dtype: None = None
) -> _DOKArray2D[np.float64] | _DOKArray2D[np.complex128]: ...
@overload  # diagonals: <complex>, format: "lil", dtype: None
def diags_array(
    diagonals: _ComplexSeq1D2D, /, *, offsets: _Offsets = 0, shape: _ToShape2D | None = None, format: _FmtLIL, dtype: None = None
) -> _LILArray[np.float64] | _LILArray[np.complex128]: ...

#
@overload  # diagonals: <known>, format: "dia" | None, dtype: None
def diags_array(
    diagonals: _ToArray1D2D[_SCT],
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDIA | None = None,
    dtype: None = None,
) -> _DIAArray[_SCT]: ...
@overload  # diagonals: <known>, format: "bsr", dtype: None
def diags_array(
    diagonals: _ToArray1D2D[_SCT],
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtBSR,
    dtype: None = None,
) -> _BSRArray[_SCT]: ...
@overload  # diagonals: <known>, format: "coo", dtype: None
def diags_array(
    diagonals: _ToArray1D2D[_SCT],
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCOO,
    dtype: None = None,
) -> _COOArray2D[_SCT]: ...
@overload  # diagonals: <known>, format: "csc", dtype: None
def diags_array(
    diagonals: _ToArray1D2D[_SCT],
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCSC,
    dtype: None = None,
) -> _CSCArray[_SCT]: ...
@overload  # diagonals: <known>, format: "csr", dtype: None
def diags_array(
    diagonals: _ToArray1D2D[_SCT],
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCSR,
    dtype: None = None,
) -> _CSRArray2D[_SCT]: ...
@overload  # diagonals: <known>, format: "dok", dtype: None
def diags_array(
    diagonals: _ToArray1D2D[_SCT],
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDOK,
    dtype: None = None,
) -> _DOKArray2D[_SCT]: ...
@overload  # diagonals: <known>, format: "lil", dtype: None
def diags_array(
    diagonals: _ToArray1D2D[_SCT],
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtLIL,
    dtype: None = None,
) -> _LILArray[_SCT]: ...

#
@overload  # diagonals: <unknown>, format: "dia" | None, dtype: bool-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDIA | None = None,
    dtype: onp.AnyBoolDType,
) -> _DIAArray[np.bool_]: ...
@overload  # diagonals: <unknown>, format: "bsr", dtype: bool-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtBSR,
    dtype: onp.AnyBoolDType,
) -> _BSRArray[np.bool_]: ...
@overload  # diagonals: <unknown>, format: "coo", dtype: bool-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCOO,
    dtype: onp.AnyBoolDType,
) -> _COOArray2D[np.bool_]: ...
@overload  # diagonals: <unknown>, format: "csc", dtype: bool-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCSC,
    dtype: onp.AnyBoolDType,
) -> _CSCArray[np.bool_]: ...
@overload  # diagonals: <unknown>, format: "csr", dtype: bool-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCSR,
    dtype: onp.AnyBoolDType,
) -> _CSRArray2D[np.bool_]: ...
@overload  # diagonals: <unknown>, format: "dok", dtype: bool-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDOK,
    dtype: onp.AnyBoolDType,
) -> _DOKArray2D[np.bool_]: ...
@overload  # diagonals: <unknown>, format: "lil", dtype: bool-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtLIL,
    dtype: onp.AnyBoolDType,
) -> _LILArray[np.bool_]: ...

#
@overload  # diagonals: <unknown>, format: "dia" | None, dtype: int-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDIA | None = None,
    dtype: onp.AnyIntDType,
) -> _DIAArray[np.int_]: ...
@overload  # diagonals: <unknown>, format: "bsr", dtype: int-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtBSR,
    dtype: onp.AnyIntDType,
) -> _BSRArray[np.int_]: ...
@overload  # diagonals: <unknown>, format: "coo", dtype: int-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCOO,
    dtype: onp.AnyIntDType,
) -> _COOArray2D[np.int_]: ...
@overload  # diagonals: <unknown>, format: "csc", dtype: int-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCSC,
    dtype: onp.AnyIntDType,
) -> _CSCArray[np.int_]: ...
@overload  # diagonals: <unknown>, format: "csr", dtype: int-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCSR,
    dtype: onp.AnyIntDType,
) -> _CSRArray2D[np.int_]: ...
@overload  # diagonals: <unknown>, format: "dok", dtype: int-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDOK,
    dtype: onp.AnyIntDType,
) -> _DOKArray2D[np.int_]: ...
@overload  # diagonals: <unknown>, format: "lil", dtype: int-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtLIL,
    dtype: onp.AnyIntDType,
) -> _LILArray[np.int_]: ...

#
@overload  # diagonals: <unknown>, format: "dia" | None, dtype: float64-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDIA | None = None,
    dtype: onp.AnyFloat64DType,
) -> _DIAArray[np.float64]: ...
@overload  # diagonals: <unknown>, format: "bsr", dtype: float64-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtBSR,
    dtype: onp.AnyFloat64DType,
) -> _BSRArray[np.float64]: ...
@overload  # diagonals: <unknown>, format: "coo", dtype: float64-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCOO,
    dtype: onp.AnyFloat64DType,
) -> _COOArray2D[np.float64]: ...
@overload  # diagonals: <unknown>, format: "csc", dtype: float64-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCSC,
    dtype: onp.AnyFloat64DType,
) -> _CSCArray[np.float64]: ...
@overload  # diagonals: <unknown>, format: "csr", dtype: float64-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCSR,
    dtype: onp.AnyFloat64DType,
) -> _CSRArray2D[np.float64]: ...
@overload  # diagonals: <unknown>, format: "dok", dtype: float64-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDOK,
    dtype: onp.AnyFloat64DType,
) -> _DOKArray2D[np.float64]: ...
@overload  # diagonals: <unknown>, format: "lil", dtype: float64-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtLIL,
    dtype: onp.AnyFloat64DType,
) -> _LILArray[np.float64]: ...

#
@overload  # diagonals: <unknown>, format: "dia" | None, dtype: complex128-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDIA | None = None,
    dtype: onp.AnyComplex128DType,
) -> _DIAArray[np.complex128]: ...
@overload  # diagonals: <unknown>, format: "bsr", dtype: complex128-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtBSR,
    dtype: onp.AnyComplex128DType,
) -> _BSRArray[np.complex128]: ...
@overload  # diagonals: <unknown>, format: "coo", dtype: complex128-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCOO,
    dtype: onp.AnyComplex128DType,
) -> _COOArray2D[np.complex128]: ...
@overload  # diagonals: <unknown>, format: "csc", dtype: complex128-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCSC,
    dtype: onp.AnyComplex128DType,
) -> _CSCArray[np.complex128]: ...
@overload  # diagonals: <unknown>, format: "csr", dtype: complex128-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCSR,
    dtype: onp.AnyComplex128DType,
) -> _CSRArray2D[np.complex128]: ...
@overload  # diagonals: <unknown>, format: "dok", dtype: complex128-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDOK,
    dtype: onp.AnyComplex128DType,
) -> _DOKArray2D[np.complex128]: ...
@overload  # diagonals: <unknown>, format: "lil", dtype: complex128-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtLIL,
    dtype: onp.AnyComplex128DType,
) -> _LILArray[np.complex128]: ...

#
@overload  # diagonals: <unknown>, format: "dia" | None, dtype: <known>
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDIA | None = None,
    dtype: onp.ToDType[_SCT],
) -> _DIAArray[_SCT]: ...
@overload  # diagonals: <unknown>, format: "bsr", dtype: <known>
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtBSR,
    dtype: onp.ToDType[_SCT],
) -> _BSRArray[_SCT]: ...
@overload  # diagonals: <unknown>, format: "coo", dtype: <known>
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCOO,
    dtype: onp.ToDType[_SCT],
) -> _COOArray2D[_SCT]: ...
@overload  # diagonals: <unknown>, format: "csc", dtype: <known>
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCSC,
    dtype: onp.ToDType[_SCT],
) -> _CSCArray[_SCT]: ...
@overload  # diagonals: <unknown>, format: "csr", dtype: <known>
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtCSR,
    dtype: onp.ToDType[_SCT],
) -> _CSRArray2D[_SCT]: ...
@overload  # diagonals: <unknown>, format: "dok", dtype: <known>
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDOK,
    dtype: onp.ToDType[_SCT],
) -> _DOKArray2D[_SCT]: ...
@overload  # diagonals: <unknown>, format: "lil", dtype: <known>
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtLIL,
    dtype: onp.ToDType[_SCT],
) -> _LILArray[_SCT]: ...

#
@overload  # diagonals: <unknown>, format: "dia" | None, dtype: <unknown>
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDIA | None = None,
    dtype: _ToDType,
) -> _DIAArray: ...
@overload  # diagonals: <unknown>, format: "bsr", dtype: <unknown>
def diags_array(
    diagonals: _ToComplex1D2D, /, *, offsets: _Offsets = 0, shape: _ToShape2D | None = None, format: _FmtBSR, dtype: _ToDType
) -> _BSRArray: ...
@overload  # diagonals: <unknown>, format: "coo", dtype: <unknown>
def diags_array(
    diagonals: _ToComplex1D2D, /, *, offsets: _Offsets = 0, shape: _ToShape2D | None = None, format: _FmtCOO, dtype: _ToDType
) -> _COOArray2D: ...
@overload  # diagonals: <unknown>, format: "csc", dtype: <unknown>
def diags_array(
    diagonals: _ToComplex1D2D, /, *, offsets: _Offsets = 0, shape: _ToShape2D | None = None, format: _FmtCSC, dtype: _ToDType
) -> _CSCArray: ...
@overload  # diagonals: <unknown>, format: "csr", dtype: <unknown>
def diags_array(
    diagonals: _ToComplex1D2D, /, *, offsets: _Offsets = 0, shape: _ToShape2D | None = None, format: _FmtCSR, dtype: _ToDType
) -> _CSRArray2D: ...
@overload  # diagonals: <unknown>, format: "dok", dtype: <unknown>
def diags_array(
    diagonals: _ToComplex1D2D, /, *, offsets: _Offsets = 0, shape: _ToShape2D | None = None, format: _FmtDOK, dtype: _ToDType
) -> _DOKArray2D: ...
@overload  # diagonals: <unknown>, format: "lil", dtype: <unknown>
def diags_array(
    diagonals: _ToComplex1D2D, /, *, offsets: _Offsets = 0, shape: _ToShape2D | None = None, format: _FmtLIL, dtype: _ToDType
) -> _LILArray: ...

###
# NOTE: `diags_array` should be prefered over `diags`
@overload  # diagonals: <known>, format: "dia" | None, dtype: None
def diags(
    diagonals: _ToArray1D2D[_SCT],
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDIA | None = None,
    dtype: None = None,
) -> dia_matrix[_SCT]: ...
@overload  # diagonals: <known>, format: "bsr", dtype: None
def diags(
    diagonals: _ToArray1D2D[_SCT], offsets: _Offsets = 0, shape: _ToShape2D | None = None, *, format: _FmtBSR, dtype: None = None
) -> bsr_matrix[_SCT]: ...
@overload  # diagonals: <known>, format: "coo", dtype: None
def diags(
    diagonals: _ToArray1D2D[_SCT], offsets: _Offsets = 0, shape: _ToShape2D | None = None, *, format: _FmtCOO, dtype: None = None
) -> coo_matrix[_SCT]: ...
@overload  # diagonals: <known>, format: "csc", dtype: None
def diags(
    diagonals: _ToArray1D2D[_SCT], offsets: _Offsets = 0, shape: _ToShape2D | None = None, *, format: _FmtCSC, dtype: None = None
) -> csc_matrix[_SCT]: ...
@overload  # diagonals: <known>, format: "csr", dtype: None
def diags(
    diagonals: _ToArray1D2D[_SCT], offsets: _Offsets = 0, shape: _ToShape2D | None = None, *, format: _FmtCSR, dtype: None = None
) -> csr_matrix[_SCT]: ...
@overload  # diagonals: <known>, format: "dok", dtype: None
def diags(
    diagonals: _ToArray1D2D[_SCT], offsets: _Offsets = 0, shape: _ToShape2D | None = None, *, format: _FmtDOK, dtype: None = None
) -> dok_matrix[_SCT]: ...
@overload  # diagonals: <known>, format: "lil", dtype: None
def diags(
    diagonals: _ToArray1D2D[_SCT], offsets: _Offsets = 0, shape: _ToShape2D | None = None, *, format: _FmtLIL, dtype: None = None
) -> lil_matrix[_SCT]: ...

#
@overload  # diagonals: <known>, format: "dia" | None, dtype: <known>
def diags(
    diagonals: _ToArray1D2D,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDIA | None = None,
    *,
    dtype: onp.ToDType[_SCT],
) -> dia_matrix[_SCT]: ...
@overload  # diagonals: <known>, format: "bsr", dtype: <known>
def diags(
    diagonals: _ToArray1D2D, offsets: _Offsets = 0, shape: _ToShape2D | None = None, *, format: _FmtBSR, dtype: onp.ToDType[_SCT]
) -> bsr_matrix[_SCT]: ...
@overload  # diagonals: <known>, format: "coo", dtype: <known>
def diags(
    diagonals: _ToArray1D2D, offsets: _Offsets = 0, shape: _ToShape2D | None = None, *, format: _FmtCOO, dtype: onp.ToDType[_SCT]
) -> coo_matrix[_SCT]: ...
@overload  # diagonals: <known>, format: "csr", dtype: <known>
def diags(
    diagonals: _ToArray1D2D, offsets: _Offsets = 0, shape: _ToShape2D | None = None, *, format: _FmtCSR, dtype: onp.ToDType[_SCT]
) -> csr_matrix[_SCT]: ...
@overload  # diagonals: <known>, format: "csc", dtype: <known>
def diags(
    diagonals: _ToArray1D2D, offsets: _Offsets = 0, shape: _ToShape2D | None = None, *, format: _FmtCSC, dtype: onp.ToDType[_SCT]
) -> csc_matrix[_SCT]: ...
@overload  # diagonals: <known>, format: "dok", dtype: <known>
def diags(
    diagonals: _ToArray1D2D, offsets: _Offsets = 0, shape: _ToShape2D | None = None, *, format: _FmtDOK, dtype: onp.ToDType[_SCT]
) -> dok_matrix[_SCT]: ...
@overload  # diagonals: <known>, format: "lil", dtype: <known>
def diags(
    diagonals: _ToArray1D2D, offsets: _Offsets = 0, shape: _ToShape2D | None = None, *, format: _FmtLIL, dtype: onp.ToDType[_SCT]
) -> lil_matrix[_SCT]: ...

#
@overload  # diagonals: <unknown>, format: "dia" | None, dtype: <unknown>
def diags(
    diagonals: _ToComplex1D2D,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    format: _FmtDIA | None = None,
    dtype: _ToDType | None = None,
) -> dia_matrix: ...
@overload  # diagonals: <unknown>, format: "bsr", dtype: <unknown>
def diags(
    diagonals: _ToComplex1D2D,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    *,
    format: _FmtBSR,
    dtype: _ToDType | None = None,
) -> bsr_matrix: ...
@overload  # diagonals: <unknown>, format: "coo", dtype: <unknown>
def diags(
    diagonals: _ToComplex1D2D,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    *,
    format: _FmtCOO,
    dtype: _ToDType | None = None,
) -> coo_matrix: ...
@overload  # diagonals: <unknown>, format: "csr", dtype: <unknown>
def diags(
    diagonals: _ToComplex1D2D,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    *,
    format: _FmtCSR,
    dtype: _ToDType | None = None,
) -> csr_matrix: ...
@overload  # diagonals: <unknown>, format: "csc", dtype: <unknown>
def diags(
    diagonals: _ToComplex1D2D,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    *,
    format: _FmtCSC,
    dtype: _ToDType | None = None,
) -> csc_matrix: ...
@overload  # diagonals: <unknown>, format: "dok", dtype: <unknown>
def diags(
    diagonals: _ToComplex1D2D,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    *,
    format: _FmtDOK,
    dtype: _ToDType | None = None,
) -> dok_matrix: ...
@overload  # diagonals: <unknown>, format: "lil", dtype: <unknown>
def diags(
    diagonals: _ToComplex1D2D,
    offsets: _Offsets = 0,
    shape: _ToShape2D | None = None,
    *,
    format: _FmtLIL,
    dtype: _ToDType | None = None,
) -> lil_matrix: ...

###
# NOTE: `diags_array` should be prefered over `spdiags`
@overload
def spdiags(
    data: _ToArray1D2D[_SCT], diags: _Offsets, m: onp.ToJustInt, n: onp.ToJustInt, format: _FmtDIA | None = None
) -> dia_matrix[_SCT]: ...
@overload
def spdiags(
    data: _ToArray1D2D[_SCT], diags: _Offsets, m: onp.ToJustInt, n: onp.ToJustInt, format: _FmtBSR
) -> bsr_matrix[_SCT]: ...
@overload
def spdiags(
    data: _ToArray1D2D[_SCT], diags: _Offsets, m: onp.ToJustInt, n: onp.ToJustInt, format: _FmtCOO
) -> coo_matrix[_SCT]: ...
@overload
def spdiags(
    data: _ToArray1D2D[_SCT], diags: _Offsets, m: onp.ToJustInt, n: onp.ToJustInt, format: _FmtCSC
) -> csc_matrix[_SCT]: ...
@overload
def spdiags(
    data: _ToArray1D2D[_SCT], diags: _Offsets, m: onp.ToJustInt, n: onp.ToJustInt, format: _FmtCSR
) -> csr_matrix[_SCT]: ...
@overload
def spdiags(
    data: _ToArray1D2D[_SCT], diags: _Offsets, m: onp.ToJustInt, n: onp.ToJustInt, format: _FmtDOK
) -> dok_matrix[_SCT]: ...
@overload
def spdiags(
    data: _ToArray1D2D[_SCT], diags: _Offsets, m: onp.ToJustInt, n: onp.ToJustInt, format: _FmtLIL
) -> lil_matrix[_SCT]: ...

#
@overload
def spdiags(
    data: _ToArray1D2D[_SCT],
    diags: _Offsets,
    m: tuple[onp.ToJustInt, onp.ToJustInt] | None = None,
    n: None = None,
    format: _FmtDIA | None = None,
) -> dia_matrix[_SCT]: ...
@overload
def spdiags(
    data: _ToArray1D2D[_SCT],
    diags: _Offsets,
    m: tuple[onp.ToJustInt, onp.ToJustInt] | None = None,
    n: None = None,
    *,
    format: _FmtBSR,
) -> bsr_matrix[_SCT]: ...
@overload
def spdiags(
    data: _ToArray1D2D[_SCT],
    diags: _Offsets,
    m: tuple[onp.ToJustInt, onp.ToJustInt] | None = None,
    n: None = None,
    *,
    format: _FmtCOO,
) -> coo_matrix[_SCT]: ...
@overload
def spdiags(
    data: _ToArray1D2D[_SCT],
    diags: _Offsets,
    m: tuple[onp.ToJustInt, onp.ToJustInt] | None = None,
    n: None = None,
    *,
    format: _FmtCSC,
) -> csc_matrix[_SCT]: ...
@overload
def spdiags(
    data: _ToArray1D2D[_SCT],
    diags: _Offsets,
    m: tuple[onp.ToJustInt, onp.ToJustInt] | None = None,
    n: None = None,
    *,
    format: _FmtCSR,
) -> csr_matrix[_SCT]: ...
@overload
def spdiags(
    data: _ToArray1D2D[_SCT],
    diags: _Offsets,
    m: tuple[onp.ToJustInt, onp.ToJustInt] | None = None,
    n: None = None,
    *,
    format: _FmtDOK,
) -> dok_matrix[_SCT]: ...
@overload
def spdiags(
    data: _ToArray1D2D[_SCT],
    diags: _Offsets,
    m: tuple[onp.ToJustInt, onp.ToJustInt] | None = None,
    n: None = None,
    *,
    format: _FmtLIL,
) -> lil_matrix[_SCT]: ...

###
@overload  # dtype: float64-like (default), format: "dia" | None
def identity(n: opt.AnyInt, dtype: onp.AnyFloat64DType = "d", format: _FmtDIA | None = None) -> dia_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "bsr"
def identity(n: opt.AnyInt, dtype: onp.AnyFloat64DType = "d", *, format: _FmtBSR) -> bsr_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "coo"
def identity(n: opt.AnyInt, dtype: onp.AnyFloat64DType = "d", *, format: _FmtCOO) -> coo_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "csc"
def identity(n: opt.AnyInt, dtype: onp.AnyFloat64DType = "d", *, format: _FmtCSC) -> csc_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "csr"
def identity(n: opt.AnyInt, dtype: onp.AnyFloat64DType = "d", *, format: _FmtCSR) -> csr_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "dok"
def identity(n: opt.AnyInt, dtype: onp.AnyFloat64DType = "d", *, format: _FmtDOK) -> dok_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "lil"
def identity(n: opt.AnyInt, dtype: onp.AnyFloat64DType = "d", *, format: _FmtLIL) -> lil_matrix[np.float64]: ...

#
@overload  # dtype: bool-like, format: "dia" | None
def identity(n: opt.AnyInt, dtype: onp.AnyBoolDType, format: _FmtDIA | None = None) -> dia_matrix[np.bool_]: ...
@overload  # dtype: bool-like, format: "bsr"
def identity(n: opt.AnyInt, dtype: onp.AnyBoolDType, format: _FmtBSR) -> bsr_matrix[np.bool_]: ...
@overload  # dtype: bool-like, format: "coo"
def identity(n: opt.AnyInt, dtype: onp.AnyBoolDType, format: _FmtCOO) -> coo_matrix[np.bool_]: ...
@overload  # dtype: bool-like, format: "csc"
def identity(n: opt.AnyInt, dtype: onp.AnyBoolDType, format: _FmtCSC) -> csc_matrix[np.bool_]: ...
@overload  # dtype: bool-like, format: "csr"
def identity(n: opt.AnyInt, dtype: onp.AnyBoolDType, format: _FmtCSR) -> csr_matrix[np.bool_]: ...
@overload  # dtype: bool-like, format: "dok"
def identity(n: opt.AnyInt, dtype: onp.AnyBoolDType, format: _FmtDOK) -> dok_matrix[np.bool_]: ...
@overload  # dtype: bool-like, format: "lil"
def identity(n: opt.AnyInt, dtype: onp.AnyBoolDType, format: _FmtLIL) -> lil_matrix[np.bool_]: ...

#
@overload  # dtype: int-like, format: "dia" | None
def identity(n: opt.AnyInt, dtype: onp.AnyIntDType, format: _FmtDIA | None = None) -> dia_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "bsr"
def identity(n: opt.AnyInt, dtype: onp.AnyIntDType, format: _FmtBSR) -> bsr_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "coo"
def identity(n: opt.AnyInt, dtype: onp.AnyIntDType, format: _FmtCOO) -> coo_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "csc"
def identity(n: opt.AnyInt, dtype: onp.AnyIntDType, format: _FmtCSC) -> csc_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "csr"
def identity(n: opt.AnyInt, dtype: onp.AnyIntDType, format: _FmtCSR) -> csr_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "dok"
def identity(n: opt.AnyInt, dtype: onp.AnyIntDType, format: _FmtDOK) -> dok_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "lil"
def identity(n: opt.AnyInt, dtype: onp.AnyIntDType, format: _FmtLIL) -> lil_matrix[np.int_]: ...

#
@overload  # dtype: complex128-like, format: "dia" | None
def identity(n: opt.AnyInt, dtype: onp.AnyComplex128DType, format: _FmtDIA | None = None) -> dia_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "bsr"
def identity(n: opt.AnyInt, dtype: onp.AnyComplex128DType, format: _FmtBSR) -> bsr_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "coo"
def identity(n: opt.AnyInt, dtype: onp.AnyComplex128DType, format: _FmtCOO) -> coo_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "csc"
def identity(n: opt.AnyInt, dtype: onp.AnyComplex128DType, format: _FmtCSC) -> csc_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "csr"
def identity(n: opt.AnyInt, dtype: onp.AnyComplex128DType, format: _FmtCSR) -> csr_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "dok"
def identity(n: opt.AnyInt, dtype: onp.AnyComplex128DType, format: _FmtDOK) -> dok_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "lil"
def identity(n: opt.AnyInt, dtype: onp.AnyComplex128DType, format: _FmtLIL) -> lil_matrix[np.complex128]: ...

#
@overload  # dtype: <known>, format: "dia" | None
def identity(n: opt.AnyInt, dtype: onp.ToDType[_SCT], format: _FmtDIA | None = None) -> dia_matrix[_SCT]: ...
@overload  # dtype: <known>, format: "bsr"
def identity(n: opt.AnyInt, dtype: onp.ToDType[_SCT], format: _FmtBSR) -> bsr_matrix[_SCT]: ...
@overload  # dtype: <known>, format: "coo"
def identity(n: opt.AnyInt, dtype: onp.ToDType[_SCT], format: _FmtCOO) -> coo_matrix[_SCT]: ...
@overload  # dtype: <known>, format: "csc"
def identity(n: opt.AnyInt, dtype: onp.ToDType[_SCT], format: _FmtCSC) -> csc_matrix[_SCT]: ...
@overload  # dtype: <known>, format: "csr"
def identity(n: opt.AnyInt, dtype: onp.ToDType[_SCT], format: _FmtCSR) -> csr_matrix[_SCT]: ...
@overload  # dtype: <known>, format: "dok"
def identity(n: opt.AnyInt, dtype: onp.ToDType[_SCT], format: _FmtDOK) -> dok_matrix[_SCT]: ...
@overload  # dtype: <known>, format: "lil"
def identity(n: opt.AnyInt, dtype: onp.ToDType[_SCT], format: _FmtLIL) -> lil_matrix[_SCT]: ...

#
@overload  # dtype: <unknown>, format: "dia" | None
def identity(n: opt.AnyInt, dtype: _ToDType, format: _FmtDIA | None = None) -> dia_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "bsr"
def identity(n: opt.AnyInt, dtype: _ToDType, format: _FmtBSR) -> bsr_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "coo"
def identity(n: opt.AnyInt, dtype: _ToDType, format: _FmtCOO) -> coo_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "csc"
def identity(n: opt.AnyInt, dtype: _ToDType, format: _FmtCSC) -> csc_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "csr"
def identity(n: opt.AnyInt, dtype: _ToDType, format: _FmtCSR) -> csr_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "dok"
def identity(n: opt.AnyInt, dtype: _ToDType, format: _FmtDOK) -> dok_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "lil"
def identity(n: opt.AnyInt, dtype: _ToDType, format: _FmtLIL) -> lil_matrix[Incomplete]: ...

###
@overload  # dtype: float64-like (default), format: "dia" | None
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtDIA | None = None
) -> _DIAArray[np.float64]: ...
@overload  # dtype: float64-like (default), format: "bsr"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtBSR
) -> _BSRArray[np.float64]: ...
@overload  # dtype: float64-like (default), format: "coo"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtCOO
) -> _COOArray2D[np.float64]: ...
@overload  # dtype: float64-like (default), format: "csc"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtCSC
) -> _CSCArray[np.float64]: ...
@overload  # dtype: float64-like (default), format: "csr"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtCSR
) -> _CSRArray2D[np.float64]: ...
@overload  # dtype: float64-like (default), format: "dok"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtDOK
) -> _DOKArray2D[np.float64]: ...
@overload  # dtype: float64-like (default), format: "lil"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtLIL
) -> _LILArray[np.float64]: ...

#
@overload  # dtype: bool-like, format: "dia" | None
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyBoolDType, format: _FmtDIA | None = None
) -> _DIAArray[np.bool_]: ...
@overload  # dtype: bool-like, format: "bsr"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyBoolDType, format: _FmtBSR
) -> _BSRArray[np.bool_]: ...
@overload  # dtype: bool-like, format: "coo"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyBoolDType, format: _FmtCOO
) -> _COOArray2D[np.bool_]: ...
@overload  # dtype: bool-like, format: "csc"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyBoolDType, format: _FmtCSC
) -> _CSCArray[np.bool_]: ...
@overload  # dtype: bool-like, format: "csr"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyBoolDType, format: _FmtCSR
) -> _CSRArray2D[np.bool_]: ...
@overload  # dtype: bool-like, format: "dok"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyBoolDType, format: _FmtDOK
) -> _DOKArray2D[np.bool_]: ...
@overload  # dtype: bool-like, format: "lil"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyBoolDType, format: _FmtLIL
) -> _LILArray[np.bool_]: ...

#
@overload  # dtype: int-like, format: "dia" | None
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyIntDType, format: _FmtDIA | None = None
) -> _DIAArray[np.int_]: ...
@overload  # dtype: int-like, format: "bsr"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyIntDType, format: _FmtBSR
) -> _BSRArray[np.int_]: ...
@overload  # dtype: int-like, format: "coo"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyIntDType, format: _FmtCOO
) -> _COOArray2D[np.int_]: ...
@overload  # dtype: int-like, format: "csc"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyIntDType, format: _FmtCSC
) -> _CSCArray[np.int_]: ...
@overload  # dtype: int-like, format: "csr"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyIntDType, format: _FmtCSR
) -> _CSRArray2D[np.int_]: ...
@overload  # dtype: int-like, format: "dok"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyIntDType, format: _FmtDOK
) -> _DOKArray2D[np.int_]: ...
@overload  # dtype: int-like, format: "lil"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyIntDType, format: _FmtLIL
) -> _LILArray[np.int_]: ...

#
@overload  # dtype: complex128-like, format: "dia" | None
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyComplex128DType, format: _FmtDIA | None = None
) -> _DIAArray[np.complex128]: ...
@overload  # dtype: complex128-like, format: "bsr"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyComplex128DType, format: _FmtBSR
) -> _BSRArray[np.complex128]: ...
@overload  # dtype: complex128-like, format: "coo"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyComplex128DType, format: _FmtCOO
) -> _COOArray2D[np.complex128]: ...
@overload  # dtype: complex128-like, format: "csc"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyComplex128DType, format: _FmtCSC
) -> _CSCArray[np.complex128]: ...
@overload  # dtype: complex128-like, format: "csr"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyComplex128DType, format: _FmtCSR
) -> _CSRArray2D[np.complex128]: ...
@overload  # dtype: complex128-like, format: "dok"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyComplex128DType, format: _FmtDOK
) -> _DOKArray2D[np.complex128]: ...
@overload  # dtype: complex128-like, format: "lil"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.AnyComplex128DType, format: _FmtLIL
) -> _LILArray[np.complex128]: ...

#
@overload  # dtype: <known>, format: "dia" | None
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.ToDType[_SCT], format: _FmtDIA | None = None
) -> _DIAArray[_SCT]: ...
@overload  # dtype: <known>, format: "bsr"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.ToDType[_SCT], format: _FmtBSR
) -> _BSRArray[_SCT]: ...
@overload  # dtype: <known>, format: "coo"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.ToDType[_SCT], format: _FmtCOO
) -> _COOArray2D[_SCT]: ...
@overload  # dtype: <known>, format: "csc"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.ToDType[_SCT], format: _FmtCSC
) -> _CSCArray[_SCT]: ...
@overload  # dtype: <known>, format: "csr"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.ToDType[_SCT], format: _FmtCSR
) -> _CSRArray2D[_SCT]: ...
@overload  # dtype: <known>, format: "dok"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.ToDType[_SCT], format: _FmtDOK
) -> _DOKArray2D[_SCT]: ...
@overload  # dtype: <known>, format: "lil"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: onp.ToDType[_SCT], format: _FmtLIL
) -> _LILArray[_SCT]: ...

#
@overload  # dtype: <unknown>, format: "dia" | None
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: _ToDType, format: _FmtDIA | None = None
) -> _DIAArray[Incomplete]: ...
@overload  # dtype: <unknown>, format: "bsr"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: _ToDType, format: _FmtBSR
) -> _BSRArray[Incomplete]: ...
@overload  # dtype: <unknown>, format: "coo"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: _ToDType, format: _FmtCOO
) -> _COOArray2D[Incomplete]: ...
@overload  # dtype: <unknown>, format: "csc"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: _ToDType, format: _FmtCSC
) -> _CSCArray[Incomplete]: ...
@overload  # dtype: <unknown>, format: "csr"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: _ToDType, format: _FmtCSR
) -> _CSRArray2D[Incomplete]: ...
@overload  # dtype: <unknown>, format: "dok"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: _ToDType, format: _FmtDOK
) -> _DOKArray2D[Incomplete]: ...
@overload  # dtype: <unknown>, format: "lil"
def eye_array(
    m: opt.AnyInt, n: opt.AnyInt | None = None, *, k: int = 0, dtype: _ToDType, format: _FmtLIL
) -> _LILArray[Incomplete]: ...

###
# NOTE: `eye_array` should be prefered over `eye`
@overload  # dtype: float64-like (default), format: "dia" | None
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtDIA | None = None
) -> dia_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "bsr"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, dtype: onp.AnyFloat64DType = ..., *, format: _FmtBSR
) -> bsr_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "coo"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, dtype: onp.AnyFloat64DType = ..., *, format: _FmtCOO
) -> coo_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "csc"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, dtype: onp.AnyFloat64DType = ..., *, format: _FmtCSC
) -> csc_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "csr"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, dtype: onp.AnyFloat64DType = ..., *, format: _FmtCSR
) -> csr_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "dok"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, dtype: onp.AnyFloat64DType = ..., *, format: _FmtDOK
) -> dok_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "lil"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, dtype: onp.AnyFloat64DType = ..., *, format: _FmtLIL
) -> lil_matrix[np.float64]: ...

#
@overload  # dtype: bool-like, format: "dia" | None
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyBoolDType, format: _FmtDIA | None = None
) -> dia_matrix[np.bool_]: ...
@overload  # dtype: bool-like, format: "bsr"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyBoolDType, format: _FmtBSR
) -> bsr_matrix[np.bool_]: ...
@overload  # dtype: bool-like, format: "coo"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyBoolDType, format: _FmtCOO
) -> coo_matrix[np.bool_]: ...
@overload  # dtype: bool-like, format: "csc"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyBoolDType, format: _FmtCSC
) -> csc_matrix[np.bool_]: ...
@overload  # dtype: bool-like, format: "csr"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyBoolDType, format: _FmtCSR
) -> csr_matrix[np.bool_]: ...
@overload  # dtype: bool-like, format: "dok"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyBoolDType, format: _FmtDOK
) -> dok_matrix[np.bool_]: ...
@overload  # dtype: bool-like, format: "lil"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyBoolDType, format: _FmtLIL
) -> lil_matrix[np.bool_]: ...

#
@overload  # dtype: int-like, format: "dia" | None
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyIntDType, format: _FmtDIA | None = None
) -> dia_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "bsr"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyIntDType, format: _FmtBSR
) -> bsr_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "coo"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyIntDType, format: _FmtCOO
) -> coo_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "csc"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyIntDType, format: _FmtCSC
) -> csc_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "csr"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyIntDType, format: _FmtCSR
) -> csr_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "dok"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyIntDType, format: _FmtDOK
) -> dok_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "lil"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyIntDType, format: _FmtLIL
) -> lil_matrix[np.int_]: ...

#
@overload  # dtype: complex128-like, format: "dia" | None
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyComplex128DType, format: _FmtDIA | None = None
) -> dia_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "bsr"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyComplex128DType, format: _FmtBSR
) -> bsr_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "coo"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyComplex128DType, format: _FmtCOO
) -> coo_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "csc"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyComplex128DType, format: _FmtCSC
) -> csc_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "csr"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyComplex128DType, format: _FmtCSR
) -> csr_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "dok"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyComplex128DType, format: _FmtDOK
) -> dok_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "lil"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.AnyComplex128DType, format: _FmtLIL
) -> lil_matrix[np.complex128]: ...

#
@overload  # dtype: <known>, format: "dia" | None
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.ToDType[_SCT], format: _FmtDIA | None = None
) -> dia_matrix[_SCT]: ...
@overload  # dtype: <known>, format: "bsr"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.ToDType[_SCT], format: _FmtBSR
) -> bsr_matrix[_SCT]: ...
@overload  # dtype: <known>, format: "coo"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.ToDType[_SCT], format: _FmtCOO
) -> coo_matrix[_SCT]: ...
@overload  # dtype: <known>, format: "csc"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.ToDType[_SCT], format: _FmtCSC
) -> csc_matrix[_SCT]: ...
@overload  # dtype: <known>, format: "csr"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.ToDType[_SCT], format: _FmtCSR
) -> csr_matrix[_SCT]: ...
@overload  # dtype: <known>, format: "dok"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.ToDType[_SCT], format: _FmtDOK
) -> dok_matrix[_SCT]: ...
@overload  # dtype: <known>, format: "lil"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: onp.ToDType[_SCT], format: _FmtLIL
) -> lil_matrix[_SCT]: ...

#
@overload  # dtype: <unknown>, format: "dia" | None
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: _ToDType, format: _FmtDIA | None = None
) -> dia_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "bsr"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: _ToDType, format: _FmtBSR
) -> bsr_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "coo"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: _ToDType, format: _FmtCOO
) -> coo_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "csc"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: _ToDType, format: _FmtCSC
) -> csc_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "csr"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: _ToDType, format: _FmtCSR
) -> csr_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "dok"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: _ToDType, format: _FmtDOK
) -> dok_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "lil"
def eye(
    m: opt.AnyInt, n: opt.AnyInt | None = None, k: int = 0, *, dtype: _ToDType, format: _FmtLIL
) -> lil_matrix[Incomplete]: ...

###
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "bsr" | None
def kron(A: _ToSpMatrix[_SCT], B: _ToSpMatrix[_SCT], format: _FmtBSR | None = None) -> bsr_matrix[_SCT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "coo"
def kron(A: _ToSpMatrix[_SCT], B: _ToSpMatrix[_SCT], format: _FmtCOO) -> coo_matrix[_SCT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "csc"
def kron(A: _ToSpMatrix[_SCT], B: _ToSpMatrix[_SCT], format: _FmtCSC) -> csc_matrix[_SCT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "csr"
def kron(A: _ToSpMatrix[_SCT], B: _ToSpMatrix[_SCT], format: _FmtCSR) -> csr_matrix[_SCT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "dia"
def kron(A: _ToSpMatrix[_SCT], B: _ToSpMatrix[_SCT], format: _FmtDIA) -> dia_matrix[_SCT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "dok"
def kron(A: _ToSpMatrix[_SCT], B: _ToSpMatrix[_SCT], format: _FmtDOK) -> dok_matrix[_SCT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "lil"
def kron(A: _ToSpMatrix[_SCT], B: _ToSpMatrix[_SCT], format: _FmtLIL) -> lil_matrix[_SCT]: ...

#
@overload  # A: sparray, B: 2D sparse, format: "bsr" | None
def kron(A: sparray[_SCT, tuple[int, int]], B: _ToSparse2D[_SCT], format: _FmtBSR | None = None) -> _BSRArray[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "coo"
def kron(A: sparray[_SCT, tuple[int, int]], B: _ToSparse2D[_SCT], format: _FmtCOO) -> _COOArray2D[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "csc"
def kron(A: sparray[_SCT, tuple[int, int]], B: _ToSparse2D[_SCT], format: _FmtCSC) -> _CSCArray[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "csr"
def kron(A: sparray[_SCT, tuple[int, int]], B: _ToSparse2D[_SCT], format: _FmtCSR) -> _CSRArray2D[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "dia"
def kron(A: sparray[_SCT, tuple[int, int]], B: _ToSparse2D[_SCT], format: _FmtDIA) -> _DIAArray[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "dok"
def kron(A: sparray[_SCT, tuple[int, int]], B: _ToSparse2D[_SCT], format: _FmtDOK) -> _DOKArray2D[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "lil"
def kron(A: sparray[_SCT, tuple[int, int]], B: _ToSparse2D[_SCT], format: _FmtLIL) -> _LILArray[_SCT]: ...

#
@overload  # A: sparse, B: sparray, format: "bsr" | None
def kron(A: _ToSparse2D[_SCT], B: sparray[_SCT, tuple[int, int]], format: _FmtBSR | None = None) -> _BSRArray[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "coo"
def kron(A: _ToSparse2D[_SCT], B: sparray[_SCT, tuple[int, int]], format: _FmtCOO) -> _COOArray2D[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "csc"
def kron(A: _ToSparse2D[_SCT], B: sparray[_SCT, tuple[int, int]], format: _FmtCSC) -> _CSCArray[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "csr"
def kron(A: _ToSparse2D[_SCT], B: sparray[_SCT, tuple[int, int]], format: _FmtCSR) -> _CSRArray2D[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "dia"
def kron(A: _ToSparse2D[_SCT], B: sparray[_SCT, tuple[int, int]], format: _FmtDIA) -> _DIAArray[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "dok"
def kron(A: _ToSparse2D[_SCT], B: sparray[_SCT, tuple[int, int]], format: _FmtDOK) -> _DOKArray2D[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "lil"
def kron(A: _ToSparse2D[_SCT], B: sparray[_SCT, tuple[int, int]], format: _FmtLIL) -> _LILArray[_SCT]: ...
@overload  # A: unknown array-like, B: unknown array-like  (catch-all)
def kron(A: onp.ToComplex2D, B: onp.ToComplex2D, format: _Format | None = None) -> _SpBase2D[Incomplete]: ...

###
# NOTE: The `overload-overlap` mypy errors are false positives.
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "csr" | None
def kronsum(A: _ToSpMatrix[_SCT], B: _ToSpMatrix[_SCT], format: _FmtCSR | None = None) -> csr_matrix[_SCT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "bsr"
def kronsum(A: _ToSpMatrix[_SCT], B: _ToSpMatrix[_SCT], format: _FmtBSR) -> bsr_matrix[_SCT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "coo"
def kronsum(A: _ToSpMatrix[_SCT], B: _ToSpMatrix[_SCT], format: _FmtCOO) -> coo_matrix[_SCT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "csc"
def kronsum(A: _ToSpMatrix[_SCT], B: _ToSpMatrix[_SCT], format: _FmtCSC) -> csc_matrix[_SCT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "dia"
def kronsum(A: _ToSpMatrix[_SCT], B: _ToSpMatrix[_SCT], format: _FmtDIA) -> dia_matrix[_SCT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "dok"
def kronsum(A: _ToSpMatrix[_SCT], B: _ToSpMatrix[_SCT], format: _FmtDOK) -> dok_matrix[_SCT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "lil"
def kronsum(A: _ToSpMatrix[_SCT], B: _ToSpMatrix[_SCT], format: _FmtLIL) -> lil_matrix[_SCT]: ...

#
@overload  # A: sparray, B: sparse, format: "csr" | None
def kronsum(A: sparray[_SCT, tuple[int, int]], B: _ToSparse2D[_SCT], format: _FmtCSR | None = None) -> _CSRArray2D[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "bsr"
def kronsum(A: sparray[_SCT, tuple[int, int]], B: _ToSparse2D[_SCT], format: _FmtBSR) -> _BSRArray[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "coo"
def kronsum(A: sparray[_SCT, tuple[int, int]], B: _ToSparse2D[_SCT], format: _FmtCOO) -> _COOArray2D[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "csc"
def kronsum(A: sparray[_SCT, tuple[int, int]], B: _ToSparse2D[_SCT], format: _FmtCSC) -> _CSCArray[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "dia"
def kronsum(A: sparray[_SCT, tuple[int, int]], B: _ToSparse2D[_SCT], format: _FmtDIA) -> _DIAArray[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "dok"
def kronsum(A: sparray[_SCT, tuple[int, int]], B: _ToSparse2D[_SCT], format: _FmtDOK) -> _DOKArray2D[_SCT]: ...
@overload  # A: sparray, B: sparse, format: "lil"
def kronsum(A: sparray[_SCT, tuple[int, int]], B: _ToSparse2D[_SCT], format: _FmtLIL) -> _LILArray[_SCT]: ...

#
@overload  # A: sparse, B: sparray, format: "csr" | None
def kronsum(A: _ToSparse2D[_SCT], B: sparray[_SCT, tuple[int, int]], format: _FmtCSR | None = None) -> _CSRArray2D[_SCT]: ...
@overload  # A: sparse, B: sparray, format: "bsr"
def kronsum(A: _ToSparse2D[_SCT], B: sparray[_SCT, tuple[int, int]], format: _FmtBSR) -> _BSRArray[_SCT]: ...
@overload  # A: sparse, B: sparray, format: "coo"
def kronsum(A: _ToSparse2D[_SCT], B: sparray[_SCT, tuple[int, int]], format: _FmtCOO) -> _COOArray2D[_SCT]: ...
@overload  # A: sparse, B: sparray, format: "csc"
def kronsum(A: _ToSparse2D[_SCT], B: sparray[_SCT, tuple[int, int]], format: _FmtCSC) -> _CSCArray[_SCT]: ...
@overload  # A: sparse, B: sparray, format: "dia"
def kronsum(A: _ToSparse2D[_SCT], B: sparray[_SCT, tuple[int, int]], format: _FmtDIA) -> _DIAArray[_SCT]: ...
@overload  # A: sparse, B: sparray, format: "dok"
def kronsum(A: _ToSparse2D[_SCT], B: sparray[_SCT, tuple[int, int]], format: _FmtDOK) -> _DOKArray2D[_SCT]: ...
@overload  # A: sparse, B: sparray, format: "lil"
def kronsum(A: _ToSparse2D[_SCT], B: sparray[_SCT, tuple[int, int]], format: _FmtLIL) -> _LILArray[_SCT]: ...
@overload  # A: unknown array-like, B: unknown array-like  (catch-all)
def kronsum(A: onp.ToComplex2D, B: onp.ToComplex2D, format: _Format | None = None) -> _SpBase2D[Incomplete]: ...

###
# NOTE: keep in sync with `vstack`
@overload  # sparray, format: <default>, dtype: <default>
def hstack(blocks: Seq[_CanStack[_T]], format: None = None, dtype: None = None) -> _T: ...
@overload  # sparray, format: "bsr", dtype: <default>
def hstack(blocks: Seq[_SpArray2D[_SCT]], format: _FmtBSR, dtype: None = None) -> _BSRArray[_SCT]: ...
@overload  # sparray, format: "coo", dtype: <default>
def hstack(blocks: Seq[_SpArray2D[_SCT]], format: _FmtCOO, dtype: None = None) -> _COOArray2D[_SCT]: ...
@overload  # sparray, format: "csc", dtype: <default>
def hstack(blocks: Seq[_SpArray2D[_SCT]], format: _FmtCSC, dtype: None = None) -> _CSCArray[_SCT]: ...
@overload  # sparray, format: "csr", dtype: <default>
def hstack(blocks: Seq[_SpArray2D[_SCT]], format: _FmtCSR, dtype: None = None) -> _CSRArray2D[_SCT]: ...
@overload  # sparray, format: "dia", dtype: <default>
def hstack(blocks: Seq[_SpArray2D[_SCT]], format: _FmtDIA, dtype: None = None) -> _DIAArray[_SCT]: ...
@overload  # sparray, format: "dok", dtype: <default>
def hstack(blocks: Seq[_SpArray2D[_SCT]], format: _FmtDOK, dtype: None = None) -> _DOKArray2D[_SCT]: ...
@overload  # sparray, format: "lil", dtype: <default>
def hstack(blocks: Seq[_SpArray2D[_SCT]], format: _FmtLIL, dtype: None = None) -> _LILArray[_SCT]: ...

#
@overload  # sparray, format: <default>, dtype: bool-like
def hstack(blocks: Seq[_CanStackAs[np.bool_, _T]], format: None = None, *, dtype: onp.AnyBoolDType) -> _T: ...
@overload  # sparray, format: "bsr", dtype: bool-like
def hstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyBoolDType) -> _BSRArray[np.bool_]: ...
@overload  # sparray, format: "coo", dtype: bool-like
def hstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyBoolDType) -> _COOArray2D[np.bool_]: ...
@overload  # sparray, format: "csc", dtype: bool-like
def hstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyBoolDType) -> _CSCArray[np.bool_]: ...
@overload  # sparray, format: "csr", dtype: bool-like
def hstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyBoolDType) -> _CSRArray2D[np.bool_]: ...
@overload  # sparray, format: "dia", dtype: bool-like
def hstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyBoolDType) -> _DIAArray[np.bool_]: ...
@overload  # sparray, format: "dok", dtype: bool-like
def hstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyBoolDType) -> _DOKArray2D[np.bool_]: ...
@overload  # sparray, format: "lil", dtype: bool-like
def hstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyBoolDType) -> _LILArray[np.bool_]: ...

#
@overload  # sparray, format: <default>, dtype: int-like
def hstack(blocks: Seq[_CanStackAs[np.int_, _T]], format: None = None, *, dtype: onp.AnyIntDType) -> _T: ...
@overload  # sparray, format: "bsr", dtype: int-like
def hstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyIntDType) -> _BSRArray[np.int_]: ...
@overload  # sparray, format: "coo", dtype: int-like
def hstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyIntDType) -> _COOArray2D[np.int_]: ...
@overload  # sparray, format: "csc", dtype: int-like
def hstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyIntDType) -> _CSCArray[np.int_]: ...
@overload  # sparray, format: "csr", dtype: int-like
def hstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyIntDType) -> _CSRArray2D[np.int_]: ...
@overload  # sparray, format: "dia", dtype: int-like
def hstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyIntDType) -> _DIAArray[np.int_]: ...
@overload  # sparray, format: "dok", dtype: int-like
def hstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyIntDType) -> _DOKArray2D[np.int_]: ...
@overload  # sparray, format: "lil", dtype: int-like
def hstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyIntDType) -> _LILArray[np.int_]: ...

#
@overload  # sparray, format: <default>, dtype: float64-like
def hstack(blocks: Seq[_CanStackAs[np.float64, _T]], format: None = None, *, dtype: onp.AnyFloat64DType) -> _T: ...
@overload  # sparray, format: "bsr", dtype: float64-like
def hstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyFloat64DType) -> _BSRArray[np.float64]: ...
@overload  # sparray, format: "coo", dtype: float64-like
def hstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyFloat64DType) -> _COOArray2D[np.float64]: ...
@overload  # sparray, format: "csc", dtype: float64-like
def hstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyFloat64DType) -> _CSCArray[np.float64]: ...
@overload  # sparray, format: "csr", dtype: float64-like
def hstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyFloat64DType) -> _CSRArray2D[np.float64]: ...
@overload  # sparray, format: "dia", dtype: float64-like
def hstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyFloat64DType) -> _DIAArray[np.float64]: ...
@overload  # sparray, format: "dok", dtype: float64-like
def hstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyFloat64DType) -> _DOKArray2D[np.float64]: ...
@overload  # sparray, format: "lil", dtype: float64-like
def hstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyFloat64DType) -> _LILArray[np.float64]: ...

#
@overload  # sparray, format: <default>, dtype: complex128-like
def hstack(blocks: Seq[_CanStackAs[np.complex128, _T]], format: None = None, *, dtype: onp.AnyComplex128DType) -> _T: ...
@overload  # sparray, format: "bsr", dtype: complex128-like
def hstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyComplex128DType) -> _BSRArray[np.complex128]: ...
@overload  # sparray, format: "coo", dtype: complex128-like
def hstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyComplex128DType) -> _COOArray2D[np.complex128]: ...
@overload  # sparray, format: "csc", dtype: complex128-like
def hstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyComplex128DType) -> _CSCArray[np.complex128]: ...
@overload  # sparray, format: "csr", dtype: complex128-like
def hstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyComplex128DType) -> _CSRArray2D[np.complex128]: ...
@overload  # sparray, format: "dia", dtype: complex128-like
def hstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyComplex128DType) -> _DIAArray[np.complex128]: ...
@overload  # sparray, format: "dok", dtype: complex128-like
def hstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyComplex128DType) -> _DOKArray2D[np.complex128]: ...
@overload  # sparray, format: "lil", dtype: complex128-like
def hstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyComplex128DType) -> _LILArray[np.complex128]: ...

#
@overload  # sparray, format: <default>, dtype: <known>
def hstack(blocks: Seq[_CanStackAs[_SCT0, _T]], format: None = None, *, dtype: onp.ToDType[_SCT0]) -> _T: ...
@overload  # sparray, format: "bsr", dtype: <known>
def hstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.ToDType[_SCT0]) -> _BSRArray[_SCT0]: ...
@overload  # sparray, format: "coo", dtype: <known>
def hstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.ToDType[_SCT0]) -> _COOArray2D[_SCT0]: ...
@overload  # sparray, format: "csc", dtype: <known>
def hstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.ToDType[_SCT0]) -> _CSCArray[_SCT0]: ...
@overload  # sparray, format: "csr", dtype: <known>
def hstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.ToDType[_SCT0]) -> _CSRArray2D[_SCT0]: ...
@overload  # sparray, format: "dia", dtype: <known>
def hstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.ToDType[_SCT0]) -> _DIAArray[_SCT0]: ...
@overload  # sparray, format: "dok", dtype: <known>
def hstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.ToDType[_SCT0]) -> _DOKArray2D[_SCT0]: ...
@overload  # sparray, format: "lil", dtype: <known>
def hstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.ToDType[_SCT0]) -> _LILArray[_SCT0]: ...

#
@overload  # sparray, format: <default>, dtype: <unknown>
def hstack(blocks: Seq[_CanStackAs[Any, _T]], format: None = None, *, dtype: _ToDType) -> _T: ...
@overload  # sparray, format: "bsr", dtype: <unknown>
def hstack(blocks: Seq[sparray], format: _FmtBSR, dtype: _ToDType) -> _BSRArray: ...
@overload  # sparray, format: "coo", dtype: <unknown>
def hstack(blocks: Seq[sparray], format: _FmtCOO, dtype: _ToDType) -> _COOArray2D: ...
@overload  # sparray, format: "csc", dtype: <unknown>
def hstack(blocks: Seq[sparray], format: _FmtCSC, dtype: _ToDType) -> _CSCArray: ...
@overload  # sparray, format: "csr", dtype: <unknown>
def hstack(blocks: Seq[sparray], format: _FmtCSR, dtype: _ToDType) -> _CSRArray2D: ...
@overload  # sparray, format: "dia", dtype: <unknown>
def hstack(blocks: Seq[sparray], format: _FmtDIA, dtype: _ToDType) -> _DIAArray: ...
@overload  # sparray, format: "dok", dtype: <unknown>
def hstack(blocks: Seq[sparray], format: _FmtDOK, dtype: _ToDType) -> _DOKArray2D: ...
@overload  # sparray, format: "lil", dtype: <unknown>
def hstack(blocks: Seq[sparray], format: _FmtLIL, dtype: _ToDType) -> _LILArray: ...

#
@overload
def hstack(blocks: Seq[_spbase], format: _Format, dtype: _ToDType | None = None) -> Incomplete: ...

###
# NOTE: keep in sync with `hstack`
@overload  # sparray, format: <default>, dtype: <default>
def vstack(blocks: Seq[_CanStack[_T]], format: None = None, dtype: None = None) -> _T: ...
@overload  # sparray, format: "bsr", dtype: <default>
def vstack(blocks: Seq[_SpArray2D[_SCT]], format: _FmtBSR, dtype: None = None) -> _BSRArray[_SCT]: ...
@overload  # sparray, format: "coo", dtype: <default>
def vstack(blocks: Seq[_SpArray2D[_SCT]], format: _FmtCOO, dtype: None = None) -> _COOArray2D[_SCT]: ...
@overload  # sparray, format: "csc", dtype: <default>
def vstack(blocks: Seq[_SpArray2D[_SCT]], format: _FmtCSC, dtype: None = None) -> _CSCArray[_SCT]: ...
@overload  # sparray, format: "csr", dtype: <default>
def vstack(blocks: Seq[_SpArray2D[_SCT]], format: _FmtCSR, dtype: None = None) -> _CSRArray2D[_SCT]: ...
@overload  # sparray, format: "dia", dtype: <default>
def vstack(blocks: Seq[_SpArray2D[_SCT]], format: _FmtDIA, dtype: None = None) -> _DIAArray[_SCT]: ...
@overload  # sparray, format: "dok", dtype: <default>
def vstack(blocks: Seq[_SpArray2D[_SCT]], format: _FmtDOK, dtype: None = None) -> _DOKArray2D[_SCT]: ...
@overload  # sparray, format: "lil", dtype: <default>
def vstack(blocks: Seq[_SpArray2D[_SCT]], format: _FmtLIL, dtype: None = None) -> _LILArray[_SCT]: ...

#
@overload  # sparray, format: <default>, dtype: bool-like
def vstack(blocks: Seq[_CanStackAs[np.bool_, _T]], format: None = None, *, dtype: onp.AnyBoolDType) -> _T: ...
@overload  # sparray, format: "bsr", dtype: bool-like
def vstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyBoolDType) -> _BSRArray[np.bool_]: ...
@overload  # sparray, format: "coo", dtype: bool-like
def vstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyBoolDType) -> _COOArray2D[np.bool_]: ...
@overload  # sparray, format: "csc", dtype: bool-like
def vstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyBoolDType) -> _CSCArray[np.bool_]: ...
@overload  # sparray, format: "csr", dtype: bool-like
def vstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyBoolDType) -> _CSRArray2D[np.bool_]: ...
@overload  # sparray, format: "dia", dtype: bool-like
def vstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyBoolDType) -> _DIAArray[np.bool_]: ...
@overload  # sparray, format: "dok", dtype: bool-like
def vstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyBoolDType) -> _DOKArray2D[np.bool_]: ...
@overload  # sparray, format: "lil", dtype: bool-like
def vstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyBoolDType) -> _LILArray[np.bool_]: ...

#
@overload  # sparray, format: <default>, dtype: int-like
def vstack(blocks: Seq[_CanStackAs[np.int_, _T]], format: None = None, *, dtype: onp.AnyIntDType) -> _T: ...
@overload  # sparray, format: "bsr", dtype: int-like
def vstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyIntDType) -> _BSRArray[np.int_]: ...
@overload  # sparray, format: "coo", dtype: int-like
def vstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyIntDType) -> _COOArray2D[np.int_]: ...
@overload  # sparray, format: "csc", dtype: int-like
def vstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyIntDType) -> _CSCArray[np.int_]: ...
@overload  # sparray, format: "csr", dtype: int-like
def vstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyIntDType) -> _CSRArray2D[np.int_]: ...
@overload  # sparray, format: "dia", dtype: int-like
def vstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyIntDType) -> _DIAArray[np.int_]: ...
@overload  # sparray, format: "dok", dtype: int-like
def vstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyIntDType) -> _DOKArray2D[np.int_]: ...
@overload  # sparray, format: "lil", dtype: int-like
def vstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyIntDType) -> _LILArray[np.int_]: ...

#
@overload  # sparray, format: <default>, dtype: float64-like
def vstack(blocks: Seq[_CanStackAs[np.float64, _T]], format: None = None, *, dtype: onp.AnyFloat64DType) -> _T: ...
@overload  # sparray, format: "bsr", dtype: float64-like
def vstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyFloat64DType) -> _BSRArray[np.float64]: ...
@overload  # sparray, format: "coo", dtype: float64-like
def vstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyFloat64DType) -> _COOArray2D[np.float64]: ...
@overload  # sparray, format: "csc", dtype: float64-like
def vstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyFloat64DType) -> _CSCArray[np.float64]: ...
@overload  # sparray, format: "csr", dtype: float64-like
def vstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyFloat64DType) -> _CSRArray2D[np.float64]: ...
@overload  # sparray, format: "dia", dtype: float64-like
def vstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyFloat64DType) -> _DIAArray[np.float64]: ...
@overload  # sparray, format: "dok", dtype: float64-like
def vstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyFloat64DType) -> _DOKArray2D[np.float64]: ...
@overload  # sparray, format: "lil", dtype: float64-like
def vstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyFloat64DType) -> _LILArray[np.float64]: ...

#
@overload  # sparray, format: <default>, dtype: complex128-like
def vstack(blocks: Seq[_CanStackAs[np.complex128, _T]], format: None = None, *, dtype: onp.AnyComplex128DType) -> _T: ...
@overload  # sparray, format: "bsr", dtype: complex128-like
def vstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyComplex128DType) -> _BSRArray[np.complex128]: ...
@overload  # sparray, format: "coo", dtype: complex128-like
def vstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyComplex128DType) -> _COOArray2D[np.complex128]: ...
@overload  # sparray, format: "csc", dtype: complex128-like
def vstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyComplex128DType) -> _CSCArray[np.complex128]: ...
@overload  # sparray, format: "csr", dtype: complex128-like
def vstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyComplex128DType) -> _CSRArray2D[np.complex128]: ...
@overload  # sparray, format: "dia", dtype: complex128-like
def vstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyComplex128DType) -> _DIAArray[np.complex128]: ...
@overload  # sparray, format: "dok", dtype: complex128-like
def vstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyComplex128DType) -> _DOKArray2D[np.complex128]: ...
@overload  # sparray, format: "lil", dtype: complex128-like
def vstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyComplex128DType) -> _LILArray[np.complex128]: ...

#
@overload  # sparray, format: <default>, dtype: <known>
def vstack(blocks: Seq[_CanStackAs[_SCT0, _T]], format: None = None, *, dtype: onp.ToDType[_SCT0]) -> _T: ...
@overload  # sparray, format: "bsr", dtype: <known>
def vstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.ToDType[_SCT0]) -> _BSRArray[_SCT0]: ...
@overload  # sparray, format: "coo", dtype: <known>
def vstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.ToDType[_SCT0]) -> _COOArray2D[_SCT0]: ...
@overload  # sparray, format: "csc", dtype: <known>
def vstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.ToDType[_SCT0]) -> _CSCArray[_SCT0]: ...
@overload  # sparray, format: "csr", dtype: <known>
def vstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.ToDType[_SCT0]) -> _CSRArray2D[_SCT0]: ...
@overload  # sparray, format: "dia", dtype: <known>
def vstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.ToDType[_SCT0]) -> _DIAArray[_SCT0]: ...
@overload  # sparray, format: "dok", dtype: <known>
def vstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.ToDType[_SCT0]) -> _DOKArray2D[_SCT0]: ...
@overload  # sparray, format: "lil", dtype: <known>
def vstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.ToDType[_SCT0]) -> _LILArray[_SCT0]: ...

#
@overload  # sparray, format: <default>, dtype: <unknown>
def vstack(blocks: Seq[_CanStackAs[Any, _T]], format: None = None, *, dtype: _ToDType) -> _T: ...
@overload  # sparray, format: "bsr", dtype: <unknown>
def vstack(blocks: Seq[sparray], format: _FmtBSR, dtype: _ToDType) -> _BSRArray: ...
@overload  # sparray, format: "coo", dtype: <unknown>
def vstack(blocks: Seq[sparray], format: _FmtCOO, dtype: _ToDType) -> _COOArray2D: ...
@overload  # sparray, format: "csc", dtype: <unknown>
def vstack(blocks: Seq[sparray], format: _FmtCSC, dtype: _ToDType) -> _CSCArray: ...
@overload  # sparray, format: "csr", dtype: <unknown>
def vstack(blocks: Seq[sparray], format: _FmtCSR, dtype: _ToDType) -> _CSRArray2D: ...
@overload  # sparray, format: "dia", dtype: <unknown>
def vstack(blocks: Seq[sparray], format: _FmtDIA, dtype: _ToDType) -> _DIAArray: ...
@overload  # sparray, format: "dok", dtype: <unknown>
def vstack(blocks: Seq[sparray], format: _FmtDOK, dtype: _ToDType) -> _DOKArray2D: ...
@overload  # sparray, format: "lil", dtype: <unknown>
def vstack(blocks: Seq[sparray], format: _FmtLIL, dtype: _ToDType) -> _LILArray: ...

#
@overload
def vstack(blocks: Seq[_spbase], format: _Format, dtype: _ToDType | None = None) -> Incomplete: ...

###
@overload  # blocks: <known, known>, format: <default>, dtype: <default>
def block_array(blocks: _ToBlocksCanStack[_T], *, format: None = None, dtype: None = None) -> _T: ...
@overload  # blocks: <array, known>, format: "bsr", dtype: <default>
def block_array(blocks: _ToBlocksSPArray[_SCT], *, format: _FmtBSR, dtype: None = None) -> _BSRArray[_SCT]: ...
@overload  # blocks: <array, known>, format: "coo", dtype: <default>
def block_array(blocks: _ToBlocksSPArray[_SCT], *, format: _FmtCOO, dtype: None = None) -> _COOArray2D[_SCT]: ...
@overload  # blocks: <array, known>, format: "csc", dtype: <default>
def block_array(blocks: _ToBlocksSPArray[_SCT], *, format: _FmtCSC, dtype: None = None) -> _CSCArray[_SCT]: ...
@overload  # blocks: <array, known>, format: "csr", dtype: <default>
def block_array(blocks: _ToBlocksSPArray[_SCT], *, format: _FmtCSR, dtype: None = None) -> _CSRArray2D[_SCT]: ...
@overload  # blocks: <array, known>, format: "dia", dtype: <default>
def block_array(blocks: _ToBlocksSPArray[_SCT], *, format: _FmtDIA, dtype: None = None) -> _DIAArray[_SCT]: ...
@overload  # blocks: <array, known>, format: "dok", dtype: <default>
def block_array(blocks: _ToBlocksSPArray[_SCT], *, format: _FmtDOK, dtype: None = None) -> _DOKArray2D[_SCT]: ...
@overload  # blocks: <array, known>, format: "lil", dtype: <default>
def block_array(blocks: _ToBlocksSPArray[_SCT], *, format: _FmtLIL, dtype: None = None) -> _LILArray[_SCT]: ...

#
@overload  # blocks: <known, bool_>, format: <default>, dtype: bool-like
def block_array(blocks: _ToBlocksCanStackAs[np.bool_, _T], *, format: None = None, dtype: onp.AnyBoolDType) -> _T: ...
@overload  # blocks: <unknown, unknown>, format: "bsr", dtype: bool-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtBSR, dtype: onp.AnyBoolDType) -> _BSRArray[np.bool_]: ...
@overload  # blocks: <unknown, unknown>, format: "coo", dtype: bool-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCOO, dtype: onp.AnyBoolDType) -> _COOArray2D[np.bool_]: ...
@overload  # blocks: <unknown, unknown>, format: "csc", dtype: bool-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSC, dtype: onp.AnyBoolDType) -> _CSCArray[np.bool_]: ...
@overload  # blocks: <unknown, unknown>, format: "csr", dtype: bool-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSR, dtype: onp.AnyBoolDType) -> _CSRArray2D[np.bool_]: ...
@overload  # blocks: <unknown, unknown>, format: "dia", dtype: bool-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDIA, dtype: onp.AnyBoolDType) -> _DIAArray[np.bool_]: ...
@overload  # blocks: <unknown, unknown>, format: "dok", dtype: bool-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDOK, dtype: onp.AnyBoolDType) -> _DOKArray2D[np.bool_]: ...
@overload  # blocks: <unknown, unknown>, format: "lil", dtype: bool-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtLIL, dtype: onp.AnyBoolDType) -> _LILArray[np.bool_]: ...

#
@overload  # blocks: <known, int_>, format: <default>, dtype: int-like
def block_array(blocks: _ToBlocksCanStackAs[np.int64, _T], *, format: None = None, dtype: onp.AnyIntDType) -> _T: ...
@overload  # blocks: <unknown, unknown>, format: "bsr", dtype: int-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtBSR, dtype: onp.AnyIntDType) -> _BSRArray[np.int_]: ...
@overload  # blocks: <unknown, unknown>, format: "coo", dtype: int-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCOO, dtype: onp.AnyIntDType) -> _COOArray2D[np.int_]: ...
@overload  # blocks: <unknown, unknown>, format: "csc", dtype: int-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSC, dtype: onp.AnyIntDType) -> _CSCArray[np.int_]: ...
@overload  # blocks: <unknown, unknown>, format: "csr", dtype: int-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSR, dtype: onp.AnyIntDType) -> _CSRArray2D[np.int_]: ...
@overload  # blocks: <unknown, unknown>, format: "dia", dtype: int-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDIA, dtype: onp.AnyIntDType) -> _DIAArray[np.int_]: ...
@overload  # blocks: <unknown, unknown>, format: "dok", dtype: int-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDOK, dtype: onp.AnyIntDType) -> _DOKArray2D[np.int_]: ...
@overload  # blocks: <unknown, unknown>, format: "lil", dtype: int-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtLIL, dtype: onp.AnyIntDType) -> _LILArray[np.int_]: ...

#
@overload  # blocks: <known, float64>, format: <default>, dtype: float64-like
def block_array(blocks: _ToBlocksCanStackAs[np.float64, _T], *, format: None = None, dtype: onp.AnyFloat64DType) -> _T: ...
@overload  # blocks: <unknown, unknown>, format: "bsr", dtype: float64-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtBSR, dtype: onp.AnyFloat64DType) -> _BSRArray[np.float64]: ...
@overload  # blocks: <unknown, unknown>, format: "coo", dtype: float64-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCOO, dtype: onp.AnyFloat64DType) -> _COOArray2D[np.float64]: ...
@overload  # blocks: <unknown, unknown>, format: "csc", dtype: float64-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSC, dtype: onp.AnyFloat64DType) -> _CSCArray[np.float64]: ...
@overload  # blocks: <unknown, unknown>, format: "csr", dtype: float64-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSR, dtype: onp.AnyFloat64DType) -> _CSRArray2D[np.float64]: ...
@overload  # blocks: <unknown, unknown>, format: "dia", dtype: float64-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDIA, dtype: onp.AnyFloat64DType) -> _DIAArray[np.float64]: ...
@overload  # blocks: <unknown, unknown>, format: "dok", dtype: float64-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDOK, dtype: onp.AnyFloat64DType) -> _DOKArray2D[np.float64]: ...
@overload  # blocks: <unknown, unknown>, format: "lil", dtype: float64-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtLIL, dtype: onp.AnyFloat64DType) -> _LILArray[np.float64]: ...

#
@overload  # blocks: <known, complex128>, format: <default>, dtype: complex128-like
def block_array(blocks: _ToBlocksCanStackAs[np.complex128, _T], *, format: None = None, dtype: onp.AnyComplex128DType) -> _T: ...
@overload  # blocks: <unknown, unknown>, format: "bsr", dtype: complex128-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtBSR, dtype: onp.AnyComplex128DType) -> _BSRArray[np.complex128]: ...
@overload  # blocks: <unknown, unknown>, format: "coo", dtype: complex128-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCOO, dtype: onp.AnyComplex128DType) -> _COOArray2D[np.complex128]: ...
@overload  # blocks: <unknown, unknown>, format: "csc", dtype: complex128-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSC, dtype: onp.AnyComplex128DType) -> _CSCArray[np.complex128]: ...
@overload  # blocks: <unknown, unknown>, format: "csr", dtype: complex128-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSR, dtype: onp.AnyComplex128DType) -> _CSRArray2D[np.complex128]: ...
@overload  # blocks: <unknown, unknown>, format: "dia", dtype: complex128-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDIA, dtype: onp.AnyComplex128DType) -> _DIAArray[np.complex128]: ...
@overload  # blocks: <unknown, unknown>, format: "dok", dtype: complex128-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDOK, dtype: onp.AnyComplex128DType) -> _DOKArray2D[np.complex128]: ...
@overload  # blocks: <unknown, unknown>, format: "lil", dtype: complex128-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtLIL, dtype: onp.AnyComplex128DType) -> _LILArray[np.complex128]: ...

#
@overload  # blocks: <known, known>, format: <default>, dtype: <known>
def block_array(blocks: _ToBlocksCanStackAs[_SCT0, _T], *, format: None = None, dtype: onp.ToDType[_SCT0]) -> _T: ...
@overload  # blocks: <unknown, unknown>, format: "bsr", dtype: <known>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtBSR, dtype: onp.ToDType[_SCT0]) -> _BSRArray[_SCT0]: ...
@overload  # blocks: <unknown, unknown>, format: "coo", dtype: <known>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCOO, dtype: onp.ToDType[_SCT0]) -> _COOArray2D[_SCT0]: ...
@overload  # blocks: <unknown, unknown>, format: "csc", dtype: <known>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSC, dtype: onp.ToDType[_SCT0]) -> _CSCArray[_SCT0]: ...
@overload  # blocks: <unknown, unknown>, format: "csr", dtype: <known>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSR, dtype: onp.ToDType[_SCT0]) -> _CSRArray2D[_SCT0]: ...
@overload  # blocks: <unknown, unknown>, format: "dia", dtype: <known>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDIA, dtype: onp.ToDType[_SCT0]) -> _DIAArray[_SCT0]: ...
@overload  # blocks: <unknown>, format: "dok", dtype: <known>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDOK, dtype: onp.ToDType[_SCT0]) -> _DOKArray2D[_SCT0]: ...
@overload  # blocks: <unknown, unknown>, format: "lil", dtype: <known>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtLIL, dtype: onp.ToDType[_SCT0]) -> _LILArray[_SCT0]: ...

#
@overload  # blocks: <known, unknown>, format: <default>, dtype: <unknown>
def block_array(blocks: _ToBlocksCanStackAs[Any, _T], *, format: None = None, dtype: _ToDType | None = None) -> _T: ...
@overload  # blocks: <unknown, unknown>, format: "bsr", dtype: <unknown>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtBSR, dtype: _ToDType | None = None) -> _BSRArray: ...
@overload  # blocks: <unknown, unknown>, format: "coo", dtype: <unknown>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCOO, dtype: _ToDType | None = None) -> _COOArray2D: ...
@overload  # blocks: <unknown, unknown>, format: "csc", dtype: <unknown>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSC, dtype: _ToDType | None = None) -> _CSCArray: ...
@overload  # blocks: <unknown, unknown>, format: "csr", dtype: <unknown>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSR, dtype: _ToDType | None = None) -> _CSRArray2D: ...
@overload  # blocks: <unknown, unknown>, format: "dia", dtype: <unknown>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDIA, dtype: _ToDType | None = None) -> _DIAArray: ...
@overload  # blocks: <unknown, unknown>, format: "dok", dtype: <unknown>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDOK, dtype: _ToDType | None = None) -> _DOKArray2D: ...
@overload  # blocks: <unknown, unknown>, format: "lil", dtype: <unknown>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtLIL, dtype: _ToDType | None = None) -> _LILArray: ...

###
@overload  # blocks: <known, known>, format: <default>, dtype: <default>
def bmat(blocks: _ToBlocksCanStack[_T], format: None = None, dtype: None = None) -> _T: ...
@overload  # blocks: <matrix, known>, format: <otherwise>, dtype: <default>
def bmat(blocks: _ToBlocksSPMatrix[_SCT], format: _Format, dtype: None = None) -> spmatrix[_SCT]: ...

#
@overload  # blocks: <known, known>, format: <default>, dtype: <known>
def bmat(blocks: _ToBlocksCanStackAs[_SCT0, _T], format: None = None, *, dtype: onp.ToDType[_SCT0]) -> _T: ...
@overload  # blocks: <matrix, known>, format: <otherwise>, dtype: <known>
def bmat(blocks: _ToBlocksSPMatrix[_SCT0], format: _Format, dtype: onp.ToDType[_SCT0]) -> spmatrix[_SCT0]: ...
@overload  # blocks: <unknown, unknown>, format: <otherwise>, dtype: <known>
def bmat(blocks: _ToBlocksUnkown, format: _Format, dtype: onp.ToDType[_SCT0]) -> spmatrix[_SCT0] | _SpArray2D[_SCT0]: ...

#
@overload  # blocks: <known, unknown>, format: <default>, dtype: <unknown>
def bmat(blocks: _ToBlocksCanStackAs[Any, _T], format: None = None, *, dtype: _ToDType) -> _T: ...
@overload  # blocks: <matrix, unknown>, format: <otherwise>, dtype: <unknown>
def bmat(blocks: _ToBlocksSPMatrix, format: _Format, *, dtype: _ToDType) -> spmatrix: ...
@overload  # blocks: <unknown, unknown>, format: <otherwise>, dtype: <unknown>
def bmat(blocks: _ToBlocksUnkown, format: _Format, *, dtype: _ToDType) -> spmatrix | _SpArray2D: ...

###
@overload  # mats: <array, known>, format: <default>, dtype: None
def block_diag(mats: Iterable[sparray[_SCT]], format: _FmtCOO | None = None, dtype: None = None) -> _COOArray2D[_SCT]: ...
@overload  # mats: <array, known>, format: "bsr", dtype: None
def block_diag(mats: Iterable[sparray[_SCT]], format: _FmtBSR, dtype: None = None) -> _BSRArray[_SCT]: ...
@overload  # mats: <array, known>, format: "csc", dtype: None
def block_diag(mats: Iterable[sparray[_SCT]], format: _FmtCSC, dtype: None = None) -> _CSCArray[_SCT]: ...
@overload  # mats: <array, known>, format: "csr", dtype: None
def block_diag(mats: Iterable[sparray[_SCT]], format: _FmtCSR, dtype: None = None) -> _CSRArray2D[_SCT]: ...
@overload  # mats: <array, known>, format: "dia", dtype: None
def block_diag(mats: Iterable[sparray[_SCT]], format: _FmtDIA, dtype: None = None) -> _DIAArray[_SCT]: ...
@overload  # mats: <array, known>, format: "dok", dtype: None
def block_diag(mats: Iterable[sparray[_SCT]], format: _FmtDOK, dtype: None = None) -> _DOKArray2D[_SCT]: ...
@overload  # mats: <array, known>, format: "lil", dtype: None
def block_diag(mats: Iterable[sparray[_SCT]], format: _FmtLIL, dtype: None = None) -> _LILArray[_SCT]: ...

#
@overload  # mats: <array, unknown>, format: <default>, dtype: bool-like
def block_diag(mats: Iterable[sparray], format: _FmtCOO | None = None, *, dtype: onp.AnyBoolDType) -> _COOArray2D[np.bool_]: ...
@overload  # mats: <array, unknown>, format: "bsr", dtype: bool-like
def block_diag(mats: Iterable[sparray], format: _FmtBSR, dtype: onp.AnyBoolDType) -> _BSRArray[np.bool_]: ...
@overload  # mats: <array, unknown>, format: "csc", dtype: bool-like
def block_diag(mats: Iterable[sparray], format: _FmtCSC, dtype: onp.AnyBoolDType) -> _CSCArray[np.bool_]: ...
@overload  # mats: <array, unknown>, format: "csr", dtype: bool-like
def block_diag(mats: Iterable[sparray], format: _FmtCSR, dtype: onp.AnyBoolDType) -> _CSRArray2D[np.bool_]: ...
@overload  # mats: <array, unknown>, format: "dia", dtype: bool-like
def block_diag(mats: Iterable[sparray], format: _FmtDIA, dtype: onp.AnyBoolDType) -> _DIAArray[np.bool_]: ...
@overload  # mats: <array, unknown>, format: "dok", dtype: bool-like
def block_diag(mats: Iterable[sparray], format: _FmtDOK, dtype: onp.AnyBoolDType) -> _DOKArray2D[np.bool_]: ...
@overload  # mats: <array, unknown>, format: "lil", dtype: bool-like
def block_diag(mats: Iterable[sparray], format: _FmtLIL, dtype: onp.AnyBoolDType) -> _LILArray[np.bool_]: ...

#
@overload  # mats: <array, unknown>, format: <default>, dtype: int-like
def block_diag(mats: Iterable[sparray], format: _FmtCOO | None = None, *, dtype: onp.AnyIntDType) -> _COOArray2D[np.int_]: ...
@overload  # mats: <array, unknown>, format: "bsr", dtype: int-like
def block_diag(mats: Iterable[sparray], format: _FmtBSR, dtype: onp.AnyIntDType) -> _BSRArray[np.int_]: ...
@overload  # mats: <array, unknown>, format: "csc", dtype: int-like
def block_diag(mats: Iterable[sparray], format: _FmtCSC, dtype: onp.AnyIntDType) -> _CSCArray[np.int_]: ...
@overload  # mats: <array, unknown>, format: "csr", dtype: int-like
def block_diag(mats: Iterable[sparray], format: _FmtCSR, dtype: onp.AnyIntDType) -> _CSRArray2D[np.int_]: ...
@overload  # mats: <array, unknown>, format: "dia", dtype: int-like
def block_diag(mats: Iterable[sparray], format: _FmtDIA, dtype: onp.AnyIntDType) -> _DIAArray[np.int_]: ...
@overload  # mats: <array, unknown>, format: "dok", dtype: int-like
def block_diag(mats: Iterable[sparray], format: _FmtDOK, dtype: onp.AnyIntDType) -> _DOKArray2D[np.int_]: ...
@overload  # mats: <array, unknown>, format: "lil", dtype: int-like
def block_diag(mats: Iterable[sparray], format: _FmtLIL, dtype: onp.AnyIntDType) -> _LILArray[np.int_]: ...

#
@overload  # mats: <array, unknown>, format: <default>, dtype: float64-like
def block_diag(
    mats: Iterable[sparray], format: _FmtCOO | None = None, *, dtype: onp.AnyFloat64DType
) -> _COOArray2D[np.float64]: ...
@overload  # mats: <array, unknown>, format: "bsr", dtype: float64-like
def block_diag(mats: Iterable[sparray], format: _FmtBSR, dtype: onp.AnyFloat64DType) -> _BSRArray[np.float64]: ...
@overload  # mats: <array, unknown>, format: "csc", dtype: float64-like
def block_diag(mats: Iterable[sparray], format: _FmtCSC, dtype: onp.AnyFloat64DType) -> _CSCArray[np.float64]: ...
@overload  # mats: <array, unknown>, format: "csr", dtype: float64-like
def block_diag(mats: Iterable[sparray], format: _FmtCSR, dtype: onp.AnyFloat64DType) -> _CSRArray2D[np.float64]: ...
@overload  # mats: <array, unknown>, format: "dia", dtype: float64-like
def block_diag(mats: Iterable[sparray], format: _FmtDIA, dtype: onp.AnyFloat64DType) -> _DIAArray[np.float64]: ...
@overload  # mats: <array, unknown>, format: "dok", dtype: float64-like
def block_diag(mats: Iterable[sparray], format: _FmtDOK, dtype: onp.AnyFloat64DType) -> _DOKArray2D[np.float64]: ...
@overload  # mats: <array, unknown>, format: "lil", dtype: float64-like
def block_diag(mats: Iterable[sparray], format: _FmtLIL, dtype: onp.AnyFloat64DType) -> _LILArray[np.float64]: ...

#
@overload  # mats: <array, unknown>, format: <default>, dtype: complex128-like
def block_diag(
    mats: Iterable[sparray], format: _FmtCOO | None = None, *, dtype: onp.AnyComplex128DType
) -> _COOArray2D[np.complex128]: ...
@overload  # mats: <array, unknown>, format: "bsr", dtype: complex128-like
def block_diag(mats: Iterable[sparray], format: _FmtBSR, dtype: onp.AnyComplex128DType) -> _BSRArray[np.complex128]: ...
@overload  # mats: <array, unknown>, format: "csc", dtype: complex128-like
def block_diag(mats: Iterable[sparray], format: _FmtCSC, dtype: onp.AnyComplex128DType) -> _CSCArray[np.complex128]: ...
@overload  # mats: <array, unknown>, format: "csr", dtype: complex128-like
def block_diag(mats: Iterable[sparray], format: _FmtCSR, dtype: onp.AnyComplex128DType) -> _CSRArray2D[np.complex128]: ...
@overload  # mats: <array, unknown>, format: "dia", dtype: complex128-like
def block_diag(mats: Iterable[sparray], format: _FmtDIA, dtype: onp.AnyComplex128DType) -> _DIAArray[np.complex128]: ...
@overload  # mats: <array, unknown>, format: "dok", dtype: complex128-like
def block_diag(mats: Iterable[sparray], format: _FmtDOK, dtype: onp.AnyComplex128DType) -> _DOKArray2D[np.complex128]: ...
@overload  # mats: <array, unknown>, format: "lil", dtype: complex128-like
def block_diag(mats: Iterable[sparray], format: _FmtLIL, dtype: onp.AnyComplex128DType) -> _LILArray[np.complex128]: ...

#
@overload  # mats: <array, unknown>, format: <default>, dtype: <known>
def block_diag(mats: Iterable[sparray], format: _FmtCOO | None = None, *, dtype: onp.ToDType[_SCT]) -> _COOArray2D[_SCT]: ...
@overload  # mats: <array, unknown>, format: "bsr", dtype: <known>
def block_diag(mats: Iterable[sparray], format: _FmtBSR, dtype: onp.ToDType[_SCT]) -> _BSRArray[_SCT]: ...
@overload  # mats: <array, unknown>, format: "csc", dtype: <known>
def block_diag(mats: Iterable[sparray], format: _FmtCSC, dtype: onp.ToDType[_SCT]) -> _CSCArray[_SCT]: ...
@overload  # mats: <array, unknown>, format: "csr", dtype: <known>
def block_diag(mats: Iterable[sparray], format: _FmtCSR, dtype: onp.ToDType[_SCT]) -> _CSRArray2D[_SCT]: ...
@overload  # mats: <array, unknown>, format: "dia", dtype: <known>
def block_diag(mats: Iterable[sparray], format: _FmtDIA, dtype: onp.ToDType[_SCT]) -> _DIAArray[_SCT]: ...
@overload  # mats: <array, unknown>, format: "dok", dtype: <known>
def block_diag(mats: Iterable[sparray], format: _FmtDOK, dtype: onp.ToDType[_SCT]) -> _DOKArray2D[_SCT]: ...
@overload  # mats: <array, unknown>, format: "lil", dtype: <known>
def block_diag(mats: Iterable[sparray], format: _FmtLIL, dtype: onp.ToDType[_SCT]) -> _LILArray[_SCT]: ...

#
@overload  # mats: <unknown, known>, format: <default>, dtype: None
def block_diag(
    mats: _ToMatsDiag[_SCT], format: _FmtCOO | None = None, dtype: None = None
) -> _COOArray2D[_SCT] | coo_matrix[_SCT]: ...
@overload  # mats: <unknown, known>, format: "bsr", dtype: None
def block_diag(mats: _ToMatsDiag[_SCT], format: _FmtBSR, dtype: None = None) -> _BSRArray[_SCT] | bsr_matrix[_SCT]: ...
@overload  # mats: <unknown, known>, format: "csc", dtype: None
def block_diag(mats: _ToMatsDiag[_SCT], format: _FmtCSC, dtype: None = None) -> _CSCArray[_SCT] | csc_matrix[_SCT]: ...
@overload  # mats: <unknown, known>, format: "csr", dtype: None
def block_diag(mats: _ToMatsDiag[_SCT], format: _FmtCSR, dtype: None = None) -> _CSRArray2D[_SCT] | csr_matrix[_SCT]: ...
@overload  # mats: <unknown, known>, format: "dia", dtype: None
def block_diag(mats: _ToMatsDiag[_SCT], format: _FmtDIA, dtype: None = None) -> _DIAArray[_SCT] | dia_matrix[_SCT]: ...
@overload  # mats: <unknown, known>, format: "dok", dtype: None
def block_diag(mats: _ToMatsDiag[_SCT], format: _FmtDOK, dtype: None = None) -> _DOKArray2D[_SCT] | dok_matrix[_SCT]: ...
@overload  # mats: <unknown, known>, format: "lil", dtype: None
def block_diag(mats: _ToMatsDiag[_SCT], format: _FmtLIL, dtype: None = None) -> _LILArray[_SCT] | lil_matrix[_SCT]: ...

#
@overload  # mats: <unknown, unknown>, format: <default>, dtype: <known>
def block_diag(
    mats: _ToMatsDiagUnknown, format: _FmtCOO | None = None, *, dtype: onp.ToDType[_SCT]
) -> _COOArray2D[_SCT] | coo_matrix[_SCT]: ...
@overload  # mats: <unknown, unknown>, format: "bsr", dtype: <known>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtBSR, dtype: onp.ToDType[_SCT]) -> _BSRArray[_SCT] | bsr_matrix[_SCT]: ...
@overload  # mats: <unknown, unknown>, format: "csc", dtype: <known>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtCSC, dtype: onp.ToDType[_SCT]) -> _CSCArray[_SCT] | csc_matrix[_SCT]: ...
@overload  # mats: <unknown, unknown>, format: "csr", dtype: <known>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtCSR, dtype: onp.ToDType[_SCT]) -> _CSRArray2D[_SCT] | csr_matrix[_SCT]: ...
@overload  # mats: <unknown, unknown>, format: "dia", dtype: <known>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtDIA, dtype: onp.ToDType[_SCT]) -> _DIAArray[_SCT] | dia_matrix[_SCT]: ...
@overload  # mats: <unknown, unknown>, format: "dok", dtype: <known>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtDOK, dtype: onp.ToDType[_SCT]) -> _DOKArray2D[_SCT] | dok_matrix[_SCT]: ...
@overload  # mats: <unknown, unknown>, format: "lil", dtype: <known>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtLIL, dtype: onp.ToDType[_SCT]) -> _LILArray[_SCT] | lil_matrix[_SCT]: ...

#
@overload  # mats: <unknown, unknown>, format: <default>, dtype: <unknown>
def block_diag(
    mats: _ToMatsDiagUnknown, format: _FmtCOO | None = None, dtype: _ToDType | None = None
) -> _COOArray2D | coo_matrix: ...
@overload  # mats: <unknown, unknown>, format: "bsr", dtype: <unknown>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtBSR, dtype: _ToDType | None = None) -> _BSRArray | bsr_matrix: ...
@overload  # mats: <unknown, unknown>, format: "csc", dtype: <unknown>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtCSC, dtype: _ToDType | None = None) -> _CSCArray | csc_matrix: ...
@overload  # mats: <unknown, unknown>, format: "csr", dtype: <unknown>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtCSR, dtype: _ToDType | None = None) -> _CSRArray2D | csr_matrix: ...
@overload  # mats: <unknown, unknown>, format: "dia", dtype: <unknown>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtDIA, dtype: _ToDType | None = None) -> _DIAArray | dia_matrix: ...
@overload  # mats: <unknown, unknown>, format: "dok", dtype: <unknown>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtDOK, dtype: _ToDType | None = None) -> _DOKArray2D | dok_matrix: ...
@overload  # mats: <unknown, unknown>, format: "lil", dtype: <unknown>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtLIL, dtype: _ToDType | None = None) -> _LILArray | lil_matrix: ...

###
@overload  # shape: T, format: <default>, dtype: <default>
def random_array(
    shape: _ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    dtype: onp.AnyFloat64DType | None = None,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> coo_array[np.float64, _ShapeT]: ...
@overload  # shape: T, format: <otherwise>, dtype: <default>
def random_array(
    shape: _ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtNonCOO,
    dtype: onp.AnyFloat64DType | None = None,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> sparray[np.float64, _ShapeT]: ...

#
@overload  # shape: T, format: <default>, dtype: <known>
def random_array(
    shape: _ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    dtype: onp.ToDType[_SCT],
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> coo_array[_SCT, _ShapeT]: ...
@overload  # shape: T, format: <otherwise>, dtype: <known>
def random_array(
    shape: _ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtNonCOO,
    dtype: onp.ToDType[_SCT],
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> sparray[_SCT, _ShapeT]: ...

#
@overload  # shape: T, format: <default>, dtype: complex
def random_array(
    shape: _ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    dtype: onp.AnyComplex128DType,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> coo_array[np.complex128, _ShapeT]: ...
@overload  # shape: T, format: <otherwise>, dtype: complex
def random_array(
    shape: _ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtNonCOO,
    dtype: onp.AnyComplex128DType,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> sparray[np.complex128, _ShapeT]: ...

#
@overload  # shape: T, format: <default>, dtype: <unknown>
def random_array(
    shape: _ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    dtype: _ToDType,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> coo_array[Any, _ShapeT]: ...
@overload  # shape: T, format: <otherwise>, dtype: <unknown>
def random_array(
    shape: _ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtNonCOO,
    dtype: _ToDType,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> sparray[Any, _ShapeT]: ...

###
# NOTE: `random_array` should be prefered over `random`
@overload  # format: <default>, dtype: <default>
def random(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    dtype: onp.AnyFloat64DType | None = None,
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[np.float64]: ...
@overload  # format: <otherwise>, dtype: <default>
def random(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    *,
    format: _FmtNonCOO,
    dtype: onp.AnyFloat64DType | None = None,
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[np.float64]: ...

#
@overload  # format: <default>, dtype: <known> (keyword)
def random(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    *,
    dtype: onp.ToDType[_SCT],
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[_SCT]: ...
@overload  # format: <otherwise>, dtype: <known> (keyword)
def random(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    *,
    format: _FmtNonCOO,
    dtype: onp.ToDType[_SCT],
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[_SCT]: ...

#
@overload  # format: <default>, dtype: <known> (positional)
def random(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating,
    format: _FmtCOO,
    dtype: onp.ToDType[_SCT],
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[_SCT]: ...
@overload  # format: <otherwise>, dtype: <known> (positional)
def random(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating,
    format: _FmtNonCOO,
    dtype: onp.ToDType[_SCT],
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[_SCT]: ...

#
@overload  # format: <default>, dtype: complex (keyword)
def random(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    *,
    dtype: onp.AnyComplex128DType,
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[np.complex128]: ...
@overload  # format: <otherwise>, dtype: complex (keyword)
def random(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    *,
    format: _FmtNonCOO,
    dtype: onp.AnyComplex128DType,
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[np.complex128]: ...

#
@overload  # format: <default>, dtype: complex (positional)
def random(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating,
    format: _FmtCOO,
    dtype: onp.AnyComplex128DType,
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[np.complex128]: ...
@overload  # format: <otherwise>, dtype: complex (positional)
def random(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating,
    format: _FmtNonCOO,
    dtype: onp.AnyComplex128DType,
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[np.complex128]: ...

#
@overload  # format: <default>, dtype: <unknown>
def random(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    dtype: _ToDType | None = None,
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix: ...
@overload  # format: <otherwise> (keyword), dtype: <unknown>
def random(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    *,
    format: _FmtNonCOO,
    dtype: _ToDType | None = None,
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix: ...
@overload  # format: <otherwise> (positional), dtype: <unknown>
def random(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating,
    format: _FmtNonCOO,
    dtype: _ToDType | None = None,
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix: ...

###
# NOTE: `random_array` should be prefered over `rand`
@overload  # format: <default>, dtype: <default>
def rand(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    dtype: onp.AnyFloat64DType | None = None,
    rng: onp.random.ToRNG | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[np.float64]: ...
@overload  # format: <otherwise>, dtype: <default>
def rand(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    *,
    format: _FmtNonCOO,
    dtype: onp.AnyFloat64DType | None = None,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[np.float64]: ...

#
@overload  # format: <default>, dtype: <known> (keyword)
def rand(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    *,
    dtype: onp.ToDType[_SCT],
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[_SCT]: ...
@overload  # format: <otherwise>, dtype: <known> (keyword)
def rand(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    *,
    format: _FmtNonCOO,
    dtype: onp.ToDType[_SCT],
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[_SCT]: ...

#
@overload  # format: <default>, dtype: <known> (positional)
def rand(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating,
    format: _FmtCOO,
    dtype: onp.ToDType[_SCT],
    rng: onp.random.ToRNG | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[_SCT]: ...
@overload  # format: <otherwise>, dtype: <known> (positional)
def rand(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating,
    format: _FmtNonCOO,
    dtype: onp.ToDType[_SCT],
    rng: onp.random.ToRNG | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[_SCT]: ...

#
@overload  # format: <default>, dtype: complex (keyword)
def rand(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    *,
    dtype: onp.AnyComplex128DType,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[np.complex128]: ...
@overload  # format: <otherwise>, dtype: complex (keyword)
def rand(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    *,
    format: _FmtNonCOO,
    dtype: onp.AnyComplex128DType,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[np.complex128]: ...

#
@overload  # format: <default>, dtype: complex (positional)
def rand(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating,
    format: _FmtCOO,
    dtype: onp.AnyComplex128DType,
    rng: onp.random.ToRNG | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[np.complex128]: ...
@overload  # format: <otherwise>, dtype: complex (positional)
def rand(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating,
    format: _FmtNonCOO,
    dtype: onp.AnyComplex128DType,
    rng: onp.random.ToRNG | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[np.complex128]: ...

#
@overload  # format: <default>, dtype: <unknown>
def rand(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    dtype: _ToDType | None = None,
    rng: onp.random.ToRNG | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix: ...
@overload  # format: <otherwise> (keyword), dtype: <unknown>
def rand(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating = 0.01,
    *,
    format: _FmtNonCOO,
    dtype: _ToDType | None = None,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix: ...
@overload  # format: <otherwise> (positional), dtype: <unknown>
def rand(
    m: opt.AnyInt,
    n: opt.AnyInt,
    density: float | npc.floating,
    format: _FmtNonCOO,
    dtype: _ToDType | None = None,
    rng: onp.random.ToRNG | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix: ...
