from _typeshed import Incomplete
from collections.abc import Callable, Iterable, Sequence as Seq
from typing import Any, Literal, Never, Protocol, overload, type_check_only

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._base import _spbase, sparray
from ._bsr import bsr_array, bsr_matrix
from ._coo import coo_array, coo_matrix
from ._csc import csc_array, csc_matrix
from ._csr import csr_array, csr_matrix
from ._dia import dia_array, dia_matrix
from ._dok import dok_array, dok_matrix
from ._lil import lil_array, lil_matrix
from ._matrix import spmatrix
from ._typing import _CanStack, _CanStackAs, _Format

__all__ = [
    "block_array",
    "block_diag",
    "bmat",
    "diags",
    "diags_array",
    "expand_dims",
    "eye",
    "eye_array",
    "hstack",
    "identity",
    "kron",
    "kronsum",
    "matrix_transpose",
    "permute_dims",
    "rand",
    "random",
    "random_array",
    "spdiags",
    "swapaxes",
    "vstack",
]

###

type _Numeric = npc.number | np.bool

type _ToDType = type[complex | _Numeric] | np.dtype[_Numeric] | str

type _SpArray2D[ScalarT: _Numeric] = sparray[ScalarT, tuple[int, int]]
type _COOArray2D[ScalarT: _Numeric] = coo_array[ScalarT, tuple[int, int]]
type _CSRArray2D[ScalarT: _Numeric] = csr_array[ScalarT, tuple[int, int]]
type _DOKArray2D[ScalarT: _Numeric] = dok_array[ScalarT, tuple[int, int]]

type _ToArray1D2D[ScalarT: _Numeric] = onp.CanArray[tuple[int] | tuple[int, int], np.dtype[ScalarT]] | Seq[ScalarT | Seq[ScalarT]]
type _ToSpMatrix[ScalarT: _Numeric] = spmatrix[ScalarT] | onp.ToArray2D[Never, ScalarT]
type _ToSparse2D[ScalarT: _Numeric] = _spbase[ScalarT, tuple[int, int]] | onp.ToArray2D[Never, ScalarT]

type _FmtBSR = Literal["bsr"]
type _FmtCOO = Literal["coo"]
type _FmtCSC = Literal["csc"]
type _FmtCSR = Literal["csr"]
type _FmtDIA = Literal["dia"]
type _FmtDOK = Literal["dok"]
type _FmtLIL = Literal["lil"]
type _FmtNonCOO = Literal["bsr", "csc", "csr", "dia", "dok", "lil"]

# TODO(julvandenbroeck): find a way to separate float and complex
type _ComplexSeq1D2D = Seq[Seq[complex] | complex]
type _ToComplex1D2D = onp.CanArray[tuple[int] | tuple[int, int], np.dtype[_Numeric]] | _ComplexSeq1D2D
type _Offsets = int | Seq[int] | onp.Array1D[npc.integer]

type _DataRVS = Callable[[int], onp.ArrayND[_Numeric]]

type _ToBlocks[T] = Seq[Seq[T | None]]
type _ToBlocksSpArray[ScalarT: _Numeric] = _ToBlocks[_SpArray2D[ScalarT]]
type _ToBlocksCanStackAs[ScalarT: _Numeric, T] = _ToBlocks[_CanStackAs[ScalarT, T]]
type _ToBlocksUnkown = _ToBlocksSpArray[_Numeric] | onp.ArrayND[np.object_]

type _ToMatsDiag[ScalarT: _Numeric] = Iterable[_spbase[ScalarT] | onp.ToArrayND[ScalarT]]
type _ToMatsDiagUnknown = Iterable[_spbase | onp.ArrayND[_Numeric] | complex | Seq[onp.ToComplex] | Seq[onp.ToComplex1D]]

@type_check_only
class _DataSampler(Protocol):
    def __call__(self, /, *, size: int) -> onp.ArrayND[_Numeric]: ...

###
#

@overload  # nasty workaround for https://github.com/microsoft/pyright/issues/10232
def expand_dims[ScalarT: _Numeric](  # type: ignore[overload-overlap]
    A: sparray[ScalarT, tuple[Never] | tuple[Never, Never]], /, *, axis: int = 0
) -> coo_array[ScalarT, tuple[int, int, *tuple[Any, ...]]]: ...
@overload
def expand_dims[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int]], /, *, axis: int = 0
) -> coo_array[ScalarT, tuple[int, int]]: ...
@overload
def expand_dims[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], /, *, axis: int = 0
) -> coo_array[ScalarT, tuple[int, int, int]]: ...
@overload
def expand_dims[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int, int]], /, *, axis: int = 0
) -> coo_array[ScalarT, tuple[int, int, int, int]]: ...
@overload  # TODO(@jorenham): shape-typing for coo_matrix
def expand_dims[ScalarT: _Numeric](A: spmatrix[ScalarT], /, *, axis: int = 0) -> coo_matrix[ScalarT]: ...

#
@overload
def swapaxes[ScalarT: _Numeric, ShapeT: tuple[int, *tuple[int, ...]]](
    A: sparray[ScalarT, ShapeT], axis1: int, axis2: int
) -> coo_array[ScalarT, ShapeT]: ...
@overload  # TODO(@jorenham): shape-typing for coo_matrix
def swapaxes[ScalarT: _Numeric](A: spmatrix[ScalarT], axis1: int, axis2: int) -> coo_matrix[ScalarT]: ...

#
@overload
def permute_dims[ScalarT: _Numeric, ShapeT: tuple[int, *tuple[int, ...]]](
    A: sparray[ScalarT, ShapeT], axes: Seq[int] | None = None, copy: bool = False
) -> coo_array[ScalarT, ShapeT]: ...
@overload  # TODO(@jorenham): shape-typing for coo_matrix
def permute_dims[ScalarT: _Numeric](
    A: spmatrix[ScalarT], axes: Seq[int] | None = None, copy: bool = False
) -> coo_matrix[ScalarT]: ...

###
@overload  # diagonals: <complex>, format: "dia" | None, dtype: None
def diags_array(
    diagonals: _ComplexSeq1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDIA | None = None,
    dtype: op.JustObject | None = ...,
) -> dia_array[np.float64] | dia_array[np.complex128]: ...
@overload  # diagonals: <complex>, format: "bsr", dtype: None
def diags_array(
    diagonals: _ComplexSeq1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtBSR,
    dtype: op.JustObject | None = ...,
) -> bsr_array[np.float64] | bsr_array[np.complex128]: ...
@overload  # diagonals: <complex>, format: "coo", dtype: None
def diags_array(
    diagonals: _ComplexSeq1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCOO,
    dtype: op.JustObject | None = ...,
) -> _COOArray2D[np.float64] | _COOArray2D[np.complex128]: ...
@overload  # diagonals: <complex>, format: "csc", dtype: None
def diags_array(
    diagonals: _ComplexSeq1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCSC,
    dtype: op.JustObject | None = ...,
) -> csc_array[np.float64] | csc_array[np.complex128]: ...
@overload  # diagonals: <complex>, format: "csr", dtype: None
def diags_array(
    diagonals: _ComplexSeq1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCSR,
    dtype: op.JustObject | None = ...,
) -> _CSRArray2D[np.float64] | _CSRArray2D[np.complex128]: ...
@overload  # diagonals: <complex>, format: "dok", dtype: None
def diags_array(
    diagonals: _ComplexSeq1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDOK,
    dtype: op.JustObject | None = ...,
) -> _DOKArray2D[np.float64] | _DOKArray2D[np.complex128]: ...
@overload  # diagonals: <complex>, format: "lil", dtype: None
def diags_array(
    diagonals: _ComplexSeq1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtLIL,
    dtype: op.JustObject | None = ...,
) -> lil_array[np.float64] | lil_array[np.complex128]: ...

#
@overload  # diagonals: <known>, format: "dia" | None, dtype: None
def diags_array[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[ScalarT],
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDIA | None = None,
    dtype: op.JustObject | None = ...,
) -> dia_array[ScalarT]: ...
@overload  # diagonals: <known>, format: "bsr", dtype: None
def diags_array[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[ScalarT],
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtBSR,
    dtype: op.JustObject | None = ...,
) -> bsr_array[ScalarT]: ...
@overload  # diagonals: <known>, format: "coo", dtype: None
def diags_array[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[ScalarT],
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCOO,
    dtype: op.JustObject | None = ...,
) -> _COOArray2D[ScalarT]: ...
@overload  # diagonals: <known>, format: "csc", dtype: None
def diags_array[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[ScalarT],
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCSC,
    dtype: op.JustObject | None = ...,
) -> csc_array[ScalarT]: ...
@overload  # diagonals: <known>, format: "csr", dtype: None
def diags_array[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[ScalarT],
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCSR,
    dtype: op.JustObject | None = ...,
) -> _CSRArray2D[ScalarT]: ...
@overload  # diagonals: <known>, format: "dok", dtype: None
def diags_array[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[ScalarT],
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDOK,
    dtype: op.JustObject | None = ...,
) -> _DOKArray2D[ScalarT]: ...
@overload  # diagonals: <known>, format: "lil", dtype: None
def diags_array[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[ScalarT],
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtLIL,
    dtype: op.JustObject | None = ...,
) -> lil_array[ScalarT]: ...

#
@overload  # diagonals: <unknown>, format: "dia" | None, dtype: bool-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDIA | None = None,
    dtype: onp.AnyBoolDType,
) -> dia_array[np.bool]: ...
@overload  # diagonals: <unknown>, format: "bsr", dtype: bool-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtBSR,
    dtype: onp.AnyBoolDType,
) -> bsr_array[np.bool]: ...
@overload  # diagonals: <unknown>, format: "coo", dtype: bool-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCOO,
    dtype: onp.AnyBoolDType,
) -> _COOArray2D[np.bool]: ...
@overload  # diagonals: <unknown>, format: "csc", dtype: bool-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCSC,
    dtype: onp.AnyBoolDType,
) -> csc_array[np.bool]: ...
@overload  # diagonals: <unknown>, format: "csr", dtype: bool-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCSR,
    dtype: onp.AnyBoolDType,
) -> _CSRArray2D[np.bool]: ...
@overload  # diagonals: <unknown>, format: "dok", dtype: bool-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDOK,
    dtype: onp.AnyBoolDType,
) -> _DOKArray2D[np.bool]: ...
@overload  # diagonals: <unknown>, format: "lil", dtype: bool-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtLIL,
    dtype: onp.AnyBoolDType,
) -> lil_array[np.bool]: ...

#
@overload  # diagonals: <unknown>, format: "dia" | None, dtype: int-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDIA | None = None,
    dtype: onp.AnyIntDType,
) -> dia_array[np.int_]: ...
@overload  # diagonals: <unknown>, format: "bsr", dtype: int-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtBSR,
    dtype: onp.AnyIntDType,
) -> bsr_array[np.int_]: ...
@overload  # diagonals: <unknown>, format: "coo", dtype: int-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCOO,
    dtype: onp.AnyIntDType,
) -> _COOArray2D[np.int_]: ...
@overload  # diagonals: <unknown>, format: "csc", dtype: int-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCSC,
    dtype: onp.AnyIntDType,
) -> csc_array[np.int_]: ...
@overload  # diagonals: <unknown>, format: "csr", dtype: int-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCSR,
    dtype: onp.AnyIntDType,
) -> _CSRArray2D[np.int_]: ...
@overload  # diagonals: <unknown>, format: "dok", dtype: int-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDOK,
    dtype: onp.AnyIntDType,
) -> _DOKArray2D[np.int_]: ...
@overload  # diagonals: <unknown>, format: "lil", dtype: int-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtLIL,
    dtype: onp.AnyIntDType,
) -> lil_array[np.int_]: ...

#
@overload  # diagonals: <unknown>, format: "dia" | None, dtype: float64-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDIA | None = None,
    dtype: onp.AnyFloat64DType,
) -> dia_array[np.float64]: ...
@overload  # diagonals: <unknown>, format: "bsr", dtype: float64-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtBSR,
    dtype: onp.AnyFloat64DType,
) -> bsr_array[np.float64]: ...
@overload  # diagonals: <unknown>, format: "coo", dtype: float64-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCOO,
    dtype: onp.AnyFloat64DType,
) -> _COOArray2D[np.float64]: ...
@overload  # diagonals: <unknown>, format: "csc", dtype: float64-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCSC,
    dtype: onp.AnyFloat64DType,
) -> csc_array[np.float64]: ...
@overload  # diagonals: <unknown>, format: "csr", dtype: float64-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCSR,
    dtype: onp.AnyFloat64DType,
) -> _CSRArray2D[np.float64]: ...
@overload  # diagonals: <unknown>, format: "dok", dtype: float64-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDOK,
    dtype: onp.AnyFloat64DType,
) -> _DOKArray2D[np.float64]: ...
@overload  # diagonals: <unknown>, format: "lil", dtype: float64-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtLIL,
    dtype: onp.AnyFloat64DType,
) -> lil_array[np.float64]: ...

#
@overload  # diagonals: <unknown>, format: "dia" | None, dtype: complex128-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDIA | None = None,
    dtype: onp.AnyComplex128DType,
) -> dia_array[np.complex128]: ...
@overload  # diagonals: <unknown>, format: "bsr", dtype: complex128-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtBSR,
    dtype: onp.AnyComplex128DType,
) -> bsr_array[np.complex128]: ...
@overload  # diagonals: <unknown>, format: "coo", dtype: complex128-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCOO,
    dtype: onp.AnyComplex128DType,
) -> _COOArray2D[np.complex128]: ...
@overload  # diagonals: <unknown>, format: "csc", dtype: complex128-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCSC,
    dtype: onp.AnyComplex128DType,
) -> csc_array[np.complex128]: ...
@overload  # diagonals: <unknown>, format: "csr", dtype: complex128-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCSR,
    dtype: onp.AnyComplex128DType,
) -> _CSRArray2D[np.complex128]: ...
@overload  # diagonals: <unknown>, format: "dok", dtype: complex128-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDOK,
    dtype: onp.AnyComplex128DType,
) -> _DOKArray2D[np.complex128]: ...
@overload  # diagonals: <unknown>, format: "lil", dtype: complex128-like
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtLIL,
    dtype: onp.AnyComplex128DType,
) -> lil_array[np.complex128]: ...

#
@overload  # diagonals: <unknown>, format: "dia" | None, dtype: <known>
def diags_array[ScalarT: _Numeric](
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDIA | None = None,
    dtype: onp.ToDType[ScalarT],
) -> dia_array[ScalarT]: ...
@overload  # diagonals: <unknown>, format: "bsr", dtype: <known>
def diags_array[ScalarT: _Numeric](
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtBSR,
    dtype: onp.ToDType[ScalarT],
) -> bsr_array[ScalarT]: ...
@overload  # diagonals: <unknown>, format: "coo", dtype: <known>
def diags_array[ScalarT: _Numeric](
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCOO,
    dtype: onp.ToDType[ScalarT],
) -> _COOArray2D[ScalarT]: ...
@overload  # diagonals: <unknown>, format: "csc", dtype: <known>
def diags_array[ScalarT: _Numeric](
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCSC,
    dtype: onp.ToDType[ScalarT],
) -> csc_array[ScalarT]: ...
@overload  # diagonals: <unknown>, format: "csr", dtype: <known>
def diags_array[ScalarT: _Numeric](
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtCSR,
    dtype: onp.ToDType[ScalarT],
) -> _CSRArray2D[ScalarT]: ...
@overload  # diagonals: <unknown>, format: "dok", dtype: <known>
def diags_array[ScalarT: _Numeric](
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDOK,
    dtype: onp.ToDType[ScalarT],
) -> _DOKArray2D[ScalarT]: ...
@overload  # diagonals: <unknown>, format: "lil", dtype: <known>
def diags_array[ScalarT: _Numeric](
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtLIL,
    dtype: onp.ToDType[ScalarT],
) -> lil_array[ScalarT]: ...

#
@overload  # diagonals: <unknown>, format: "dia" | None, dtype: <unknown>
def diags_array(
    diagonals: _ToComplex1D2D,
    /,
    *,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDIA | None = None,
    dtype: _ToDType,
) -> dia_array: ...
@overload  # diagonals: <unknown>, format: "bsr", dtype: <unknown>
def diags_array(
    diagonals: _ToComplex1D2D, /, *, offsets: _Offsets = 0, shape: tuple[int, int] | None = None, format: _FmtBSR, dtype: _ToDType
) -> bsr_array: ...
@overload  # diagonals: <unknown>, format: "coo", dtype: <unknown>
def diags_array(
    diagonals: _ToComplex1D2D, /, *, offsets: _Offsets = 0, shape: tuple[int, int] | None = None, format: _FmtCOO, dtype: _ToDType
) -> _COOArray2D[Any]: ...
@overload  # diagonals: <unknown>, format: "csc", dtype: <unknown>
def diags_array(
    diagonals: _ToComplex1D2D, /, *, offsets: _Offsets = 0, shape: tuple[int, int] | None = None, format: _FmtCSC, dtype: _ToDType
) -> csc_array: ...
@overload  # diagonals: <unknown>, format: "csr", dtype: <unknown>
def diags_array(
    diagonals: _ToComplex1D2D, /, *, offsets: _Offsets = 0, shape: tuple[int, int] | None = None, format: _FmtCSR, dtype: _ToDType
) -> _CSRArray2D[Any]: ...
@overload  # diagonals: <unknown>, format: "dok", dtype: <unknown>
def diags_array(
    diagonals: _ToComplex1D2D, /, *, offsets: _Offsets = 0, shape: tuple[int, int] | None = None, format: _FmtDOK, dtype: _ToDType
) -> _DOKArray2D[Any]: ...
@overload  # diagonals: <unknown>, format: "lil", dtype: <unknown>
def diags_array(
    diagonals: _ToComplex1D2D, /, *, offsets: _Offsets = 0, shape: tuple[int, int] | None = None, format: _FmtLIL, dtype: _ToDType
) -> lil_array: ...

###
# NOTE: `diags_array` should be prefered over `diags`
@overload  # diagonals: <known>, format: "dia" | None, dtype: None
def diags[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[ScalarT],
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDIA | None = None,
    dtype: op.JustObject | None = ...,
) -> dia_matrix[ScalarT]: ...
@overload  # diagonals: <known>, format: "bsr", dtype: None
def diags[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[ScalarT],
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtBSR,
    dtype: op.JustObject | None = ...,
) -> bsr_matrix[ScalarT]: ...
@overload  # diagonals: <known>, format: "coo", dtype: None
def diags[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[ScalarT],
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtCOO,
    dtype: op.JustObject | None = ...,
) -> coo_matrix[ScalarT]: ...
@overload  # diagonals: <known>, format: "csc", dtype: None
def diags[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[ScalarT],
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtCSC,
    dtype: op.JustObject | None = ...,
) -> csc_matrix[ScalarT]: ...
@overload  # diagonals: <known>, format: "csr", dtype: None
def diags[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[ScalarT],
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtCSR,
    dtype: op.JustObject | None = ...,
) -> csr_matrix[ScalarT]: ...
@overload  # diagonals: <known>, format: "dok", dtype: None
def diags[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[ScalarT],
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtDOK,
    dtype: op.JustObject | None = ...,
) -> dok_matrix[ScalarT]: ...
@overload  # diagonals: <known>, format: "lil", dtype: None
def diags[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[ScalarT],
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtLIL,
    dtype: op.JustObject | None = ...,
) -> lil_matrix[ScalarT]: ...

#
@overload  # diagonals: <known>, format: "dia" | None, dtype: <known>
def diags[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[_Numeric],
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDIA | None = None,
    *,
    dtype: onp.ToDType[ScalarT],
) -> dia_matrix[ScalarT]: ...
@overload  # diagonals: <known>, format: "bsr", dtype: <known>
def diags[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[_Numeric],
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtBSR,
    dtype: onp.ToDType[ScalarT],
) -> bsr_matrix[ScalarT]: ...
@overload  # diagonals: <known>, format: "coo", dtype: <known>
def diags[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[_Numeric],
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtCOO,
    dtype: onp.ToDType[ScalarT],
) -> coo_matrix[ScalarT]: ...
@overload  # diagonals: <known>, format: "csr", dtype: <known>
def diags[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[_Numeric],
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtCSR,
    dtype: onp.ToDType[ScalarT],
) -> csr_matrix[ScalarT]: ...
@overload  # diagonals: <known>, format: "csc", dtype: <known>
def diags[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[_Numeric],
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtCSC,
    dtype: onp.ToDType[ScalarT],
) -> csc_matrix[ScalarT]: ...
@overload  # diagonals: <known>, format: "dok", dtype: <known>
def diags[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[_Numeric],
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtDOK,
    dtype: onp.ToDType[ScalarT],
) -> dok_matrix[ScalarT]: ...
@overload  # diagonals: <known>, format: "lil", dtype: <known>
def diags[ScalarT: _Numeric](
    diagonals: _ToArray1D2D[_Numeric],
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtLIL,
    dtype: onp.ToDType[ScalarT],
) -> lil_matrix[ScalarT]: ...

#
@overload  # diagonals: <unknown>, format: "dia" | None, dtype: <unknown>
def diags(
    diagonals: _ToComplex1D2D,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    format: _FmtDIA | None = None,
    dtype: _ToDType | op.JustObject | None = ...,
) -> dia_matrix: ...
@overload  # diagonals: <unknown>, format: "bsr", dtype: <unknown>
def diags(
    diagonals: _ToComplex1D2D,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtBSR,
    dtype: _ToDType | op.JustObject | None = ...,
) -> bsr_matrix: ...
@overload  # diagonals: <unknown>, format: "coo", dtype: <unknown>
def diags(
    diagonals: _ToComplex1D2D,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtCOO,
    dtype: _ToDType | op.JustObject | None = ...,
) -> coo_matrix: ...
@overload  # diagonals: <unknown>, format: "csr", dtype: <unknown>
def diags(
    diagonals: _ToComplex1D2D,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtCSR,
    dtype: _ToDType | op.JustObject | None = ...,
) -> csr_matrix: ...
@overload  # diagonals: <unknown>, format: "csc", dtype: <unknown>
def diags(
    diagonals: _ToComplex1D2D,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtCSC,
    dtype: _ToDType | op.JustObject | None = ...,
) -> csc_matrix: ...
@overload  # diagonals: <unknown>, format: "dok", dtype: <unknown>
def diags(
    diagonals: _ToComplex1D2D,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtDOK,
    dtype: _ToDType | op.JustObject | None = ...,
) -> dok_matrix: ...
@overload  # diagonals: <unknown>, format: "lil", dtype: <unknown>
def diags(
    diagonals: _ToComplex1D2D,
    offsets: _Offsets = 0,
    shape: tuple[int, int] | None = None,
    *,
    format: _FmtLIL,
    dtype: _ToDType | op.JustObject | None = ...,
) -> lil_matrix: ...

###
# NOTE: `diags_array` should be prefered over `spdiags`
@overload
def spdiags[ScalarT: _Numeric](
    data: _ToArray1D2D[ScalarT], diags: _Offsets, m: int, n: int, format: _FmtDIA | None = None
) -> dia_matrix[ScalarT]: ...
@overload
def spdiags[ScalarT: _Numeric](
    data: _ToArray1D2D[ScalarT], diags: _Offsets, m: int, n: int, format: _FmtBSR
) -> bsr_matrix[ScalarT]: ...
@overload
def spdiags[ScalarT: _Numeric](
    data: _ToArray1D2D[ScalarT], diags: _Offsets, m: int, n: int, format: _FmtCOO
) -> coo_matrix[ScalarT]: ...
@overload
def spdiags[ScalarT: _Numeric](
    data: _ToArray1D2D[ScalarT], diags: _Offsets, m: int, n: int, format: _FmtCSC
) -> csc_matrix[ScalarT]: ...
@overload
def spdiags[ScalarT: _Numeric](
    data: _ToArray1D2D[ScalarT], diags: _Offsets, m: int, n: int, format: _FmtCSR
) -> csr_matrix[ScalarT]: ...
@overload
def spdiags[ScalarT: _Numeric](
    data: _ToArray1D2D[ScalarT], diags: _Offsets, m: int, n: int, format: _FmtDOK
) -> dok_matrix[ScalarT]: ...
@overload
def spdiags[ScalarT: _Numeric](
    data: _ToArray1D2D[ScalarT], diags: _Offsets, m: int, n: int, format: _FmtLIL
) -> lil_matrix[ScalarT]: ...

#
@overload
def spdiags[ScalarT: _Numeric](
    data: _ToArray1D2D[ScalarT], diags: _Offsets, m: tuple[int, int] | None = None, n: None = None, format: _FmtDIA | None = None
) -> dia_matrix[ScalarT]: ...
@overload
def spdiags[ScalarT: _Numeric](
    data: _ToArray1D2D[ScalarT], diags: _Offsets, m: tuple[int, int] | None = None, n: None = None, *, format: _FmtBSR
) -> bsr_matrix[ScalarT]: ...
@overload
def spdiags[ScalarT: _Numeric](
    data: _ToArray1D2D[ScalarT], diags: _Offsets, m: tuple[int, int] | None = None, n: None = None, *, format: _FmtCOO
) -> coo_matrix[ScalarT]: ...
@overload
def spdiags[ScalarT: _Numeric](
    data: _ToArray1D2D[ScalarT], diags: _Offsets, m: tuple[int, int] | None = None, n: None = None, *, format: _FmtCSC
) -> csc_matrix[ScalarT]: ...
@overload
def spdiags[ScalarT: _Numeric](
    data: _ToArray1D2D[ScalarT], diags: _Offsets, m: tuple[int, int] | None = None, n: None = None, *, format: _FmtCSR
) -> csr_matrix[ScalarT]: ...
@overload
def spdiags[ScalarT: _Numeric](
    data: _ToArray1D2D[ScalarT], diags: _Offsets, m: tuple[int, int] | None = None, n: None = None, *, format: _FmtDOK
) -> dok_matrix[ScalarT]: ...
@overload
def spdiags[ScalarT: _Numeric](
    data: _ToArray1D2D[ScalarT], diags: _Offsets, m: tuple[int, int] | None = None, n: None = None, *, format: _FmtLIL
) -> lil_matrix[ScalarT]: ...

###
@overload  # dtype: float64-like (default), format: "dia" | None
def identity(n: int, dtype: onp.AnyFloat64DType = "d", format: _FmtDIA | None = None) -> dia_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "bsr"
def identity(n: int, dtype: onp.AnyFloat64DType = "d", *, format: _FmtBSR) -> bsr_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "coo"
def identity(n: int, dtype: onp.AnyFloat64DType = "d", *, format: _FmtCOO) -> coo_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "csc"
def identity(n: int, dtype: onp.AnyFloat64DType = "d", *, format: _FmtCSC) -> csc_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "csr"
def identity(n: int, dtype: onp.AnyFloat64DType = "d", *, format: _FmtCSR) -> csr_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "dok"
def identity(n: int, dtype: onp.AnyFloat64DType = "d", *, format: _FmtDOK) -> dok_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "lil"
def identity(n: int, dtype: onp.AnyFloat64DType = "d", *, format: _FmtLIL) -> lil_matrix[np.float64]: ...

#
@overload  # dtype: bool-like, format: "dia" | None
def identity(n: int, dtype: onp.AnyBoolDType, format: _FmtDIA | None = None) -> dia_matrix[np.bool]: ...
@overload  # dtype: bool-like, format: "bsr"
def identity(n: int, dtype: onp.AnyBoolDType, format: _FmtBSR) -> bsr_matrix[np.bool]: ...
@overload  # dtype: bool-like, format: "coo"
def identity(n: int, dtype: onp.AnyBoolDType, format: _FmtCOO) -> coo_matrix[np.bool]: ...
@overload  # dtype: bool-like, format: "csc"
def identity(n: int, dtype: onp.AnyBoolDType, format: _FmtCSC) -> csc_matrix[np.bool]: ...
@overload  # dtype: bool-like, format: "csr"
def identity(n: int, dtype: onp.AnyBoolDType, format: _FmtCSR) -> csr_matrix[np.bool]: ...
@overload  # dtype: bool-like, format: "dok"
def identity(n: int, dtype: onp.AnyBoolDType, format: _FmtDOK) -> dok_matrix[np.bool]: ...
@overload  # dtype: bool-like, format: "lil"
def identity(n: int, dtype: onp.AnyBoolDType, format: _FmtLIL) -> lil_matrix[np.bool]: ...

#
@overload  # dtype: int-like, format: "dia" | None
def identity(n: int, dtype: onp.AnyIntDType, format: _FmtDIA | None = None) -> dia_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "bsr"
def identity(n: int, dtype: onp.AnyIntDType, format: _FmtBSR) -> bsr_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "coo"
def identity(n: int, dtype: onp.AnyIntDType, format: _FmtCOO) -> coo_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "csc"
def identity(n: int, dtype: onp.AnyIntDType, format: _FmtCSC) -> csc_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "csr"
def identity(n: int, dtype: onp.AnyIntDType, format: _FmtCSR) -> csr_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "dok"
def identity(n: int, dtype: onp.AnyIntDType, format: _FmtDOK) -> dok_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "lil"
def identity(n: int, dtype: onp.AnyIntDType, format: _FmtLIL) -> lil_matrix[np.int_]: ...

#
@overload  # dtype: complex128-like, format: "dia" | None
def identity(n: int, dtype: onp.AnyComplex128DType, format: _FmtDIA | None = None) -> dia_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "bsr"
def identity(n: int, dtype: onp.AnyComplex128DType, format: _FmtBSR) -> bsr_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "coo"
def identity(n: int, dtype: onp.AnyComplex128DType, format: _FmtCOO) -> coo_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "csc"
def identity(n: int, dtype: onp.AnyComplex128DType, format: _FmtCSC) -> csc_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "csr"
def identity(n: int, dtype: onp.AnyComplex128DType, format: _FmtCSR) -> csr_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "dok"
def identity(n: int, dtype: onp.AnyComplex128DType, format: _FmtDOK) -> dok_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "lil"
def identity(n: int, dtype: onp.AnyComplex128DType, format: _FmtLIL) -> lil_matrix[np.complex128]: ...

#
@overload  # dtype: <known>, format: "dia" | None
def identity[ScalarT: _Numeric](n: int, dtype: onp.ToDType[ScalarT], format: _FmtDIA | None = None) -> dia_matrix[ScalarT]: ...
@overload  # dtype: <known>, format: "bsr"
def identity[ScalarT: _Numeric](n: int, dtype: onp.ToDType[ScalarT], format: _FmtBSR) -> bsr_matrix[ScalarT]: ...
@overload  # dtype: <known>, format: "coo"
def identity[ScalarT: _Numeric](n: int, dtype: onp.ToDType[ScalarT], format: _FmtCOO) -> coo_matrix[ScalarT]: ...
@overload  # dtype: <known>, format: "csc"
def identity[ScalarT: _Numeric](n: int, dtype: onp.ToDType[ScalarT], format: _FmtCSC) -> csc_matrix[ScalarT]: ...
@overload  # dtype: <known>, format: "csr"
def identity[ScalarT: _Numeric](n: int, dtype: onp.ToDType[ScalarT], format: _FmtCSR) -> csr_matrix[ScalarT]: ...
@overload  # dtype: <known>, format: "dok"
def identity[ScalarT: _Numeric](n: int, dtype: onp.ToDType[ScalarT], format: _FmtDOK) -> dok_matrix[ScalarT]: ...
@overload  # dtype: <known>, format: "lil"
def identity[ScalarT: _Numeric](n: int, dtype: onp.ToDType[ScalarT], format: _FmtLIL) -> lil_matrix[ScalarT]: ...

#
@overload  # dtype: <unknown>, format: "dia" | None
def identity(n: int, dtype: _ToDType, format: _FmtDIA | None = None) -> dia_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "bsr"
def identity(n: int, dtype: _ToDType, format: _FmtBSR) -> bsr_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "coo"
def identity(n: int, dtype: _ToDType, format: _FmtCOO) -> coo_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "csc"
def identity(n: int, dtype: _ToDType, format: _FmtCSC) -> csc_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "csr"
def identity(n: int, dtype: _ToDType, format: _FmtCSR) -> csr_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "dok"
def identity(n: int, dtype: _ToDType, format: _FmtDOK) -> dok_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "lil"
def identity(n: int, dtype: _ToDType, format: _FmtLIL) -> lil_matrix[Incomplete]: ...

###
@overload  # dtype: float64-like (default), format: "dia" | None
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtDIA | None = None
) -> dia_array[np.float64]: ...
@overload  # dtype: float64-like (default), format: "bsr"
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtBSR
) -> bsr_array[np.float64]: ...
@overload  # dtype: float64-like (default), format: "coo"
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtCOO
) -> _COOArray2D[np.float64]: ...
@overload  # dtype: float64-like (default), format: "csc"
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtCSC
) -> csc_array[np.float64]: ...
@overload  # dtype: float64-like (default), format: "csr"
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtCSR
) -> _CSRArray2D[np.float64]: ...
@overload  # dtype: float64-like (default), format: "dok"
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtDOK
) -> _DOKArray2D[np.float64]: ...
@overload  # dtype: float64-like (default), format: "lil"
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtLIL
) -> lil_array[np.float64]: ...

#
@overload  # dtype: bool-like, format: "dia" | None
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyBoolDType, format: _FmtDIA | None = None
) -> dia_array[np.bool]: ...
@overload  # dtype: bool-like, format: "bsr"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyBoolDType, format: _FmtBSR) -> bsr_array[np.bool]: ...
@overload  # dtype: bool-like, format: "coo"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyBoolDType, format: _FmtCOO) -> _COOArray2D[np.bool]: ...
@overload  # dtype: bool-like, format: "csc"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyBoolDType, format: _FmtCSC) -> csc_array[np.bool]: ...
@overload  # dtype: bool-like, format: "csr"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyBoolDType, format: _FmtCSR) -> _CSRArray2D[np.bool]: ...
@overload  # dtype: bool-like, format: "dok"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyBoolDType, format: _FmtDOK) -> _DOKArray2D[np.bool]: ...
@overload  # dtype: bool-like, format: "lil"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyBoolDType, format: _FmtLIL) -> lil_array[np.bool]: ...

#
@overload  # dtype: int-like, format: "dia" | None
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyIntDType, format: _FmtDIA | None = None
) -> dia_array[np.int_]: ...
@overload  # dtype: int-like, format: "bsr"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyIntDType, format: _FmtBSR) -> bsr_array[np.int_]: ...
@overload  # dtype: int-like, format: "coo"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyIntDType, format: _FmtCOO) -> _COOArray2D[np.int_]: ...
@overload  # dtype: int-like, format: "csc"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyIntDType, format: _FmtCSC) -> csc_array[np.int_]: ...
@overload  # dtype: int-like, format: "csr"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyIntDType, format: _FmtCSR) -> _CSRArray2D[np.int_]: ...
@overload  # dtype: int-like, format: "dok"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyIntDType, format: _FmtDOK) -> _DOKArray2D[np.int_]: ...
@overload  # dtype: int-like, format: "lil"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyIntDType, format: _FmtLIL) -> lil_array[np.int_]: ...

#
@overload  # dtype: complex128-like, format: "dia" | None
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyComplex128DType, format: _FmtDIA | None = None
) -> dia_array[np.complex128]: ...
@overload  # dtype: complex128-like, format: "bsr"
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyComplex128DType, format: _FmtBSR
) -> bsr_array[np.complex128]: ...
@overload  # dtype: complex128-like, format: "coo"
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyComplex128DType, format: _FmtCOO
) -> _COOArray2D[np.complex128]: ...
@overload  # dtype: complex128-like, format: "csc"
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyComplex128DType, format: _FmtCSC
) -> csc_array[np.complex128]: ...
@overload  # dtype: complex128-like, format: "csr"
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyComplex128DType, format: _FmtCSR
) -> _CSRArray2D[np.complex128]: ...
@overload  # dtype: complex128-like, format: "dok"
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyComplex128DType, format: _FmtDOK
) -> _DOKArray2D[np.complex128]: ...
@overload  # dtype: complex128-like, format: "lil"
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.AnyComplex128DType, format: _FmtLIL
) -> lil_array[np.complex128]: ...

#
@overload  # dtype: <known>, format: "dia" | None
def eye_array[ScalarT: _Numeric](
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.ToDType[ScalarT], format: _FmtDIA | None = None
) -> dia_array[ScalarT]: ...
@overload  # dtype: <known>, format: "bsr"
def eye_array[ScalarT: _Numeric](
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.ToDType[ScalarT], format: _FmtBSR
) -> bsr_array[ScalarT]: ...
@overload  # dtype: <known>, format: "coo"
def eye_array[ScalarT: _Numeric](
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.ToDType[ScalarT], format: _FmtCOO
) -> _COOArray2D[ScalarT]: ...
@overload  # dtype: <known>, format: "csc"
def eye_array[ScalarT: _Numeric](
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.ToDType[ScalarT], format: _FmtCSC
) -> csc_array[ScalarT]: ...
@overload  # dtype: <known>, format: "csr"
def eye_array[ScalarT: _Numeric](
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.ToDType[ScalarT], format: _FmtCSR
) -> _CSRArray2D[ScalarT]: ...
@overload  # dtype: <known>, format: "dok"
def eye_array[ScalarT: _Numeric](
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.ToDType[ScalarT], format: _FmtDOK
) -> _DOKArray2D[ScalarT]: ...
@overload  # dtype: <known>, format: "lil"
def eye_array[ScalarT: _Numeric](
    m: int, n: int | None = None, *, k: int = 0, dtype: onp.ToDType[ScalarT], format: _FmtLIL
) -> lil_array[ScalarT]: ...

#
@overload  # dtype: <unknown>, format: "dia" | None
def eye_array(
    m: int, n: int | None = None, *, k: int = 0, dtype: _ToDType, format: _FmtDIA | None = None
) -> dia_array[Incomplete]: ...
@overload  # dtype: <unknown>, format: "bsr"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: _ToDType, format: _FmtBSR) -> bsr_array[Incomplete]: ...
@overload  # dtype: <unknown>, format: "coo"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: _ToDType, format: _FmtCOO) -> _COOArray2D[Incomplete]: ...
@overload  # dtype: <unknown>, format: "csc"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: _ToDType, format: _FmtCSC) -> csc_array[Incomplete]: ...
@overload  # dtype: <unknown>, format: "csr"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: _ToDType, format: _FmtCSR) -> _CSRArray2D[Incomplete]: ...
@overload  # dtype: <unknown>, format: "dok"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: _ToDType, format: _FmtDOK) -> _DOKArray2D[Incomplete]: ...
@overload  # dtype: <unknown>, format: "lil"
def eye_array(m: int, n: int | None = None, *, k: int = 0, dtype: _ToDType, format: _FmtLIL) -> lil_array[Incomplete]: ...

###
# NOTE: `eye_array` should be prefered over `eye`
@overload  # dtype: float64-like (default), format: "dia" | None
def eye(
    m: int, n: int | None = None, k: int = 0, dtype: onp.AnyFloat64DType = ..., format: _FmtDIA | None = None
) -> dia_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "bsr"
def eye(
    m: int, n: int | None = None, k: int = 0, dtype: onp.AnyFloat64DType = ..., *, format: _FmtBSR
) -> bsr_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "coo"
def eye(
    m: int, n: int | None = None, k: int = 0, dtype: onp.AnyFloat64DType = ..., *, format: _FmtCOO
) -> coo_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "csc"
def eye(
    m: int, n: int | None = None, k: int = 0, dtype: onp.AnyFloat64DType = ..., *, format: _FmtCSC
) -> csc_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "csr"
def eye(
    m: int, n: int | None = None, k: int = 0, dtype: onp.AnyFloat64DType = ..., *, format: _FmtCSR
) -> csr_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "dok"
def eye(
    m: int, n: int | None = None, k: int = 0, dtype: onp.AnyFloat64DType = ..., *, format: _FmtDOK
) -> dok_matrix[np.float64]: ...
@overload  # dtype: float64-like (default), format: "lil"
def eye(
    m: int, n: int | None = None, k: int = 0, dtype: onp.AnyFloat64DType = ..., *, format: _FmtLIL
) -> lil_matrix[np.float64]: ...

#
@overload  # dtype: bool-like, format: "dia" | None
def eye(
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyBoolDType, format: _FmtDIA | None = None
) -> dia_matrix[np.bool]: ...
@overload  # dtype: bool-like, format: "bsr"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyBoolDType, format: _FmtBSR) -> bsr_matrix[np.bool]: ...
@overload  # dtype: bool-like, format: "coo"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyBoolDType, format: _FmtCOO) -> coo_matrix[np.bool]: ...
@overload  # dtype: bool-like, format: "csc"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyBoolDType, format: _FmtCSC) -> csc_matrix[np.bool]: ...
@overload  # dtype: bool-like, format: "csr"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyBoolDType, format: _FmtCSR) -> csr_matrix[np.bool]: ...
@overload  # dtype: bool-like, format: "dok"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyBoolDType, format: _FmtDOK) -> dok_matrix[np.bool]: ...
@overload  # dtype: bool-like, format: "lil"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyBoolDType, format: _FmtLIL) -> lil_matrix[np.bool]: ...

#
@overload  # dtype: int-like, format: "dia" | None
def eye(
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyIntDType, format: _FmtDIA | None = None
) -> dia_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "bsr"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyIntDType, format: _FmtBSR) -> bsr_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "coo"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyIntDType, format: _FmtCOO) -> coo_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "csc"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyIntDType, format: _FmtCSC) -> csc_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "csr"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyIntDType, format: _FmtCSR) -> csr_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "dok"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyIntDType, format: _FmtDOK) -> dok_matrix[np.int_]: ...
@overload  # dtype: int-like, format: "lil"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyIntDType, format: _FmtLIL) -> lil_matrix[np.int_]: ...

#
@overload  # dtype: complex128-like, format: "dia" | None
def eye(
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyComplex128DType, format: _FmtDIA | None = None
) -> dia_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "bsr"
def eye(
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyComplex128DType, format: _FmtBSR
) -> bsr_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "coo"
def eye(
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyComplex128DType, format: _FmtCOO
) -> coo_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "csc"
def eye(
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyComplex128DType, format: _FmtCSC
) -> csc_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "csr"
def eye(
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyComplex128DType, format: _FmtCSR
) -> csr_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "dok"
def eye(
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyComplex128DType, format: _FmtDOK
) -> dok_matrix[np.complex128]: ...
@overload  # dtype: complex128-like, format: "lil"
def eye(
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.AnyComplex128DType, format: _FmtLIL
) -> lil_matrix[np.complex128]: ...

#
@overload  # dtype: <known>, format: "dia" | None
def eye[ScalarT: _Numeric](
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.ToDType[ScalarT], format: _FmtDIA | None = None
) -> dia_matrix[ScalarT]: ...
@overload  # dtype: <known>, format: "bsr"
def eye[ScalarT: _Numeric](
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.ToDType[ScalarT], format: _FmtBSR
) -> bsr_matrix[ScalarT]: ...
@overload  # dtype: <known>, format: "coo"
def eye[ScalarT: _Numeric](
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.ToDType[ScalarT], format: _FmtCOO
) -> coo_matrix[ScalarT]: ...
@overload  # dtype: <known>, format: "csc"
def eye[ScalarT: _Numeric](
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.ToDType[ScalarT], format: _FmtCSC
) -> csc_matrix[ScalarT]: ...
@overload  # dtype: <known>, format: "csr"
def eye[ScalarT: _Numeric](
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.ToDType[ScalarT], format: _FmtCSR
) -> csr_matrix[ScalarT]: ...
@overload  # dtype: <known>, format: "dok"
def eye[ScalarT: _Numeric](
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.ToDType[ScalarT], format: _FmtDOK
) -> dok_matrix[ScalarT]: ...
@overload  # dtype: <known>, format: "lil"
def eye[ScalarT: _Numeric](
    m: int, n: int | None = None, k: int = 0, *, dtype: onp.ToDType[ScalarT], format: _FmtLIL
) -> lil_matrix[ScalarT]: ...

#
@overload  # dtype: <unknown>, format: "dia" | None
def eye(
    m: int, n: int | None = None, k: int = 0, *, dtype: _ToDType, format: _FmtDIA | None = None
) -> dia_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "bsr"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: _ToDType, format: _FmtBSR) -> bsr_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "coo"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: _ToDType, format: _FmtCOO) -> coo_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "csc"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: _ToDType, format: _FmtCSC) -> csc_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "csr"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: _ToDType, format: _FmtCSR) -> csr_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "dok"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: _ToDType, format: _FmtDOK) -> dok_matrix[Incomplete]: ...
@overload  # dtype: <unknown>, format: "lil"
def eye(m: int, n: int | None = None, k: int = 0, *, dtype: _ToDType, format: _FmtLIL) -> lil_matrix[Incomplete]: ...

###
@overload  # A: spmatrix or 2d array-like, B: spmatrix, format: None = None
def kron[ScalarT: _Numeric](
    A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: None = None
) -> bsr_matrix[ScalarT] | coo_matrix[ScalarT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "bsr"
def kron[ScalarT: _Numeric](A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: _FmtBSR) -> bsr_matrix[ScalarT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "coo"
def kron[ScalarT: _Numeric](A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: _FmtCOO) -> coo_matrix[ScalarT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "csc"
def kron[ScalarT: _Numeric](A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: _FmtCSC) -> csc_matrix[ScalarT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "csr"
def kron[ScalarT: _Numeric](A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: _FmtCSR) -> csr_matrix[ScalarT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "dia"
def kron[ScalarT: _Numeric](A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: _FmtDIA) -> dia_matrix[ScalarT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "dok"
def kron[ScalarT: _Numeric](A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: _FmtDOK) -> dok_matrix[ScalarT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "lil"
def kron[ScalarT: _Numeric](A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: _FmtLIL) -> lil_matrix[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: None = ...
def kron[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: None = None
) -> bsr_array[ScalarT] | _COOArray2D[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "bsr"
def kron[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: _FmtBSR
) -> bsr_array[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "coo"
def kron[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: _FmtCOO
) -> _COOArray2D[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "csc"
def kron[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: _FmtCSC
) -> csc_array[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "csr"
def kron[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: _FmtCSR
) -> _CSRArray2D[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "dia"
def kron[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: _FmtDIA
) -> dia_array[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "dok"
def kron[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: _FmtDOK
) -> _DOKArray2D[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "lil"
def kron[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: _FmtLIL
) -> lil_array[ScalarT]: ...
@overload  # A: sparse, B: sparray, format: None = ...
def kron[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: None = None
) -> bsr_array[ScalarT] | _COOArray2D[ScalarT]: ...
@overload  # A: sparse, B: sparray, format: "bsr"
def kron[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: _FmtBSR
) -> bsr_array[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "coo"
def kron[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: _FmtCOO
) -> _COOArray2D[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "csc"
def kron[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: _FmtCSC
) -> csc_array[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "csr"
def kron[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: _FmtCSR
) -> _CSRArray2D[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "dia"
def kron[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: _FmtDIA
) -> dia_array[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "dok"
def kron[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: _FmtDOK
) -> _DOKArray2D[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "lil"
def kron[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: _FmtLIL
) -> lil_array[ScalarT]: ...
@overload  # A: unknown array-like, B: unknown array-like  (catch-all)
def kron(A: onp.ToComplex2D, B: onp.ToComplex2D, format: _Format | None = None) -> Incomplete: ...

###
# NOTE: The `overload-overlap` mypy errors are false positives.
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "csr" | None
def kronsum[ScalarT: _Numeric](
    A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: _FmtCSR | None = None
) -> csr_matrix[ScalarT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "bsr"
def kronsum[ScalarT: _Numeric](A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: _FmtBSR) -> bsr_matrix[ScalarT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "coo"
def kronsum[ScalarT: _Numeric](A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: _FmtCOO) -> coo_matrix[ScalarT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "csc"
def kronsum[ScalarT: _Numeric](A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: _FmtCSC) -> csc_matrix[ScalarT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "dia"
def kronsum[ScalarT: _Numeric](A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: _FmtDIA) -> dia_matrix[ScalarT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "dok"
def kronsum[ScalarT: _Numeric](A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: _FmtDOK) -> dok_matrix[ScalarT]: ...
@overload  # A: spmatrix or 2d array-like, B: spmatrix or 2d array-like, format: "lil"
def kronsum[ScalarT: _Numeric](A: _ToSpMatrix[ScalarT], B: _ToSpMatrix[ScalarT], format: _FmtLIL) -> lil_matrix[ScalarT]: ...

#
@overload  # A: sparray, B: sparse, format: "csr" | None
def kronsum[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: _FmtCSR | None = None
) -> _CSRArray2D[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "bsr"
def kronsum[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: _FmtBSR
) -> bsr_array[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "coo"
def kronsum[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: _FmtCOO
) -> _COOArray2D[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "csc"
def kronsum[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: _FmtCSC
) -> csc_array[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "dia"
def kronsum[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: _FmtDIA
) -> dia_array[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "dok"
def kronsum[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: _FmtDOK
) -> _DOKArray2D[ScalarT]: ...
@overload  # A: sparray, B: sparse, format: "lil"
def kronsum[ScalarT: _Numeric](
    A: sparray[ScalarT, tuple[int, int]], B: _ToSparse2D[ScalarT], format: _FmtLIL
) -> lil_array[ScalarT]: ...

#
@overload  # A: sparse, B: sparray, format: "csr" | None
def kronsum[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: _FmtCSR | None = None
) -> _CSRArray2D[ScalarT]: ...
@overload  # A: sparse, B: sparray, format: "bsr"
def kronsum[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: _FmtBSR
) -> bsr_array[ScalarT]: ...
@overload  # A: sparse, B: sparray, format: "coo"
def kronsum[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: _FmtCOO
) -> _COOArray2D[ScalarT]: ...
@overload  # A: sparse, B: sparray, format: "csc"
def kronsum[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: _FmtCSC
) -> csc_array[ScalarT]: ...
@overload  # A: sparse, B: sparray, format: "dia"
def kronsum[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: _FmtDIA
) -> dia_array[ScalarT]: ...
@overload  # A: sparse, B: sparray, format: "dok"
def kronsum[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: _FmtDOK
) -> _DOKArray2D[ScalarT]: ...
@overload  # A: sparse, B: sparray, format: "lil"
def kronsum[ScalarT: _Numeric](
    A: _ToSparse2D[ScalarT], B: sparray[ScalarT, tuple[int, int]], format: _FmtLIL
) -> lil_array[ScalarT]: ...
@overload  # A: unknown array-like, B: unknown array-like  (catch-all)
def kronsum(A: onp.ToComplex2D, B: onp.ToComplex2D, format: _Format | None = None) -> Incomplete: ...

###
# NOTE: keep in sync with `vstack`
@overload  # sparray, format: <default>, dtype: <default>
def hstack[T](blocks: Seq[_CanStack[T]], format: None = None, dtype: None = None) -> T: ...
@overload  # sparray, format: "bsr", dtype: <default>
def hstack[ScalarT: _Numeric](blocks: Seq[_SpArray2D[ScalarT]], format: _FmtBSR, dtype: None = None) -> bsr_array[ScalarT]: ...
@overload  # sparray, format: "coo", dtype: <default>
def hstack[ScalarT: _Numeric](blocks: Seq[_SpArray2D[ScalarT]], format: _FmtCOO, dtype: None = None) -> _COOArray2D[ScalarT]: ...
@overload  # sparray, format: "csc", dtype: <default>
def hstack[ScalarT: _Numeric](blocks: Seq[_SpArray2D[ScalarT]], format: _FmtCSC, dtype: None = None) -> csc_array[ScalarT]: ...
@overload  # sparray, format: "csr", dtype: <default>
def hstack[ScalarT: _Numeric](blocks: Seq[_SpArray2D[ScalarT]], format: _FmtCSR, dtype: None = None) -> _CSRArray2D[ScalarT]: ...
@overload  # sparray, format: "dia", dtype: <default>
def hstack[ScalarT: _Numeric](blocks: Seq[_SpArray2D[ScalarT]], format: _FmtDIA, dtype: None = None) -> dia_array[ScalarT]: ...
@overload  # sparray, format: "dok", dtype: <default>
def hstack[ScalarT: _Numeric](blocks: Seq[_SpArray2D[ScalarT]], format: _FmtDOK, dtype: None = None) -> _DOKArray2D[ScalarT]: ...
@overload  # sparray, format: "lil", dtype: <default>
def hstack[ScalarT: _Numeric](blocks: Seq[_SpArray2D[ScalarT]], format: _FmtLIL, dtype: None = None) -> lil_array[ScalarT]: ...

#
@overload  # sparray, format: <default>, dtype: bool-like
def hstack[T](blocks: Seq[_CanStackAs[np.bool, T]], format: None = None, *, dtype: onp.AnyBoolDType) -> T: ...
@overload  # sparray, format: "bsr", dtype: bool-like
def hstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyBoolDType) -> bsr_array[np.bool]: ...
@overload  # sparray, format: "coo", dtype: bool-like
def hstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyBoolDType) -> _COOArray2D[np.bool]: ...
@overload  # sparray, format: "csc", dtype: bool-like
def hstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyBoolDType) -> csc_array[np.bool]: ...
@overload  # sparray, format: "csr", dtype: bool-like
def hstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyBoolDType) -> _CSRArray2D[np.bool]: ...
@overload  # sparray, format: "dia", dtype: bool-like
def hstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyBoolDType) -> dia_array[np.bool]: ...
@overload  # sparray, format: "dok", dtype: bool-like
def hstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyBoolDType) -> _DOKArray2D[np.bool]: ...
@overload  # sparray, format: "lil", dtype: bool-like
def hstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyBoolDType) -> lil_array[np.bool]: ...

#
@overload  # sparray, format: <default>, dtype: int-like
def hstack[T](blocks: Seq[_CanStackAs[np.int_, T]], format: None = None, *, dtype: onp.AnyIntDType) -> T: ...
@overload  # sparray, format: "bsr", dtype: int-like
def hstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyIntDType) -> bsr_array[np.int_]: ...
@overload  # sparray, format: "coo", dtype: int-like
def hstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyIntDType) -> _COOArray2D[np.int_]: ...
@overload  # sparray, format: "csc", dtype: int-like
def hstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyIntDType) -> csc_array[np.int_]: ...
@overload  # sparray, format: "csr", dtype: int-like
def hstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyIntDType) -> _CSRArray2D[np.int_]: ...
@overload  # sparray, format: "dia", dtype: int-like
def hstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyIntDType) -> dia_array[np.int_]: ...
@overload  # sparray, format: "dok", dtype: int-like
def hstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyIntDType) -> _DOKArray2D[np.int_]: ...
@overload  # sparray, format: "lil", dtype: int-like
def hstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyIntDType) -> lil_array[np.int_]: ...

#
@overload  # sparray, format: <default>, dtype: float64-like
def hstack[T](blocks: Seq[_CanStackAs[np.float64, T]], format: None = None, *, dtype: onp.AnyFloat64DType) -> T: ...
@overload  # sparray, format: "bsr", dtype: float64-like
def hstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyFloat64DType) -> bsr_array[np.float64]: ...
@overload  # sparray, format: "coo", dtype: float64-like
def hstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyFloat64DType) -> _COOArray2D[np.float64]: ...
@overload  # sparray, format: "csc", dtype: float64-like
def hstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyFloat64DType) -> csc_array[np.float64]: ...
@overload  # sparray, format: "csr", dtype: float64-like
def hstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyFloat64DType) -> _CSRArray2D[np.float64]: ...
@overload  # sparray, format: "dia", dtype: float64-like
def hstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyFloat64DType) -> dia_array[np.float64]: ...
@overload  # sparray, format: "dok", dtype: float64-like
def hstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyFloat64DType) -> _DOKArray2D[np.float64]: ...
@overload  # sparray, format: "lil", dtype: float64-like
def hstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyFloat64DType) -> lil_array[np.float64]: ...

#
@overload  # sparray, format: <default>, dtype: complex128-like
def hstack[T](blocks: Seq[_CanStackAs[np.complex128, T]], format: None = None, *, dtype: onp.AnyComplex128DType) -> T: ...
@overload  # sparray, format: "bsr", dtype: complex128-like
def hstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyComplex128DType) -> bsr_array[np.complex128]: ...
@overload  # sparray, format: "coo", dtype: complex128-like
def hstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyComplex128DType) -> _COOArray2D[np.complex128]: ...
@overload  # sparray, format: "csc", dtype: complex128-like
def hstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyComplex128DType) -> csc_array[np.complex128]: ...
@overload  # sparray, format: "csr", dtype: complex128-like
def hstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyComplex128DType) -> _CSRArray2D[np.complex128]: ...
@overload  # sparray, format: "dia", dtype: complex128-like
def hstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyComplex128DType) -> dia_array[np.complex128]: ...
@overload  # sparray, format: "dok", dtype: complex128-like
def hstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyComplex128DType) -> _DOKArray2D[np.complex128]: ...
@overload  # sparray, format: "lil", dtype: complex128-like
def hstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyComplex128DType) -> lil_array[np.complex128]: ...

#
@overload  # sparray, format: <default>, dtype: <known>
def hstack[ScalarT: _Numeric, T](
    blocks: Seq[_CanStackAs[ScalarT, T]], format: None = None, *, dtype: onp.ToDType[ScalarT]
) -> T: ...
@overload  # sparray, format: "bsr", dtype: <known>
def hstack[ScalarT: _Numeric](blocks: Seq[sparray], format: _FmtBSR, dtype: onp.ToDType[ScalarT]) -> bsr_array[ScalarT]: ...
@overload  # sparray, format: "coo", dtype: <known>
def hstack[ScalarT: _Numeric](blocks: Seq[sparray], format: _FmtCOO, dtype: onp.ToDType[ScalarT]) -> _COOArray2D[ScalarT]: ...
@overload  # sparray, format: "csc", dtype: <known>
def hstack[ScalarT: _Numeric](blocks: Seq[sparray], format: _FmtCSC, dtype: onp.ToDType[ScalarT]) -> csc_array[ScalarT]: ...
@overload  # sparray, format: "csr", dtype: <known>
def hstack[ScalarT: _Numeric](blocks: Seq[sparray], format: _FmtCSR, dtype: onp.ToDType[ScalarT]) -> _CSRArray2D[ScalarT]: ...
@overload  # sparray, format: "dia", dtype: <known>
def hstack[ScalarT: _Numeric](blocks: Seq[sparray], format: _FmtDIA, dtype: onp.ToDType[ScalarT]) -> dia_array[ScalarT]: ...
@overload  # sparray, format: "dok", dtype: <known>
def hstack[ScalarT: _Numeric](blocks: Seq[sparray], format: _FmtDOK, dtype: onp.ToDType[ScalarT]) -> _DOKArray2D[ScalarT]: ...
@overload  # sparray, format: "lil", dtype: <known>
def hstack[ScalarT: _Numeric](blocks: Seq[sparray], format: _FmtLIL, dtype: onp.ToDType[ScalarT]) -> lil_array[ScalarT]: ...

#
@overload  # sparray, format: <default>, dtype: <unknown>
def hstack[T](blocks: Seq[_CanStackAs[Any, T]], format: None = None, *, dtype: _ToDType) -> T: ...
@overload  # sparray, format: "bsr", dtype: <unknown>
def hstack(blocks: Seq[sparray], format: _FmtBSR, dtype: _ToDType) -> bsr_array: ...
@overload  # sparray, format: "coo", dtype: <unknown>
def hstack(blocks: Seq[sparray], format: _FmtCOO, dtype: _ToDType) -> _COOArray2D[Any]: ...
@overload  # sparray, format: "csc", dtype: <unknown>
def hstack(blocks: Seq[sparray], format: _FmtCSC, dtype: _ToDType) -> csc_array: ...
@overload  # sparray, format: "csr", dtype: <unknown>
def hstack(blocks: Seq[sparray], format: _FmtCSR, dtype: _ToDType) -> _CSRArray2D[Any]: ...
@overload  # sparray, format: "dia", dtype: <unknown>
def hstack(blocks: Seq[sparray], format: _FmtDIA, dtype: _ToDType) -> dia_array: ...
@overload  # sparray, format: "dok", dtype: <unknown>
def hstack(blocks: Seq[sparray], format: _FmtDOK, dtype: _ToDType) -> _DOKArray2D[Any]: ...
@overload  # sparray, format: "lil", dtype: <unknown>
def hstack(blocks: Seq[sparray], format: _FmtLIL, dtype: _ToDType) -> lil_array: ...

#
@overload
def hstack(blocks: Seq[_spbase], format: _Format, dtype: _ToDType | None = None) -> Incomplete: ...

###
# NOTE: keep in sync with `hstack`
@overload  # sparray, format: <default>, dtype: <default>
def vstack[T](blocks: Seq[_CanStack[T]], format: None = None, dtype: None = None) -> T: ...
@overload  # sparray, format: "bsr", dtype: <default>
def vstack[ScalarT: _Numeric](blocks: Seq[_SpArray2D[ScalarT]], format: _FmtBSR, dtype: None = None) -> bsr_array[ScalarT]: ...
@overload  # sparray, format: "coo", dtype: <default>
def vstack[ScalarT: _Numeric](blocks: Seq[_SpArray2D[ScalarT]], format: _FmtCOO, dtype: None = None) -> _COOArray2D[ScalarT]: ...
@overload  # sparray, format: "csc", dtype: <default>
def vstack[ScalarT: _Numeric](blocks: Seq[_SpArray2D[ScalarT]], format: _FmtCSC, dtype: None = None) -> csc_array[ScalarT]: ...
@overload  # sparray, format: "csr", dtype: <default>
def vstack[ScalarT: _Numeric](blocks: Seq[_SpArray2D[ScalarT]], format: _FmtCSR, dtype: None = None) -> _CSRArray2D[ScalarT]: ...
@overload  # sparray, format: "dia", dtype: <default>
def vstack[ScalarT: _Numeric](blocks: Seq[_SpArray2D[ScalarT]], format: _FmtDIA, dtype: None = None) -> dia_array[ScalarT]: ...
@overload  # sparray, format: "dok", dtype: <default>
def vstack[ScalarT: _Numeric](blocks: Seq[_SpArray2D[ScalarT]], format: _FmtDOK, dtype: None = None) -> _DOKArray2D[ScalarT]: ...
@overload  # sparray, format: "lil", dtype: <default>
def vstack[ScalarT: _Numeric](blocks: Seq[_SpArray2D[ScalarT]], format: _FmtLIL, dtype: None = None) -> lil_array[ScalarT]: ...

#
@overload  # sparray, format: <default>, dtype: bool-like
def vstack[T](blocks: Seq[_CanStackAs[np.bool, T]], format: None = None, *, dtype: onp.AnyBoolDType) -> T: ...
@overload  # sparray, format: "bsr", dtype: bool-like
def vstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyBoolDType) -> bsr_array[np.bool]: ...
@overload  # sparray, format: "coo", dtype: bool-like
def vstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyBoolDType) -> _COOArray2D[np.bool]: ...
@overload  # sparray, format: "csc", dtype: bool-like
def vstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyBoolDType) -> csc_array[np.bool]: ...
@overload  # sparray, format: "csr", dtype: bool-like
def vstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyBoolDType) -> _CSRArray2D[np.bool]: ...
@overload  # sparray, format: "dia", dtype: bool-like
def vstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyBoolDType) -> dia_array[np.bool]: ...
@overload  # sparray, format: "dok", dtype: bool-like
def vstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyBoolDType) -> _DOKArray2D[np.bool]: ...
@overload  # sparray, format: "lil", dtype: bool-like
def vstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyBoolDType) -> lil_array[np.bool]: ...

#
@overload  # sparray, format: <default>, dtype: int-like
def vstack[T](blocks: Seq[_CanStackAs[np.int_, T]], format: None = None, *, dtype: onp.AnyIntDType) -> T: ...
@overload  # sparray, format: "bsr", dtype: int-like
def vstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyIntDType) -> bsr_array[np.int_]: ...
@overload  # sparray, format: "coo", dtype: int-like
def vstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyIntDType) -> _COOArray2D[np.int_]: ...
@overload  # sparray, format: "csc", dtype: int-like
def vstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyIntDType) -> csc_array[np.int_]: ...
@overload  # sparray, format: "csr", dtype: int-like
def vstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyIntDType) -> _CSRArray2D[np.int_]: ...
@overload  # sparray, format: "dia", dtype: int-like
def vstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyIntDType) -> dia_array[np.int_]: ...
@overload  # sparray, format: "dok", dtype: int-like
def vstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyIntDType) -> _DOKArray2D[np.int_]: ...
@overload  # sparray, format: "lil", dtype: int-like
def vstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyIntDType) -> lil_array[np.int_]: ...

#
@overload  # sparray, format: <default>, dtype: float64-like
def vstack[T](blocks: Seq[_CanStackAs[np.float64, T]], format: None = None, *, dtype: onp.AnyFloat64DType) -> T: ...
@overload  # sparray, format: "bsr", dtype: float64-like
def vstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyFloat64DType) -> bsr_array[np.float64]: ...
@overload  # sparray, format: "coo", dtype: float64-like
def vstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyFloat64DType) -> _COOArray2D[np.float64]: ...
@overload  # sparray, format: "csc", dtype: float64-like
def vstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyFloat64DType) -> csc_array[np.float64]: ...
@overload  # sparray, format: "csr", dtype: float64-like
def vstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyFloat64DType) -> _CSRArray2D[np.float64]: ...
@overload  # sparray, format: "dia", dtype: float64-like
def vstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyFloat64DType) -> dia_array[np.float64]: ...
@overload  # sparray, format: "dok", dtype: float64-like
def vstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyFloat64DType) -> _DOKArray2D[np.float64]: ...
@overload  # sparray, format: "lil", dtype: float64-like
def vstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyFloat64DType) -> lil_array[np.float64]: ...

#
@overload  # sparray, format: <default>, dtype: complex128-like
def vstack[T](blocks: Seq[_CanStackAs[np.complex128, T]], format: None = None, *, dtype: onp.AnyComplex128DType) -> T: ...
@overload  # sparray, format: "bsr", dtype: complex128-like
def vstack(blocks: Seq[sparray], format: _FmtBSR, dtype: onp.AnyComplex128DType) -> bsr_array[np.complex128]: ...
@overload  # sparray, format: "coo", dtype: complex128-like
def vstack(blocks: Seq[sparray], format: _FmtCOO, dtype: onp.AnyComplex128DType) -> _COOArray2D[np.complex128]: ...
@overload  # sparray, format: "csc", dtype: complex128-like
def vstack(blocks: Seq[sparray], format: _FmtCSC, dtype: onp.AnyComplex128DType) -> csc_array[np.complex128]: ...
@overload  # sparray, format: "csr", dtype: complex128-like
def vstack(blocks: Seq[sparray], format: _FmtCSR, dtype: onp.AnyComplex128DType) -> _CSRArray2D[np.complex128]: ...
@overload  # sparray, format: "dia", dtype: complex128-like
def vstack(blocks: Seq[sparray], format: _FmtDIA, dtype: onp.AnyComplex128DType) -> dia_array[np.complex128]: ...
@overload  # sparray, format: "dok", dtype: complex128-like
def vstack(blocks: Seq[sparray], format: _FmtDOK, dtype: onp.AnyComplex128DType) -> _DOKArray2D[np.complex128]: ...
@overload  # sparray, format: "lil", dtype: complex128-like
def vstack(blocks: Seq[sparray], format: _FmtLIL, dtype: onp.AnyComplex128DType) -> lil_array[np.complex128]: ...

#
@overload  # sparray, format: <default>, dtype: <known>
def vstack[ScalarT: _Numeric, T](
    blocks: Seq[_CanStackAs[ScalarT, T]], format: None = None, *, dtype: onp.ToDType[ScalarT]
) -> T: ...
@overload  # sparray, format: "bsr", dtype: <known>
def vstack[ScalarT: _Numeric](blocks: Seq[sparray], format: _FmtBSR, dtype: onp.ToDType[ScalarT]) -> bsr_array[ScalarT]: ...
@overload  # sparray, format: "coo", dtype: <known>
def vstack[ScalarT: _Numeric](blocks: Seq[sparray], format: _FmtCOO, dtype: onp.ToDType[ScalarT]) -> _COOArray2D[ScalarT]: ...
@overload  # sparray, format: "csc", dtype: <known>
def vstack[ScalarT: _Numeric](blocks: Seq[sparray], format: _FmtCSC, dtype: onp.ToDType[ScalarT]) -> csc_array[ScalarT]: ...
@overload  # sparray, format: "csr", dtype: <known>
def vstack[ScalarT: _Numeric](blocks: Seq[sparray], format: _FmtCSR, dtype: onp.ToDType[ScalarT]) -> _CSRArray2D[ScalarT]: ...
@overload  # sparray, format: "dia", dtype: <known>
def vstack[ScalarT: _Numeric](blocks: Seq[sparray], format: _FmtDIA, dtype: onp.ToDType[ScalarT]) -> dia_array[ScalarT]: ...
@overload  # sparray, format: "dok", dtype: <known>
def vstack[ScalarT: _Numeric](blocks: Seq[sparray], format: _FmtDOK, dtype: onp.ToDType[ScalarT]) -> _DOKArray2D[ScalarT]: ...
@overload  # sparray, format: "lil", dtype: <known>
def vstack[ScalarT: _Numeric](blocks: Seq[sparray], format: _FmtLIL, dtype: onp.ToDType[ScalarT]) -> lil_array[ScalarT]: ...

#
@overload  # sparray, format: <default>, dtype: <unknown>
def vstack[T](blocks: Seq[_CanStackAs[Any, T]], format: None = None, *, dtype: _ToDType) -> T: ...
@overload  # sparray, format: "bsr", dtype: <unknown>
def vstack(blocks: Seq[sparray], format: _FmtBSR, dtype: _ToDType) -> bsr_array: ...
@overload  # sparray, format: "coo", dtype: <unknown>
def vstack(blocks: Seq[sparray], format: _FmtCOO, dtype: _ToDType) -> _COOArray2D[Any]: ...
@overload  # sparray, format: "csc", dtype: <unknown>
def vstack(blocks: Seq[sparray], format: _FmtCSC, dtype: _ToDType) -> csc_array: ...
@overload  # sparray, format: "csr", dtype: <unknown>
def vstack(blocks: Seq[sparray], format: _FmtCSR, dtype: _ToDType) -> _CSRArray2D[Any]: ...
@overload  # sparray, format: "dia", dtype: <unknown>
def vstack(blocks: Seq[sparray], format: _FmtDIA, dtype: _ToDType) -> dia_array: ...
@overload  # sparray, format: "dok", dtype: <unknown>
def vstack(blocks: Seq[sparray], format: _FmtDOK, dtype: _ToDType) -> _DOKArray2D[Any]: ...
@overload  # sparray, format: "lil", dtype: <unknown>
def vstack(blocks: Seq[sparray], format: _FmtLIL, dtype: _ToDType) -> lil_array: ...

#
@overload
def vstack(blocks: Seq[_spbase], format: _Format, dtype: _ToDType | None = None) -> Incomplete: ...

###
@overload  # blocks: <known, known>, format: <default>, dtype: <default>
def block_array[T](blocks: _ToBlocks[_CanStack[T]], *, format: None = None, dtype: None = None) -> T: ...
@overload  # blocks: <array, known>, format: "bsr", dtype: <default>
def block_array[ScalarT: _Numeric](
    blocks: _ToBlocksSpArray[ScalarT], *, format: _FmtBSR, dtype: None = None
) -> bsr_array[ScalarT]: ...
@overload  # blocks: <array, known>, format: "coo", dtype: <default>
def block_array[ScalarT: _Numeric](
    blocks: _ToBlocksSpArray[ScalarT], *, format: _FmtCOO, dtype: None = None
) -> _COOArray2D[ScalarT]: ...
@overload  # blocks: <array, known>, format: "csc", dtype: <default>
def block_array[ScalarT: _Numeric](
    blocks: _ToBlocksSpArray[ScalarT], *, format: _FmtCSC, dtype: None = None
) -> csc_array[ScalarT]: ...
@overload  # blocks: <array, known>, format: "csr", dtype: <default>
def block_array[ScalarT: _Numeric](
    blocks: _ToBlocksSpArray[ScalarT], *, format: _FmtCSR, dtype: None = None
) -> _CSRArray2D[ScalarT]: ...
@overload  # blocks: <array, known>, format: "dia", dtype: <default>
def block_array[ScalarT: _Numeric](
    blocks: _ToBlocksSpArray[ScalarT], *, format: _FmtDIA, dtype: None = None
) -> dia_array[ScalarT]: ...
@overload  # blocks: <array, known>, format: "dok", dtype: <default>
def block_array[ScalarT: _Numeric](
    blocks: _ToBlocksSpArray[ScalarT], *, format: _FmtDOK, dtype: None = None
) -> _DOKArray2D[ScalarT]: ...
@overload  # blocks: <array, known>, format: "lil", dtype: <default>
def block_array[ScalarT: _Numeric](
    blocks: _ToBlocksSpArray[ScalarT], *, format: _FmtLIL, dtype: None = None
) -> lil_array[ScalarT]: ...

#
@overload  # blocks: <known, bool_>, format: <default>, dtype: bool-like
def block_array[T](blocks: _ToBlocksCanStackAs[np.bool, T], *, format: None = None, dtype: onp.AnyBoolDType) -> T: ...
@overload  # blocks: <unknown, unknown>, format: "bsr", dtype: bool-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtBSR, dtype: onp.AnyBoolDType) -> bsr_array[np.bool]: ...
@overload  # blocks: <unknown, unknown>, format: "coo", dtype: bool-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCOO, dtype: onp.AnyBoolDType) -> _COOArray2D[np.bool]: ...
@overload  # blocks: <unknown, unknown>, format: "csc", dtype: bool-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSC, dtype: onp.AnyBoolDType) -> csc_array[np.bool]: ...
@overload  # blocks: <unknown, unknown>, format: "csr", dtype: bool-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSR, dtype: onp.AnyBoolDType) -> _CSRArray2D[np.bool]: ...
@overload  # blocks: <unknown, unknown>, format: "dia", dtype: bool-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDIA, dtype: onp.AnyBoolDType) -> dia_array[np.bool]: ...
@overload  # blocks: <unknown, unknown>, format: "dok", dtype: bool-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDOK, dtype: onp.AnyBoolDType) -> _DOKArray2D[np.bool]: ...
@overload  # blocks: <unknown, unknown>, format: "lil", dtype: bool-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtLIL, dtype: onp.AnyBoolDType) -> lil_array[np.bool]: ...

#
@overload  # blocks: <known, int_>, format: <default>, dtype: int-like
def block_array[T](blocks: _ToBlocksCanStackAs[np.int64, T], *, format: None = None, dtype: onp.AnyIntDType) -> T: ...
@overload  # blocks: <unknown, unknown>, format: "bsr", dtype: int-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtBSR, dtype: onp.AnyIntDType) -> bsr_array[np.int_]: ...
@overload  # blocks: <unknown, unknown>, format: "coo", dtype: int-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCOO, dtype: onp.AnyIntDType) -> _COOArray2D[np.int_]: ...
@overload  # blocks: <unknown, unknown>, format: "csc", dtype: int-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSC, dtype: onp.AnyIntDType) -> csc_array[np.int_]: ...
@overload  # blocks: <unknown, unknown>, format: "csr", dtype: int-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSR, dtype: onp.AnyIntDType) -> _CSRArray2D[np.int_]: ...
@overload  # blocks: <unknown, unknown>, format: "dia", dtype: int-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDIA, dtype: onp.AnyIntDType) -> dia_array[np.int_]: ...
@overload  # blocks: <unknown, unknown>, format: "dok", dtype: int-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDOK, dtype: onp.AnyIntDType) -> _DOKArray2D[np.int_]: ...
@overload  # blocks: <unknown, unknown>, format: "lil", dtype: int-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtLIL, dtype: onp.AnyIntDType) -> lil_array[np.int_]: ...

#
@overload  # blocks: <known, float64>, format: <default>, dtype: float64-like
def block_array[T](blocks: _ToBlocksCanStackAs[np.float64, T], *, format: None = None, dtype: onp.AnyFloat64DType) -> T: ...
@overload  # blocks: <unknown, unknown>, format: "bsr", dtype: float64-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtBSR, dtype: onp.AnyFloat64DType) -> bsr_array[np.float64]: ...
@overload  # blocks: <unknown, unknown>, format: "coo", dtype: float64-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCOO, dtype: onp.AnyFloat64DType) -> _COOArray2D[np.float64]: ...
@overload  # blocks: <unknown, unknown>, format: "csc", dtype: float64-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSC, dtype: onp.AnyFloat64DType) -> csc_array[np.float64]: ...
@overload  # blocks: <unknown, unknown>, format: "csr", dtype: float64-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSR, dtype: onp.AnyFloat64DType) -> _CSRArray2D[np.float64]: ...
@overload  # blocks: <unknown, unknown>, format: "dia", dtype: float64-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDIA, dtype: onp.AnyFloat64DType) -> dia_array[np.float64]: ...
@overload  # blocks: <unknown, unknown>, format: "dok", dtype: float64-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDOK, dtype: onp.AnyFloat64DType) -> _DOKArray2D[np.float64]: ...
@overload  # blocks: <unknown, unknown>, format: "lil", dtype: float64-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtLIL, dtype: onp.AnyFloat64DType) -> lil_array[np.float64]: ...

#
@overload  # blocks: <known, complex128>, format: <default>, dtype: complex128-like
def block_array[T](blocks: _ToBlocksCanStackAs[np.complex128, T], *, format: None = None, dtype: onp.AnyComplex128DType) -> T: ...
@overload  # blocks: <unknown, unknown>, format: "bsr", dtype: complex128-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtBSR, dtype: onp.AnyComplex128DType) -> bsr_array[np.complex128]: ...
@overload  # blocks: <unknown, unknown>, format: "coo", dtype: complex128-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCOO, dtype: onp.AnyComplex128DType) -> _COOArray2D[np.complex128]: ...
@overload  # blocks: <unknown, unknown>, format: "csc", dtype: complex128-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSC, dtype: onp.AnyComplex128DType) -> csc_array[np.complex128]: ...
@overload  # blocks: <unknown, unknown>, format: "csr", dtype: complex128-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSR, dtype: onp.AnyComplex128DType) -> _CSRArray2D[np.complex128]: ...
@overload  # blocks: <unknown, unknown>, format: "dia", dtype: complex128-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDIA, dtype: onp.AnyComplex128DType) -> dia_array[np.complex128]: ...
@overload  # blocks: <unknown, unknown>, format: "dok", dtype: complex128-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDOK, dtype: onp.AnyComplex128DType) -> _DOKArray2D[np.complex128]: ...
@overload  # blocks: <unknown, unknown>, format: "lil", dtype: complex128-like
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtLIL, dtype: onp.AnyComplex128DType) -> lil_array[np.complex128]: ...

#
@overload  # blocks: <known, known>, format: <default>, dtype: <known>
def block_array[ScalarT: _Numeric, T](
    blocks: _ToBlocksCanStackAs[ScalarT, T], *, format: None = None, dtype: onp.ToDType[ScalarT]
) -> T: ...
@overload  # blocks: <unknown, unknown>, format: "bsr", dtype: <known>
def block_array[ScalarT: _Numeric](
    blocks: _ToBlocksUnkown, *, format: _FmtBSR, dtype: onp.ToDType[ScalarT]
) -> bsr_array[ScalarT]: ...
@overload  # blocks: <unknown, unknown>, format: "coo", dtype: <known>
def block_array[ScalarT: _Numeric](
    blocks: _ToBlocksUnkown, *, format: _FmtCOO, dtype: onp.ToDType[ScalarT]
) -> _COOArray2D[ScalarT]: ...
@overload  # blocks: <unknown, unknown>, format: "csc", dtype: <known>
def block_array[ScalarT: _Numeric](
    blocks: _ToBlocksUnkown, *, format: _FmtCSC, dtype: onp.ToDType[ScalarT]
) -> csc_array[ScalarT]: ...
@overload  # blocks: <unknown, unknown>, format: "csr", dtype: <known>
def block_array[ScalarT: _Numeric](
    blocks: _ToBlocksUnkown, *, format: _FmtCSR, dtype: onp.ToDType[ScalarT]
) -> _CSRArray2D[ScalarT]: ...
@overload  # blocks: <unknown, unknown>, format: "dia", dtype: <known>
def block_array[ScalarT: _Numeric](
    blocks: _ToBlocksUnkown, *, format: _FmtDIA, dtype: onp.ToDType[ScalarT]
) -> dia_array[ScalarT]: ...
@overload  # blocks: <unknown>, format: "dok", dtype: <known>
def block_array[ScalarT: _Numeric](
    blocks: _ToBlocksUnkown, *, format: _FmtDOK, dtype: onp.ToDType[ScalarT]
) -> _DOKArray2D[ScalarT]: ...
@overload  # blocks: <unknown, unknown>, format: "lil", dtype: <known>
def block_array[ScalarT: _Numeric](
    blocks: _ToBlocksUnkown, *, format: _FmtLIL, dtype: onp.ToDType[ScalarT]
) -> lil_array[ScalarT]: ...

#
@overload  # blocks: <known, unknown>, format: <default>, dtype: <unknown>
def block_array[T](blocks: _ToBlocksCanStackAs[Any, T], *, format: None = None, dtype: _ToDType | None = None) -> T: ...
@overload  # blocks: <unknown, unknown>, format: "bsr", dtype: <unknown>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtBSR, dtype: _ToDType | None = None) -> bsr_array: ...
@overload  # blocks: <unknown, unknown>, format: "coo", dtype: <unknown>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCOO, dtype: _ToDType | None = None) -> _COOArray2D[Any]: ...
@overload  # blocks: <unknown, unknown>, format: "csc", dtype: <unknown>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSC, dtype: _ToDType | None = None) -> csc_array: ...
@overload  # blocks: <unknown, unknown>, format: "csr", dtype: <unknown>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtCSR, dtype: _ToDType | None = None) -> _CSRArray2D[Any]: ...
@overload  # blocks: <unknown, unknown>, format: "dia", dtype: <unknown>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDIA, dtype: _ToDType | None = None) -> dia_array: ...
@overload  # blocks: <unknown, unknown>, format: "dok", dtype: <unknown>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtDOK, dtype: _ToDType | None = None) -> _DOKArray2D[Any]: ...
@overload  # blocks: <unknown, unknown>, format: "lil", dtype: <unknown>
def block_array(blocks: _ToBlocksUnkown, *, format: _FmtLIL, dtype: _ToDType | None = None) -> lil_array: ...

###
@overload  # blocks: <known, known>, format: <default>, dtype: <default>
def bmat[T](blocks: _ToBlocks[_CanStack[T]], format: None = None, dtype: None = None) -> T: ...
@overload  # blocks: <matrix, known>, format: <otherwise>, dtype: <default>
def bmat[ScalarT: _Numeric](blocks: _ToBlocks[spmatrix[ScalarT]], format: _Format, dtype: None = None) -> spmatrix[ScalarT]: ...

#
@overload  # blocks: <known, known>, format: <default>, dtype: <known>
def bmat[ScalarT: _Numeric, T](
    blocks: _ToBlocksCanStackAs[ScalarT, T], format: None = None, *, dtype: onp.ToDType[ScalarT]
) -> T: ...
@overload  # blocks: <matrix, known>, format: <otherwise>, dtype: <known>
def bmat[ScalarT: _Numeric](
    blocks: _ToBlocks[spmatrix[ScalarT]], format: _Format, dtype: onp.ToDType[ScalarT]
) -> spmatrix[ScalarT]: ...
@overload  # blocks: <unknown, unknown>, format: <otherwise>, dtype: <known>
def bmat[ScalarT: _Numeric](
    blocks: _ToBlocksUnkown, format: _Format, dtype: onp.ToDType[ScalarT]
) -> spmatrix[ScalarT] | _SpArray2D[ScalarT]: ...

#
@overload  # blocks: <known, unknown>, format: <default>, dtype: <unknown>
def bmat[T](blocks: _ToBlocksCanStackAs[Any, T], format: None = None, *, dtype: _ToDType) -> T: ...
@overload  # blocks: <matrix, unknown>, format: <otherwise>, dtype: <unknown>
def bmat(blocks: _ToBlocks[spmatrix[_Numeric]], format: _Format, *, dtype: _ToDType) -> spmatrix: ...
@overload  # blocks: <unknown, unknown>, format: <otherwise>, dtype: <unknown>
def bmat(blocks: _ToBlocksUnkown, format: _Format, *, dtype: _ToDType) -> spmatrix | _SpArray2D[Any]: ...

###
@overload  # mats: <array, known>, format: <default>, dtype: None
def block_diag[ScalarT: _Numeric](
    mats: Iterable[sparray[ScalarT]], format: _FmtCOO | None = None, dtype: None = None
) -> _COOArray2D[ScalarT]: ...
@overload  # mats: <array, known>, format: "bsr", dtype: None
def block_diag[ScalarT: _Numeric](
    mats: Iterable[sparray[ScalarT]], format: _FmtBSR, dtype: None = None
) -> bsr_array[ScalarT]: ...
@overload  # mats: <array, known>, format: "csc", dtype: None
def block_diag[ScalarT: _Numeric](
    mats: Iterable[sparray[ScalarT]], format: _FmtCSC, dtype: None = None
) -> csc_array[ScalarT]: ...
@overload  # mats: <array, known>, format: "csr", dtype: None
def block_diag[ScalarT: _Numeric](
    mats: Iterable[sparray[ScalarT]], format: _FmtCSR, dtype: None = None
) -> _CSRArray2D[ScalarT]: ...
@overload  # mats: <array, known>, format: "dia", dtype: None
def block_diag[ScalarT: _Numeric](
    mats: Iterable[sparray[ScalarT]], format: _FmtDIA, dtype: None = None
) -> dia_array[ScalarT]: ...
@overload  # mats: <array, known>, format: "dok", dtype: None
def block_diag[ScalarT: _Numeric](
    mats: Iterable[sparray[ScalarT]], format: _FmtDOK, dtype: None = None
) -> _DOKArray2D[ScalarT]: ...
@overload  # mats: <array, known>, format: "lil", dtype: None
def block_diag[ScalarT: _Numeric](
    mats: Iterable[sparray[ScalarT]], format: _FmtLIL, dtype: None = None
) -> lil_array[ScalarT]: ...

#
@overload  # mats: <array, unknown>, format: <default>, dtype: bool-like
def block_diag(mats: Iterable[sparray], format: _FmtCOO | None = None, *, dtype: onp.AnyBoolDType) -> _COOArray2D[np.bool]: ...
@overload  # mats: <array, unknown>, format: "bsr", dtype: bool-like
def block_diag(mats: Iterable[sparray], format: _FmtBSR, dtype: onp.AnyBoolDType) -> bsr_array[np.bool]: ...
@overload  # mats: <array, unknown>, format: "csc", dtype: bool-like
def block_diag(mats: Iterable[sparray], format: _FmtCSC, dtype: onp.AnyBoolDType) -> csc_array[np.bool]: ...
@overload  # mats: <array, unknown>, format: "csr", dtype: bool-like
def block_diag(mats: Iterable[sparray], format: _FmtCSR, dtype: onp.AnyBoolDType) -> _CSRArray2D[np.bool]: ...
@overload  # mats: <array, unknown>, format: "dia", dtype: bool-like
def block_diag(mats: Iterable[sparray], format: _FmtDIA, dtype: onp.AnyBoolDType) -> dia_array[np.bool]: ...
@overload  # mats: <array, unknown>, format: "dok", dtype: bool-like
def block_diag(mats: Iterable[sparray], format: _FmtDOK, dtype: onp.AnyBoolDType) -> _DOKArray2D[np.bool]: ...
@overload  # mats: <array, unknown>, format: "lil", dtype: bool-like
def block_diag(mats: Iterable[sparray], format: _FmtLIL, dtype: onp.AnyBoolDType) -> lil_array[np.bool]: ...

#
@overload  # mats: <array, unknown>, format: <default>, dtype: int-like
def block_diag(mats: Iterable[sparray], format: _FmtCOO | None = None, *, dtype: onp.AnyIntDType) -> _COOArray2D[np.int_]: ...
@overload  # mats: <array, unknown>, format: "bsr", dtype: int-like
def block_diag(mats: Iterable[sparray], format: _FmtBSR, dtype: onp.AnyIntDType) -> bsr_array[np.int_]: ...
@overload  # mats: <array, unknown>, format: "csc", dtype: int-like
def block_diag(mats: Iterable[sparray], format: _FmtCSC, dtype: onp.AnyIntDType) -> csc_array[np.int_]: ...
@overload  # mats: <array, unknown>, format: "csr", dtype: int-like
def block_diag(mats: Iterable[sparray], format: _FmtCSR, dtype: onp.AnyIntDType) -> _CSRArray2D[np.int_]: ...
@overload  # mats: <array, unknown>, format: "dia", dtype: int-like
def block_diag(mats: Iterable[sparray], format: _FmtDIA, dtype: onp.AnyIntDType) -> dia_array[np.int_]: ...
@overload  # mats: <array, unknown>, format: "dok", dtype: int-like
def block_diag(mats: Iterable[sparray], format: _FmtDOK, dtype: onp.AnyIntDType) -> _DOKArray2D[np.int_]: ...
@overload  # mats: <array, unknown>, format: "lil", dtype: int-like
def block_diag(mats: Iterable[sparray], format: _FmtLIL, dtype: onp.AnyIntDType) -> lil_array[np.int_]: ...

#
@overload  # mats: <array, unknown>, format: <default>, dtype: float64-like
def block_diag(
    mats: Iterable[sparray], format: _FmtCOO | None = None, *, dtype: onp.AnyFloat64DType
) -> _COOArray2D[np.float64]: ...
@overload  # mats: <array, unknown>, format: "bsr", dtype: float64-like
def block_diag(mats: Iterable[sparray], format: _FmtBSR, dtype: onp.AnyFloat64DType) -> bsr_array[np.float64]: ...
@overload  # mats: <array, unknown>, format: "csc", dtype: float64-like
def block_diag(mats: Iterable[sparray], format: _FmtCSC, dtype: onp.AnyFloat64DType) -> csc_array[np.float64]: ...
@overload  # mats: <array, unknown>, format: "csr", dtype: float64-like
def block_diag(mats: Iterable[sparray], format: _FmtCSR, dtype: onp.AnyFloat64DType) -> _CSRArray2D[np.float64]: ...
@overload  # mats: <array, unknown>, format: "dia", dtype: float64-like
def block_diag(mats: Iterable[sparray], format: _FmtDIA, dtype: onp.AnyFloat64DType) -> dia_array[np.float64]: ...
@overload  # mats: <array, unknown>, format: "dok", dtype: float64-like
def block_diag(mats: Iterable[sparray], format: _FmtDOK, dtype: onp.AnyFloat64DType) -> _DOKArray2D[np.float64]: ...
@overload  # mats: <array, unknown>, format: "lil", dtype: float64-like
def block_diag(mats: Iterable[sparray], format: _FmtLIL, dtype: onp.AnyFloat64DType) -> lil_array[np.float64]: ...

#
@overload  # mats: <array, unknown>, format: <default>, dtype: complex128-like
def block_diag(
    mats: Iterable[sparray], format: _FmtCOO | None = None, *, dtype: onp.AnyComplex128DType
) -> _COOArray2D[np.complex128]: ...
@overload  # mats: <array, unknown>, format: "bsr", dtype: complex128-like
def block_diag(mats: Iterable[sparray], format: _FmtBSR, dtype: onp.AnyComplex128DType) -> bsr_array[np.complex128]: ...
@overload  # mats: <array, unknown>, format: "csc", dtype: complex128-like
def block_diag(mats: Iterable[sparray], format: _FmtCSC, dtype: onp.AnyComplex128DType) -> csc_array[np.complex128]: ...
@overload  # mats: <array, unknown>, format: "csr", dtype: complex128-like
def block_diag(mats: Iterable[sparray], format: _FmtCSR, dtype: onp.AnyComplex128DType) -> _CSRArray2D[np.complex128]: ...
@overload  # mats: <array, unknown>, format: "dia", dtype: complex128-like
def block_diag(mats: Iterable[sparray], format: _FmtDIA, dtype: onp.AnyComplex128DType) -> dia_array[np.complex128]: ...
@overload  # mats: <array, unknown>, format: "dok", dtype: complex128-like
def block_diag(mats: Iterable[sparray], format: _FmtDOK, dtype: onp.AnyComplex128DType) -> _DOKArray2D[np.complex128]: ...
@overload  # mats: <array, unknown>, format: "lil", dtype: complex128-like
def block_diag(mats: Iterable[sparray], format: _FmtLIL, dtype: onp.AnyComplex128DType) -> lil_array[np.complex128]: ...

#
@overload  # mats: <array, unknown>, format: <default>, dtype: <known>
def block_diag[ScalarT: _Numeric](
    mats: Iterable[sparray], format: _FmtCOO | None = None, *, dtype: onp.ToDType[ScalarT]
) -> _COOArray2D[ScalarT]: ...
@overload  # mats: <array, unknown>, format: "bsr", dtype: <known>
def block_diag[ScalarT: _Numeric](
    mats: Iterable[sparray], format: _FmtBSR, dtype: onp.ToDType[ScalarT]
) -> bsr_array[ScalarT]: ...
@overload  # mats: <array, unknown>, format: "csc", dtype: <known>
def block_diag[ScalarT: _Numeric](
    mats: Iterable[sparray], format: _FmtCSC, dtype: onp.ToDType[ScalarT]
) -> csc_array[ScalarT]: ...
@overload  # mats: <array, unknown>, format: "csr", dtype: <known>
def block_diag[ScalarT: _Numeric](
    mats: Iterable[sparray], format: _FmtCSR, dtype: onp.ToDType[ScalarT]
) -> _CSRArray2D[ScalarT]: ...
@overload  # mats: <array, unknown>, format: "dia", dtype: <known>
def block_diag[ScalarT: _Numeric](
    mats: Iterable[sparray], format: _FmtDIA, dtype: onp.ToDType[ScalarT]
) -> dia_array[ScalarT]: ...
@overload  # mats: <array, unknown>, format: "dok", dtype: <known>
def block_diag[ScalarT: _Numeric](
    mats: Iterable[sparray], format: _FmtDOK, dtype: onp.ToDType[ScalarT]
) -> _DOKArray2D[ScalarT]: ...
@overload  # mats: <array, unknown>, format: "lil", dtype: <known>
def block_diag[ScalarT: _Numeric](
    mats: Iterable[sparray], format: _FmtLIL, dtype: onp.ToDType[ScalarT]
) -> lil_array[ScalarT]: ...

#
@overload  # mats: <unknown, known>, format: <default>, dtype: None
def block_diag[ScalarT: _Numeric](
    mats: _ToMatsDiag[ScalarT], format: _FmtCOO | None = None, dtype: None = None
) -> _COOArray2D[ScalarT] | coo_matrix[ScalarT]: ...
@overload  # mats: <unknown, known>, format: "bsr", dtype: None
def block_diag[ScalarT: _Numeric](
    mats: _ToMatsDiag[ScalarT], format: _FmtBSR, dtype: None = None
) -> bsr_array[ScalarT] | bsr_matrix[ScalarT]: ...
@overload  # mats: <unknown, known>, format: "csc", dtype: None
def block_diag[ScalarT: _Numeric](
    mats: _ToMatsDiag[ScalarT], format: _FmtCSC, dtype: None = None
) -> csc_array[ScalarT] | csc_matrix[ScalarT]: ...
@overload  # mats: <unknown, known>, format: "csr", dtype: None
def block_diag[ScalarT: _Numeric](
    mats: _ToMatsDiag[ScalarT], format: _FmtCSR, dtype: None = None
) -> _CSRArray2D[ScalarT] | csr_matrix[ScalarT]: ...
@overload  # mats: <unknown, known>, format: "dia", dtype: None
def block_diag[ScalarT: _Numeric](
    mats: _ToMatsDiag[ScalarT], format: _FmtDIA, dtype: None = None
) -> dia_array[ScalarT] | dia_matrix[ScalarT]: ...
@overload  # mats: <unknown, known>, format: "dok", dtype: None
def block_diag[ScalarT: _Numeric](
    mats: _ToMatsDiag[ScalarT], format: _FmtDOK, dtype: None = None
) -> _DOKArray2D[ScalarT] | dok_matrix[ScalarT]: ...
@overload  # mats: <unknown, known>, format: "lil", dtype: None
def block_diag[ScalarT: _Numeric](
    mats: _ToMatsDiag[ScalarT], format: _FmtLIL, dtype: None = None
) -> lil_array[ScalarT] | lil_matrix[ScalarT]: ...

#
@overload  # mats: <unknown, unknown>, format: <default>, dtype: <known>
def block_diag[ScalarT: _Numeric](
    mats: _ToMatsDiagUnknown, format: _FmtCOO | None = None, *, dtype: onp.ToDType[ScalarT]
) -> _COOArray2D[ScalarT] | coo_matrix[ScalarT]: ...
@overload  # mats: <unknown, unknown>, format: "bsr", dtype: <known>
def block_diag[ScalarT: _Numeric](
    mats: _ToMatsDiagUnknown, format: _FmtBSR, dtype: onp.ToDType[ScalarT]
) -> bsr_array[ScalarT] | bsr_matrix[ScalarT]: ...
@overload  # mats: <unknown, unknown>, format: "csc", dtype: <known>
def block_diag[ScalarT: _Numeric](
    mats: _ToMatsDiagUnknown, format: _FmtCSC, dtype: onp.ToDType[ScalarT]
) -> csc_array[ScalarT] | csc_matrix[ScalarT]: ...
@overload  # mats: <unknown, unknown>, format: "csr", dtype: <known>
def block_diag[ScalarT: _Numeric](
    mats: _ToMatsDiagUnknown, format: _FmtCSR, dtype: onp.ToDType[ScalarT]
) -> _CSRArray2D[ScalarT] | csr_matrix[ScalarT]: ...
@overload  # mats: <unknown, unknown>, format: "dia", dtype: <known>
def block_diag[ScalarT: _Numeric](
    mats: _ToMatsDiagUnknown, format: _FmtDIA, dtype: onp.ToDType[ScalarT]
) -> dia_array[ScalarT] | dia_matrix[ScalarT]: ...
@overload  # mats: <unknown, unknown>, format: "dok", dtype: <known>
def block_diag[ScalarT: _Numeric](
    mats: _ToMatsDiagUnknown, format: _FmtDOK, dtype: onp.ToDType[ScalarT]
) -> _DOKArray2D[ScalarT] | dok_matrix[ScalarT]: ...
@overload  # mats: <unknown, unknown>, format: "lil", dtype: <known>
def block_diag[ScalarT: _Numeric](
    mats: _ToMatsDiagUnknown, format: _FmtLIL, dtype: onp.ToDType[ScalarT]
) -> lil_array[ScalarT] | lil_matrix[ScalarT]: ...

#
@overload  # mats: <unknown, unknown>, format: <default>, dtype: <unknown>
def block_diag(
    mats: _ToMatsDiagUnknown, format: _FmtCOO | None = None, dtype: _ToDType | None = None
) -> _COOArray2D[Any] | coo_matrix: ...
@overload  # mats: <unknown, unknown>, format: "bsr", dtype: <unknown>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtBSR, dtype: _ToDType | None = None) -> bsr_array | bsr_matrix: ...
@overload  # mats: <unknown, unknown>, format: "csc", dtype: <unknown>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtCSC, dtype: _ToDType | None = None) -> csc_array | csc_matrix: ...
@overload  # mats: <unknown, unknown>, format: "csr", dtype: <unknown>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtCSR, dtype: _ToDType | None = None) -> _CSRArray2D[Any] | csr_matrix: ...
@overload  # mats: <unknown, unknown>, format: "dia", dtype: <unknown>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtDIA, dtype: _ToDType | None = None) -> dia_array | dia_matrix: ...
@overload  # mats: <unknown, unknown>, format: "dok", dtype: <unknown>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtDOK, dtype: _ToDType | None = None) -> _DOKArray2D[Any] | dok_matrix: ...
@overload  # mats: <unknown, unknown>, format: "lil", dtype: <unknown>
def block_diag(mats: _ToMatsDiagUnknown, format: _FmtLIL, dtype: _ToDType | None = None) -> lil_array | lil_matrix: ...

###
@overload  # shape: T, format: <default>, dtype: <default>
def random_array[ShapeT: tuple[int, *tuple[int, ...]]](
    shape: ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    dtype: onp.AnyFloat64DType | None = None,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> coo_array[np.float64, ShapeT]: ...
@overload  # shape: T, format: <otherwise>, dtype: <default>
def random_array[ShapeT: tuple[int, *tuple[int, ...]]](
    shape: ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtNonCOO,
    dtype: onp.AnyFloat64DType | None = None,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> sparray[np.float64, ShapeT]: ...

#
@overload  # shape: T, format: <default>, dtype: <known>
def random_array[ShapeT: tuple[int, *tuple[int, ...]], ScalarT: _Numeric](
    shape: ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    dtype: onp.ToDType[ScalarT],
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> coo_array[ScalarT, ShapeT]: ...
@overload  # shape: T, format: <otherwise>, dtype: <known>
def random_array[ShapeT: tuple[int, *tuple[int, ...]], ScalarT: _Numeric](
    shape: ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtNonCOO,
    dtype: onp.ToDType[ScalarT],
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> sparray[ScalarT, ShapeT]: ...

#
@overload  # shape: T, format: <default>, dtype: complex
def random_array[ShapeT: tuple[int, *tuple[int, ...]]](
    shape: ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    dtype: onp.AnyComplex128DType,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> coo_array[np.complex128, ShapeT]: ...
@overload  # shape: T, format: <otherwise>, dtype: complex
def random_array[ShapeT: tuple[int, *tuple[int, ...]]](
    shape: ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtNonCOO,
    dtype: onp.AnyComplex128DType,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> sparray[np.complex128, ShapeT]: ...

#
@overload  # shape: T, format: <default>, dtype: <unknown>
def random_array[ShapeT: tuple[int, *tuple[int, ...]]](
    shape: ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    dtype: _ToDType,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> coo_array[Any, ShapeT]: ...
@overload  # shape: T, format: <otherwise>, dtype: <unknown>
def random_array[ShapeT: tuple[int, *tuple[int, ...]]](
    shape: ShapeT,
    *,
    density: float | npc.floating = 0.01,
    format: _FmtNonCOO,
    dtype: _ToDType,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
    data_sampler: _DataSampler | None = None,
) -> sparray[Any, ShapeT]: ...

###
# NOTE: `random_array` should be prefered over `random`
@overload  # format: <default>, dtype: <default>
def random(
    m: int,
    n: int,
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
    m: int,
    n: int,
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
def random[ScalarT: _Numeric](
    m: int,
    n: int,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    *,
    dtype: onp.ToDType[ScalarT],
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[ScalarT]: ...
@overload  # format: <otherwise>, dtype: <known> (keyword)
def random[ScalarT: _Numeric](
    m: int,
    n: int,
    density: float | npc.floating = 0.01,
    *,
    format: _FmtNonCOO,
    dtype: onp.ToDType[ScalarT],
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[ScalarT]: ...

#
@overload  # format: <default>, dtype: <known> (positional)
def random[ScalarT: _Numeric](
    m: int,
    n: int,
    density: float | npc.floating,
    format: _FmtCOO,
    dtype: onp.ToDType[ScalarT],
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[ScalarT]: ...
@overload  # format: <otherwise>, dtype: <known> (positional)
def random[ScalarT: _Numeric](
    m: int,
    n: int,
    density: float | npc.floating,
    format: _FmtNonCOO,
    dtype: onp.ToDType[ScalarT],
    rng: onp.random.ToRNG | None = None,
    data_rvs: _DataRVS | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[ScalarT]: ...

#
@overload  # format: <default>, dtype: complex (keyword)
def random(
    m: int,
    n: int,
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
    m: int,
    n: int,
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
    m: int,
    n: int,
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
    m: int,
    n: int,
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
    m: int,
    n: int,
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
    m: int,
    n: int,
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
    m: int,
    n: int,
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
    m: int,
    n: int,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    dtype: onp.AnyFloat64DType | None = None,
    rng: onp.random.ToRNG | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[np.float64]: ...
@overload  # format: <otherwise>, dtype: <default>
def rand(
    m: int,
    n: int,
    density: float | npc.floating = 0.01,
    *,
    format: _FmtNonCOO,
    dtype: onp.AnyFloat64DType | None = None,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[np.float64]: ...

#
@overload  # format: <default>, dtype: <known> (keyword)
def rand[ScalarT: _Numeric](
    m: int,
    n: int,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    *,
    dtype: onp.ToDType[ScalarT],
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[ScalarT]: ...
@overload  # format: <otherwise>, dtype: <known> (keyword)
def rand[ScalarT: _Numeric](
    m: int,
    n: int,
    density: float | npc.floating = 0.01,
    *,
    format: _FmtNonCOO,
    dtype: onp.ToDType[ScalarT],
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[ScalarT]: ...

#
@overload  # format: <default>, dtype: <known> (positional)
def rand[ScalarT: _Numeric](
    m: int,
    n: int,
    density: float | npc.floating,
    format: _FmtCOO,
    dtype: onp.ToDType[ScalarT],
    rng: onp.random.ToRNG | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[ScalarT]: ...
@overload  # format: <otherwise>, dtype: <known> (positional)
def rand[ScalarT: _Numeric](
    m: int,
    n: int,
    density: float | npc.floating,
    format: _FmtNonCOO,
    dtype: onp.ToDType[ScalarT],
    rng: onp.random.ToRNG | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix[ScalarT]: ...

#
@overload  # format: <default>, dtype: complex (keyword)
def rand(
    m: int,
    n: int,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    *,
    dtype: onp.AnyComplex128DType,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[np.complex128]: ...
@overload  # format: <otherwise>, dtype: complex (keyword)
def rand(
    m: int,
    n: int,
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
    m: int,
    n: int,
    density: float | npc.floating,
    format: _FmtCOO,
    dtype: onp.AnyComplex128DType,
    rng: onp.random.ToRNG | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix[np.complex128]: ...
@overload  # format: <otherwise>, dtype: complex (positional)
def rand(
    m: int,
    n: int,
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
    m: int,
    n: int,
    density: float | npc.floating = 0.01,
    format: _FmtCOO = "coo",
    dtype: _ToDType | None = None,
    rng: onp.random.ToRNG | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> coo_matrix: ...
@overload  # format: <otherwise> (keyword), dtype: <unknown>
def rand(
    m: int,
    n: int,
    density: float | npc.floating = 0.01,
    *,
    format: _FmtNonCOO,
    dtype: _ToDType | None = None,
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix: ...
@overload  # format: <otherwise> (positional), dtype: <unknown>
def rand(
    m: int,
    n: int,
    density: float | npc.floating,
    format: _FmtNonCOO,
    dtype: _ToDType | None = None,
    rng: onp.random.ToRNG | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> spmatrix: ...

@type_check_only
class _HasMT[MT](Protocol):
    @property
    def mT(self) -> MT: ...

#
def matrix_transpose[MT](A: _HasMT[MT]) -> MT: ...
