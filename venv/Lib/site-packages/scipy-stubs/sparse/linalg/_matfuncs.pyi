from typing import Any, Final, Generic, Literal, Self, SupportsIndex, overload, override
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._interface import LinearOperator
from scipy.sparse import (
    bsr_array,
    bsr_matrix,
    coo_array,
    coo_matrix,
    csc_array,
    csc_matrix,
    csr_array,
    csr_matrix,
    dia_array,
    dia_matrix,
    dok_array,
    dok_matrix,
    lil_array,
    lil_matrix,
)
from scipy.sparse._base import _spbase

__all__ = ["expm", "inv", "matrix_power"]

###

type _Structure = Literal["upper_triangular"]

type _CS[ST: npc.number | np.bool] = csc_array[ST] | csc_matrix[ST] | csr_array[ST] | csr_matrix[ST]
type _NonCS[ST: npc.number | np.bool] = (
    bsr_array[ST]
    | bsr_matrix[ST]
    | coo_array[ST]
    | coo_matrix[ST]
    | dok_array[ST]
    | dok_matrix[ST]
    | dia_array[ST]
    | dia_matrix[ST]
    | lil_array[ST]
    | lil_matrix[ST]
)
type _ToCSCArray[ST: npc.number | np.bool] = bsr_array[ST] | bsr_matrix[ST] | csc_array[ST] | dia_array[ST] | dia_matrix[ST]
type _ToCSRArray[ST: npc.number | np.bool] = coo_array[ST] | csr_array[ST] | dok_array[ST] | lil_array[ST]
type _ToCSRMatrix[ST: npc.number | np.bool] = coo_matrix[ST] | csr_matrix[ST] | dok_matrix[ST] | lil_matrix[ST]
type _ToSelf[ST: npc.number | np.bool] = (
    bsr_array[ST] | bsr_matrix[ST] | csc_array[ST] | csc_matrix[ST] | dia_array[ST] | dia_matrix[ST]
)

type _SubF64 = npc.integer64 | npc.integer32
type _SubF32 = npc.integer16 | npc.integer8 | np.bool
type _SubF = npc.integer | np.bool

_SCT_co = TypeVar("_SCT_co", bound=npc.number | np.bool, default=Any, covariant=True)

###

UPPER_TRIANGULAR: Final[_Structure] = "upper_triangular"

class MatrixPowerOperator(LinearOperator[_SCT_co], Generic[_SCT_co]):
    @property
    @override
    # pyrefly: ignore [bad-override]
    def T(self, /) -> Self: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]
    def __init__(self, /, A: onp.Array2D[_SCT_co] | _spbase, p: int, structure: _Structure | None = None) -> None: ...

class ProductOperator(LinearOperator[_SCT_co], Generic[_SCT_co]):
    @property
    @override
    # pyrefly: ignore [bad-override]
    def T(self, /) -> Self: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]
    def __init__(self, /, *args: onp.Array2D[_SCT_co] | _spbase, structure: _Structure | None = None) -> None: ...

#
@overload
def inv[ScalarT: npc.inexact](A: _NonCS[ScalarT]) -> csc_array[ScalarT]: ...
@overload
def inv[SparseT: _CS[npc.inexact]](A: SparseT) -> SparseT: ...
@overload
def inv(A: _NonCS[_SubF64] | csc_array[_SubF64]) -> csc_array[np.float64]: ...  # type: ignore[overload-overlap]
@overload
def inv(A: _NonCS[_SubF32] | csc_array[_SubF32]) -> csc_array[np.float32]: ...
@overload
def inv(A: csr_array[_SubF64]) -> csr_array[np.float64]: ...  # type: ignore[overload-overlap]
@overload
def inv(A: csr_array[_SubF32]) -> csr_array[np.float32]: ...
@overload
def inv(A: csr_matrix[_SubF64]) -> csr_matrix[np.float64]: ...  # type: ignore[overload-overlap]
@overload
def inv(A: csr_matrix[_SubF32]) -> csr_matrix[np.float32]: ...
@overload
def inv(A: csc_matrix[_SubF64]) -> csc_matrix[np.float64]: ...  # type: ignore[overload-overlap]
@overload
def inv(A: csc_matrix[_SubF32]) -> csc_matrix[np.float32]: ...

#
@overload
def expm[ScalarT: npc.inexact](A: _ToCSCArray[ScalarT]) -> csc_array[ScalarT]: ...
@overload
def expm(A: _ToCSCArray[_SubF]) -> csc_array[np.float64]: ...
@overload
def expm[ScalarT: npc.inexact](A: _ToCSRArray[ScalarT]) -> csr_array[ScalarT]: ...
@overload
def expm(A: _ToCSRArray[_SubF]) -> csr_array[np.float64]: ...
@overload
def expm[ScalarT: npc.inexact](A: _ToCSRMatrix[ScalarT]) -> csr_matrix[ScalarT]: ...
@overload
def expm(A: _ToCSRMatrix[_SubF]) -> csr_matrix[np.float64]: ...
@overload
def expm[ScalarT: npc.inexact](A: csc_matrix[ScalarT]) -> csc_matrix[ScalarT]: ...
@overload
def expm(A: csc_matrix[_SubF]) -> csc_matrix[np.float64]: ...

#
@overload
def matrix_power[ScalarT: npc.number | np.bool](A: _ToCSRArray[ScalarT], power: SupportsIndex) -> csr_array[ScalarT]: ...
@overload
def matrix_power[ScalarT: npc.number | np.bool](A: _ToCSRMatrix[ScalarT], power: SupportsIndex) -> csr_matrix[ScalarT]: ...
@overload
def matrix_power[SparseT: _ToSelf[npc.number | np.bool]](A: SparseT, power: SupportsIndex) -> SparseT: ...
