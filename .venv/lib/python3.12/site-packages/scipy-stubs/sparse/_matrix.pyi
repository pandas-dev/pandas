# needed (once) for `numpy>=2.2.0`
# mypy: disable-error-code="overload-overlap"

from collections.abc import Sequence
from types import GenericAlias
from typing import Any, Generic, Literal as L, Self, SupportsIndex, TypeAlias, overload, type_check_only
from typing_extensions import TypeAliasType, TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._base import _spbase
from ._bsr import bsr_matrix
from ._coo import coo_matrix
from ._csc import csc_matrix
from ._csr import csr_matrix
from ._dia import dia_matrix
from ._dok import dok_matrix
from ._lil import lil_matrix
from ._typing import _Format

_T = TypeVar("_T")
_ScalarT = TypeVar("_ScalarT", bound=npc.number | np.bool_)
_ScalarT_co = TypeVar("_ScalarT_co", bound=npc.number | np.bool_, default=Any, covariant=True)

_SpMatrixT = TypeVar("_SpMatrixT", bound=spmatrix)

_SpFromInT = TypeVar("_SpFromInT", bound=spmatrix[npc.number])
_SpFromFloatT = TypeVar("_SpFromFloatT", bound=spmatrix[npc.inexact])
_SpFromComplexT = TypeVar("_SpFromComplexT", bound=spmatrix[npc.complexfloating])

_ToInt8: TypeAlias = np.int8 | np.bool_
_ToInt: TypeAlias = npc.integer | np.bool_
_ToFloat32: TypeAlias = np.float32 | _ToInt
_ToFloat: TypeAlias = npc.floating | _ToInt
_ToComplex64: TypeAlias = np.complex64 | _ToFloat

_DualMatrixLike = TypeAliasType("_DualMatrixLike", _T | _ScalarT | _spbase[_ScalarT], type_params=(_T, _ScalarT))
_DualArrayLike = TypeAliasType(
    "_DualArrayLike",
    Sequence[Sequence[_T | _ScalarT] | onp.CanArrayND[_ScalarT]] | onp.CanArrayND[_ScalarT],
    type_params=(_T, _ScalarT),
)

_SpMatrixOut: TypeAlias = bsr_matrix[_ScalarT] | csc_matrix[_ScalarT] | csr_matrix[_ScalarT]

_StackedSparseMatrix: TypeAlias = coo_matrix[_ScalarT] | csc_matrix[_ScalarT] | csr_matrix[_ScalarT]

###

class spmatrix(Generic[_ScalarT_co]):
    # NOTE: These two methods do not exist at runtime.
    # See the relevant comment in `sparse._base._spbase` for more information.
    @type_check_only
    def __assoc_stacked__(self, /) -> _StackedSparseMatrix[_ScalarT_co]: ...
    @type_check_only
    def __assoc_stacked_as__(self, sctype: _ScalarT, /) -> _StackedSparseMatrix[_ScalarT]: ...

    #
    @property
    def _bsr_container(self, /) -> bsr_matrix[_ScalarT_co]: ...
    @property
    def _coo_container(self, /) -> coo_matrix[_ScalarT_co]: ...
    @property
    def _csc_container(self, /) -> csc_matrix[_ScalarT_co]: ...
    @property
    def _csr_container(self, /) -> csr_matrix[_ScalarT_co]: ...
    @property
    def _dia_container(self, /) -> dia_matrix[_ScalarT_co]: ...
    @property
    def _dok_container(self, /) -> dok_matrix[_ScalarT_co]: ...
    @property
    def _lil_container(self, /) -> lil_matrix[_ScalarT_co]: ...

    #
    @property
    def shape(self, /) -> tuple[int, int]: ...
    def get_shape(self, /) -> tuple[int, int]: ...
    def set_shape(self, /, shape: tuple[SupportsIndex, SupportsIndex]) -> None: ...

    #
    @overload  # Self[-Bool], other: scalar-like +Bool
    def __mul__(self, other: bool | np.bool_, /) -> Self: ...
    @overload  # Self[-Int], other: scalar-like +Int
    def __mul__(self: _SpFromInT, other: onp.ToInt, /) -> _SpFromInT: ...
    @overload  # Self[-Float], other: scalar-like +Float
    def __mul__(self: _SpFromFloatT, other: onp.ToFloat, /) -> _SpFromFloatT: ...
    @overload  # Self[-Complex], other: scalar-like +Complex
    def __mul__(self: _SpFromComplexT, other: onp.ToComplex, /) -> _SpFromComplexT: ...
    @overload  # spmatrix, other: spmatrix
    def __mul__(self: _SpMatrixT, other: _SpMatrixT, /) -> _SpMatrixT: ...
    @overload  # spmatrix[-Bool], other: sparse +Bool
    def __mul__(self: spmatrix, other: _spbase[np.bool_], /) -> _SpMatrixOut[_ScalarT_co]: ...
    @overload  # spmatrix[-Bool], other: array-like +Bool
    def __mul__(self: spmatrix, other: _DualArrayLike[bool, np.bool_], /) -> onp.Array2D[_ScalarT_co]: ...
    @overload  # spmatrix[-Int], other: sparse +Int
    def __mul__(self: spmatrix[npc.number], other: _spbase[_ToInt8], /) -> _SpMatrixOut[_ScalarT_co]: ...
    @overload  # spmatrix[-Int], other: array-like +Int
    def __mul__(self: spmatrix[npc.number], other: _DualArrayLike[bool, _ToInt8], /) -> onp.Array2D[_ScalarT_co]: ...
    @overload  # spmatrix[-Float], other: sparse +Float
    def __mul__(self: spmatrix[npc.inexact], other: _spbase[_ToFloat32 | _ScalarT_co], /) -> _SpMatrixOut[_ScalarT_co]: ...
    @overload  # spmatrix[-Float], other: array-like +Float
    def __mul__(self: spmatrix[npc.inexact], other: _DualArrayLike[int, _ToFloat32], /) -> onp.Array2D[_ScalarT_co]: ...
    @overload  # spmatrix[-Complex], other: sparse +Complex
    def __mul__(
        self: spmatrix[npc.complexfloating], other: _spbase[_ToComplex64 | _ScalarT_co], /
    ) -> _SpMatrixOut[_ScalarT_co]: ...
    @overload  # spmatrix[-Complex], other: array-like +Complex
    def __mul__(
        self: spmatrix[npc.complexfloating], other: _DualArrayLike[float, _ToComplex64], /
    ) -> onp.Array2D[_ScalarT_co]: ...
    @overload  # spmatrix[+Bool], other: scalar- or matrix-like ~Int
    def __mul__(self: spmatrix[np.bool_], other: _DualMatrixLike[op.JustInt, npc.integer], /) -> spmatrix[npc.integer]: ...
    @overload  # spmatrix[+Bool], other: array-like ~Int
    def __mul__(self: spmatrix[np.bool_], other: _DualArrayLike[op.JustInt, npc.integer], /) -> onp.Array2D[npc.integer]: ...
    @overload  # spmatrix[+Int], other: scalar- or matrix-like ~Float
    def __mul__(self: spmatrix[_ToInt], other: _DualMatrixLike[op.JustFloat, npc.floating], /) -> spmatrix[npc.floating]: ...
    @overload  # spmatrix[+Int], other: array-like ~Float
    def __mul__(self: spmatrix[_ToInt], other: _DualArrayLike[op.JustFloat, npc.floating], /) -> onp.Array2D[npc.floating]: ...
    @overload  # spmatrix[+Float], other: scalar- or matrix-like ~Complex
    def __mul__(
        self: spmatrix[_ToFloat], other: _DualMatrixLike[op.JustComplex, npc.complexfloating], /
    ) -> spmatrix[npc.complexfloating]: ...
    @overload  # spmatrix[+Float], other: array-like ~Complex
    def __mul__(
        self: spmatrix[_ToFloat], other: _DualArrayLike[op.JustComplex, npc.complexfloating], /
    ) -> onp.Array2D[npc.complexfloating]: ...
    @overload  # catch-all
    def __mul__(
        self, other: _DualArrayLike[complex, npc.number | np.bool_] | _spbase, /
    ) -> _spbase[Any, Any] | onp.Array[Any, Any]: ...
    __rmul__ = __mul__

    #
    def __pow__(self, rhs: op.CanIndex, /) -> Self: ...

    #
    def getmaxprint(self, /) -> int: ...
    def getformat(self, /) -> _Format: ...
    # NOTE: `axis` is only supported by `{coo,csc,csr,lil}_matrix`
    def getnnz(self, /, axis: None = None) -> int: ...
    def getH(self, /) -> Self: ...
    def getcol(self, /, j: onp.ToJustInt) -> csc_matrix[_ScalarT_co]: ...
    def getrow(self, /, i: onp.ToJustInt) -> csr_matrix[_ScalarT_co]: ...

    # NOTE: mypy reports a false positive for overlapping overloads
    @overload
    def asfptype(self: spmatrix[np.bool_ | npc.integer8 | npc.integer16], /) -> spmatrix[np.float32]: ...
    @overload
    def asfptype(self: spmatrix[npc.integer32 | npc.integer64], /) -> spmatrix[np.float64]: ...
    @overload
    def asfptype(self, /) -> Self: ...

    #
    @overload
    def todense(self, /, order: L["C", "F"] | None = None, out: None = None) -> onp.Matrix[_ScalarT_co]: ...
    @overload
    def todense(self, /, order: L["C", "F"] | None, out: onp.ArrayND[_ScalarT]) -> onp.Matrix[_ScalarT]: ...
    @overload
    def todense(self, /, order: L["C", "F"] | None = None, *, out: onp.ArrayND[_ScalarT]) -> onp.Matrix[_ScalarT]: ...

    #
    @classmethod
    def __class_getitem__(cls, arg: type | object, /) -> GenericAlias: ...
