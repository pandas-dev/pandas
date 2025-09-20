from types import EllipsisType
from typing import Any, Generic, Self, TypeAlias, overload
from typing_extensions import Buffer, TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._base import _spbase, sparray
from ._coo import coo_array
from ._matrix import spmatrix

_ScalarT_co = TypeVar("_ScalarT_co", bound=npc.number | np.bool_, default=Any, covariant=True)
_ShapeT_co = TypeVar("_ShapeT_co", bound=onp.AtLeast1D, default=onp.AtLeast0D[Any], covariant=True)

_Self2T = TypeVar("_Self2T", bound=IndexMixin[Any, _2D])
_SelfMatrixT = TypeVar("_SelfMatrixT", bound=spmatrix)

_ToNumber: TypeAlias = op.CanComplex | op.CanFloat | op.CanInt | op.CanIndex | str | Buffer

_1D: TypeAlias = tuple[int]  # noqa: PYI042
_2D: TypeAlias = tuple[int, int]  # noqa: PYI042

_ToIndex1: TypeAlias = op.CanIndex | tuple[op.CanIndex]
_ToIndex2: TypeAlias = tuple[op.CanIndex, op.CanIndex]
_ToIndex1Of2: TypeAlias = _ToIndex1 | tuple[op.CanIndex, _ToSlice | onp.ToInt1D] | tuple[_ToSlice | onp.ToInt1D, op.CanIndex]
_ToIndex2Of2: TypeAlias = tuple[onp.ToInt1D, onp.ToInt1D]

_ToSlice: TypeAlias = slice | EllipsisType
_ToSlice1: TypeAlias = _ToSlice | onp.ToInt1D | tuple[op.CanIndex, None] | tuple[None, op.CanIndex]
_ToSlice2: TypeAlias = _ToSlice | _spbase[np.bool_, _2D] | list[np.bool_] | list[bool] | list[int] | onp.ToBool2D

###

INT_TYPES: tuple[type[int], type[npc.integer]] = ...

class IndexMixin(Generic[_ScalarT_co, _ShapeT_co]):
    @overload
    def __getitem__(self: IndexMixin[Any, _1D], ix: op.CanIndex, /) -> _ScalarT_co: ...
    @overload
    def __getitem__(self: IndexMixin[Any, _1D], ix: None, /) -> coo_array[_ScalarT_co, tuple[int, int]]: ...
    @overload
    def __getitem__(self: IndexMixin[Any, _2D], ix: _ToIndex2, /) -> _ScalarT_co: ...
    @overload
    def __getitem__(self, ixs: _ToSlice1, /) -> Self: ...
    @overload
    def __getitem__(self: _Self2T, ixs: _ToSlice2, /) -> _Self2T: ...
    @overload  # indexing one axis of a csc/csr/dok/lil matrix returns another csc/csr matrix with the first axis set to 1.
    def __getitem__(self: _SelfMatrixT, ixs: _ToIndex1Of2, /) -> _SelfMatrixT: ...  # type: ignore[misc]
    @overload  # indexing one axis of a csc/csr/dok/lil array returns a 1d coo array
    def __getitem__(self: sparray, ixs: _ToIndex1Of2, /) -> coo_array[_ScalarT_co, tuple[int]]: ...  # type: ignore[misc]
    @overload  # multiindexing two axes of a numpy matrix matrix returns a 1d numpy.matrix
    def __getitem__(self: spmatrix, ixs: _ToIndex2Of2, /) -> onp.Matrix[_ScalarT_co]: ...  # type: ignore[misc]
    @overload  # multiindexing two axes of a 2d sparse array returns a 1d ndarray
    def __getitem__(self: sparray[Any, _2D], ixs: _ToIndex2Of2, /) -> onp.Array1D[_ScalarT_co]: ...  # type: ignore[misc]

    #
    @overload
    def __setitem__(self: IndexMixin[Any, _1D], ix: _ToIndex1, x: _ToNumber, /) -> None: ...
    @overload
    def __setitem__(self: IndexMixin[Any, _2D], ix: _ToIndex1Of2 | _ToIndex2, x: _ToNumber, /) -> None: ...
