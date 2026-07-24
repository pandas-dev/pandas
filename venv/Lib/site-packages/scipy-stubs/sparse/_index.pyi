from collections.abc import Buffer
from types import EllipsisType
from typing import Any, Generic, Self, SupportsComplex, SupportsFloat, SupportsIndex, SupportsInt, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._base import _spbase, sparray
from ._coo import coo_array
from ._matrix import spmatrix

###

type _ToNumber = SupportsIndex | SupportsInt | SupportsFloat | SupportsComplex | str | Buffer

type _1D = tuple[int]  # noqa: PYI042
type _2D = tuple[int, int]  # noqa: PYI042

# 1d simple index
type _ToIndex1 = int | tuple[int]
# 2d simple index
type _ToIndex2 = tuple[int, int]
# 1d multi-index of a 2d array
type _ToIndex1Of2 = _ToIndex1 | tuple[int, _ToSlice | onp.ToInt1D] | tuple[_ToSlice | onp.ToInt1D, int]
# 2d multi-index of a 2d array
type _ToIndex2Of2 = tuple[onp.ToInt1D, onp.ToInt1D]

# single slice-like index (maintains axis)
type _ToSlice = slice[int | None, int | None, int | None] | EllipsisType
# axis-wise slice (maintains 1d or 2d shape)
type _ToSlice1 = (
    _ToSlice | onp.ToInt1D | tuple[int, None] | tuple[None, int] | tuple[_ToSlice, onp.ToInt1D] | tuple[onp.ToInt1D, _ToSlice]
)
# axis-wise slice for 2d arrays (maintains only 2d shape)
type _ToSlice2 = _ToSlice | tuple[_ToSlice, _ToSlice] | _spbase[np.bool, _2D] | list[np.bool] | list[bool] | list[int]

_ScalarT = TypeVar("_ScalarT", bound=npc.number | np.bool)
_ScalarT_co = TypeVar("_ScalarT_co", bound=npc.number | np.bool, default=Any, covariant=True)
_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, *tuple[int, ...]], default=tuple[Any, ...], covariant=True)

_Self2T = TypeVar("_Self2T", bound=IndexMixin[Any, _2D])
_SelfMatrixT = TypeVar("_SelfMatrixT", bound=spmatrix)

###

INT_TYPES: tuple[type[int], type[npc.integer]] = ...

class IndexMixin(Generic[_ScalarT_co, _ShapeT_co]):
    @overload
    def __getitem__(self: IndexMixin[Any, _1D], ix: int, /) -> _ScalarT_co: ...
    @overload
    def __getitem__(self: IndexMixin[Any, _1D], ix: None, /) -> coo_array[_ScalarT_co, tuple[int, int]]: ...
    @overload
    def __getitem__(self: IndexMixin[Any, _2D], ix: _ToIndex2, /) -> _ScalarT_co: ...
    @overload
    def __getitem__(self, ixs: _ToSlice1, /) -> Self: ...
    @overload
    def __getitem__(self: _Self2T, ixs: _ToSlice2, /) -> _Self2T: ...
    @overload
    def __getitem__(self: _SelfMatrixT, ixs: _ToIndex1Of2, /) -> _SelfMatrixT: ...  # type: ignore[misc]
    @overload
    def __getitem__(self: sparray[_ScalarT], ixs: _ToIndex1Of2, /) -> coo_array[_ScalarT, tuple[int]]: ...  # type: ignore[misc]
    @overload
    def __getitem__(self: spmatrix, ixs: _ToIndex2Of2, /) -> onp.Matrix[_ScalarT_co]: ...  # type: ignore[misc]
    @overload
    def __getitem__(self: sparray[_ScalarT, _2D], ixs: _ToIndex2Of2, /) -> onp.Array1D[_ScalarT]: ...  # type: ignore[misc]

    #
    @overload
    def __setitem__(self: IndexMixin[Any, _1D], ix: _ToIndex1, x: _ToNumber, /) -> None: ...
    @overload
    def __setitem__(
        self: IndexMixin[Any, _2D], ix: _ToIndex1Of2 | _ToIndex2Of2 | _ToIndex2 | _ToSlice1 | _ToSlice2, x: _ToNumber, /
    ) -> None: ...
