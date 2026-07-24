import sys
from typing import TYPE_CHECKING, Any, Literal, Protocol, Self, overload

if sys.version_info >= (3, 13):
    from types import CapsuleType
    from typing import TypeVar, runtime_checkable
else:
    from typing_extensions import CapsuleType, TypeVar, runtime_checkable

import numpy as np

from optype._utils import set_module

__all__ = ["Scalar"]


###

_PT_co = TypeVar("_PT_co", covariant=True)
_NB_co = TypeVar("_NB_co", bound=int, default=int, covariant=True)
_DT = TypeVar("_DT", bound=np.dtype[Any], default=np.dtype[Any])

type _L0 = Literal[0]
type _L1 = Literal[1]
type _Array0D[DT: np.dtype[Any]] = np.ndarray[tuple[()], DT]


###


@runtime_checkable
@set_module("optype.numpy")
class Scalar(Protocol[_PT_co, _NB_co]):
    """
    A lightweight `numpy.generic` interface that's actually generic, and
    doesn't require all that nasty `numpy.typing.NBitBase` stuff.
    """

    # unfortunately `| int` is required for compat with `numpy.__init__.pyi`
    @property
    def itemsize(self, /) -> _NB_co | int: ...

    @property
    def base(self, /) -> None: ...
    @property
    def data(self, /) -> memoryview: ...
    @property
    def dtype(self, /) -> np.dtype[Self]: ...  # type: ignore[type-var]  # pyright: ignore[reportInvalidTypeArguments]  # ty: ignore[invalid-type-arguments]
    @property
    def flags(self, /) -> Any: ...
    @property
    def nbytes(self, /) -> int: ...
    @property
    def ndim(self, /) -> _L0: ...
    @property
    def shape(self, /) -> tuple[()]: ...
    @property
    def size(self, /) -> _L1: ...
    @property
    def strides(self, /) -> tuple[()]: ...

    def item(self, k: _L0 | tuple[()] | tuple[_L0] = ..., /) -> _PT_co: ...

    @property
    def __array_priority__(self, /) -> float: ...  # -1000000.0
    @property
    def __array_interface__(self, /) -> dict[str, Any]: ...
    @property
    def __array_struct__(self, /) -> CapsuleType: ...

    def __bool__(self, /) -> bool: ...

    # Unlike `numpy/__init__.pyi` suggests, there exists no `__bytes__` method
    # in `np.generic`. Instead, it implements the (C) buffer protocol.

    def __buffer__(self, flags: int, /) -> memoryview: ...

    def __copy__(self, /) -> Self: ...
    def __deepcopy__(self, memo: dict[int, object] | None, /) -> Self: ...

    @overload
    def __array__(self, /) -> _Array0D[np.dtype[Any]]: ...
    @overload
    def __array__(self, dtype: _DT, /) -> _Array0D[_DT]: ...


generic = np.generic
flexible = np.flexible
character = np.character

if TYPE_CHECKING:
    from numpy._typing import (
        _8Bit,
        _16Bit,
        _32Bit,
        _64Bit,
        _96Bit,
        _128Bit,
    )

    type number = np.number[Any]  # noqa: PYI042
    type integer = np.integer[Any]  # noqa: PYI042
    type uinteger = np.unsignedinteger[Any]  # noqa: PYI042
    type sinteger = np.signedinteger[Any]  # noqa: PYI042
    type inexact = np.inexact[Any]  # noqa: PYI042
    type floating = np.floating[Any]  # noqa: PYI042
    type cfloating = np.complexfloating[Any, Any]  # noqa: PYI042

    type integer8 = np.integer[_8Bit]  # noqa: PYI042
    type integer16 = np.integer[_16Bit]  # noqa: PYI042
    type integer32 = np.integer[_32Bit]  # noqa: PYI042
    type integer64 = np.integer[_64Bit]  # noqa: PYI042

    type floating16 = np.floating[_16Bit]  # noqa: PYI042
    type floating32 = np.floating[_32Bit]  # noqa: PYI042
    type floating64 = np.floating[_64Bit]  # noqa: PYI042
    # float96, float128, and longdouble
    type floating80 = np.floating[_96Bit] | np.floating[_128Bit]  # noqa: PYI042

    type cfloating32 = np.complexfloating[_32Bit, _32Bit]  # noqa: PYI042
    type cfloating64 = np.complexfloating[_64Bit, _64Bit]  # noqa: PYI042
    # complex192, complex256, and clongdouble
    type cfloating80 = (  # noqa: PYI042
        np.complexfloating[_96Bit, _96Bit] | np.complexfloating[_128Bit, _128Bit]
    )

    type inexact32 = np.inexact[_32Bit]  # noqa: PYI042
    type inexact64 = np.inexact[_64Bit]  # noqa: PYI042
    # float96, complex192, float128, complex256, longdouble, and clongdouble
    type inexact80 = np.inexact[_96Bit] | np.inexact[_128Bit]  # noqa: PYI042

    type number8 = np.number[_8Bit]  # noqa: PYI042
    type number16 = np.number[_16Bit]  # noqa: PYI042
    type number32 = np.number[_32Bit]  # noqa: PYI042
    type number64 = np.number[_64Bit]  # noqa: PYI042
else:
    number = np.number
    integer = np.integer
    uinteger = np.unsignedinteger
    sinteger = np.signedinteger
    inexact = np.inexact
    floating = np.floating
    cfloating = np.complexfloating

    integer8 = np.int8 | np.uint8
    integer16 = np.int16 | np.uint16
    integer32 = np.int32 | np.uint32
    integer64 = np.int64 | np.uint64

    floating16 = np.float16
    floating32 = np.float32
    floating64 = np.float64
    floating80 = np.longdouble

    cfloating32 = np.complex64
    cfloating64 = np.complex128
    cfloating80 = np.clongdouble

    inexact32 = floating32 | cfloating32
    inexact64 = floating64 | cfloating64
    inexact80 = floating80 | cfloating80

    number8 = integer8
    number16 = integer16 | np.float16
    number32 = integer32 | inexact32
    number64 = integer64 | inexact64
