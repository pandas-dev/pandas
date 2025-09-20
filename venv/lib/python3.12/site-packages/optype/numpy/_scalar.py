import sys
from typing import Any, Literal, Protocol, Self, TypeAlias, overload

if sys.version_info >= (3, 13):
    from types import CapsuleType
    from typing import TypeAliasType, TypeVar, runtime_checkable
else:
    from typing_extensions import CapsuleType, TypeAliasType, TypeVar, runtime_checkable

import numpy as np
from numpy._typing import (
    _8Bit,  # noqa: PLC2701
    _16Bit,  # noqa: PLC2701
    _32Bit,  # noqa: PLC2701
    _64Bit,  # noqa: PLC2701
    _96Bit,  # noqa: PLC2701
    _128Bit,  # noqa: PLC2701
)

from optype._utils import set_module

__all__ = ["Scalar"]


###

_PT_co = TypeVar("_PT_co", covariant=True)
_NB_co = TypeVar("_NB_co", bound=int, default=int, covariant=True)
_DT = TypeVar("_DT", bound=np.dtype[Any], default=np.dtype[Any])

_L0: TypeAlias = Literal[0]
_L1: TypeAlias = Literal[1]
_Array0D: TypeAlias = np.ndarray[tuple[()], _DT]


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
    def dtype(self, /) -> np.dtype[Self]: ...  # type: ignore[type-var]  # pyright: ignore[reportInvalidTypeArguments]
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

    if sys.version_info >= (3, 12):

        def __buffer__(self, flags: int, /) -> memoryview: ...

    def __copy__(self, /) -> Self: ...
    def __deepcopy__(self, memo: dict[int, object] | None, /) -> Self: ...

    @overload
    def __array__(self, /) -> _Array0D: ...
    @overload
    def __array__(self, dtype: _DT, /) -> _Array0D[_DT]: ...


generic = np.generic
flexible = np.flexible
character = np.character

number = TypeAliasType("number", np.number[Any])
integer = TypeAliasType("integer", np.integer[Any])
uinteger = TypeAliasType("uinteger", np.unsignedinteger[Any])
sinteger = TypeAliasType("sinteger", np.signedinteger[Any])
inexact = TypeAliasType("inexact", np.inexact[Any])
floating = TypeAliasType("floating", np.floating[Any])
cfloating = TypeAliasType("cfloating", np.complexfloating[Any, Any])

integer8 = TypeAliasType("integer8", np.integer[_8Bit])
integer16 = TypeAliasType("integer16", np.integer[_16Bit])
integer32 = TypeAliasType("integer32", np.integer[_32Bit])
integer64 = TypeAliasType("integer64", np.integer[_64Bit])

floating16 = TypeAliasType("floating16", np.floating[_16Bit])
floating32 = TypeAliasType("floating32", np.floating[_32Bit])
floating64 = TypeAliasType("floating64", np.floating[_64Bit])
# float96, float128, and longdouble
floating80 = TypeAliasType("floating80", np.floating[_96Bit] | np.floating[_128Bit])

cfloating32 = TypeAliasType("cfloating32", np.complexfloating[_32Bit, _32Bit])
cfloating64 = TypeAliasType("cfloating64", np.complexfloating[_64Bit, _64Bit])
# complex192, complex256, and clongdouble
cfloating80 = TypeAliasType(
    "cfloating80",
    np.complexfloating[_96Bit, _96Bit] | np.complexfloating[_128Bit, _128Bit],
)

inexact32 = TypeAliasType("inexact32", np.inexact[_32Bit])
inexact64 = TypeAliasType("inexact64", np.inexact[_64Bit])
# float96, complex192, float128, complex256, longdouble, and clongdouble
inexact80 = TypeAliasType("inexact80", np.inexact[_96Bit] | np.inexact[_128Bit])

number8 = TypeAliasType("number8", np.number[_8Bit])
number16 = TypeAliasType("number16", np.number[_16Bit])
number32 = TypeAliasType("number32", np.number[_32Bit])
number64 = TypeAliasType("number64", np.number[_64Bit])
