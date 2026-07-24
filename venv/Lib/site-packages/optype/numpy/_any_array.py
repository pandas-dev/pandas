import sys
from collections.abc import Iterator
from typing import Any, Protocol, TypeAliasType

from optype._utils import set_module

if sys.version_info >= (3, 13):
    from typing import TypeVar
else:
    from typing_extensions import TypeVar

import numpy as np
import numpy_typing_compat

import optype.numpy._scalar as _sc
from ._shape import AnyShape
from optype._core import CanBuffer, JustComplex, JustFloat, JustInt, JustObject

# ruff: noqa: RUF022
__all__ = [
    "AnyArray",
    "AnyNumberArray",
    "AnyIntegerArray",
    "AnyUnsignedIntegerArray",
    "AnySignedIntegerArray",
    "AnyInexactArray",
    "AnyFloatingArray",
    "AnyComplexFloatingArray",
    "AnyFlexibleArray",
    "AnyCharacterArray",

    "AnyBoolArray",

    "AnyUIntArray", "AnyIntArray",
    "AnyUInt8Array", "AnyInt8Array",
    "AnyUInt16Array", "AnyInt16Array",
    "AnyUInt32Array", "AnyInt32Array",
    "AnyUInt64Array", "AnyInt64Array",
    "AnyUByteArray", "AnyByteArray",
    "AnyUShortArray", "AnyShortArray",
    "AnyUIntCArray", "AnyIntCArray",
    "AnyUIntPArray", "AnyIntPArray",
    "AnyULongArray", "AnyLongArray",
    "AnyULongLongArray", "AnyLongLongArray",

    "AnyFloat16Array",
    "AnyFloat32Array", "AnyComplex64Array",
    "AnyFloat64Array", "AnyComplex128Array",
    "AnyLongDoubleArray", "AnyCLongDoubleArray",

    "AnyDateTime64Array",
    "AnyTimeDelta64Array",

    "AnyBytesArray",
    "AnyStrArray",
    "AnyVoidArray",
    "AnyObjectArray",
    "AnyStringArray",
]  # fmt: skip


def __dir__() -> list[str]:
    return __all__


###


_ST = TypeVar("_ST", bound=np.generic, default=Any)
_VT = TypeVar("_VT", default=_ST)

_T_co = TypeVar("_T_co", covariant=True)
_ST_co = TypeVar("_ST_co", bound=np.generic, covariant=True)


# NOTE: Does not include scalar types
class _AnyArrayNP(Protocol[_ST_co]):
    def __len__(self, /) -> int: ...
    def __array__(self, /) -> np.ndarray[AnyShape, np.dtype[_ST_co]]: ...


# NOTE: does not include tuple
class _AnyArrayPY0(Protocol[_T_co]):
    def __len__(self, /) -> int: ...
    def __getitem__(self, i: int, /) -> "_T_co | _AnyArrayPY0[_T_co]": ...
    def __reversed__(self, /) -> Iterator["_T_co | _AnyArrayPY0[_T_co]"]: ...
    def index(self, x: Any, /) -> int: ...


type _AnyArrayPY[T] = tuple[T, ...] | _AnyArrayPY0[T]
_AnyArray = TypeAliasType(  # noqa: UP040
    "_AnyArray",
    _AnyArrayNP[_ST] | _AnyArrayPY[_VT] | _AnyArrayPY[_AnyArrayNP[_ST]],
    type_params=(_ST, _VT),
)

###

AnyArray = TypeAliasType(  # noqa: UP040
    "AnyArray",
    _AnyArray[_ST, object] | CanBuffer,
    type_params=(_ST,),
)

type AnyNumberArray = _AnyArray[
    _sc.number,
    _sc.number | JustInt | JustFloat | JustComplex,
]
type AnyIntegerArray = _AnyArray[_sc.integer, _sc.integer | JustInt]
type AnySignedIntegerArray = _AnyArray[_sc.sinteger, _sc.sinteger | JustInt]
type AnyUnsignedIntegerArray = _AnyArray[_sc.uinteger]
type AnyInexactArray = _AnyArray[
    _sc.inexact,
    _sc.inexact | JustFloat | JustComplex,
]

type AnyBoolArray = _AnyArray[np.bool, np.bool | bool]

type AnyUInt8Array = _AnyArray[np.uint8, np.uint8 | CanBuffer] | CanBuffer
AnyUByteArray = AnyUInt8Array
type AnyUInt16Array = _AnyArray[np.uint16]
AnyUShortArray = AnyUInt16Array
type AnyUInt32Array = _AnyArray[np.uint32]
type AnyUInt64Array = _AnyArray[np.uint64]
type AnyUIntCArray = _AnyArray[np.uintc]
type AnyULongLongArray = _AnyArray[np.ulonglong]
type AnyULongArray = _AnyArray[np.ulong]
type AnyUIntPArray = _AnyArray[np.uintp]
type AnyUIntArray = _AnyArray[np.uint]

type AnyInt8Array = _AnyArray[np.int8]
AnyByteArray = AnyInt8Array
type AnyInt16Array = _AnyArray[np.int16]
AnyShortArray = AnyInt16Array
type AnyInt32Array = _AnyArray[np.int32]
type AnyInt64Array = _AnyArray[np.int64]
type AnyIntCArray = _AnyArray[np.intc]
type AnyLongLongArray = _AnyArray[np.longlong]
type AnyLongArray = _AnyArray[np.long]  # no int (numpy<=1)
type AnyIntPArray = _AnyArray[np.intp]  # no int (numpy>=2)
type AnyIntArray = _AnyArray[np.int_, np.int_ | JustInt]

type AnyFloatingArray = _AnyArray[_sc.floating, _sc.floating | JustFloat]
type AnyFloat16Array = _AnyArray[np.float16]
type AnyFloat32Array = _AnyArray[np.float32]
type AnyFloat64Array = _AnyArray[_sc.floating64, _sc.floating64 | JustFloat]
type AnyLongDoubleArray = _AnyArray[_sc.floating80]

type AnyComplexFloatingArray = _AnyArray[
    _sc.cfloating,
    _sc.cfloating | JustComplex,
]
type AnyComplex64Array = _AnyArray[np.complex64]
type AnyComplex128Array = _AnyArray[
    _sc.cfloating64,
    _sc.cfloating64 | JustComplex,
]
type AnyCLongDoubleArray = _AnyArray[_sc.cfloating80]

type AnyCharacterArray = _AnyArray[np.character, np.character | bytes | str]
type AnyBytesArray = _AnyArray[np.bytes_, np.bytes_ | bytes]
type AnyStrArray = _AnyArray[np.str_, np.str_ | str]

type AnyFlexibleArray = _AnyArray[np.flexible, np.flexible | bytes | str]

# TODO(jorenham): structured types
# https://github.com/jorenham/optype/issues/371
type AnyVoidArray = _AnyArray[np.void]

type AnyDateTime64Array = _AnyArray[np.datetime64]
type AnyTimeDelta64Array = _AnyArray[np.timedelta64]
type AnyObjectArray = _AnyArray[np.object_, np.object_ | JustObject]


@set_module("optype.numpy")
class AnyStringArray(Protocol):
    def __array__(self, /) -> np.ndarray[Any, numpy_typing_compat.StringDType]: ...
