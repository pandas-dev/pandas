import sys
from collections.abc import Iterator
from typing import Any, Protocol, TypeAlias

from optype._utils import set_module

if sys.version_info >= (3, 13):
    from typing import TypeAliasType, TypeVar
else:
    from typing_extensions import TypeAliasType, TypeVar

import numpy as np
import numpy_typing_compat as nptc

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


_T = TypeVar("_T")
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


_AnyArrayPY: TypeAlias = tuple[_T, ...] | _AnyArrayPY0[_T]
_AnyArray = TypeAliasType(
    "_AnyArray",
    _AnyArrayNP[_ST] | _AnyArrayPY[_VT] | _AnyArrayPY[_AnyArrayNP[_ST]],
    type_params=(_ST, _VT),
)

###

AnyArray: TypeAlias = _AnyArray[_ST, object] | CanBuffer

AnyNumberArray: TypeAlias = _AnyArray[
    _sc.number,
    _sc.number | JustInt | JustFloat | JustComplex,
]
AnyIntegerArray: TypeAlias = _AnyArray[_sc.integer, _sc.integer | JustInt]
AnySignedIntegerArray: TypeAlias = _AnyArray[_sc.sinteger, _sc.sinteger | JustInt]
AnyUnsignedIntegerArray: TypeAlias = _AnyArray[_sc.uinteger]
AnyInexactArray: TypeAlias = _AnyArray[
    _sc.inexact,
    _sc.inexact | JustFloat | JustComplex,
]

AnyBoolArray: TypeAlias = _AnyArray[np.bool_, np.bool_ | bool]

AnyUInt8Array: TypeAlias = _AnyArray[np.uint8, np.uint8 | CanBuffer] | CanBuffer
AnyUByteArray = AnyUInt8Array
AnyUInt16Array: TypeAlias = _AnyArray[np.uint16]
AnyUShortArray = AnyUInt16Array
AnyUInt32Array: TypeAlias = _AnyArray[np.uint32]
AnyUInt64Array: TypeAlias = _AnyArray[np.uint64]
AnyUIntCArray: TypeAlias = _AnyArray[np.uintc]
AnyULongLongArray: TypeAlias = _AnyArray[np.ulonglong]
AnyULongArray: TypeAlias = _AnyArray[nptc.ulong]
AnyUIntPArray: TypeAlias = _AnyArray[np.uintp]
AnyUIntArray: TypeAlias = _AnyArray[np.uint]

AnyInt8Array: TypeAlias = _AnyArray[np.int8]
AnyByteArray = AnyInt8Array
AnyInt16Array: TypeAlias = _AnyArray[np.int16]
AnyShortArray = AnyInt16Array
AnyInt32Array: TypeAlias = _AnyArray[np.int32]
AnyInt64Array: TypeAlias = _AnyArray[np.int64]
AnyIntCArray: TypeAlias = _AnyArray[np.intc]
AnyLongLongArray: TypeAlias = _AnyArray[np.longlong]
AnyLongArray: TypeAlias = _AnyArray[nptc.long]  # no int (numpy<=1)
AnyIntPArray: TypeAlias = _AnyArray[np.intp]  # no int (numpy>=2)
AnyIntArray: TypeAlias = _AnyArray[np.int_, np.int_ | JustInt]

AnyFloatingArray: TypeAlias = _AnyArray[_sc.floating, _sc.floating | JustFloat]
AnyFloat16Array: TypeAlias = _AnyArray[np.float16]
AnyFloat32Array: TypeAlias = _AnyArray[np.float32]
AnyFloat64Array: TypeAlias = _AnyArray[_sc.floating64, _sc.floating64 | JustFloat]
AnyLongDoubleArray: TypeAlias = _AnyArray[_sc.floating80]

AnyComplexFloatingArray: TypeAlias = _AnyArray[
    _sc.cfloating,
    _sc.cfloating | JustComplex,
]
AnyComplex64Array: TypeAlias = _AnyArray[np.complex64]
AnyComplex128Array: TypeAlias = _AnyArray[
    _sc.cfloating64,
    _sc.cfloating64 | JustComplex,
]
AnyCLongDoubleArray: TypeAlias = _AnyArray[_sc.cfloating80]

AnyCharacterArray: TypeAlias = _AnyArray[np.character, np.character | bytes | str]
AnyBytesArray: TypeAlias = _AnyArray[np.bytes_, np.bytes_ | bytes]
AnyStrArray: TypeAlias = _AnyArray[np.str_, np.str_ | str]

AnyFlexibleArray: TypeAlias = _AnyArray[np.flexible, np.flexible | bytes | str]

# TODO(jorenham): structured types
# https://github.com/jorenham/optype/issues/371
AnyVoidArray: TypeAlias = _AnyArray[np.void]

AnyDateTime64Array: TypeAlias = _AnyArray[np.datetime64]
AnyTimeDelta64Array: TypeAlias = _AnyArray[np.timedelta64]
AnyObjectArray: TypeAlias = _AnyArray[np.object_, np.object_ | JustObject]


@set_module("optype.numpy")
class AnyStringArray(Protocol):
    def __array__(self, /) -> np.ndarray[Any, "nptc.StringDType"]: ...
