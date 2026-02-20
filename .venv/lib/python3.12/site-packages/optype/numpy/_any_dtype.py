"""
The allowed `np.dtype` arguments for specific scalar types.
The names are analogous to those in `numpy.dtypes`.
"""

import sys
from collections.abc import Sequence
from typing import Any, Never, Protocol, TypeAlias

if sys.version_info >= (3, 13):
    from typing import TypeAliasType
else:
    from typing_extensions import TypeAliasType

import numpy as np
import numpy_typing_compat

import optype.numpy._dtype_attr as a
import optype.numpy._scalar as _sc
from ._dtype import ToDType as To
from optype._core._just import Just, JustComplex, JustFloat, JustInt

# ruff: noqa: RUF022
__all__ = [
    "AnyDType",
    "AnyNumberDType",
    "AnyIntegerDType",
    "AnyInexactDType",
    "AnyFlexibleDType",
    "AnyUnsignedIntegerDType",
    "AnySignedIntegerDType",
    "AnyFloatingDType",
    "AnyComplexFloatingDType",
    "AnyCharacterDType",

    "AnyBoolDType",

    "AnyUIntDType",
    "AnyUInt8DType",
    "AnyUInt8DType",
    "AnyUInt16DType",
    "AnyUInt32DType",
    "AnyUInt64DType",
    "AnyUIntPDType",
    "AnyUByteDType",
    "AnyUShortDType",
    "AnyUIntCDType",
    "AnyULongDType",
    "AnyULongLongDType",

    "AnyIntDType",
    "AnyInt8DType",
    "AnyInt8DType",
    "AnyInt16DType",
    "AnyInt32DType",
    "AnyInt64DType",
    "AnyIntPDType",
    "AnyByteDType",
    "AnyShortDType",
    "AnyIntCDType",
    "AnyLongDType",
    "AnyLongLongDType",

    "AnyFloat16DType",
    "AnyFloat32DType",
    "AnyFloat64DType",
    "AnyLongDoubleDType",

    "AnyComplex64DType",
    "AnyComplex128DType",
    "AnyCLongDoubleDType",

    "AnyDateTime64DType",
    "AnyTimeDelta64DType",

    "AnyBytesDType",
    "AnyBytes8DType",
    "AnyStrDType",
    "AnyVoidDType",
    "AnyObjectDType",

    "AnyStringDType",
]  # fmt: skip


def __dir__() -> list[str]:
    return __all__


###

fc_cls: TypeAlias = type[JustFloat | JustComplex]  # noqa: PYI042
ifc_cls: TypeAlias = type[JustInt | JustFloat | JustComplex]  # noqa: PYI042

O_cls: TypeAlias = type[Just[object]]
S_cls: TypeAlias = type[Just[bytes]]
U_cls: TypeAlias = type[Just[str]]
V_cls: TypeAlias = type[memoryview]  # final
SU_cls: TypeAlias = S_cls | U_cls
SUV_cls: TypeAlias = SU_cls | V_cls

###

# mypy: disable-error-code="no-redef"
# pyright: reportRedeclaration=false

# b1
AnyBoolDType = TypeAliasType("AnyBoolDType", type[bool] | To[np.bool_] | a.b1_code)

# i1
AnyInt8DType = TypeAliasType("AnyInt8DType", To[np.int8] | a.i1_code)
AnyByteDType = AnyInt8DType  # deprecated
# u1
AnyUInt8DType = TypeAliasType("AnyUInt8DType", To[np.uint8] | a.u1_code)
AnyUByteDType = AnyUInt8DType  # deprecated
# i2
AnyInt16DType = TypeAliasType("AnyInt16DType", To[np.int16] | a.i2_code)
AnyShortDType = AnyInt16DType  # deprecated
# u2
AnyUInt16DType = TypeAliasType("AnyUInt16DType", To[np.uint16] | a.u2_code)
AnyUShortDType = AnyUInt16DType  # deprecated
# i4
AnyInt32DType = TypeAliasType("AnyInt32DType", To[np.int32] | a.i4_code)
AnyIntCDType = AnyInt32DType  # deprecated
# u4
AnyUInt32DType = TypeAliasType("AnyUInt32DType", To[np.uint32] | a.u4_code)
AnyUIntCDType = AnyUInt32DType  # deprecated
# i8
AnyInt64DType = TypeAliasType("AnyInt64DType", To[np.int64] | a.i8_code)
AnyLongLongDType = AnyInt64DType  # deprecated
# u8
AnyUInt64DType = TypeAliasType("AnyUInt64DType", To[np.uint64] | a.u8_code)
AnyULongLongDType = AnyUInt64DType  # deprecated
# int_ / intp / long
if numpy_typing_compat.NUMPY_GE_2_0:
    AnyIntPDType = TypeAliasType(
        "AnyIntPDType",
        type[JustInt] | To[np.int64] | a.i0_code,
    )
    AnyIntDType = TypeAliasType("AnyIntDType", AnyIntPDType)
    AnyLongDType = TypeAliasType(
        "AnyLongDType",
        To[numpy_typing_compat.long] | a.l_code,
    )
else:
    AnyIntPDType = TypeAliasType("AnyIntPDType", To[np.int64] | a.i0_code)  # type: ignore[misc]
    AnyIntDType = TypeAliasType(  # type: ignore[misc]
        "AnyIntDType",
        type[JustInt] | To[numpy_typing_compat.long] | a.l_code,
    )
    AnyLongDType = AnyIntDType  # type: ignore[misc]

# uint / uintp / ulong
AnyUIntPDType = TypeAliasType("AnyUIntPDType", To[np.uint64] | a.u0_code)
AnyULongDType = TypeAliasType("AnyULongDType", To[numpy_typing_compat.ulong] | a.L_code)
AnyUIntDType = AnyULongDType


# f2
AnyFloat16DType = TypeAliasType("AnyFloat16DType", To[_sc.floating16] | a.f2_code)
# f4
AnyFloat32DType = TypeAliasType("AnyFloat32DType", To[_sc.floating32] | a.f4_code)
# f8
AnyFloat64DType = TypeAliasType(
    "AnyFloat64DType",
    type[JustFloat] | To[_sc.floating64] | a.f8_code,
)
# f12 | f16
AnyLongDoubleDType = TypeAliasType("AnyLongDoubleDType", To[_sc.floating80] | a.g_code)

# c8
AnyComplex64DType = TypeAliasType("AnyComplex64DType", To[_sc.cfloating32] | a.c8_code)
# c16
AnyComplex128DType = TypeAliasType(
    "AnyComplex128DType",
    type[JustComplex] | To[_sc.cfloating64] | a.c16_code,
)
# c24 | c32
AnyCLongDoubleDType = TypeAliasType(
    "AnyCLongDoubleDType",
    To[_sc.cfloating80] | a.G_code,
)

# M / N8
AnyDateTime64DType = TypeAliasType("AnyDateTime64DType", To[np.datetime64] | a.M_code)
# m / m8
AnyTimeDelta64DType = TypeAliasType(
    "AnyTimeDelta64DType",
    To[np.timedelta64] | a.m_code,
)

# flexible

AnyBytesDType = TypeAliasType("AnyBytesDType", S_cls | To[np.bytes_] | a.S0_code)
AnyBytes8DType = TypeAliasType("AnyBytes8DType", a.S1_code)
AnyStrDType = TypeAliasType("AnyStrDType", S_cls | To[np.str_] | a.U0_code)

_ToShape: TypeAlias = tuple[int, ...] | int
_ToName: TypeAlias = str | tuple[str, str]
_ToDType: TypeAlias = To[Any] | str
_ToStructured: TypeAlias = (
    tuple[_ToDType, _ToShape]
    | Sequence[tuple[_ToName, _ToDType] | tuple[_ToName, _ToDType, _ToShape]]
)
AnyVoidDType = TypeAliasType(
    "AnyVoidDType",
    V_cls | To[np.void] | a.V0_code | _ToStructured,
)

# object

AnyObjectDType = TypeAliasType("AnyObjectDType", O_cls | To[np.object_] | a.O_code)
AnyDType = TypeAliasType("AnyDType", _ToDType | _ToStructured)

# abstract

AnySignedIntegerDType = TypeAliasType(
    "AnySignedIntegerDType",
    type[JustInt] | To[_sc.sinteger] | a.ix_code,
)
AnyUnsignedIntegerDType = TypeAliasType(
    "AnyUnsignedIntegerDType",
    To[_sc.uinteger] | a.ux_code,
)
AnyFloatingDType = TypeAliasType(
    "AnyFloatingDType",
    type[JustFloat] | To[_sc.floating] | a.fx_code,
)
AnyComplexFloatingDType = TypeAliasType(
    "AnyComplexFloatingDType",
    type[JustComplex] | To[_sc.cfloating] | a.cx_code,
)

AnyIntegerDType = TypeAliasType(
    "AnyIntegerDType",
    type[JustInt] | To[_sc.integer] | a.iu_code,
)
AnyInexactDType = TypeAliasType("AnyInexactDType", fc_cls | To[_sc.inexact] | a.fc_code)
AnyNumberDType = TypeAliasType("AnyNumberDType", ifc_cls | To[_sc.number] | a.iufc_code)

AnyCharacterDType = TypeAliasType(
    "AnyCharacterDType",
    SU_cls | To[np.character] | a.SU_code,
)
AnyFlexibleDType = TypeAliasType(
    "AnyFlexibleDType",
    SUV_cls | To[np.flexible] | a.SUV_code,
)

if numpy_typing_compat.NUMPY_GE_2_0:

    class _HasStringDType(Protocol):
        @property
        def dtype(self) -> numpy_typing_compat.StringDType: ...

    AnyStringDType = TypeAliasType(
        "AnyStringDType",
        numpy_typing_compat.StringDType | _HasStringDType | a.T_code,
    )
else:
    AnyStringDType = TypeAliasType("AnyStringDType", Never)  # type: ignore[misc]
