"""
The allowed `np.dtype` arguments for specific scalar types.
The names are analogous to those in `numpy.dtypes`.
"""

from collections.abc import Sequence
from typing import Any, Protocol

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

type fc_cls = type[JustFloat | JustComplex]  # noqa: PYI042
type ifc_cls = type[JustInt | JustFloat | JustComplex]  # noqa: PYI042

type O_cls = type[Just[object]]
type S_cls = type[Just[bytes]]
type U_cls = type[Just[str]]
type V_cls = type[memoryview]  # final
type SU_cls = S_cls | U_cls
type SUV_cls = SU_cls | V_cls

###

# b1
type AnyBoolDType = type[bool] | To[np.bool] | a.b1_code

# i1
type AnyInt8DType = To[np.int8] | a.i1_code
type AnyByteDType = AnyInt8DType  # deprecated
# u1
type AnyUInt8DType = To[np.uint8] | a.u1_code
type AnyUByteDType = AnyUInt8DType  # deprecated
# i2
type AnyInt16DType = To[np.int16] | a.i2_code
type AnyShortDType = AnyInt16DType  # deprecated
# u2
type AnyUInt16DType = To[np.uint16] | a.u2_code
type AnyUShortDType = AnyUInt16DType  # deprecated
# i4
type AnyInt32DType = To[np.int32] | a.i4_code
type AnyIntCDType = AnyInt32DType  # deprecated
# u4
type AnyUInt32DType = To[np.uint32] | a.u4_code
type AnyUIntCDType = AnyUInt32DType  # deprecated
# i8
type AnyInt64DType = To[np.int64] | a.i8_code
type AnyLongLongDType = AnyInt64DType  # deprecated
# u8
type AnyUInt64DType = To[np.uint64] | a.u8_code
type AnyULongLongDType = AnyUInt64DType  # deprecated
# int_ / intp / long
type AnyIntPDType = type[JustInt] | To[np.int64] | a.i0_code
type AnyIntDType = AnyIntPDType
type AnyLongDType = To[np.long] | a.l_code

# uint / uintp / ulong
type AnyUIntPDType = To[np.uint64] | a.u0_code
type AnyULongDType = To[np.ulong] | a.L_code
type AnyUIntDType = AnyULongDType


# f2
type AnyFloat16DType = To[_sc.floating16] | a.f2_code
# f4
type AnyFloat32DType = To[_sc.floating32] | a.f4_code
# f8
type AnyFloat64DType = type[JustFloat] | To[_sc.floating64] | a.f8_code
# f12 | f16
type AnyLongDoubleDType = To[_sc.floating80] | a.g_code

# c8
type AnyComplex64DType = To[_sc.cfloating32] | a.c8_code
# c16
type AnyComplex128DType = type[JustComplex] | To[_sc.cfloating64] | a.c16_code
# c24 | c32
type AnyCLongDoubleDType = To[_sc.cfloating80] | a.G_code

# M / N8
type AnyDateTime64DType = To[np.datetime64] | a.M_code
# m / m8
type AnyTimeDelta64DType = To[np.timedelta64] | a.m_code

# flexible

type AnyBytesDType = S_cls | To[np.bytes_] | a.S0_code
type AnyBytes8DType = a.S1_code
type AnyStrDType = S_cls | To[np.str_] | a.U0_code

type _ToShape = tuple[int, ...] | int
type _ToName = str | tuple[str, str]
type _ToDType = To[Any] | str
type _ToStructured = (
    tuple[_ToDType, _ToShape]
    | Sequence[tuple[_ToName, _ToDType] | tuple[_ToName, _ToDType, _ToShape]]
)
type AnyVoidDType = V_cls | To[np.void] | a.V0_code | _ToStructured

# object

type AnyObjectDType = O_cls | To[np.object_] | a.O_code
type AnyDType = _ToDType | _ToStructured

# abstract

type AnySignedIntegerDType = type[JustInt] | To[_sc.sinteger] | a.ix_code
type AnyUnsignedIntegerDType = To[_sc.uinteger] | a.ux_code
type AnyFloatingDType = type[JustFloat] | To[_sc.floating] | a.fx_code
type AnyComplexFloatingDType = type[JustComplex] | To[_sc.cfloating] | a.cx_code

type AnyIntegerDType = type[JustInt] | To[_sc.integer] | a.iu_code
type AnyInexactDType = fc_cls | To[_sc.inexact] | a.fc_code
type AnyNumberDType = ifc_cls | To[_sc.number] | a.iufc_code

type AnyCharacterDType = SU_cls | To[np.character] | a.SU_code
type AnyFlexibleDType = SUV_cls | To[np.flexible] | a.SUV_code


class _HasStringDType(Protocol):
    @property
    def dtype(self) -> numpy_typing_compat.StringDType: ...


type AnyStringDType = numpy_typing_compat.StringDType | _HasStringDType | a.T_code
