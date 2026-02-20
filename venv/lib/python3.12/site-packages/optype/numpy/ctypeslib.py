"""
A collection of `ctypes` type aliases for several numpy scalar types.

NOTE:
    `optype.numpy` assumes a `C99`-compatible C compiler, a 32 or 64-bit
    system, and a `ILP32`, `LLP64` or `LP64` data model; see
    https://en.cppreference.com/w/c/language/arithmetic_types for more info.

    If this is not the case for you, then please open an issue at:
    https://github.com/jorenham/optype/issues
"""

import ctypes as ct
import sys

# ruff: noqa: N812
from ctypes import (
    c_bool as Bool,
    c_byte as Byte,
    c_char as Bytes,
    c_double as Float64,
    c_float as Float32,
    c_int as IntC,
    c_int8 as Int8,
    c_int16 as Int16,
    c_int32 as Int32,
    c_int64 as Int64,
    c_long as Long,
    # NOTE: `longdouble` only works as type, not as value!
    c_longdouble as LongDouble,
    c_longlong as LongLong,
    c_short as Short,
    c_size_t as UIntP,  # `void_p` on numpy<2, but almost always the same
    c_ssize_t as IntP,
    c_ubyte as UByte,
    c_uint as UIntC,
    c_uint8 as UInt8,
    c_uint16 as UInt16,
    c_uint32 as UInt32,
    c_uint64 as UInt64,
    c_ulong as ULong,
    c_ulonglong as ULongLong,
    c_ushort as UShort,
    py_object as Object,
)
from typing import TYPE_CHECKING, Final, Literal, Never, TypeAlias, cast

if sys.version_info >= (3, 13):
    from typing import TypeAliasType, TypeVar
else:
    from typing_extensions import TypeAliasType, TypeVar

if sys.version_info >= (3, 14) and sys.platform != "win32":
    from ctypes import (  # noqa: I001
        c_float_complex as Complex64,
        c_double_complex as Complex128,
        c_longdouble_complex as CLongDouble,
    )
else:
    Complex64 = Never
    Complex128 = Never
    CLongDouble = Never

from ._ctypeslib import CScalar, CType

# ruff: noqa: RUF022
__all__ = [
    "CType",
    "CScalar",
    "Array",

    "Generic",
    "Number",
    "Integer",
    "Inexact",
    "UnsignedInteger",
    "SignedInteger",
    "Floating",
    "ComplexFloating",
    "Flexible",

    "Bool",

    "UInt8", "Int8",
    "UInt16", "Int16",
    "UInt32", "Int32",
    "UInt64", "Int64",
    "UByte", "Byte",
    "UShort", "Short",
    "UIntC", "IntC",
    "UIntP", "IntP",
    "ULong", "Long",
    "ULongLong", "LongLong",

    "Float32",
    "Float64",
    "LongDouble",

    "Complex64",
    "Complex128",
    "CLongDouble",

    "Bytes",
    "Void",
    "Object",
]  # fmt: skip


def __dir__() -> list[str]:
    return __all__


###


SIZE_BYTE: Final = cast("Literal[1]", ct.sizeof(ct.c_byte))
SIZE_SHORT: Final = cast("Literal[2]", ct.sizeof(ct.c_short))
SIZE_INTC: Final = cast("Literal[4]", ct.sizeof(ct.c_int))
SIZE_INTP: Final = cast("Literal[4, 8]", ct.sizeof(ct.c_ssize_t))
SIZE_LONG: Final = cast("Literal[4, 8]", ct.sizeof(ct.c_long))
SIZE_LONGLONG: Final = cast("Literal[8]", ct.sizeof(ct.c_longlong))

SIZE_SINGLE: Final = cast("Literal[4]", ct.sizeof(ct.c_float))
SIZE_DOUBLE: Final = cast("Literal[8]", ct.sizeof(ct.c_double))
SIZE_LONGDOUBLE: Final = cast("Literal[8, 10, 12, 16]", ct.sizeof(ct.c_longdouble))


def __is_dev() -> bool:
    from importlib import metadata  # noqa: PLC0415

    return "dev" in metadata.version((__package__ or "optype").removesuffix(".numpy"))


if __is_dev():
    assert SIZE_BYTE == 1, f"`sizeof(byte) = {SIZE_BYTE}`, expected 1"
    assert SIZE_SHORT == 2, f"`sizeof(short) = {SIZE_SHORT}`, expected 2"
    # If you run a 16-bit system and this assertion fails, please open an issue, or even
    # better, upgrade to something that's younger than Joe fucking Biden.
    assert SIZE_INTC == 4, f"`sizeof(int) = {SIZE_INTC}`, expected 4"
    assert SIZE_INTP in {4, 8}, f"`sizeof(ssize_t) = {SIZE_INTP}`, expected 4 or 8"
    assert SIZE_LONG in {4, 8}, f"`sizeof(long int) = {SIZE_LONG}`, expected 4 or 8"
    assert SIZE_LONGLONG == 8, f"`sizeof(long long int) = {SIZE_LONGLONG}`, expected 8"

    assert SIZE_SINGLE == 4, f"`sizeof(float) = {SIZE_SINGLE}`, expected 4"
    assert SIZE_DOUBLE == 8, f"`sizeof(double) = {SIZE_DOUBLE}`, expected 8"
    assert SIZE_LONGDOUBLE in {8, 10, 12, 16}, (
        f"`sizeof(long double) = {SIZE_LONGDOUBLE}`, expected 8, 10, 12 or 16",
    )  # fmt: skip


CT = TypeVar("CT", bound="CType")
_Array: TypeAlias = ct.Array[CT] | ct.Array["_Array[CT]"]
Array = TypeAliasType("Array", _Array[CT], type_params=(CT,))

# `c_(u)byte` is an alias for `c_(u)int8`
_SignedInteger: TypeAlias = (
    Int8 | Int16 | Int32 | Int64 | Short | IntC | IntP | Long | LongLong
)
_UnsignedInteger: TypeAlias = (
    UInt8 | UInt16 | UInt32 | UInt64 | UShort | UIntC | UIntP | ULong | ULongLong
)

if TYPE_CHECKING:
    # exploit the fact that `ct._SimpleCData` is invariant in `T`
    _Integer: TypeAlias = CScalar[int]
    _Floating: TypeAlias = CScalar[float]
    _ComplexFloating: TypeAlias = CScalar[complex]
else:
    _Integer = CScalar
    _Floating = CScalar
    _ComplexFloating = CScalar

SignedInteger = TypeAliasType("SignedInteger", _SignedInteger)
UnsignedInteger = TypeAliasType("UnsignedInteger", _UnsignedInteger)
Integer = TypeAliasType("Integer", _Integer)
Floating = TypeAliasType("Floating", _Floating)
ComplexFloating = TypeAliasType("ComplexFloating", _ComplexFloating)
Inexact = TypeAliasType("Inexact", _Floating | _ComplexFloating)
Number = TypeAliasType("Number", _Integer | _Floating | _ComplexFloating)

Void = TypeAliasType("Void", ct.Structure | ct.Union)
Flexible = TypeAliasType("Flexible", Bytes | Void)

Generic = TypeAliasType("Generic", Bool | Number | Flexible | Object)
