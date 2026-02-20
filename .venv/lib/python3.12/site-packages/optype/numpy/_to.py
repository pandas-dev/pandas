import sys
from collections.abc import Sequence as Seq
from typing import Any, Literal, Protocol, TypeAlias

if sys.version_info >= (3, 13):
    from typing import TypeAliasType, TypeVar
else:
    from typing_extensions import TypeAliasType, TypeVar

import numpy as np
import numpy_typing_compat as nptc

import optype.numpy.compat as npc
from ._sequence_nd import SequenceND as SeqND
from optype._core._just import JustComplex, JustFloat, JustInt

__all__ = [  # noqa: RUF022
    "ToScalar",
    "ToArray1D", "ToArrayStrict1D",
    "ToArray2D", "ToArrayStrict2D",
    "ToArray3D", "ToArrayStrict3D",
    "ToArrayND",

    "ToFalse", "ToTrue",
    "ToJustFalse", "ToJustTrue",

    "ToBool", "ToJustBool",
    "ToBool1D", "ToJustBool1D", "ToBoolStrict1D", "ToJustBoolStrict1D",
    "ToBool2D", "ToJustBool2D", "ToBoolStrict2D", "ToJustBoolStrict2D",
    "ToBool3D", "ToJustBool3D", "ToBoolStrict3D", "ToJustBoolStrict3D",
    "ToBoolND", "ToJustBoolND",

    "ToInt", "ToJustInt",
    "ToInt1D", "ToJustInt1D", "ToIntStrict1D", "ToJustIntStrict1D",
    "ToInt2D", "ToJustInt2D", "ToIntStrict2D", "ToJustIntStrict2D",
    "ToInt3D", "ToJustInt3D", "ToIntStrict3D", "ToJustIntStrict3D",
    "ToIntND", "ToJustIntND",

    "ToJustInt64",
    "ToJustInt64_1D", "ToJustInt64Strict1D",
    "ToJustInt64_2D", "ToJustInt64Strict2D",
    "ToJustInt64_3D", "ToJustInt64Strict3D",
    "ToJustInt64_ND",

    "ToFloat16", "ToJustFloat16",
    "ToFloat16_1D", "ToJustFloat16_1D", "ToFloat16Strict1D", "ToJustFloat16Strict1D",
    "ToFloat16_2D", "ToJustFloat16_2D", "ToFloat16Strict2D", "ToJustFloat16Strict2D",
    "ToFloat16_3D", "ToJustFloat16_3D", "ToFloat16Strict3D", "ToJustFloat16Strict3D",
    "ToFloat16_ND", "ToJustFloat16_ND",

    "ToFloat32", "ToJustFloat32",
    "ToFloat32_1D", "ToJustFloat32_1D", "ToFloat32Strict1D", "ToJustFloat32Strict1D",
    "ToFloat32_2D", "ToJustFloat32_2D", "ToFloat32Strict2D", "ToJustFloat32Strict2D",
    "ToFloat32_3D", "ToJustFloat32_3D", "ToFloat32Strict3D", "ToJustFloat32Strict3D",
    "ToFloat32_ND", "ToJustFloat32_ND",

    "ToFloat64", "ToJustFloat64",
    "ToFloat64_1D", "ToJustFloat64_1D", "ToFloat64Strict1D", "ToJustFloat64Strict1D",
    "ToFloat64_2D", "ToJustFloat64_2D", "ToFloat64Strict2D", "ToJustFloat64Strict2D",
    "ToFloat64_3D", "ToJustFloat64_3D", "ToFloat64Strict3D", "ToJustFloat64Strict3D",
    "ToFloat64_ND", "ToJustFloat64_ND",

    "ToJustLongDouble",
    "ToJustLongDouble1D", "ToJustLongDoubleStrict1D",
    "ToJustLongDouble2D", "ToJustLongDoubleStrict2D",
    "ToJustLongDouble3D", "ToJustLongDoubleStrict3D",
    "ToJustLongDoubleND",

    "ToFloat", "ToJustFloat",
    "ToFloat1D", "ToJustFloat1D", "ToFloatStrict1D", "ToJustFloatStrict1D",
    "ToFloat2D", "ToJustFloat2D", "ToFloatStrict2D", "ToJustFloatStrict2D",
    "ToFloat3D", "ToJustFloat3D", "ToFloatStrict3D", "ToJustFloatStrict3D",
    "ToFloatND", "ToJustFloatND",

    "ToComplex64", "ToJustComplex64",
    "ToComplex64_1D", "ToJustComplex64_1D",
    "ToComplex64_2D", "ToJustComplex64_2D",
    "ToComplex64_3D", "ToJustComplex64_3D",
    "ToComplex64_ND", "ToJustComplex64_ND",
    "ToComplex64Strict1D", "ToJustComplex64Strict1D",
    "ToComplex64Strict2D", "ToJustComplex64Strict2D",
    "ToComplex64Strict3D", "ToJustComplex64Strict3D",

    "ToComplex128", "ToJustComplex128",
    "ToComplex128_1D", "ToJustComplex128_1D",
    "ToComplex128_2D", "ToJustComplex128_2D",
    "ToComplex128_3D", "ToJustComplex128_3D",
    "ToComplex128_ND", "ToJustComplex128_ND",
    "ToComplex128Strict1D", "ToJustComplex128Strict1D",
    "ToComplex128Strict2D", "ToJustComplex128Strict2D",
    "ToComplex128Strict3D", "ToJustComplex128Strict3D",

    "ToJustCLongDouble",
    "ToJustCLongDouble1D", "ToJustCLongDoubleStrict1D",
    "ToJustCLongDouble2D", "ToJustCLongDoubleStrict2D",
    "ToJustCLongDouble3D", "ToJustCLongDoubleStrict3D",
    "ToJustCLongDoubleND",

    "ToComplex", "ToJustComplex",
    "ToComplex1D", "ToJustComplex1D", "ToComplexStrict1D", "ToJustComplexStrict1D",
    "ToComplex2D", "ToJustComplex2D", "ToComplexStrict2D", "ToJustComplexStrict2D",
    "ToComplex3D", "ToJustComplex3D", "ToComplexStrict3D", "ToJustComplexStrict3D",
    "ToComplexND", "ToJustComplexND",
]  # fmt: skip


def __dir__() -> list[str]:
    return __all__


###


_PyBool: TypeAlias = bool | Literal[0, 1]  # 0 and 1 are sometimes used as bool values
_PyScalar: TypeAlias = complex | bytes | str  # `complex` equivs `complex | float | int`


T = TypeVar("T", default=_PyScalar)
SCT = TypeVar("SCT", bound=np.generic, default=Any)
SCT_co = TypeVar("SCT_co", bound=np.generic, covariant=True)


# unlike `optype.numpy.CanArray0D` and `CanArrayND`, these one also accepts scalar types
# (and aren't runtime checkable)


_CanArrayStrict1D = TypeAliasType(
    "_CanArrayStrict1D",
    nptc.CanArray[tuple[int], np.dtype[SCT]],
    type_params=(SCT,),
)
_CanArrayStrict2D = TypeAliasType(
    "_CanArrayStrict2D",
    nptc.CanArray[tuple[int, int], np.dtype[SCT]],
    type_params=(SCT,),
)
_CanArrayStrict3D = TypeAliasType(
    "_CanArrayStrict3D",
    nptc.CanArray[tuple[int, int, int], np.dtype[SCT]],
    type_params=(SCT,),
)


class _CanArrayND(Protocol[SCT_co]):
    def __len__(self, /) -> int: ...
    def __array__(self, /) -> np.ndarray[Any, np.dtype[SCT_co]]: ...


class _CanArray(Protocol[SCT_co]):
    def __array__(self, /) -> np.ndarray[Any, np.dtype[SCT_co]]: ...


_To1D1 = TypeAliasType("_To1D1", _CanArrayND[SCT] | Seq[SCT], type_params=(SCT,))
_To1D2 = TypeAliasType("_To1D2", _CanArrayND[SCT] | Seq[T | SCT], type_params=(T, SCT))

_To2D1 = TypeAliasType(
    "_To2D1",
    _CanArrayND[SCT] | Seq[_To1D1[SCT]],
    type_params=(SCT,),
)
_To2D2 = TypeAliasType(
    "_To2D2",
    _CanArrayND[SCT] | Seq[_To1D2[T, SCT]],
    type_params=(T, SCT),
)

_To3D1 = TypeAliasType(
    "_To3D1",
    _CanArrayND[SCT] | Seq[_To2D1[SCT]],
    type_params=(SCT,),
)
_To3D2 = TypeAliasType(
    "_To3D2",
    _CanArrayND[SCT] | Seq[_To2D2[T, SCT]],
    type_params=(T, SCT),
)

_ToND1 = TypeAliasType(
    "_ToND1",
    _CanArrayND[SCT] | SeqND[_CanArray[SCT]],
    type_params=(SCT,),
)
_ToND2 = TypeAliasType(
    "_ToND2",
    _CanArrayND[SCT] | SeqND[T | _CanArray[SCT]],
    type_params=(T, SCT),
)

_ToStrict1D1 = TypeAliasType(
    "_ToStrict1D1",
    _CanArrayStrict1D[SCT] | Seq[SCT],
    type_params=(SCT,),
)
_ToStrict1D2 = TypeAliasType(
    "_ToStrict1D2",
    _CanArrayStrict1D[SCT] | Seq[T | SCT],
    type_params=(T, SCT),
)

_ToStrict2D1 = TypeAliasType(
    "_ToStrict2D1",
    _CanArrayStrict2D[SCT] | Seq[_ToStrict1D1[SCT]],
    type_params=(SCT,),
)
_ToStrict2D2 = TypeAliasType(
    "_ToStrict2D2",
    _CanArrayStrict2D[SCT] | Seq[_ToStrict1D2[T, SCT]],
    type_params=(T, SCT),
)

_ToStrict3D1 = TypeAliasType(
    "_ToStrict3D1",
    _CanArrayStrict3D[SCT] | Seq[_ToStrict2D1[SCT]],
    type_params=(SCT,),
)
_ToStrict3D2 = TypeAliasType(
    "_ToStrict3D2",
    _CanArrayStrict3D[SCT] | Seq[_ToStrict2D2[T, SCT]],
    type_params=(T, SCT),
)


###

# TODO(jorenham): export & document
# https://github.com/jorenham/optype/issues/373

integer_co = TypeAliasType("integer_co", npc.integer | np.bool_)
floating_co = TypeAliasType("floating_co", npc.floating | npc.integer | np.bool_)
complexfloating_co = TypeAliasType("complexfloating_co", npc.number | np.bool_)

# promotion rules with safe casting mode
f16_co = TypeAliasType("f16_co", npc.floating16 | npc.integer8 | np.bool_)
f32_co = TypeAliasType(
    "f32_co",
    npc.floating32 | npc.floating16 | npc.integer16 | npc.integer8 | np.bool_,
)
c64_co = TypeAliasType(
    "c64_co",
    npc.inexact32 | npc.number16 | npc.integer8 | np.bool_,
)
f64_co = TypeAliasType(
    "f64_co",
    npc.floating64 | npc.floating32 | npc.floating16 | npc.integer | np.bool_,
)
c128_co = TypeAliasType(
    "c128_co",
    npc.number64 | npc.number32 | npc.number16 | npc.integer | np.bool_,
)


###

# scalar- and array-likes, with "coercible" shape-types

ToScalar = TypeAliasType("ToScalar", _PyScalar | np.generic)
ToArray1D = TypeAliasType("ToArray1D", _To1D2[T, SCT], type_params=(T, SCT))
ToArray2D = TypeAliasType("ToArray2D", _To2D2[T, SCT], type_params=(T, SCT))
ToArray3D = TypeAliasType("ToArray3D", _To3D2[T, SCT], type_params=(T, SCT))
ToArrayND = TypeAliasType("ToArrayND", _ToND2[T, SCT], type_params=(T, SCT))

ToFalse = TypeAliasType("ToFalse", nptc.LiteralFalse | Literal[0])
ToTrue = TypeAliasType("ToTrue", nptc.LiteralTrue | Literal[1])

ToJustFalse = TypeAliasType("ToJustFalse", nptc.LiteralFalse)
ToJustTrue = TypeAliasType("ToJustTrue", nptc.LiteralTrue)

ToBool = TypeAliasType("ToBool", _PyBool | np.bool_)
ToBool1D = TypeAliasType("ToBool1D", _To1D2[_PyBool, np.bool_])
ToBool2D = TypeAliasType("ToBool2D", _To2D2[_PyBool, np.bool_])
ToBool3D = TypeAliasType("ToBool3D", _To3D2[_PyBool, np.bool_])
ToBoolND = TypeAliasType("ToBoolND", _ToND2[_PyBool, np.bool_])

ToInt = TypeAliasType("ToInt", int | integer_co)
ToInt1D = TypeAliasType("ToInt1D", _To1D2[int, integer_co])
ToInt2D = TypeAliasType("ToInt2D", _To2D2[int, integer_co])
ToInt3D = TypeAliasType("ToInt3D", _To3D2[int, integer_co])
ToIntND = TypeAliasType("ToIntND", _ToND2[int, integer_co])

ToFloat16 = TypeAliasType("ToFloat16", f16_co)
ToFloat16_1D = TypeAliasType("ToFloat16_1D", _To1D1[f16_co])
ToFloat16_2D = TypeAliasType("ToFloat16_2D", _To2D1[f16_co])
ToFloat16_3D = TypeAliasType("ToFloat16_3D", _To3D1[f16_co])
ToFloat16_ND = TypeAliasType("ToFloat16_ND", _ToND1[f16_co])

ToFloat32 = TypeAliasType("ToFloat32", f32_co)
ToFloat32_1D = TypeAliasType("ToFloat32_1D", _To1D1[f32_co])
ToFloat32_2D = TypeAliasType("ToFloat32_2D", _To2D1[f32_co])
ToFloat32_3D = TypeAliasType("ToFloat32_3D", _To3D1[f32_co])
ToFloat32_ND = TypeAliasType("ToFloat32_ND", _ToND1[f32_co])

ToFloat64 = TypeAliasType("ToFloat64", float | f64_co)
ToFloat64_1D = TypeAliasType("ToFloat64_1D", _To1D2[float, f64_co])
ToFloat64_2D = TypeAliasType("ToFloat64_2D", _To2D2[float, f64_co])
ToFloat64_3D = TypeAliasType("ToFloat64_3D", _To3D2[float, f64_co])
ToFloat64_ND = TypeAliasType("ToFloat64_ND", _ToND2[float, f64_co])

ToFloat = TypeAliasType("ToFloat", float | floating_co)
ToFloat1D = TypeAliasType("ToFloat1D", _To1D2[float, floating_co])
ToFloat2D = TypeAliasType("ToFloat2D", _To2D2[float, floating_co])
ToFloat3D = TypeAliasType("ToFloat3D", _To3D2[float, floating_co])
ToFloatND = TypeAliasType("ToFloatND", _ToND2[float, floating_co])

ToComplex64 = TypeAliasType("ToComplex64", c64_co)
ToComplex64_1D = TypeAliasType("ToComplex64_1D", _To1D1[c64_co])
ToComplex64_2D = TypeAliasType("ToComplex64_2D", _To2D1[c64_co])
ToComplex64_3D = TypeAliasType("ToComplex64_3D", _To3D1[c64_co])
ToComplex64_ND = TypeAliasType("ToComplex64_ND", _ToND1[c64_co])

ToComplex128 = TypeAliasType("ToComplex128", complex | c128_co)
ToComplex128_1D = TypeAliasType("ToComplex128_1D", _To1D2[complex, c128_co])
ToComplex128_2D = TypeAliasType("ToComplex128_2D", _To2D2[complex, c128_co])
ToComplex128_3D = TypeAliasType("ToComplex128_3D", _To3D2[complex, c128_co])
ToComplex128_ND = TypeAliasType("ToComplex128_ND", _ToND2[complex, c128_co])

ToComplex = TypeAliasType("ToComplex", complex | complexfloating_co)
ToComplex1D = TypeAliasType("ToComplex1D", _To1D2[complex, complexfloating_co])
ToComplex2D = TypeAliasType("ToComplex2D", _To2D2[complex, complexfloating_co])
ToComplex3D = TypeAliasType("ToComplex3D", _To3D2[complex, complexfloating_co])
ToComplexND = TypeAliasType("ToComplexND", _ToND2[complex, complexfloating_co])

# scalar- and array-likes, with "just" that scalar type

ToJustBool = TypeAliasType("ToJustBool", bool | np.bool_)
ToJustBool1D = TypeAliasType("ToJustBool1D", _To1D2[bool, np.bool_])
ToJustBool2D = TypeAliasType("ToJustBool2D", _To2D2[bool, np.bool_])
ToJustBool3D = TypeAliasType("ToJustBool3D", _To3D2[bool, np.bool_])
ToJustBoolND = TypeAliasType("ToJustBoolND", _ToND2[bool, np.bool_])

ToJustInt64 = TypeAliasType("ToJustInt64", JustInt | np.intp)
ToJustInt64_1D = TypeAliasType("ToJustInt64_1D", _To1D2[JustInt, np.intp])
ToJustInt64_2D = TypeAliasType("ToJustInt64_2D", _To2D2[JustInt, np.intp])
ToJustInt64_3D = TypeAliasType("ToJustInt64_3D", _To3D2[JustInt, np.intp])
ToJustInt64_ND = TypeAliasType("ToJustInt64_ND", _ToND2[JustInt, np.intp])

ToJustInt = TypeAliasType("ToJustInt", JustInt | npc.integer)
ToJustInt1D = TypeAliasType("ToJustInt1D", _To1D2[JustInt, npc.integer])
ToJustInt2D = TypeAliasType("ToJustInt2D", _To2D2[JustInt, npc.integer])
ToJustInt3D = TypeAliasType("ToJustInt3D", _To3D2[JustInt, npc.integer])
ToJustIntND = TypeAliasType("ToJustIntND", _ToND2[JustInt, npc.integer])

ToJustFloat16 = TypeAliasType("ToJustFloat16", npc.floating16)
ToJustFloat16_1D = TypeAliasType("ToJustFloat16_1D", _To1D1[npc.floating16])
ToJustFloat16_2D = TypeAliasType("ToJustFloat16_2D", _To2D1[npc.floating16])
ToJustFloat16_3D = TypeAliasType("ToJustFloat16_3D", _To3D1[npc.floating16])
ToJustFloat16_ND = TypeAliasType("ToJustFloat16_ND", _ToND1[npc.floating16])

ToJustFloat32 = TypeAliasType("ToJustFloat32", npc.floating32)
ToJustFloat32_1D = TypeAliasType("ToJustFloat32_1D", _To1D1[npc.floating32])
ToJustFloat32_2D = TypeAliasType("ToJustFloat32_2D", _To2D1[npc.floating32])
ToJustFloat32_3D = TypeAliasType("ToJustFloat32_3D", _To3D1[npc.floating32])
ToJustFloat32_ND = TypeAliasType("ToJustFloat32_ND", _ToND1[npc.floating32])

ToJustFloat64 = TypeAliasType("ToJustFloat64", JustFloat | npc.floating64)
ToJustFloat64_1D = TypeAliasType("ToJustFloat64_1D", _To1D2[JustFloat, npc.floating64])
ToJustFloat64_2D = TypeAliasType("ToJustFloat64_2D", _To2D2[JustFloat, npc.floating64])
ToJustFloat64_3D = TypeAliasType("ToJustFloat64_3D", _To3D2[JustFloat, npc.floating64])
ToJustFloat64_ND = TypeAliasType("ToJustFloat64_ND", _ToND2[JustFloat, npc.floating64])

ToJustLongDouble = TypeAliasType("ToJustLongDouble", npc.floating80)
ToJustLongDouble1D = TypeAliasType("ToJustLongDouble1D", _To1D1[npc.floating80])
ToJustLongDouble2D = TypeAliasType("ToJustLongDouble2D", _To2D1[npc.floating80])
ToJustLongDouble3D = TypeAliasType("ToJustLongDouble3D", _To3D1[npc.floating80])
ToJustLongDoubleND = TypeAliasType("ToJustLongDoubleND", _ToND1[npc.floating80])

ToJustFloat = TypeAliasType("ToJustFloat", JustFloat | npc.floating)
ToJustFloat1D = TypeAliasType("ToJustFloat1D", _To1D2[JustFloat, npc.floating])
ToJustFloat2D = TypeAliasType("ToJustFloat2D", _To2D2[JustFloat, npc.floating])
ToJustFloat3D = TypeAliasType("ToJustFloat3D", _To3D2[JustFloat, npc.floating])
ToJustFloatND = TypeAliasType("ToJustFloatND", _ToND2[JustFloat, npc.floating])

ToJustComplex64 = TypeAliasType("ToJustComplex64", npc.complexfloating64)
ToJustComplex64_1D = TypeAliasType("ToJustComplex64_1D", _To1D1[npc.complexfloating64])
ToJustComplex64_2D = TypeAliasType("ToJustComplex64_2D", _To2D1[npc.complexfloating64])
ToJustComplex64_3D = TypeAliasType("ToJustComplex64_3D", _To3D1[npc.complexfloating64])
ToJustComplex64_ND = TypeAliasType("ToJustComplex64_ND", _ToND1[npc.complexfloating64])

ToJustComplex128 = TypeAliasType(
    "ToJustComplex128",
    JustComplex | npc.complexfloating128,
)
ToJustComplex128_1D = TypeAliasType(
    "ToJustComplex128_1D",
    _To1D2[JustComplex, npc.complexfloating128],
)
ToJustComplex128_2D = TypeAliasType(
    "ToJustComplex128_2D",
    _To2D2[JustComplex, npc.complexfloating128],
)
ToJustComplex128_3D = TypeAliasType(
    "ToJustComplex128_3D",
    _To3D2[JustComplex, npc.complexfloating128],
)
ToJustComplex128_ND = TypeAliasType(
    "ToJustComplex128_ND",
    _ToND2[JustComplex, npc.complexfloating128],
)

ToJustCLongDouble = TypeAliasType("ToJustCLongDouble", npc.complexfloating160)
ToJustCLongDouble1D = TypeAliasType(
    "ToJustCLongDouble1D",
    _To1D1[npc.complexfloating160],
)
ToJustCLongDouble2D = TypeAliasType(
    "ToJustCLongDouble2D",
    _To2D1[npc.complexfloating160],
)
ToJustCLongDouble3D = TypeAliasType(
    "ToJustCLongDouble3D",
    _To3D1[npc.complexfloating160],
)
ToJustCLongDoubleND = TypeAliasType(
    "ToJustCLongDoubleND",
    _ToND1[npc.complexfloating160],
)

ToJustComplex = TypeAliasType("ToJustComplex", JustComplex | npc.complexfloating)
ToJustComplex1D = TypeAliasType(
    "ToJustComplex1D",
    _To1D2[JustComplex, npc.complexfloating],
)
ToJustComplex2D = TypeAliasType(
    "ToJustComplex2D",
    _To2D2[JustComplex, npc.complexfloating],
)
ToJustComplex3D = TypeAliasType(
    "ToJustComplex3D",
    _To3D2[JustComplex, npc.complexfloating],
)
ToJustComplexND = TypeAliasType(
    "ToJustComplexND",
    _ToND2[JustComplex, npc.complexfloating],
)

# array-likes, with "coercible" shape-types, and "strict" shape-types

ToArrayStrict1D = TypeAliasType(
    "ToArrayStrict1D",
    _ToStrict1D2[T, SCT],
    type_params=(T, SCT),
)
ToArrayStrict2D = TypeAliasType(
    "ToArrayStrict2D",
    _ToStrict2D2[T, SCT],
    type_params=(T, SCT),
)
ToArrayStrict3D = TypeAliasType(
    "ToArrayStrict3D",
    _ToStrict3D2[T, SCT],
    type_params=(T, SCT),
)

ToBoolStrict1D = TypeAliasType("ToBoolStrict1D", _ToStrict1D2[_PyBool, np.bool_])
ToBoolStrict2D = TypeAliasType("ToBoolStrict2D", _ToStrict2D2[_PyBool, np.bool_])
ToBoolStrict3D = TypeAliasType("ToBoolStrict3D", _ToStrict3D2[_PyBool, np.bool_])

ToIntStrict1D = TypeAliasType("ToIntStrict1D", _ToStrict1D2[int, integer_co])
ToIntStrict2D = TypeAliasType("ToIntStrict2D", _ToStrict2D2[int, integer_co])
ToIntStrict3D = TypeAliasType("ToIntStrict3D", _ToStrict3D2[int, integer_co])

ToFloat16Strict1D = TypeAliasType("ToFloat16Strict1D", _ToStrict1D1[f16_co])
ToFloat16Strict2D = TypeAliasType("ToFloat16Strict2D", _ToStrict2D1[f16_co])
ToFloat16Strict3D = TypeAliasType("ToFloat16Strict3D", _ToStrict3D1[f16_co])

ToFloat32Strict1D = TypeAliasType("ToFloat32Strict1D", _ToStrict1D1[f32_co])
ToFloat32Strict2D = TypeAliasType("ToFloat32Strict2D", _ToStrict2D1[f32_co])
ToFloat32Strict3D = TypeAliasType("ToFloat32Strict3D", _ToStrict3D1[f32_co])

ToFloat64Strict1D = TypeAliasType("ToFloat64Strict1D", _ToStrict1D2[float, f64_co])
ToFloat64Strict2D = TypeAliasType("ToFloat64Strict2D", _ToStrict2D2[float, f64_co])
ToFloat64Strict3D = TypeAliasType("ToFloat64Strict3D", _ToStrict3D2[float, f64_co])

ToFloatStrict1D = TypeAliasType("ToFloatStrict1D", _ToStrict1D2[float, floating_co])
ToFloatStrict2D = TypeAliasType("ToFloatStrict2D", _ToStrict2D2[float, floating_co])
ToFloatStrict3D = TypeAliasType("ToFloatStrict3D", _ToStrict3D2[float, floating_co])

ToComplex64Strict1D = TypeAliasType("ToComplex64Strict1D", _ToStrict1D1[c64_co])
ToComplex64Strict2D = TypeAliasType("ToComplex64Strict2D", _ToStrict2D1[c64_co])
ToComplex64Strict3D = TypeAliasType("ToComplex64Strict3D", _ToStrict3D1[c64_co])

ToComplex128Strict1D = TypeAliasType(
    "ToComplex128Strict1D",
    _ToStrict1D2[complex, c128_co],
)
ToComplex128Strict2D = TypeAliasType(
    "ToComplex128Strict2D",
    _ToStrict2D2[complex, c128_co],
)
ToComplex128Strict3D = TypeAliasType(
    "ToComplex128Strict3D",
    _ToStrict3D2[complex, c128_co],
)

ToComplexStrict1D = TypeAliasType(
    "ToComplexStrict1D",
    _ToStrict1D2[complex, complexfloating_co],
)
ToComplexStrict2D = TypeAliasType(
    "ToComplexStrict2D",
    _ToStrict2D2[complex, complexfloating_co],
)
ToComplexStrict3D = TypeAliasType(
    "ToComplexStrict3D",
    _ToStrict3D2[complex, complexfloating_co],
)

# array-likes, with "just" that scalar type, and "strict" shape-types

ToJustBoolStrict1D = TypeAliasType("ToJustBoolStrict1D", _ToStrict1D2[bool, np.bool_])
ToJustBoolStrict2D = TypeAliasType("ToJustBoolStrict2D", _ToStrict2D2[bool, np.bool_])
ToJustBoolStrict3D = TypeAliasType("ToJustBoolStrict3D", _ToStrict3D2[bool, np.bool_])
ToJustInt64Strict1D = TypeAliasType(
    "ToJustInt64Strict1D",
    _ToStrict1D2[JustInt, np.int64],
)
ToJustInt64Strict2D = TypeAliasType(
    "ToJustInt64Strict2D",
    _ToStrict2D2[JustInt, np.int64],
)
ToJustInt64Strict3D = TypeAliasType(
    "ToJustInt64Strict3D",
    _ToStrict3D2[JustInt, np.int64],
)
ToJustIntStrict1D = TypeAliasType(
    "ToJustIntStrict1D",
    _ToStrict1D2[JustInt, npc.integer],
)
ToJustIntStrict2D = TypeAliasType(
    "ToJustIntStrict2D",
    _ToStrict2D2[JustInt, npc.integer],
)
ToJustIntStrict3D = TypeAliasType(
    "ToJustIntStrict3D",
    _ToStrict3D2[JustInt, npc.integer],
)
ToJustFloat16Strict1D = TypeAliasType(
    "ToJustFloat16Strict1D",
    _ToStrict1D1[npc.floating16],
)
ToJustFloat16Strict2D = TypeAliasType(
    "ToJustFloat16Strict2D",
    _ToStrict2D1[npc.floating16],
)
ToJustFloat16Strict3D = TypeAliasType(
    "ToJustFloat16Strict3D",
    _ToStrict3D1[npc.floating16],
)
ToJustFloat32Strict1D = TypeAliasType(
    "ToJustFloat32Strict1D",
    _ToStrict1D1[npc.floating32],
)
ToJustFloat32Strict2D = TypeAliasType(
    "ToJustFloat32Strict2D",
    _ToStrict2D1[npc.floating32],
)
ToJustFloat32Strict3D = TypeAliasType(
    "ToJustFloat32Strict3D",
    _ToStrict3D1[npc.floating32],
)
ToJustFloat64Strict1D = TypeAliasType(
    "ToJustFloat64Strict1D",
    _ToStrict1D2[JustFloat, npc.floating64],
)
ToJustFloat64Strict2D = TypeAliasType(
    "ToJustFloat64Strict2D",
    _ToStrict2D2[JustFloat, npc.floating64],
)
ToJustFloat64Strict3D = TypeAliasType(
    "ToJustFloat64Strict3D",
    _ToStrict3D2[JustFloat, npc.floating64],
)
ToJustLongDoubleStrict1D = TypeAliasType(
    "ToJustLongDoubleStrict1D",
    _ToStrict1D1[npc.floating80],
)
ToJustLongDoubleStrict2D = TypeAliasType(
    "ToJustLongDoubleStrict2D",
    _ToStrict2D1[npc.floating80],
)
ToJustLongDoubleStrict3D = TypeAliasType(
    "ToJustLongDoubleStrict3D",
    _ToStrict3D1[npc.floating80],
)
ToJustFloatStrict1D = TypeAliasType(
    "ToJustFloatStrict1D",
    _ToStrict1D2[JustFloat, npc.floating],
)
ToJustFloatStrict2D = TypeAliasType(
    "ToJustFloatStrict2D",
    _ToStrict2D2[JustFloat, npc.floating],
)
ToJustFloatStrict3D = TypeAliasType(
    "ToJustFloatStrict3D",
    _ToStrict3D2[JustFloat, npc.floating],
)
ToJustComplex64Strict1D = TypeAliasType(
    "ToJustComplex64Strict1D",
    _ToStrict1D1[npc.complexfloating64],
)
ToJustComplex64Strict2D = TypeAliasType(
    "ToJustComplex64Strict2D",
    _ToStrict2D1[npc.complexfloating64],
)
ToJustComplex64Strict3D = TypeAliasType(
    "ToJustComplex64Strict3D",
    _ToStrict3D1[npc.complexfloating64],
)
ToJustComplex128Strict1D = TypeAliasType(
    "ToJustComplex128Strict1D",
    _ToStrict1D2[JustComplex, npc.complexfloating128],
)
ToJustComplex128Strict2D = TypeAliasType(
    "ToJustComplex128Strict2D",
    _ToStrict2D2[JustComplex, npc.complexfloating128],
)
ToJustComplex128Strict3D = TypeAliasType(
    "ToJustComplex128Strict3D",
    _ToStrict3D2[JustComplex, npc.complexfloating128],
)
ToJustCLongDoubleStrict1D = TypeAliasType(
    "ToJustCLongDoubleStrict1D",
    _ToStrict1D1[npc.complexfloating160],
)
ToJustCLongDoubleStrict2D = TypeAliasType(
    "ToJustCLongDoubleStrict2D",
    _ToStrict2D1[npc.complexfloating160],
)
ToJustCLongDoubleStrict3D = TypeAliasType(
    "ToJustCLongDoubleStrict3D",
    _ToStrict3D1[npc.complexfloating160],
)
ToJustComplexStrict1D = TypeAliasType(
    "ToJustComplexStrict1D",
    _ToStrict1D2[JustComplex, npc.complexfloating],
)
ToJustComplexStrict2D = TypeAliasType(
    "ToJustComplexStrict2D",
    _ToStrict2D2[JustComplex, npc.complexfloating],
)
ToJustComplexStrict3D = TypeAliasType(
    "ToJustComplexStrict3D",
    _ToStrict3D2[JustComplex, npc.complexfloating],
)
