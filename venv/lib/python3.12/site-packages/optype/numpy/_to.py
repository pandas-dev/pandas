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
f16_co = TypeAliasType(
    "f16_co",
    npc.floating16 | npc.integer8 | np.bool_,
)
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

ToScalar: TypeAlias = _PyScalar | np.generic
ToArray1D: TypeAlias = _To1D2[T, SCT]
ToArray2D: TypeAlias = _To2D2[T, SCT]
ToArray3D: TypeAlias = _To3D2[T, SCT]
ToArrayND: TypeAlias = _ToND2[T, SCT]

ToFalse = TypeAliasType("ToFalse", nptc.LiteralFalse | Literal[0])
ToTrue = TypeAliasType("ToTrue", nptc.LiteralTrue | Literal[1])

ToJustFalse = TypeAliasType("ToJustFalse", nptc.LiteralFalse)
ToJustTrue = TypeAliasType("ToJustTrue", nptc.LiteralTrue)

ToBool: TypeAlias = _PyBool | np.bool_
ToBool1D: TypeAlias = _To1D2[_PyBool, np.bool_]
ToBool2D: TypeAlias = _To2D2[_PyBool, np.bool_]
ToBool3D: TypeAlias = _To3D2[_PyBool, np.bool_]
ToBoolND: TypeAlias = _ToND2[_PyBool, np.bool_]

ToInt: TypeAlias = int | integer_co
ToInt1D: TypeAlias = _To1D2[int, integer_co]
ToInt2D: TypeAlias = _To2D2[int, integer_co]
ToInt3D: TypeAlias = _To3D2[int, integer_co]
ToIntND: TypeAlias = _ToND2[int, integer_co]

ToFloat16: TypeAlias = f16_co
ToFloat16_1D: TypeAlias = _To1D1[f16_co]
ToFloat16_2D: TypeAlias = _To2D1[f16_co]
ToFloat16_3D: TypeAlias = _To3D1[f16_co]
ToFloat16_ND: TypeAlias = _ToND1[f16_co]

ToFloat32: TypeAlias = f32_co
ToFloat32_1D: TypeAlias = _To1D1[f32_co]
ToFloat32_2D: TypeAlias = _To2D1[f32_co]
ToFloat32_3D: TypeAlias = _To3D1[f32_co]
ToFloat32_ND: TypeAlias = _ToND1[f32_co]

ToFloat64: TypeAlias = float | f64_co
ToFloat64_1D: TypeAlias = _To1D2[float, f64_co]
ToFloat64_2D: TypeAlias = _To2D2[float, f64_co]
ToFloat64_3D: TypeAlias = _To3D2[float, f64_co]
ToFloat64_ND: TypeAlias = _ToND2[float, f64_co]

ToFloat: TypeAlias = float | floating_co
ToFloat1D: TypeAlias = _To1D2[float, floating_co]
ToFloat2D: TypeAlias = _To2D2[float, floating_co]
ToFloat3D: TypeAlias = _To3D2[float, floating_co]
ToFloatND: TypeAlias = _ToND2[float, floating_co]

ToComplex64: TypeAlias = c64_co
ToComplex64_1D: TypeAlias = _To1D1[c64_co]
ToComplex64_2D: TypeAlias = _To2D1[c64_co]
ToComplex64_3D: TypeAlias = _To3D1[c64_co]
ToComplex64_ND: TypeAlias = _ToND1[c64_co]

ToComplex128: TypeAlias = complex | c128_co
ToComplex128_1D: TypeAlias = _To1D2[complex, c128_co]
ToComplex128_2D: TypeAlias = _To2D2[complex, c128_co]
ToComplex128_3D: TypeAlias = _To3D2[complex, c128_co]
ToComplex128_ND: TypeAlias = _ToND2[complex, c128_co]

ToComplex: TypeAlias = complex | complexfloating_co
ToComplex1D: TypeAlias = _To1D2[complex, complexfloating_co]
ToComplex2D: TypeAlias = _To2D2[complex, complexfloating_co]
ToComplex3D: TypeAlias = _To3D2[complex, complexfloating_co]
ToComplexND: TypeAlias = _ToND2[complex, complexfloating_co]

# scalar- and array-likes, with "just" that scalar type

ToJustBool: TypeAlias = bool | np.bool_
ToJustBool1D: TypeAlias = _To1D2[bool, np.bool_]
ToJustBool2D: TypeAlias = _To2D2[bool, np.bool_]
ToJustBool3D: TypeAlias = _To3D2[bool, np.bool_]
ToJustBoolND: TypeAlias = _ToND2[bool, np.bool_]

ToJustInt64: TypeAlias = JustInt | np.intp
ToJustInt64_1D: TypeAlias = _To1D2[JustInt, np.intp]
ToJustInt64_2D: TypeAlias = _To2D2[JustInt, np.intp]
ToJustInt64_3D: TypeAlias = _To3D2[JustInt, np.intp]
ToJustInt64_ND: TypeAlias = _ToND2[JustInt, np.intp]

ToJustInt: TypeAlias = JustInt | npc.integer
ToJustInt1D: TypeAlias = _To1D2[JustInt, npc.integer]
ToJustInt2D: TypeAlias = _To2D2[JustInt, npc.integer]
ToJustInt3D: TypeAlias = _To3D2[JustInt, npc.integer]
ToJustIntND: TypeAlias = _ToND2[JustInt, npc.integer]

ToJustFloat16: TypeAlias = npc.floating16
ToJustFloat16_1D: TypeAlias = _To1D1[npc.floating16]
ToJustFloat16_2D: TypeAlias = _To2D1[npc.floating16]
ToJustFloat16_3D: TypeAlias = _To3D1[npc.floating16]
ToJustFloat16_ND: TypeAlias = _ToND1[npc.floating16]

ToJustFloat32: TypeAlias = npc.floating32
ToJustFloat32_1D: TypeAlias = _To1D1[npc.floating32]
ToJustFloat32_2D: TypeAlias = _To2D1[npc.floating32]
ToJustFloat32_3D: TypeAlias = _To3D1[npc.floating32]
ToJustFloat32_ND: TypeAlias = _ToND1[npc.floating32]

ToJustFloat64: TypeAlias = JustFloat | npc.floating64
ToJustFloat64_1D: TypeAlias = _To1D2[JustFloat, npc.floating64]
ToJustFloat64_2D: TypeAlias = _To2D2[JustFloat, npc.floating64]
ToJustFloat64_3D: TypeAlias = _To3D2[JustFloat, npc.floating64]
ToJustFloat64_ND: TypeAlias = _ToND2[JustFloat, npc.floating64]

ToJustLongDouble: TypeAlias = npc.floating80
ToJustLongDouble1D: TypeAlias = _To1D1[npc.floating80]
ToJustLongDouble2D: TypeAlias = _To2D1[npc.floating80]
ToJustLongDouble3D: TypeAlias = _To3D1[npc.floating80]
ToJustLongDoubleND: TypeAlias = _ToND1[npc.floating80]

ToJustFloat: TypeAlias = JustFloat | npc.floating
ToJustFloat1D: TypeAlias = _To1D2[JustFloat, npc.floating]
ToJustFloat2D: TypeAlias = _To2D2[JustFloat, npc.floating]
ToJustFloat3D: TypeAlias = _To3D2[JustFloat, npc.floating]
ToJustFloatND: TypeAlias = _ToND2[JustFloat, npc.floating]

ToJustComplex64: TypeAlias = npc.complexfloating64
ToJustComplex64_1D: TypeAlias = _To1D1[npc.complexfloating64]
ToJustComplex64_2D: TypeAlias = _To2D1[npc.complexfloating64]
ToJustComplex64_3D: TypeAlias = _To3D1[npc.complexfloating64]
ToJustComplex64_ND: TypeAlias = _ToND1[npc.complexfloating64]

ToJustComplex128: TypeAlias = JustComplex | npc.complexfloating128
ToJustComplex128_1D: TypeAlias = _To1D2[JustComplex, npc.complexfloating128]
ToJustComplex128_2D: TypeAlias = _To2D2[JustComplex, npc.complexfloating128]
ToJustComplex128_3D: TypeAlias = _To3D2[JustComplex, npc.complexfloating128]
ToJustComplex128_ND: TypeAlias = _ToND2[JustComplex, npc.complexfloating128]

ToJustCLongDouble: TypeAlias = npc.complexfloating160
ToJustCLongDouble1D: TypeAlias = _To1D1[npc.complexfloating160]
ToJustCLongDouble2D: TypeAlias = _To2D1[npc.complexfloating160]
ToJustCLongDouble3D: TypeAlias = _To3D1[npc.complexfloating160]
ToJustCLongDoubleND: TypeAlias = _ToND1[npc.complexfloating160]

ToJustComplex: TypeAlias = JustComplex | npc.complexfloating
ToJustComplex1D: TypeAlias = _To1D2[JustComplex, npc.complexfloating]
ToJustComplex2D: TypeAlias = _To2D2[JustComplex, npc.complexfloating]
ToJustComplex3D: TypeAlias = _To3D2[JustComplex, npc.complexfloating]
ToJustComplexND: TypeAlias = _ToND2[JustComplex, npc.complexfloating]

# array-likes, with "coercible" shape-types, and "strict" shape-types

ToArrayStrict1D: TypeAlias = _ToStrict1D2[T, SCT]
ToArrayStrict2D: TypeAlias = _ToStrict2D2[T, SCT]
ToArrayStrict3D: TypeAlias = _ToStrict3D2[T, SCT]

ToBoolStrict1D: TypeAlias = _ToStrict1D2[_PyBool, np.bool_]
ToBoolStrict2D: TypeAlias = _ToStrict2D2[_PyBool, np.bool_]
ToBoolStrict3D: TypeAlias = _ToStrict3D2[_PyBool, np.bool_]

ToIntStrict1D: TypeAlias = _ToStrict1D2[int, integer_co]
ToIntStrict2D: TypeAlias = _ToStrict2D2[int, integer_co]
ToIntStrict3D: TypeAlias = _ToStrict3D2[int, integer_co]

ToFloat16Strict1D: TypeAlias = _ToStrict1D1[f16_co]
ToFloat16Strict2D: TypeAlias = _ToStrict2D1[f16_co]
ToFloat16Strict3D: TypeAlias = _ToStrict3D1[f16_co]

ToFloat32Strict1D: TypeAlias = _ToStrict1D1[f32_co]
ToFloat32Strict2D: TypeAlias = _ToStrict2D1[f32_co]
ToFloat32Strict3D: TypeAlias = _ToStrict3D1[f32_co]

ToFloat64Strict1D: TypeAlias = _ToStrict1D2[float, f64_co]
ToFloat64Strict2D: TypeAlias = _ToStrict2D2[float, f64_co]
ToFloat64Strict3D: TypeAlias = _ToStrict3D2[float, f64_co]

ToFloatStrict1D: TypeAlias = _ToStrict1D2[float, floating_co]
ToFloatStrict2D: TypeAlias = _ToStrict2D2[float, floating_co]
ToFloatStrict3D: TypeAlias = _ToStrict3D2[float, floating_co]

ToComplex64Strict1D: TypeAlias = _ToStrict1D1[c64_co]
ToComplex64Strict2D: TypeAlias = _ToStrict2D1[c64_co]
ToComplex64Strict3D: TypeAlias = _ToStrict3D1[c64_co]

ToComplex128Strict1D: TypeAlias = _ToStrict1D2[complex, c128_co]
ToComplex128Strict2D: TypeAlias = _ToStrict2D2[complex, c128_co]
ToComplex128Strict3D: TypeAlias = _ToStrict3D2[complex, c128_co]

ToComplexStrict1D: TypeAlias = _ToStrict1D2[complex, complexfloating_co]
ToComplexStrict2D: TypeAlias = _ToStrict2D2[complex, complexfloating_co]
ToComplexStrict3D: TypeAlias = _ToStrict3D2[complex, complexfloating_co]

# array-likes, with "just" that scalar type, and "strict" shape-types

ToJustBoolStrict1D: TypeAlias = _ToStrict1D2[bool, np.bool_]
ToJustBoolStrict2D: TypeAlias = _ToStrict2D2[bool, np.bool_]
ToJustBoolStrict3D: TypeAlias = _ToStrict3D2[bool, np.bool_]

ToJustInt64Strict1D: TypeAlias = _ToStrict1D2[JustInt, np.int64]
ToJustInt64Strict2D: TypeAlias = _ToStrict2D2[JustInt, np.int64]
ToJustInt64Strict3D: TypeAlias = _ToStrict3D2[JustInt, np.int64]

ToJustIntStrict1D: TypeAlias = _ToStrict1D2[JustInt, npc.integer]
ToJustIntStrict2D: TypeAlias = _ToStrict2D2[JustInt, npc.integer]
ToJustIntStrict3D: TypeAlias = _ToStrict3D2[JustInt, npc.integer]

ToJustFloat16Strict1D: TypeAlias = _ToStrict1D1[npc.floating16]
ToJustFloat16Strict2D: TypeAlias = _ToStrict2D1[npc.floating16]
ToJustFloat16Strict3D: TypeAlias = _ToStrict3D1[npc.floating16]

ToJustFloat32Strict1D: TypeAlias = _ToStrict1D1[npc.floating32]
ToJustFloat32Strict2D: TypeAlias = _ToStrict2D1[npc.floating32]
ToJustFloat32Strict3D: TypeAlias = _ToStrict3D1[npc.floating32]

ToJustFloat64Strict1D: TypeAlias = _ToStrict1D2[JustFloat, npc.floating64]
ToJustFloat64Strict2D: TypeAlias = _ToStrict2D2[JustFloat, npc.floating64]
ToJustFloat64Strict3D: TypeAlias = _ToStrict3D2[JustFloat, npc.floating64]

ToJustLongDoubleStrict1D: TypeAlias = _ToStrict1D1[npc.floating80]
ToJustLongDoubleStrict2D: TypeAlias = _ToStrict2D1[npc.floating80]
ToJustLongDoubleStrict3D: TypeAlias = _ToStrict3D1[npc.floating80]

ToJustFloatStrict1D: TypeAlias = _ToStrict1D2[JustFloat, npc.floating]
ToJustFloatStrict2D: TypeAlias = _ToStrict2D2[JustFloat, npc.floating]
ToJustFloatStrict3D: TypeAlias = _ToStrict3D2[JustFloat, npc.floating]

ToJustComplex64Strict1D: TypeAlias = _ToStrict1D1[npc.complexfloating64]
ToJustComplex64Strict2D: TypeAlias = _ToStrict2D1[npc.complexfloating64]
ToJustComplex64Strict3D: TypeAlias = _ToStrict3D1[npc.complexfloating64]

ToJustComplex128Strict1D: TypeAlias = _ToStrict1D2[JustComplex, npc.complexfloating128]
ToJustComplex128Strict2D: TypeAlias = _ToStrict2D2[JustComplex, npc.complexfloating128]
ToJustComplex128Strict3D: TypeAlias = _ToStrict3D2[JustComplex, npc.complexfloating128]

ToJustCLongDoubleStrict1D: TypeAlias = _ToStrict1D1[npc.complexfloating160]
ToJustCLongDoubleStrict2D: TypeAlias = _ToStrict2D1[npc.complexfloating160]
ToJustCLongDoubleStrict3D: TypeAlias = _ToStrict3D1[npc.complexfloating160]

ToJustComplexStrict1D: TypeAlias = _ToStrict1D2[JustComplex, npc.complexfloating]
ToJustComplexStrict2D: TypeAlias = _ToStrict2D2[JustComplex, npc.complexfloating]
ToJustComplexStrict3D: TypeAlias = _ToStrict3D2[JustComplex, npc.complexfloating]
