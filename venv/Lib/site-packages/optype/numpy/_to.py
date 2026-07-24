import sys
from collections.abc import Callable, Sequence as Seq
from typing import (
    TYPE_CHECKING,
    Annotated as Ann,
    Any,
    Literal,
    Protocol,
    TypeAliasType,
    runtime_checkable,
)

if sys.version_info >= (3, 13):
    from typing import TypeVar
else:
    from typing_extensions import TypeVar

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


type _PyBool = bool | Literal[0, 1]  # 0 and 1 are sometimes used as bool values
type _PyScalar = complex | bytes | str  # `complex` equivs `complex | float | int`


T = TypeVar("T", default=_PyScalar)
SCT = TypeVar("SCT", bound=np.generic, default=Any)
SCT_co = TypeVar("SCT_co", bound=np.generic, covariant=True)


# unlike `optype.numpy.CanArray0D` and `CanArrayND`, these one also accepts scalar types
# (and aren't runtime checkable)


@runtime_checkable
class _CanArrayND(Protocol[SCT_co]):
    def __len__(self, /) -> int: ...
    def __array__(self, /) -> np.ndarray[Any, np.dtype[SCT_co]]: ...


@runtime_checkable
class _CanArray(Protocol[SCT_co]):
    def __array__(self, /) -> np.ndarray[Any, np.dtype[SCT_co]]: ...


if TYPE_CHECKING:
    _CanArrayStrict1D = TypeAliasType(  # noqa: UP040
        "_CanArrayStrict1D",
        nptc.CanArray[tuple[int], np.dtype[SCT]],
        type_params=(SCT,),
    )
    _CanArrayStrict2D = TypeAliasType(  # noqa: UP040
        "_CanArrayStrict2D",
        nptc.CanArray[tuple[int, int], np.dtype[SCT]],
        type_params=(SCT,),
    )
    _CanArrayStrict3D = TypeAliasType(  # noqa: UP040
        "_CanArrayStrict3D",
        nptc.CanArray[tuple[int, int, int], np.dtype[SCT]],
        type_params=(SCT,),
    )

else:

    @runtime_checkable
    class _CanArrayStrict1D(Protocol[SCT_co]):
        def __array__(self) -> np.ndarray[tuple[int], np.dtype[SCT_co]]: ...

    @runtime_checkable
    class _CanArrayStrict2D(Protocol[SCT_co]):
        def __array__(self) -> np.ndarray[tuple[int, int], np.dtype[SCT_co]]: ...

    @runtime_checkable
    class _CanArrayStrict3D(Protocol[SCT_co]):
        def __array__(self) -> np.ndarray[tuple[int, int, int], np.dtype[SCT_co]]: ...


type _To1D1[SCT: np.generic] = _CanArrayND[SCT] | Seq[SCT]
type _To1D2[T, SCT: np.generic] = _CanArrayND[SCT] | Seq[T | SCT]

type _To2D1[SCT: np.generic] = _CanArrayND[SCT] | Seq[_To1D1[SCT]]
type _To2D2[T, SCT: np.generic] = _CanArrayND[SCT] | Seq[_To1D2[T, SCT]]

type _To3D1[SCT: np.generic] = _CanArrayND[SCT] | Seq[_To2D1[SCT]]
type _To3D2[T, SCT: np.generic] = _CanArrayND[SCT] | Seq[_To2D2[T, SCT]]

type _ToND1[SCT: np.generic] = _CanArrayND[SCT] | SeqND[_CanArray[SCT]]
type _ToND2[T, SCT: np.generic] = _CanArrayND[SCT] | SeqND[T | _CanArray[SCT]]

type _ToStrict1D1[SCT: np.generic] = _CanArrayStrict1D[SCT] | Seq[SCT]
type _ToStrict1D2[T, SCT: np.generic] = _CanArrayStrict1D[SCT] | Seq[T | SCT]

type _ToStrict2D1[SCT: np.generic] = _CanArrayStrict2D[SCT] | Seq[_ToStrict1D1[SCT]]
type _ToStrict2D2[T, SCT: np.generic] = (
    _CanArrayStrict2D[SCT] | Seq[_ToStrict1D2[T, SCT]]
)

type _ToStrict3D1[SCT: np.generic] = _CanArrayStrict3D[SCT] | Seq[_ToStrict2D1[SCT]]
type _ToStrict3D2[T, SCT: np.generic] = (
    _CanArrayStrict3D[SCT] | Seq[_ToStrict2D2[T, SCT]]
)


###

# TODO(jorenham): export & document
# https://github.com/jorenham/optype/issues/373

type integer_co = npc.integer | np.bool  # noqa: PYI042
type floating_co = npc.floating | npc.integer | np.bool  # noqa: PYI042
type complexfloating_co = npc.number | np.bool  # noqa: PYI042

# promotion rules with safe casting mode
type f16_co = npc.floating16 | npc.integer8 | np.bool  # noqa: PYI042
type f32_co = npc.floating32 | npc.floating16 | npc.integer16 | npc.integer8 | np.bool  # noqa: PYI042
type c64_co = npc.inexact32 | npc.number16 | npc.integer8 | np.bool  # noqa: PYI042
type f64_co = npc.floating64 | npc.floating32 | npc.floating16 | npc.integer | np.bool  # noqa: PYI042
type c128_co = npc.number64 | npc.number32 | npc.number16 | npc.integer | np.bool  # noqa: PYI042

###

# runtime type-checking compat

if TYPE_CHECKING:
    # we don't use float | int to avoid clutter in type-checker error output
    py_float = float
    py_complex = complex
else:
    py_float = float | int
    py_complex = complex | py_float

###

# scalar- and array-likes, with "coercible" shape-types

type ToScalar = _PyScalar | np.generic
ToArray1D = TypeAliasType("ToArray1D", _To1D2[T, SCT], type_params=(T, SCT))  # noqa: UP040
ToArray2D = TypeAliasType("ToArray2D", _To2D2[T, SCT], type_params=(T, SCT))  # noqa: UP040
ToArray3D = TypeAliasType("ToArray3D", _To3D2[T, SCT], type_params=(T, SCT))  # noqa: UP040
ToArrayND = TypeAliasType("ToArrayND", _ToND2[T, SCT], type_params=(T, SCT))  # noqa: UP040

type ToFalse = nptc.LiteralFalse | Literal[0]
type ToTrue = nptc.LiteralTrue | Literal[1]

type ToJustFalse = nptc.LiteralFalse
type ToJustTrue = nptc.LiteralTrue

type ToBool = _PyBool | np.bool
type ToBool1D = _To1D2[_PyBool, np.bool]
type ToBool2D = _To2D2[_PyBool, np.bool]
type ToBool3D = _To3D2[_PyBool, np.bool]
type ToBoolND = _ToND2[_PyBool, np.bool]

type ToInt = int | integer_co
type ToInt1D = _To1D2[int, integer_co]
type ToInt2D = _To2D2[int, integer_co]
type ToInt3D = _To3D2[int, integer_co]
type ToIntND = _ToND2[int, integer_co]

type ToFloat16 = f16_co
type ToFloat16_1D = _To1D1[f16_co]
type ToFloat16_2D = _To2D1[f16_co]
type ToFloat16_3D = _To3D1[f16_co]
type ToFloat16_ND = _ToND1[f16_co]

type ToFloat32 = f32_co
type ToFloat32_1D = _To1D1[f32_co]
type ToFloat32_2D = _To2D1[f32_co]
type ToFloat32_3D = _To3D1[f32_co]
type ToFloat32_ND = _ToND1[f32_co]

type ToFloat64 = py_float | f64_co
type ToFloat64_1D = _To1D2[py_float, f64_co]
type ToFloat64_2D = _To2D2[py_float, f64_co]
type ToFloat64_3D = _To3D2[py_float, f64_co]
type ToFloat64_ND = _ToND2[py_float, f64_co]

type ToFloat = py_float | floating_co
type ToFloat1D = _To1D2[py_float, floating_co]
type ToFloat2D = _To2D2[py_float, floating_co]
type ToFloat3D = _To3D2[py_float, floating_co]
type ToFloatND = _ToND2[py_float, floating_co]

type ToComplex64 = c64_co
type ToComplex64_1D = _To1D1[c64_co]
type ToComplex64_2D = _To2D1[c64_co]
type ToComplex64_3D = _To3D1[c64_co]
type ToComplex64_ND = _ToND1[c64_co]

type ToComplex128 = py_complex | c128_co
type ToComplex128_1D = _To1D2[py_complex, c128_co]
type ToComplex128_2D = _To2D2[py_complex, c128_co]
type ToComplex128_3D = _To3D2[py_complex, c128_co]
type ToComplex128_ND = _ToND2[py_complex, c128_co]

type ToComplex = py_complex | complexfloating_co
type ToComplex1D = _To1D2[py_complex, complexfloating_co]
type ToComplex2D = _To2D2[py_complex, complexfloating_co]
type ToComplex3D = _To3D2[py_complex, complexfloating_co]
type ToComplexND = _ToND2[py_complex, complexfloating_co]

# scalar- and array-likes, with "just" that scalar type

type ToJustBool = bool | np.bool
type ToJustBool1D = _To1D2[bool, np.bool]
type ToJustBool2D = _To2D2[bool, np.bool]
type ToJustBool3D = _To3D2[bool, np.bool]
type ToJustBoolND = _ToND2[bool, np.bool]

type ToJustInt64 = JustInt | np.intp
type ToJustInt64_1D = _To1D2[JustInt, np.intp]
type ToJustInt64_2D = _To2D2[JustInt, np.intp]
type ToJustInt64_3D = _To3D2[JustInt, np.intp]
type ToJustInt64_ND = _ToND2[JustInt, np.intp]

type ToJustInt = JustInt | npc.integer
type ToJustInt1D = _To1D2[JustInt, npc.integer]
type ToJustInt2D = _To2D2[JustInt, npc.integer]
type ToJustInt3D = _To3D2[JustInt, npc.integer]
type ToJustIntND = _ToND2[JustInt, npc.integer]

type ToJustFloat16 = npc.floating16
type ToJustFloat16_1D = _To1D1[npc.floating16]
type ToJustFloat16_2D = _To2D1[npc.floating16]
type ToJustFloat16_3D = _To3D1[npc.floating16]
type ToJustFloat16_ND = _ToND1[npc.floating16]

type ToJustFloat32 = npc.floating32
type ToJustFloat32_1D = _To1D1[npc.floating32]
type ToJustFloat32_2D = _To2D1[npc.floating32]
type ToJustFloat32_3D = _To3D1[npc.floating32]
type ToJustFloat32_ND = _ToND1[npc.floating32]

type ToJustFloat64 = JustFloat | npc.floating64
type ToJustFloat64_1D = _To1D2[JustFloat, npc.floating64]
type ToJustFloat64_2D = _To2D2[JustFloat, npc.floating64]
type ToJustFloat64_3D = _To3D2[JustFloat, npc.floating64]
type ToJustFloat64_ND = _ToND2[JustFloat, npc.floating64]

type ToJustLongDouble = npc.floating80
type ToJustLongDouble1D = _To1D1[npc.floating80]
type ToJustLongDouble2D = _To2D1[npc.floating80]
type ToJustLongDouble3D = _To3D1[npc.floating80]
type ToJustLongDoubleND = _ToND1[npc.floating80]

type ToJustFloat = JustFloat | npc.floating
type ToJustFloat1D = _To1D2[JustFloat, npc.floating]
type ToJustFloat2D = _To2D2[JustFloat, npc.floating]
type ToJustFloat3D = _To3D2[JustFloat, npc.floating]
type ToJustFloatND = _ToND2[JustFloat, npc.floating]

type ToJustComplex64 = npc.complexfloating64
type ToJustComplex64_1D = _To1D1[npc.complexfloating64]
type ToJustComplex64_2D = _To2D1[npc.complexfloating64]
type ToJustComplex64_3D = _To3D1[npc.complexfloating64]
type ToJustComplex64_ND = _ToND1[npc.complexfloating64]

type ToJustComplex128 = JustComplex | npc.complexfloating128
type ToJustComplex128_1D = _To1D2[JustComplex, npc.complexfloating128]
type ToJustComplex128_2D = _To2D2[JustComplex, npc.complexfloating128]
type ToJustComplex128_3D = _To3D2[JustComplex, npc.complexfloating128]
type ToJustComplex128_ND = _ToND2[JustComplex, npc.complexfloating128]

type ToJustCLongDouble = npc.complexfloating160
type ToJustCLongDouble1D = _To1D1[npc.complexfloating160]
type ToJustCLongDouble2D = _To2D1[npc.complexfloating160]
type ToJustCLongDouble3D = _To3D1[npc.complexfloating160]
type ToJustCLongDoubleND = _ToND1[npc.complexfloating160]

type ToJustComplex = JustComplex | npc.complexfloating
type ToJustComplex1D = _To1D2[JustComplex, npc.complexfloating]
type ToJustComplex2D = _To2D2[JustComplex, npc.complexfloating]
type ToJustComplex3D = _To3D2[JustComplex, npc.complexfloating]
type ToJustComplexND = _ToND2[JustComplex, npc.complexfloating]

# array-likes, with "coercible" shape-types, and "strict" shape-types

ToArrayStrict1D = TypeAliasType(  # noqa: UP040
    "ToArrayStrict1D",
    _ToStrict1D2[T, SCT],
    type_params=(T, SCT),
)
ToArrayStrict2D = TypeAliasType(  # noqa: UP040
    "ToArrayStrict2D",
    _ToStrict2D2[T, SCT],
    type_params=(T, SCT),
)
ToArrayStrict3D = TypeAliasType(  # noqa: UP040
    "ToArrayStrict3D",
    _ToStrict3D2[T, SCT],
    type_params=(T, SCT),
)

type ToBoolStrict1D = _ToStrict1D2[_PyBool, np.bool]
type ToBoolStrict2D = _ToStrict2D2[_PyBool, np.bool]
type ToBoolStrict3D = _ToStrict3D2[_PyBool, np.bool]

type ToIntStrict1D = _ToStrict1D2[int, integer_co]
type ToIntStrict2D = _ToStrict2D2[int, integer_co]
type ToIntStrict3D = _ToStrict3D2[int, integer_co]

type ToFloat16Strict1D = _ToStrict1D1[f16_co]
type ToFloat16Strict2D = _ToStrict2D1[f16_co]
type ToFloat16Strict3D = _ToStrict3D1[f16_co]

type ToFloat32Strict1D = _ToStrict1D1[f32_co]
type ToFloat32Strict2D = _ToStrict2D1[f32_co]
type ToFloat32Strict3D = _ToStrict3D1[f32_co]

type ToFloat64Strict1D = _ToStrict1D2[py_float, f64_co]
type ToFloat64Strict2D = _ToStrict2D2[py_float, f64_co]
type ToFloat64Strict3D = _ToStrict3D2[py_float, f64_co]

type ToFloatStrict1D = _ToStrict1D2[py_float, floating_co]
type ToFloatStrict2D = _ToStrict2D2[py_float, floating_co]
type ToFloatStrict3D = _ToStrict3D2[py_float, floating_co]

type ToComplex64Strict1D = _ToStrict1D1[c64_co]
type ToComplex64Strict2D = _ToStrict2D1[c64_co]
type ToComplex64Strict3D = _ToStrict3D1[c64_co]

type ToComplex128Strict1D = _ToStrict1D2[py_complex, c128_co]
type ToComplex128Strict2D = _ToStrict2D2[py_complex, c128_co]
type ToComplex128Strict3D = _ToStrict3D2[py_complex, c128_co]

type ToComplexStrict1D = _ToStrict1D2[py_complex, complexfloating_co]
type ToComplexStrict2D = _ToStrict2D2[py_complex, complexfloating_co]
type ToComplexStrict3D = _ToStrict3D2[py_complex, complexfloating_co]

# array-likes, with "just" that scalar type, and "strict" shape-types

type ToJustBoolStrict1D = _ToStrict1D2[bool, np.bool]
type ToJustBoolStrict2D = _ToStrict2D2[bool, np.bool]
type ToJustBoolStrict3D = _ToStrict3D2[bool, np.bool]
type ToJustInt64Strict1D = _ToStrict1D2[JustInt, np.int64]
type ToJustInt64Strict2D = _ToStrict2D2[JustInt, np.int64]
type ToJustInt64Strict3D = _ToStrict3D2[JustInt, np.int64]
type ToJustIntStrict1D = _ToStrict1D2[JustInt, npc.integer]
type ToJustIntStrict2D = _ToStrict2D2[JustInt, npc.integer]
type ToJustIntStrict3D = _ToStrict3D2[JustInt, npc.integer]
type ToJustFloat16Strict1D = _ToStrict1D1[npc.floating16]
type ToJustFloat16Strict2D = _ToStrict2D1[npc.floating16]
type ToJustFloat16Strict3D = _ToStrict3D1[npc.floating16]
type ToJustFloat32Strict1D = _ToStrict1D1[npc.floating32]
type ToJustFloat32Strict2D = _ToStrict2D1[npc.floating32]
type ToJustFloat32Strict3D = _ToStrict3D1[npc.floating32]
type ToJustFloat64Strict1D = _ToStrict1D2[JustFloat, npc.floating64]
type ToJustFloat64Strict2D = _ToStrict2D2[JustFloat, npc.floating64]
type ToJustFloat64Strict3D = _ToStrict3D2[JustFloat, npc.floating64]
type ToJustLongDoubleStrict1D = _ToStrict1D1[npc.floating80]
type ToJustLongDoubleStrict2D = _ToStrict2D1[npc.floating80]
type ToJustLongDoubleStrict3D = _ToStrict3D1[npc.floating80]
type ToJustFloatStrict1D = _ToStrict1D2[JustFloat, npc.floating]
type ToJustFloatStrict2D = _ToStrict2D2[JustFloat, npc.floating]
type ToJustFloatStrict3D = _ToStrict3D2[JustFloat, npc.floating]
type ToJustComplex64Strict1D = _ToStrict1D1[npc.complexfloating64]
type ToJustComplex64Strict2D = _ToStrict2D1[npc.complexfloating64]
type ToJustComplex64Strict3D = _ToStrict3D1[npc.complexfloating64]
type ToJustComplex128Strict1D = _ToStrict1D2[JustComplex, npc.complexfloating128]
type ToJustComplex128Strict2D = _ToStrict2D2[JustComplex, npc.complexfloating128]
type ToJustComplex128Strict3D = _ToStrict3D2[JustComplex, npc.complexfloating128]
type ToJustCLongDoubleStrict1D = _ToStrict1D1[npc.complexfloating160]
type ToJustCLongDoubleStrict2D = _ToStrict2D1[npc.complexfloating160]
type ToJustCLongDoubleStrict3D = _ToStrict3D1[npc.complexfloating160]
type ToJustComplexStrict1D = _ToStrict1D2[JustComplex, npc.complexfloating]
type ToJustComplexStrict2D = _ToStrict2D2[JustComplex, npc.complexfloating]
type ToJustComplexStrict3D = _ToStrict3D2[JustComplex, npc.complexfloating]


# if not TYPE_CHECKING:
try:
    from beartype.vale import Is as _BeartypeIs
except ImportError:
    pass
else:
    # `beartype.vale._core._valecore.BeartypeValidator` is private, unfortunately
    type _BeartypeValidator = Any

    def _install() -> None:
        def build(
            check: Callable[[np.dtype[np.generic], Any], bool],
        ) -> Callable[[int | None, object], _BeartypeValidator]:
            def factory(nd: int | None, target: object) -> _BeartypeValidator:
                def predicate(x: object, /) -> bool:
                    try:
                        x = np.asanyarray(x)
                    except (ValueError, TypeError):
                        return False
                    return (nd is None or x.ndim == nd) and check(x.dtype, target)

                return _BeartypeIs[predicate]

            return factory

        is_cast = build(np.can_cast)
        is_sub = build(lambda dt, t: issubclass(dt.type, t))

        def wrap(key: str, ann: _BeartypeValidator) -> None:
            if key in globals():
                globals()[key] = Ann[globals()[key], ann]  # ty:ignore[invalid-type-form]

        for name, sct in [
            ("Bool", np.bool),
            ("Int64", np.int64),
            ("Float16", np.float16),
            ("Float32", np.float32),
            ("Float64", np.float64),
            ("LongDouble", np.longdouble),
            ("Complex64", np.complex64),
            ("Complex128", np.complex128),
            ("CLongDouble", np.clongdouble),
        ]:
            name_ = f"{name}_" if name[-1].isdigit() else name

            wrap(f"To{name}", is_cast(0, sct))
            wrap(f"To{name_}1D", is_cast(1, sct))
            wrap(f"To{name_}2D", is_cast(2, sct))
            wrap(f"To{name_}3D", is_cast(3, sct))
            wrap(f"To{name_}ND", is_cast(None, sct))

            wrap(f"ToJust{name}", is_sub(0, sct))
            wrap(f"ToJust{name_}1D", is_sub(1, sct))
            wrap(f"ToJust{name_}2D", is_sub(2, sct))
            wrap(f"ToJust{name_}3D", is_sub(3, sct))
            wrap(f"ToJust{name_}ND", is_sub(None, sct))

            wrap(f"To{name}Strict1D", is_cast(1, sct))
            wrap(f"To{name}Strict2D", is_cast(2, sct))
            wrap(f"To{name}Strict3D", is_cast(3, sct))

            wrap(f"ToJust{name}Strict1D", is_sub(1, sct))
            wrap(f"ToJust{name}Strict2D", is_sub(2, sct))
            wrap(f"ToJust{name}Strict3D", is_sub(3, sct))

        for name, sct, co_sct in [
            ("Int", np.integer, (np.bool,)),
            ("Float", np.floating, (np.integer, np.bool)),
            ("Complex", np.complexfloating, (np.integer, np.floating, np.bool)),
        ]:
            co_classes = (sct, *co_sct)

            wrap(f"To{name}", is_sub(0, co_classes))
            wrap(f"To{name}1D", is_sub(1, co_classes))
            wrap(f"To{name}2D", is_sub(2, co_classes))
            wrap(f"To{name}3D", is_sub(3, co_classes))
            wrap(f"To{name}ND", is_sub(None, co_classes))

            wrap(f"ToJust{name}", is_sub(0, sct))
            wrap(f"ToJust{name}1D", is_sub(1, sct))
            wrap(f"ToJust{name}2D", is_sub(2, sct))
            wrap(f"ToJust{name}3D", is_sub(3, sct))
            wrap(f"ToJust{name}ND", is_sub(None, sct))

            wrap(f"To{name}Strict1D", is_sub(1, co_classes))
            wrap(f"To{name}Strict2D", is_sub(2, co_classes))
            wrap(f"To{name}Strict3D", is_sub(3, co_classes))

            wrap(f"ToJust{name}Strict1D", is_sub(1, sct))
            wrap(f"ToJust{name}Strict2D", is_sub(2, sct))
            wrap(f"ToJust{name}Strict3D", is_sub(3, sct))

    _install()
    del _install
