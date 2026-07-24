import abc
import types
from typing import Any, ClassVar, Final, Generic, Literal, Self, final, overload, override, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._lti_conversion import _DiscretizeMethod

__all__ = [
    "StateSpace",
    "TransferFunction",
    "ZerosPolesGain",
    "bode",
    "dbode",
    "dfreqresp",
    "dimpulse",
    "dlsim",
    "dlti",
    "dstep",
    "freqresp",
    "impulse",
    "lsim",
    "lti",
    "place_poles",
    "step",
]

###

_ZerosT_co = TypeVar("_ZerosT_co", bound=npc.inexact32 | npc.inexact64, default=Any, covariant=True)
_PolesT_co = TypeVar("_PolesT_co", bound=_Float, default=np.float64 | Any, covariant=True)
_DTT_co = TypeVar("_DTT_co", bound=onp.ToComplex | None, default=Any, covariant=True)
_RequestedPolesT_co = TypeVar("_RequestedPolesT_co", bound=npc.inexact64, default=np.float64 | np.complex128, covariant=True)

type _Tuple3[T] = tuple[T, T, T]
type _Tuple4[T] = tuple[T, T, T, T]

type _Float = np.float32 | np.float64
type _Complex = np.complex64 | np.complex128
type _Inexact = _Float | _Complex
type _Number = npc.integer | _Inexact

type _ToNumber = complex | _Number

type _Array12D[ScalarT: np.generic] = onp.ArrayND[ScalarT, tuple[int] | tuple[int, int]]

type _Float1D = onp.Array1D[_Float]
type _Float64_1D = onp.Array1D[np.float64]
type _Float64_2D = onp.Array2D[np.float64]
type _Complex1D = onp.Array1D[_Complex]

type _ToFloat12D = onp.ToFloat1D | onp.ToFloat2D
type _ToFloat012D = onp.ToFloat | _ToFloat12D
type _ToComplex12D = onp.ToComplex1D | onp.ToComplex2D
type _ToComplex012D = onp.ToComplex | _ToComplex12D

type _ToInexact32_1D = onp.ToJustFloat32_1D | onp.ToJustComplex64_1D
type _ToInexact32_2D = onp.ToJustFloat32_2D | onp.ToJustComplex64_2D
type _ToInexact64_1D = onp.ToArray1D[complex, npc.inexact64 | npc.integer]
type _ToInexact64_2D = onp.ToArray2D[complex, npc.inexact64 | npc.integer]

###

# numerator, denominator
type _ToTFContFloat = tuple[_ToFloat12D, onp.ToComplex1D]
type _ToTFContInexact = tuple[_ToComplex12D, onp.ToComplex1D]
type _ToTFContInexact32 = tuple[_ToInexact32_1D | _ToInexact32_2D, _ToInexact32_1D]
type _ToTFContInexact64 = tuple[_ToInexact64_1D | _ToInexact64_2D, _ToInexact64_1D]
# numerator, denominator, dt
type _ToTFDisc = tuple[_ToFloat12D, onp.ToComplex1D, onp.ToFloat]

# zeros, poles, gain
type _ToZPKContFloat = tuple[onp.ToFloat1D, onp.ToFloat1D, onp.ToFloat]
type _ToZPKContInexact = tuple[onp.ToComplex1D, onp.ToComplex1D, onp.ToFloat]
type _ToZPKContInexact32 = tuple[_ToInexact32_1D, _ToInexact32_1D, onp.ToFloat]
type _ToZPKContInexact64 = tuple[_ToInexact64_1D, _ToInexact64_1D, onp.ToFloat]
# zeros, poles, gain, dt
type _ToZPKDisc = tuple[onp.ToFloat1D, onp.ToFloat1D, onp.ToFloat, onp.ToFloat]

# A, B, C, D
type _ToSSContFloat = _Tuple4[onp.ToFloat2D]
type _ToSSContInexact = _Tuple4[onp.ToComplex2D]
type _ToSSContInexact32 = _Tuple4[_ToInexact32_2D]
type _ToSSContInexact64 = _Tuple4[_ToInexact64_2D]
# A, B, C, D, dt
type _ToSSDisc = tuple[onp.ToFloat2D, onp.ToFloat2D, onp.ToFloat2D, onp.ToFloat2D, onp.ToFloat]

type _ToLTIFloat = _ToTFContFloat | _ToZPKContFloat | _ToSSContFloat
type _ToLTIInexact = _ToTFContInexact | _ToZPKContInexact | _ToSSContInexact
type _ToLTIInexact32 = _ToTFContInexact32 | _ToZPKContInexact32 | _ToSSContInexact32
type _ToLTIInexact64 = _ToTFContInexact64 | _ToZPKContInexact64 | _ToSSContInexact64
type _ToDLTI = _ToTFDisc | _ToZPKDisc | _ToSSDisc

###

class LinearTimeInvariant(Generic[_ZerosT_co, _PolesT_co, _DTT_co]):
    inputs: Final[int]
    outputs: Final[int]

    @classmethod
    def __class_getitem__(cls, args: object | tuple[object, ...], /) -> types.GenericAlias: ...

    #
    def __new__(cls, /, *system: *tuple[()], dt: None = None) -> Self: ...

    #
    @abc.abstractmethod
    @type_check_only
    def to_tf(self, /) -> TransferFunction[_PolesT_co, _DTT_co]: ...
    @abc.abstractmethod
    @type_check_only
    def to_zpk(self, /) -> ZerosPolesGain[_ZerosT_co, _PolesT_co, _DTT_co]: ...
    @abc.abstractmethod
    @type_check_only
    def to_ss(self, /) -> StateSpace[_ZerosT_co, _PolesT_co, _DTT_co]: ...

    #
    @property
    def zeros(self, /) -> _Array12D[_ZerosT_co]: ...
    @property
    def poles(self, /) -> onp.Array1D[_PolesT_co]: ...
    @property
    def dt(self, /) -> _DTT_co: ...

class lti(LinearTimeInvariant[_ZerosT_co, _PolesT_co, None], Generic[_ZerosT_co, _PolesT_co], metaclass=abc.ABCMeta):
    @override
    @overload
    def __new__(cls, num: _ToFloat12D, den: onp.ToFloat1D, /) -> TransferFunctionContinuous[_Float]: ...  # pyrefly:ignore[bad-override]
    @overload
    def __new__(cls, zeros: onp.ToFloat1D, poles: onp.ToFloat1D, gain: onp.ToFloat, /) -> ZerosPolesGainContinuous[_Float]: ...
    @overload
    def __new__(cls, zeros: onp.ToComplex1D, poles: onp.ToComplex1D, gain: onp.ToFloat, /) -> ZerosPolesGainContinuous: ...
    @overload
    def __new__(cls, A: _ToFloat012D, B: _ToFloat012D, C: _ToFloat012D, D: _ToFloat012D, /) -> StateSpaceContinuous[_Float]: ...
    @overload
    def __new__(cls, A: _ToComplex012D, B: _ToComplex012D, C: _ToComplex012D, D: _ToComplex012D, /) -> StateSpaceContinuous: ...

    #
    def __init__(self, /, *system: *tuple[()]) -> None: ...

    #
    @overload
    def impulse(
        self: lti[np.float32 | np.float64],
        /,
        X0: onp.ToFloat1D | None = None,
        T: onp.ToFloat1D | None = None,
        N: int | None = None,
    ) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.float64]]: ...
    @overload
    def impulse(
        self: lti[np.complex64 | np.complex128],
        /,
        X0: onp.ToComplex1D | None = None,
        T: onp.ToFloat1D | None = None,
        N: int | None = None,
    ) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.complex128]]: ...

    #
    @overload
    def step(
        self: lti[np.float32 | np.float64],
        /,
        X0: onp.ToFloat1D | None = None,
        T: onp.ToFloat1D | None = None,
        N: int | None = None,
    ) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.float64]]: ...
    @overload
    def step(
        self: lti[np.complex64 | np.complex128],
        /,
        X0: onp.ToComplex1D | None = None,
        T: onp.ToFloat1D | None = None,
        N: int | None = None,
    ) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.complex128]]: ...

    #
    @overload
    def output(
        self: lti[np.float32 | np.float64],
        /,
        U: _ToFloat012D | None,
        T: onp.ToFloat | onp.ToFloat1D,
        X0: onp.ToComplex1D | None = None,
    ) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
    @overload
    def output(
        self: lti[np.complex64 | np.complex128],
        /,
        U: _ToFloat012D | None,
        T: onp.ToFloat | onp.ToFloat1D,
        X0: onp.ToComplex1D | None = None,
    ) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.complex128], onp.ArrayND[np.complex128]]: ...

    #
    @overload
    def bode(self: lti[npc.inexact64], /, w: onp.ToFloat1D | None = None, n: int = 100) -> _Tuple3[onp.Array1D[np.float64]]: ...
    @overload
    def bode(self: lti[npc.inexact32], /, w: onp.ToFloat1D | None = None, n: int = 100) -> _Tuple3[onp.Array1D[np.float32]]: ...

    #
    @overload
    def freqresp(
        self: lti[npc.inexact64], /, w: onp.ToFloat1D | None = None, n: int = 10_000
    ) -> tuple[onp.Array1D[np.float64], onp.Array1D[np.complex128]]: ...
    @overload
    def freqresp(
        self: lti[npc.inexact32], /, w: onp.ToFloat1D | None = None, n: int = 10_000
    ) -> tuple[onp.Array1D[np.float32], onp.Array1D[np.complex64]]: ...

    #
    @abc.abstractmethod
    def to_discrete[DTT: onp.ToComplex | None](
        self, /, dt: DTT, method: _DiscretizeMethod = "zoh", alpha: float | None = None
    ) -> dlti[_ZerosT_co, _PolesT_co, DTT]: ...

#
class dlti(LinearTimeInvariant[_ZerosT_co, _PolesT_co, _DTT_co], Generic[_ZerosT_co, _PolesT_co, _DTT_co], metaclass=abc.ABCMeta):
    @override
    @overload
    def __new__(  # pyrefly:ignore[bad-override]
        cls, num: _ToFloat12D, den: onp.ToFloat1D, /, *, dt: _DTT_co = ...
    ) -> TransferFunctionDiscrete[_Float, _DTT_co]: ...
    @overload
    def __new__(
        cls, zeros: onp.ToFloat1D, poles: onp.ToFloat1D, gain: onp.ToFloat, /, *, dt: _DTT_co = ...
    ) -> ZerosPolesGainDiscrete[_Float, _Float, _DTT_co]: ...
    @overload
    def __new__(
        cls, zeros: onp.ToComplex1D, poles: onp.ToComplex1D, gain: onp.ToFloat, /, *, dt: _DTT_co = ...
    ) -> ZerosPolesGainDiscrete[Any, _Float, _DTT_co]: ...
    @overload
    def __new__(
        cls, A: _ToFloat012D, B: _ToFloat012D, C: _ToFloat012D, D: _ToFloat012D, /, *, dt: _DTT_co = ...
    ) -> StateSpaceDiscrete[_Float, _Float, _DTT_co]: ...
    @overload
    def __new__(
        cls, A: _ToComplex012D, B: _ToComplex012D, C: _ToComplex012D, D: _ToComplex012D, /, *, dt: _DTT_co = ...
    ) -> StateSpaceDiscrete[Any, _Float, _DTT_co]: ...

    #
    def __init__(self, /, *system: *tuple[()], dt: _DTT_co) -> None: ...

    #
    def output(
        self, /, u: _ToFloat012D | None, t: onp.ToFloat1D, x0: onp.ToFloat1D | None = None
    ) -> tuple[_Float64_1D, _Float64_2D] | tuple[_Float64_1D, _Float64_2D, _Float64_2D]: ...
    def impulse(
        self, /, x0: onp.ToFloat1D | None = None, t: onp.ToFloat1D | None = None, n: int | None = None
    ) -> tuple[onp.Array1D[np.float64], tuple[onp.Array2D[np.float64], ...]]: ...
    #
    def step(
        self, /, x0: onp.ToFloat1D | None = None, t: onp.ToFloat1D | None = None, n: int | None = None
    ) -> tuple[onp.Array1D[np.float64], tuple[onp.Array2D[np.float64], ...]]: ...
    def bode(self, /, w: onp.ToFloat1D | None = None, n: int = 100) -> _Tuple3[onp.Array1D[np.float64]]: ...
    def freqresp(
        self, /, w: onp.ToFloat1D | None = None, n: int = 10_000, whole: bool = False
    ) -> tuple[onp.Array1D[np.float64], onp.Array1D[np.complex128]]: ...

#
class TransferFunction(LinearTimeInvariant[_PolesT_co, _PolesT_co, _DTT_co], Generic[_PolesT_co, _DTT_co], metaclass=abc.ABCMeta):
    @override
    @overload
    def __new__[PolesT: _Float](
        cls, system: lti[PolesT, PolesT], /, *, dt: None = None
    ) -> TransferFunctionContinuous[PolesT]: ...
    @overload
    def __new__(
        cls,
        num: onp.ToArray1D[float, npc.integer | npc.floating64] | onp.ToArray2D[float, npc.integer | npc.floating64],
        den: onp.ToFloat1D,
        /,
        *,
        dt: None = None,
    ) -> TransferFunctionContinuous[np.float64]: ...
    @overload
    def __new__(
        cls, num: onp.ToFloat1D | onp.ToFloat2D, den: onp.ToArray1D[float, npc.integer | npc.floating64], /, *, dt: None = None
    ) -> TransferFunctionContinuous[np.float64]: ...
    @overload
    def __new__(
        cls, num: onp.ToJustFloat32_1D | onp.ToJustFloat32_2D, den: onp.ToFloat32_1D, /, *, dt: None = None
    ) -> TransferFunctionContinuous[np.float32]: ...
    @overload
    def __new__(
        cls, num: onp.ToFloat32_1D | onp.ToFloat32_2D, den: onp.ToJustFloat32_1D, /, *, dt: None = None
    ) -> TransferFunctionContinuous[np.float32]: ...
    @overload
    def __new__(cls, num: _ToFloat12D, den: onp.ToFloat1D, /, *, dt: None = None) -> TransferFunctionContinuous[_Float]: ...
    @overload
    def __new__[PolesT: _Float, DTT: onp.ToComplex | None](
        cls, system: dlti[PolesT, PolesT, DTT], /, *, dt: None = None
    ) -> TransferFunctionDiscrete[PolesT, DTT]: ...
    @overload
    def __new__[DTT: onp.ToComplex](
        cls,
        num: onp.ToArray1D[float, npc.integer | npc.floating64] | onp.ToArray2D[float, npc.integer | npc.floating64],
        den: onp.ToFloat1D,
        /,
        *,
        dt: DTT,
    ) -> TransferFunctionDiscrete[np.float64, DTT]: ...
    @overload
    def __new__[DTT: onp.ToComplex](
        cls, num: onp.ToFloat1D | onp.ToFloat2D, den: onp.ToArray1D[float, npc.integer | npc.floating64], /, *, dt: DTT
    ) -> TransferFunctionDiscrete[np.float64, DTT]: ...
    @overload
    def __new__[DTT: onp.ToComplex](
        cls, num: onp.ToJustFloat32_1D | onp.ToJustFloat32_2D, den: onp.ToFloat32_1D, /, *, dt: DTT
    ) -> TransferFunctionDiscrete[np.float32, DTT]: ...
    @overload
    def __new__[DTT: onp.ToComplex](
        cls, num: onp.ToFloat32_1D | onp.ToFloat32_2D, den: onp.ToJustFloat32_1D, /, *, dt: DTT
    ) -> TransferFunctionDiscrete[np.float32, DTT]: ...
    @overload
    def __new__[DTT: onp.ToComplex](
        cls, num: _ToFloat12D, den: onp.ToFloat1D, /, *, dt: DTT
    ) -> TransferFunctionDiscrete[_Float, DTT]: ...

    #
    @overload
    def __init__(self, system: LinearTimeInvariant[_PolesT_co, _PolesT_co, _DTT_co], /, *, dt: None = None) -> None: ...
    @overload
    def __init__(self, num: _ToFloat12D, den: onp.ToFloat1D, /, *, dt: _DTT_co = ...) -> None: ...

    #
    @property
    def num(self, /) -> _Array12D[_PolesT_co]: ...
    @num.setter
    def num(self, num: _ToFloat12D, /) -> None: ...

    #
    @property
    def den(self, /) -> onp.Array1D[_PolesT_co]: ...
    @den.setter
    def den(self, den: onp.ToFloat1D, /) -> None: ...

    #
    @override
    def to_tf(self, /) -> Self: ...
    @override
    def to_zpk(self, /) -> ZerosPolesGain[_PolesT_co, _PolesT_co, _DTT_co]: ...
    @override
    def to_ss(self, /) -> StateSpace[_PolesT_co, _PolesT_co, _DTT_co]: ...

@final
class TransferFunctionContinuous(TransferFunction[_PolesT_co, None], lti[_PolesT_co, _PolesT_co], Generic[_PolesT_co]):
    @override
    def to_tf(self, /) -> Self: ...
    @override
    def to_discrete[DTT: onp.ToComplex | None](
        self, /, dt: DTT, method: _DiscretizeMethod = "zoh", alpha: float | None = None
    ) -> TransferFunctionDiscrete[_PolesT_co, DTT]: ...

@final
class TransferFunctionDiscrete(
    TransferFunction[_PolesT_co, _DTT_co], dlti[_PolesT_co, _PolesT_co, _DTT_co], Generic[_PolesT_co, _DTT_co]
):
    @overload
    def __init__(self, system: LinearTimeInvariant[_PolesT_co, _PolesT_co, _DTT_co], /, *, dt: None = None) -> None: ...
    @overload
    def __init__(self, numerator: _ToFloat12D, denominator: onp.ToFloat1D, /, *, dt: _DTT_co = ...) -> None: ...

    #
    @override
    def output(
        self, /, u: _ToFloat012D | None, t: onp.ToFloat1D, x0: onp.ToFloat1D | None = None
    ) -> tuple[_Float64_1D, _Float64_2D]: ...

#
class ZerosPolesGain(LinearTimeInvariant[_ZerosT_co, _PolesT_co, _DTT_co], Generic[_ZerosT_co, _PolesT_co, _DTT_co]):
    @override
    @overload
    def __new__(
        cls, system: lti[_ZerosT_co, _PolesT_co], /, *, dt: None = None
    ) -> ZerosPolesGainContinuous[_ZerosT_co, _PolesT_co]: ...
    @overload
    def __new__(
        cls, system: dlti[_ZerosT_co, _PolesT_co, _DTT_co], /, *, dt: None = None
    ) -> ZerosPolesGainDiscrete[_ZerosT_co, _PolesT_co, _DTT_co]: ...
    @overload
    def __new__(
        cls,
        zeros: onp.ToArray1D[float, npc.integer | npc.floating64] | onp.ToArray2D[float, npc.integer | npc.floating64],
        poles: onp.ToArray1D[float, npc.integer | npc.floating64],
        gain: onp.ToFloat,
        /,
        *,
        dt: None = None,
    ) -> ZerosPolesGainContinuous[np.float64, np.float64]: ...
    @overload
    def __new__[DTT: onp.ToComplex](
        cls,
        zeros: onp.ToArray1D[float, npc.integer | npc.floating64] | onp.ToArray2D[float, npc.integer | npc.floating64],
        poles: onp.ToArray1D[float, npc.integer | npc.floating64],
        gain: onp.ToFloat,
        /,
        *,
        dt: DTT,
    ) -> ZerosPolesGainDiscrete[np.float64, np.float64, DTT]: ...
    @overload
    def __new__(
        cls, zeros: _ToFloat12D, poles: onp.ToFloat1D, gain: onp.ToFloat, /, *, dt: None = None
    ) -> ZerosPolesGainContinuous[_Float, _Float]: ...
    @overload
    def __new__[DTT: onp.ToComplex](
        cls, zeros: _ToFloat12D, poles: onp.ToFloat1D, gain: onp.ToFloat, /, *, dt: DTT
    ) -> ZerosPolesGainDiscrete[_Float, _Float, DTT]: ...
    @overload
    def __new__(
        cls,
        zeros: onp.ToJustComplex128_1D | onp.ToJustComplex128_2D,
        poles: onp.ToArray1D[float, npc.integer | npc.floating64],
        gain: onp.ToFloat,
        /,
        *,
        dt: None = None,
    ) -> ZerosPolesGainContinuous[np.complex128, np.float64]: ...
    @overload
    def __new__[DTT: onp.ToComplex](
        cls,
        zeros: onp.ToJustComplex128_1D | onp.ToJustComplex128_2D,
        poles: onp.ToArray1D[float, npc.integer | npc.floating64],
        gain: onp.ToFloat,
        /,
        *,
        dt: DTT,
    ) -> ZerosPolesGainDiscrete[np.complex128, np.float64, DTT]: ...
    @overload
    def __new__(
        cls, zeros: _ToComplex12D, poles: onp.ToFloat1D, gain: onp.ToFloat, /, *, dt: None = None
    ) -> ZerosPolesGainContinuous[Any, _Float]: ...
    @overload
    def __new__[DTT: onp.ToComplex](
        cls, zeros: _ToComplex12D, poles: onp.ToFloat1D, gain: onp.ToFloat, /, *, dt: DTT
    ) -> ZerosPolesGainDiscrete[Any, _Float, DTT]: ...

    #
    @overload
    def __init__(self, system: LinearTimeInvariant[_ZerosT_co, _PolesT_co, _DTT_co], /, *, dt: None = None) -> None: ...
    @overload
    def __init__(self, zeros: _ToComplex12D, poles: onp.ToFloat1D, gain: onp.ToFloat, /, *, dt: _DTT_co = ...) -> None: ...

    #
    @property
    @override
    def zeros(self, /) -> _Array12D[_ZerosT_co]: ...
    @zeros.setter
    @override
    def zeros(self, zeros: _ToComplex12D, /) -> None: ...

    #
    @property
    @override
    def poles(self, /) -> onp.Array1D[_PolesT_co]: ...
    @poles.setter
    @override
    def poles(self, gain: onp.ToFloat1D, /) -> None: ...

    #
    @property
    def gain(self, /) -> float: ...
    @gain.setter
    def gain(self, gain: float, /) -> None: ...

    #
    @override
    def to_tf(self, /) -> TransferFunction[_PolesT_co, _DTT_co]: ...
    @override
    def to_zpk(self, /) -> Self: ...
    @override
    def to_ss(self, /) -> StateSpace[_ZerosT_co, _PolesT_co, _DTT_co]: ...

@final
class ZerosPolesGainContinuous(
    ZerosPolesGain[_ZerosT_co, _PolesT_co, None], lti[_ZerosT_co, _PolesT_co], Generic[_ZerosT_co, _PolesT_co]
):
    @override
    def to_zpk(self, /) -> Self: ...
    @override
    def to_discrete[DTT: onp.ToComplex | None](
        self, /, dt: DTT, method: _DiscretizeMethod = "zoh", alpha: float | None = None
    ) -> ZerosPolesGainDiscrete[_ZerosT_co, _PolesT_co, DTT]: ...

@final
class ZerosPolesGainDiscrete(
    ZerosPolesGain[_ZerosT_co, _PolesT_co, _DTT_co],
    dlti[_ZerosT_co, _PolesT_co, _DTT_co],
    Generic[_ZerosT_co, _PolesT_co, _DTT_co],
):
    @overload
    def __init__(self, system: ZerosPolesGain[_ZerosT_co, _PolesT_co, _DTT_co], /) -> None: ...
    @overload
    def __init__[DTT: onp.ToComplex | None](
        self: ZerosPolesGainDiscrete[np.float64, np.float64, DTT],
        zeros: onp.ToArray1D[float, npc.integer | npc.floating64] | onp.ToArray2D[float, npc.integer | npc.floating64],
        poles: onp.ToArray1D[float, npc.integer | npc.floating64],
        gain: onp.ToFloat,
        /,
        *,
        dt: DTT = ...,
    ) -> None: ...
    @overload
    def __init__[DTT: onp.ToComplex | None](
        self: ZerosPolesGainDiscrete[_Float, _Float, DTT],
        zeros: _ToFloat12D,
        poles: onp.ToFloat1D,
        gain: onp.ToFloat,
        /,
        *,
        dt: DTT = ...,
    ) -> None: ...
    @overload
    def __init__[DTT: onp.ToComplex | None](
        self: ZerosPolesGainDiscrete[np.complex128, np.float64, DTT],
        zeros: onp.ToJustComplex128_1D | onp.ToJustComplex128_2D,
        poles: onp.ToArray1D[float, npc.integer | npc.floating64],
        gain: onp.ToFloat,
        /,
        *,
        dt: DTT = ...,
    ) -> None: ...
    @overload
    def __init__[DTT: onp.ToComplex | None](
        self: ZerosPolesGainDiscrete[Any, _Float, DTT],
        zeros: _ToComplex12D,
        poles: onp.ToFloat1D,
        gain: onp.ToFloat,
        /,
        *,
        dt: DTT = ...,
    ) -> None: ...

    #
    @override
    def output(
        self, /, u: _ToFloat012D | None, t: onp.ToFloat1D, x0: onp.ToFloat1D | None = None
    ) -> tuple[_Float64_1D, _Float64_2D]: ...

class StateSpace(LinearTimeInvariant[_ZerosT_co, _PolesT_co, _DTT_co], Generic[_ZerosT_co, _PolesT_co, _DTT_co]):
    __array_priority__: ClassVar[float] = 100.0
    __array_ufunc__: ClassVar[None] = None

    @override
    @overload
    def __new__(
        cls, system: lti[_ZerosT_co, _PolesT_co], /, *, dt: None = None
    ) -> StateSpaceContinuous[_ZerosT_co, _PolesT_co]: ...
    @overload
    def __new__(
        cls, system: dlti[_ZerosT_co, _PolesT_co, _DTT_co], /, *, dt: None = None
    ) -> StateSpaceDiscrete[_ZerosT_co, _PolesT_co, _DTT_co]: ...
    @overload
    def __new__(
        cls,
        A: onp.ToArrayND[float, npc.integer | npc.floating64],
        B: _ToFloat012D,
        C: _ToFloat012D,
        D: _ToFloat012D,
        /,
        *,
        dt: None = None,
    ) -> StateSpaceContinuous[np.float64, np.float64]: ...
    @overload
    def __new__[DTT: onp.ToComplex](
        cls,
        A: onp.ToArrayND[float, npc.integer | npc.floating64],
        B: _ToFloat012D,
        C: _ToFloat012D,
        D: _ToFloat012D,
        /,
        *,
        dt: DTT,
    ) -> StateSpaceDiscrete[np.float64, np.float64, DTT]: ...
    @overload
    def __new__(
        cls, A: _ToFloat012D, B: _ToFloat012D, C: _ToFloat012D, D: _ToFloat012D, /, *, dt: None = None
    ) -> StateSpaceContinuous[_Float, _Float]: ...
    @overload
    def __new__[DTT: onp.ToComplex](
        cls, A: _ToFloat012D, B: _ToFloat012D, C: _ToFloat012D, D: _ToFloat012D, /, *, dt: DTT
    ) -> StateSpaceDiscrete[_Float, _Float, DTT]: ...
    @overload
    def __new__(
        cls, A: onp.ToJustComplex128_ND, B: _ToComplex012D, C: _ToComplex012D, D: _ToComplex012D, /, *, dt: None = None
    ) -> StateSpaceContinuous[np.complex128, np.float64]: ...
    @overload
    def __new__[DTT: onp.ToComplex](
        cls, A: onp.ToJustComplex128_ND, B: _ToComplex012D, C: _ToComplex012D, D: _ToComplex012D, /, *, dt: DTT
    ) -> StateSpaceDiscrete[np.complex128, np.float64, DTT]: ...
    @overload
    def __new__(
        cls, A: _ToComplex012D, B: _ToComplex012D, C: _ToComplex012D, D: _ToComplex012D, /, *, dt: None = None
    ) -> StateSpaceContinuous[Any, _Float]: ...
    @overload
    def __new__[DTT: onp.ToComplex](
        cls, A: _ToComplex012D, B: _ToComplex012D, C: _ToComplex012D, D: _ToComplex012D, /, *, dt: DTT
    ) -> StateSpaceDiscrete[Any, _Float, DTT]: ...

    #
    @overload
    def __init__(self, system: StateSpace[_ZerosT_co, _PolesT_co, _DTT_co], /) -> None: ...
    @overload
    def __init__(
        self: StateSpace[np.float64, np.float64],
        A: onp.ToArrayND[float, npc.integer | npc.floating64],
        B: _ToFloat012D,
        C: _ToFloat012D,
        D: _ToFloat012D,
        /,
        *,
        dt: _DTT_co = ...,
    ) -> None: ...
    @overload
    def __init__(
        self: StateSpace[_Float, _Float],
        A: _ToFloat012D,
        B: _ToFloat012D,
        C: _ToFloat012D,
        D: _ToFloat012D,
        /,
        *,
        dt: _DTT_co = ...,
    ) -> None: ...
    @overload
    def __init__(
        self: StateSpace[np.complex128, np.float64],
        A: onp.ToJustComplex128_ND,
        B: _ToComplex012D,
        C: _ToComplex012D,
        D: _ToComplex012D,
        /,
        *,
        dt: _DTT_co = ...,
    ) -> None: ...
    @overload
    def __init__(
        self: StateSpace[Any, _Float],
        A: _ToComplex012D,
        B: _ToComplex012D,
        C: _ToComplex012D,
        D: _ToComplex012D,
        /,
        *,
        dt: _DTT_co = ...,
    ) -> None: ...

    #
    def __neg__(self, /) -> Self: ...
    def __add__(self, other: Self | npc.number | onp.ArrayND[npc.number], /) -> Self: ...
    def __sub__(self, other: Self | npc.number | onp.ArrayND[npc.number], /) -> Self: ...
    def __mul__(self, other: Self | npc.number | onp.ArrayND[npc.number], /) -> Self: ...
    def __truediv__(self, other: _ToNumber, /) -> Self: ...
    def __radd__(self, other: Self | npc.number | onp.ArrayND[npc.number], /) -> Self: ...
    def __rsub__(self, other: Self | npc.number | onp.ArrayND[npc.number], /) -> Self: ...
    def __rmul__(self, other: Self | npc.number | onp.ArrayND[npc.number], /) -> Self: ...

    #
    @property
    def A(self, /) -> onp.Array2D[_Number]: ...
    @A.setter
    def A(self, /, A: onp.ToComplex2D) -> None: ...
    #
    @property
    def B(self, /) -> onp.Array2D[_Number]: ...
    @B.setter
    def B(self, /, B: onp.ToComplex2D) -> None: ...
    #
    @property
    def C(self, /) -> onp.Array2D[_Number]: ...
    @C.setter
    def C(self, /, C: onp.ToComplex2D) -> None: ...
    #
    @property
    def D(self, /) -> onp.Array2D[_Number]: ...
    @D.setter
    def D(self, /, D: onp.ToComplex2D) -> None: ...
    #
    @override
    def to_tf(self, /, *, input: onp.ToInt = 0) -> TransferFunction[_PolesT_co, _DTT_co]: ...
    @override
    def to_zpk(self, /, *, input: onp.ToInt = 0) -> ZerosPolesGain[_ZerosT_co, _PolesT_co, _DTT_co]: ...
    @override
    def to_ss(self, /) -> Self: ...

@final
class StateSpaceContinuous(
    StateSpace[_ZerosT_co, _PolesT_co, None], lti[_ZerosT_co, _PolesT_co], Generic[_ZerosT_co, _PolesT_co]
):
    @override
    def to_ss(self, /) -> Self: ...
    @override
    def to_discrete[DTT: onp.ToComplex | None](
        self, /, dt: DTT, method: _DiscretizeMethod = "zoh", alpha: float | None = None
    ) -> StateSpaceDiscrete[_ZerosT_co, _PolesT_co, DTT]: ...

@final
class StateSpaceDiscrete(
    StateSpace[_ZerosT_co, _PolesT_co, _DTT_co], dlti[_ZerosT_co, _PolesT_co, _DTT_co], Generic[_ZerosT_co, _PolesT_co, _DTT_co]
):
    @overload
    def __init__(self, system: StateSpace[_ZerosT_co, _PolesT_co, _DTT_co], /) -> None: ...
    @overload
    def __init__[DTT: onp.ToComplex | None](
        self: StateSpaceDiscrete[np.float64, np.float64, DTT],
        A: onp.ToArrayND[float, npc.integer | npc.floating64],
        B: _ToFloat012D,
        C: _ToFloat012D,
        D: _ToFloat012D,
        /,
        *,
        dt: DTT = ...,
    ) -> None: ...
    @overload
    def __init__[DTT: onp.ToComplex | None](
        self: StateSpaceDiscrete[_Float, _Float, DTT],
        A: _ToFloat012D,
        B: _ToFloat012D,
        C: _ToFloat012D,
        D: _ToFloat012D,
        /,
        *,
        dt: DTT = ...,
    ) -> None: ...
    @overload
    def __init__[DTT: onp.ToComplex | None](
        self: StateSpaceDiscrete[np.complex128, np.float64, DTT],
        A: onp.ToJustComplex128_ND,
        B: _ToComplex012D,
        C: _ToComplex012D,
        D: _ToComplex012D,
        /,
        *,
        dt: DTT = ...,
    ) -> None: ...
    @overload
    def __init__[DTT: onp.ToComplex | None](
        self: StateSpaceDiscrete[Any, _Float, DTT],
        A: _ToComplex012D,
        B: _ToComplex012D,
        C: _ToComplex012D,
        D: _ToComplex012D,
        /,
        *,
        dt: DTT = ...,
    ) -> None: ...

    #
    @override
    def output(
        self, /, u: _ToFloat012D | None, t: onp.ToFloat1D, x0: onp.ToFloat1D | None = None
    ) -> tuple[_Float64_1D, _Float64_2D, _Float64_2D]: ...

# NOTE: Only used as return type of `place_poles`.
@final
class Bunch(Generic[_RequestedPolesT_co]):
    gain_matrix: _Float64_2D
    computed_poles: onp.Array1D[np.complex128]
    requested_poles: onp.Array1D[_RequestedPolesT_co]
    X: onp.Array2D[np.float64 | np.complex128]  # float64 iff `B` is square and full-rank
    rtol: float
    nb_iter: float  # `nan` if `B` is square and full-rank, `int` otherwise

    def __init__(self, /) -> None: ...

#
@overload  # +f64
def place_poles(
    A: onp.ToFloat2D,
    B: onp.ToFloat2D,
    poles: onp.ToFloat64_1D,
    method: Literal["YT", "KNV0"] = "YT",
    rtol: float = 0.001,
    maxiter: int = 30,
) -> Bunch[np.float64]: ...
@overload  # +floating
def place_poles(
    A: onp.ToFloat2D,
    B: onp.ToFloat2D,
    poles: onp.ToFloat1D,
    method: Literal["YT", "KNV0"] = "YT",
    rtol: float = 0.001,
    maxiter: int = 30,
) -> Bunch[np.float64 | Any]: ...
@overload  # ~c128
def place_poles(
    A: onp.ToFloat2D,
    B: onp.ToFloat2D,
    poles: onp.ToJustComplex128_1D,
    method: Literal["YT", "KNV0"] = "YT",
    rtol: float = 0.001,
    maxiter: int = 30,
) -> Bunch[np.complex128]: ...
@overload  # ~complexfloating
def place_poles(
    A: onp.ToFloat2D,
    B: onp.ToFloat2D,
    poles: onp.ToJustComplex1D,
    method: Literal["YT", "KNV0"] = "YT",
    rtol: float = 0.001,
    maxiter: int = 30,
) -> Bunch[np.complex128 | Any]: ...
@overload  # +complexfloating
def place_poles(
    A: onp.ToFloat2D,
    B: onp.ToFloat2D,
    poles: onp.ToComplex1D,
    method: Literal["YT", "KNV0"] = "YT",
    rtol: float = 0.001,
    maxiter: int = 30,
) -> Bunch[Any]: ...

#

# keep in sync with impulse and step
@overload
def lsim(
    system: lti[np.float32 | np.float64] | _ToLTIFloat,
    U: _ToFloat012D | None,
    T: onp.ToFloat1D,
    X0: onp.ToFloat1D | None = None,
    interp: bool = True,
) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
@overload
def lsim(
    system: lti[np.complex64 | np.complex128],
    U: _ToFloat012D | None,
    T: onp.ToFloat1D,
    X0: onp.ToComplex1D | None = None,
    interp: bool = True,
) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.complex128], onp.ArrayND[np.complex128]]: ...
@overload
def lsim(
    system: _ToLTIInexact, U: _ToFloat012D | None, T: onp.ToFloat1D, X0: onp.ToComplex1D | None = None, interp: bool = True
) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.complex128 | Any], onp.ArrayND[np.complex128 | Any]]: ...

# keep in sync with lsim and step
@overload
def impulse(
    system: lti[np.float32 | np.float64] | _ToLTIFloat,
    X0: onp.ToFloat1D | None = None,
    T: onp.ToFloat1D | None = None,
    N: int | None = None,
) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.float64]]: ...
@overload
def impulse(
    system: lti[np.complex64 | np.complex128],
    X0: onp.ToComplex1D | None = None,
    T: onp.ToFloat1D | None = None,
    N: int | None = None,
) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.complex128]]: ...
@overload
def impulse(
    system: _ToLTIInexact, X0: onp.ToComplex1D | None = None, T: onp.ToFloat1D | None = None, N: int | None = None
) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.complex128 | Any]]: ...

# keep in sync with lsim and impulse
@overload
def step(
    system: lti[np.float32 | np.float64] | _ToLTIFloat,
    X0: onp.ToFloat1D | None = None,
    T: onp.ToFloat1D | None = None,
    N: int | None = None,
) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.float64]]: ...
@overload
def step(
    system: lti[np.complex64 | np.complex128],
    X0: onp.ToComplex1D | None = None,
    T: onp.ToFloat1D | None = None,
    N: int | None = None,
) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.complex128]]: ...
@overload
def step(
    system: _ToLTIInexact, X0: onp.ToComplex1D | None = None, T: onp.ToFloat1D | None = None, N: int | None = None
) -> tuple[onp.Array1D[np.float64], onp.ArrayND[np.complex128 | Any]]: ...

#
@overload
def bode(
    system: lti[npc.inexact64] | _ToLTIInexact64, w: onp.ToFloat1D | None = None, n: int = 100
) -> _Tuple3[onp.Array1D[np.float64]]: ...
@overload
def bode(
    system: lti[npc.inexact32] | _ToLTIInexact32, w: onp.ToFloat1D | None = None, n: int = 100
) -> _Tuple3[onp.Array1D[np.float32]]: ...
@overload
def bode(system: lti | _ToLTIInexact, w: onp.ToFloat1D | None = None, n: int = 100) -> _Tuple3[_Float1D]: ...

#
@overload
def freqresp(
    system: lti[npc.inexact64] | _ToLTIInexact64, w: onp.ToFloat1D | None = None, n: int = 10_000
) -> tuple[onp.Array1D[np.float64], onp.Array1D[np.complex128]]: ...
@overload
def freqresp(
    system: lti[npc.inexact32] | _ToLTIInexact32, w: onp.ToFloat1D | None = None, n: int = 10_000
) -> tuple[onp.Array1D[np.float32], onp.Array1D[np.complex64]]: ...
@overload
def freqresp(system: _ToLTIInexact, w: onp.ToFloat1D | None = None, n: int = 10_000) -> tuple[_Float1D, _Complex1D]: ...

#

#
@overload
def dlsim(
    system: StateSpaceDiscrete | _ToSSDisc,
    u: _ToFloat012D | None,
    t: onp.ToFloat1D | None = None,
    x0: onp.ToFloat1D | None = None,
) -> tuple[_Float64_1D, _Float64_2D, _Float64_2D]: ...
@overload
def dlsim(
    system: TransferFunctionDiscrete | ZerosPolesGainDiscrete | _ToTFDisc | _ToZPKDisc,
    u: _ToFloat012D | None,
    t: onp.ToFloat1D | None = None,
    x0: onp.ToFloat1D | None = None,
) -> tuple[_Float64_1D, _Float64_2D]: ...

#
def dimpulse(
    system: dlti | _ToDLTI, x0: onp.ToFloat1D | None = None, t: onp.ToFloat1D | None = None, n: int | None = None
) -> tuple[onp.Array1D[np.float64], tuple[onp.Array2D[np.float64], ...]]: ...

#
def dstep(
    system: dlti | _ToDLTI, x0: onp.ToFloat1D | None = None, t: onp.ToFloat1D | None = None, n: int | None = None
) -> tuple[onp.Array1D[np.float64], tuple[onp.Array2D[np.float64], ...]]: ...

#
def dbode(system: dlti | _ToDLTI, w: onp.ToFloat1D | None = None, n: int = 100) -> _Tuple3[onp.Array1D[np.float64]]: ...

#
def dfreqresp(
    system: dlti | _ToDLTI, w: onp.ToFloat1D | None = None, n: int = 10_000, whole: bool = False
) -> tuple[onp.Array1D[np.float64], onp.Array1D[np.complex128]]: ...
