# mypy: disable-error-code="explicit-override"
import abc
from typing import Any, ClassVar, Final, Generic, Literal, Never, Self, TypeAlias, final, overload, type_check_only
from typing_extensions import TypeVar, override

import numpy as np
import optype as op
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

_Float: TypeAlias = np.float32 | np.float64
_Complex: TypeAlias = np.complex64 | np.complex128
_Inexact: TypeAlias = _Float | _Complex
_Number: TypeAlias = npc.integer | _Inexact

_ToNumber: TypeAlias = complex | _Number
_ToNumberOrND: TypeAlias = _ToNumber | onp.ArrayND[_Number]

_SCT = TypeVar("_SCT", bound=np.generic)
_Array12D: TypeAlias = onp.ArrayND[_SCT, tuple[int] | tuple[int, int]]
_Array012D: TypeAlias = onp.ArrayND[_SCT, onp.AtMost2D]

_Float1D: TypeAlias = onp.Array1D[_Float]
_Float64_1D: TypeAlias = onp.Array1D[np.float64]
_Float64_2D: TypeAlias = onp.Array2D[np.float64]
_Complex1D: TypeAlias = onp.Array1D[_Complex]
_Inexact1D: TypeAlias = onp.Array1D[_Inexact]

_ToFloat12D: TypeAlias = onp.ToFloat1D | onp.ToFloat2D
_ToFloat012D: TypeAlias = onp.ToFloat | _ToFloat12D
_ToComplex12D: TypeAlias = onp.ToComplex1D | onp.ToComplex2D
_ToComplex012D: TypeAlias = onp.ToComplex | _ToComplex12D

###

# numerator, denominator
_ToTFContReal: TypeAlias = tuple[_ToFloat12D, onp.ToComplex1D]
_ToTFContComplex: TypeAlias = tuple[_ToComplex12D, onp.ToComplex1D]
# numerator, denominator, dt
_ToTFDiscReal: TypeAlias = tuple[_ToFloat12D, onp.ToComplex1D, onp.ToFloat]

# zeros, poles, gain
_ToZPKContReal: TypeAlias = tuple[onp.ToFloat1D, onp.ToFloat1D, onp.ToFloat]
_ToZPKContComplex: TypeAlias = tuple[onp.ToComplex1D, onp.ToComplex1D, onp.ToFloat]
# zeros, poles, gain, dt
_ToZPKDiscReal: TypeAlias = tuple[onp.ToFloat1D, onp.ToFloat1D, onp.ToFloat, onp.ToFloat]

# A, B, C, D
_ToSSContReal: TypeAlias = tuple[onp.ToFloat2D, onp.ToFloat2D, onp.ToFloat2D, onp.ToFloat2D]
_ToSSContComplex: TypeAlias = tuple[onp.ToComplex2D, onp.ToComplex2D, onp.ToComplex2D, onp.ToComplex2D]
# A, B, C, D, dt
_ToSSDiscReal: TypeAlias = tuple[onp.ToFloat2D, onp.ToFloat2D, onp.ToFloat2D, onp.ToFloat2D, onp.ToFloat]

_ToLTIReal: TypeAlias = _ToTFContReal | _ToZPKContReal | _ToSSContReal
_ToLTIComplex: TypeAlias = _ToTFContComplex | _ToZPKContComplex | _ToSSContComplex
_ToDLTI: TypeAlias = _ToTFDiscReal | _ToZPKDiscReal | _ToSSDiscReal

_ZerosT = TypeVar("_ZerosT", bound=_Inexact)
_ZerosT_co = TypeVar("_ZerosT_co", bound=_Inexact, default=_Inexact, covariant=True)
_PolesT = TypeVar("_PolesT", bound=_Float)
_PolesT_co = TypeVar("_PolesT_co", bound=_Float, default=_Float, covariant=True)
_DTT = TypeVar("_DTT", bound=onp.ToComplex | None)
_DTT_co = TypeVar("_DTT_co", bound=onp.ToComplex | None, default=Any, covariant=True)

###

class LinearTimeInvariant(Generic[_ZerosT_co, _PolesT_co, _DTT_co]):
    inputs: Final[int]
    outputs: Final[int]

    def __new__(cls, *system: Never, **kwargs: Never) -> Self: ...

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
    @overload
    def __new__(cls, *system: *tuple[_ToFloat12D, onp.ToFloat1D]) -> TransferFunctionContinuous[_Float]: ...
    @overload
    def __new__(cls, *system: *tuple[onp.ToFloat1D, onp.ToFloat1D, onp.ToFloat]) -> ZerosPolesGainContinuous[_Float]: ...
    @overload
    def __new__(cls, *system: *tuple[onp.ToComplex1D, onp.ToComplex1D, onp.ToFloat]) -> ZerosPolesGainContinuous: ...
    @overload
    def __new__(cls, *system: *tuple[_ToFloat012D, _ToFloat012D, _ToFloat012D, _ToFloat012D]) -> StateSpaceContinuous[_Float]: ...
    @overload
    def __new__(cls, *system: *tuple[_ToComplex012D, _ToComplex012D, _ToComplex012D, _ToComplex012D]) -> StateSpaceContinuous: ...

    #
    def __init__(self, /, *system: Never) -> None: ...

    #
    def impulse(
        self, /, X0: onp.ToFloat1D | None = None, T: onp.ToFloat1D | None = None, N: onp.ToJustInt | None = None
    ) -> tuple[_Float1D, _Float1D]: ...
    def step(
        self, /, X0: onp.ToComplex1D | None = None, T: onp.ToFloat1D | None = None, N: onp.ToJustInt | None = None
    ) -> tuple[_Float1D, _Inexact1D]: ...
    def output(
        self, /, U: _ToFloat012D | None, T: onp.ToFloat1D, X0: onp.ToComplex1D | None = None
    ) -> tuple[_Array12D[_Float], _Inexact1D, _Array12D[_Inexact]]: ...
    def bode(self, /, w: onp.ToFloat1D | None = None, n: onp.ToJustInt = 100) -> tuple[_Float1D, _Float1D, _Float1D]: ...
    def freqresp(self, /, w: onp.ToFloat1D | None = None, n: onp.ToJustInt = 10_000) -> tuple[_Float1D, _Complex1D]: ...

    #
    @abc.abstractmethod
    def to_discrete(
        self, /, dt: _DTT, method: _DiscretizeMethod = "zoh", alpha: onp.ToJustFloat | None = None
    ) -> dlti[_ZerosT_co, _PolesT_co, _DTT]: ...

#
class dlti(LinearTimeInvariant[_ZerosT_co, _PolesT_co, _DTT_co], Generic[_ZerosT_co, _PolesT_co, _DTT_co], metaclass=abc.ABCMeta):
    @overload
    def __new__(
        cls, *system: *tuple[_ToFloat12D, onp.ToFloat1D], dt: _DTT_co = ...
    ) -> TransferFunctionDiscrete[_Float, _DTT_co]: ...
    @overload
    def __new__(
        cls, *system: *tuple[onp.ToFloat1D, onp.ToFloat2D, onp.ToFloat1D], dt: _DTT_co = ...
    ) -> ZerosPolesGainDiscrete[_Float, _Float, _DTT_co]: ...
    @overload
    def __new__(
        cls, *system: *tuple[onp.ToComplex1D, onp.ToComplex1D, onp.ToFloat], dt: _DTT_co = ...
    ) -> ZerosPolesGainDiscrete[_Inexact, _Float, _DTT_co]: ...
    @overload
    def __new__(
        cls, *system: *tuple[_ToFloat012D, _ToFloat012D, _ToFloat012D, _ToFloat012D], dt: _DTT_co = ...
    ) -> StateSpaceDiscrete[_Float, _Float, _DTT_co]: ...
    @overload
    def __new__(
        cls, *system: *tuple[_ToComplex012D, _ToComplex012D, _ToComplex012D, _ToComplex012D], dt: _DTT_co = ...
    ) -> StateSpaceDiscrete[_Inexact, _Float, _DTT_co]: ...

    #
    def __init__(self, /, *system: Never, dt: onp.ToFloat, **kwargs: Never) -> None: ...

    #
    def impulse(
        self, /, x0: onp.ToFloat1D | None = None, t: onp.ToFloat1D | None = None, n: onp.ToJustInt | None = None
    ) -> tuple[_Float1D, _Float1D]: ...
    def step(
        self, /, x0: onp.ToFloat1D | None = None, t: onp.ToFloat1D | None = None, n: onp.ToJustInt | None = None
    ) -> tuple[_Float1D, _Float1D]: ...
    def output(
        self, /, u: _ToFloat12D | onp.ToFloat | None, t: onp.ToFloat1D, x0: onp.ToFloat1D | None = None
    ) -> tuple[_Float1D, _Float1D]: ...
    def bode(self, /, w: onp.ToFloat1D | None = None, n: onp.ToJustInt = 100) -> tuple[_Float1D, _Float1D, _Float1D]: ...
    def freqresp(
        self, /, w: onp.ToFloat1D | None = None, n: onp.ToJustInt = 10_000, whole: op.CanBool = False
    ) -> tuple[_Float1D, _Complex1D]: ...

#
class TransferFunction(LinearTimeInvariant[_PolesT_co, _PolesT_co, _DTT_co], Generic[_PolesT_co, _DTT_co], metaclass=abc.ABCMeta):
    @overload
    def __new__(cls, *system: *tuple[lti[_PolesT, _PolesT]]) -> TransferFunctionContinuous[_PolesT]: ...
    @overload
    def __new__(cls, *system: *tuple[dlti[_PolesT, _PolesT, _DTT]]) -> TransferFunctionDiscrete[_PolesT, _DTT]: ...
    @overload
    def __new__(cls, *system: *tuple[_ToFloat12D, onp.ToFloat1D]) -> TransferFunctionContinuous[_Float]: ...
    @overload
    def __new__(cls, *system: *tuple[_ToFloat12D, onp.ToFloat1D], dt: _DTT) -> TransferFunctionDiscrete[_Float, _DTT]: ...

    #
    @overload
    def __init__(self, system: LinearTimeInvariant[_PolesT_co, _PolesT_co, _DTT_co], /) -> None: ...
    @overload
    def __init__(self, numerator: _ToFloat12D, denominator: onp.ToFloat1D, /) -> None: ...

    #
    @property
    def num(self, /) -> _Array12D[_PolesT_co]: ...
    @num.setter
    def num(self, /, num: _ToFloat12D) -> None: ...

    #
    @property
    def den(self, /) -> onp.Array1D[_PolesT_co]: ...
    @den.setter
    def den(self, /, den: onp.ToFloat1D) -> None: ...

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
    def to_discrete(
        self, /, dt: _DTT, method: _DiscretizeMethod = "zoh", alpha: onp.ToJustFloat | None = None
    ) -> TransferFunctionDiscrete[_PolesT_co, _DTT]: ...

@final
class TransferFunctionDiscrete(
    TransferFunction[_PolesT_co, _DTT_co], dlti[_PolesT_co, _PolesT_co, _DTT_co], Generic[_PolesT_co, _DTT_co]
):
    @overload
    def __init__(self, system: LinearTimeInvariant[_PolesT_co, _PolesT_co, _DTT_co], /) -> None: ...
    @overload
    def __init__(self, numerator: _ToFloat12D, denominator: onp.ToFloat1D, /, *, dt: _DTT_co = ...) -> None: ...

#
class ZerosPolesGain(LinearTimeInvariant[_ZerosT_co, _PolesT_co, _DTT_co], Generic[_ZerosT_co, _PolesT_co, _DTT_co]):
    @overload
    def __new__(cls, *system: *tuple[lti[_ZerosT_co, _PolesT_co]]) -> ZerosPolesGainContinuous[_ZerosT_co, _PolesT_co]: ...
    @overload
    def __new__(
        cls, *system: *tuple[dlti[_ZerosT_co, _PolesT_co, _DTT_co]]
    ) -> ZerosPolesGainDiscrete[_ZerosT_co, _PolesT_co, _DTT_co]: ...
    @overload
    def __new__(cls, *system: *tuple[_ToFloat12D, onp.ToFloat1D, onp.ToFloat]) -> ZerosPolesGainContinuous[_Float, _Float]: ...
    @overload
    def __new__(
        cls, *system: *tuple[_ToComplex12D, onp.ToFloat1D, onp.ToFloat]
    ) -> ZerosPolesGainContinuous[_Inexact, _Float]: ...
    @overload
    def __new__(
        cls, *system: *tuple[_ToFloat12D, onp.ToFloat1D, onp.ToFloat], dt: _DTT
    ) -> ZerosPolesGainDiscrete[_Float, _Float, _DTT]: ...
    @overload
    def __new__(
        cls, *system: *tuple[_ToComplex12D, onp.ToFloat1D, onp.ToFloat], dt: _DTT
    ) -> ZerosPolesGainDiscrete[_Inexact, _Float, _DTT]: ...

    #
    @overload
    def __init__(self, system: LinearTimeInvariant[_ZerosT_co, _PolesT_co, _DTT_co], /) -> None: ...
    @overload
    def __init__(self, zeros: _ToComplex12D, poles: onp.ToFloat1D, gain: onp.ToFloat, /) -> None: ...

    #
    @property
    @override
    def zeros(self, /) -> _Array12D[_ZerosT_co]: ...
    @zeros.setter
    def zeros(self, zeros: _ToComplex12D, /) -> None: ...

    #
    @property
    @override
    def poles(self, /) -> onp.Array1D[_PolesT_co]: ...
    @poles.setter
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
    def to_discrete(
        self, /, dt: _DTT, method: _DiscretizeMethod = "zoh", alpha: onp.ToJustFloat | None = None
    ) -> ZerosPolesGainDiscrete[_ZerosT_co, _PolesT_co, _DTT]: ...

@final
class ZerosPolesGainDiscrete(
    ZerosPolesGain[_ZerosT_co, _PolesT_co, _DTT_co],
    dlti[_ZerosT_co, _PolesT_co, _DTT_co],
    Generic[_ZerosT_co, _PolesT_co, _DTT_co],
):
    @overload
    def __init__(self, system: ZerosPolesGain[_ZerosT_co, _PolesT_co, _DTT_co], /) -> None: ...
    @overload
    def __init__(
        self: ZerosPolesGainDiscrete[_Float, _Float, _DTT],
        zeros: _ToFloat12D,
        poles: onp.ToFloat1D,
        gain: onp.ToFloat,
        /,
        *,
        dt: _DTT = ...,
    ) -> None: ...
    @overload
    def __init__(
        self: ZerosPolesGainDiscrete[_Inexact, _Float, _DTT],
        zeros: _ToComplex12D,
        poles: onp.ToFloat1D,
        gain: onp.ToFloat,
        /,
        *,
        dt: _DTT = ...,
    ) -> None: ...

class StateSpace(LinearTimeInvariant[_ZerosT_co, _PolesT_co, _DTT_co], Generic[_ZerosT_co, _PolesT_co, _DTT_co]):
    __array_priority__: ClassVar[float] = 100.0
    __array_ufunc__: ClassVar[None] = None

    @overload
    def __new__(cls, *system: *tuple[lti[_ZerosT_co, _PolesT_co]]) -> StateSpaceContinuous[_ZerosT_co, _PolesT_co]: ...
    @overload
    def __new__(
        cls, *system: *tuple[dlti[_ZerosT_co, _PolesT_co, _DTT_co]]
    ) -> StateSpaceDiscrete[_ZerosT_co, _PolesT_co, _DTT_co]: ...
    @overload
    def __new__(
        cls, *system: *tuple[_ToFloat012D, _ToFloat012D, _ToFloat012D, _ToFloat012D]
    ) -> StateSpaceContinuous[_Float, _Float]: ...
    @overload
    def __new__(
        cls, *system: *tuple[_ToComplex012D, _ToComplex012D, _ToComplex012D, _ToComplex012D]
    ) -> StateSpaceContinuous[_Inexact, _Float]: ...
    @overload
    def __new__(
        cls, *system: *tuple[_ToFloat012D, _ToFloat012D, _ToFloat012D, _ToFloat012D], dt: _DTT
    ) -> StateSpaceDiscrete[_Float, _Float, _DTT]: ...
    @overload
    def __new__(
        cls, *system: *tuple[_ToComplex012D, _ToComplex012D, _ToComplex012D, _ToComplex012D], dt: _DTT
    ) -> StateSpaceDiscrete[_Inexact, _Float, _DTT]: ...

    #
    @overload
    def __init__(self, system: StateSpace[_ZerosT_co, _PolesT_co, _DTT_co], /) -> None: ...
    @overload
    def __init__(
        self: StateSpace[_Float, _Float], A: _ToFloat012D, B: _ToFloat012D, C: _ToFloat012D, D: _ToFloat012D, /
    ) -> None: ...
    @overload
    def __init__(
        self: StateSpace[_Inexact, _Float], A: _ToComplex012D, B: _ToComplex012D, C: _ToComplex012D, D: _ToComplex012D, /
    ) -> None: ...

    #
    def __neg__(self, /) -> Self: ...
    def __add__(self, other: Self | _ToNumberOrND, /) -> Self: ...
    def __sub__(self, other: Self | _ToNumberOrND, /) -> Self: ...
    def __mul__(self, other: Self | _ToNumberOrND, /) -> Self: ...
    def __truediv__(self, other: _ToNumber, /) -> Self: ...
    # ehh mypy, u ok?
    def __radd__(self, other: _ToNumberOrND, /) -> Self: ...  # type: ignore[misc]
    def __rsub__(self, other: _ToNumberOrND, /) -> Self: ...  # type: ignore[misc]
    def __rmul__(self, other: _ToNumberOrND, /) -> Self: ...  # type: ignore[misc]

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
    def to_discrete(
        self, /, dt: _DTT, method: _DiscretizeMethod = "zoh", alpha: onp.ToJustFloat | None = None
    ) -> StateSpaceDiscrete[_ZerosT_co, _PolesT_co, _DTT]: ...

@final
class StateSpaceDiscrete(
    StateSpace[_ZerosT_co, _PolesT_co, _DTT_co], dlti[_ZerosT_co, _PolesT_co, _DTT_co], Generic[_ZerosT_co, _PolesT_co, _DTT_co]
):
    @overload
    def __init__(self, system: StateSpace[_ZerosT_co, _PolesT_co, _DTT_co], /) -> None: ...
    @overload
    def __init__(
        self: StateSpaceDiscrete[_Float, _Float, _DTT],
        A: onp.ToFloat2D,
        B: onp.ToFloat2D,
        C: onp.ToFloat2D,
        D: onp.ToFloat2D,
        /,
        *,
        dt: _DTT = ...,
    ) -> None: ...
    @overload
    def __init__(
        self: StateSpaceDiscrete[_Inexact, _Float, _DTT],
        A: onp.ToComplex2D,
        B: onp.ToComplex2D,
        C: onp.ToComplex2D,
        D: onp.ToComplex2D,
        /,
        *,
        dt: _DTT = ...,
    ) -> None: ...

# NOTE: Only used as return type of `place_poles`.
class Bunch:
    gain_matrix: _Float64_2D
    computed_poles: _Float64_1D
    requested_poles: _Float64_1D
    X: onp.Array2D[np.complex128]
    rtol: float
    nb_iter: int

    def __init__(self, /, **kwds: Never) -> None: ...

#
def place_poles(
    A: onp.ToFloat2D,
    B: onp.ToFloat2D,
    poles: onp.ToComplex1D,
    method: Literal["YT", "KNV0"] = "YT",
    rtol: float = 0.001,
    maxiter: int = 30,
) -> Bunch: ...

#
@overload
def lsim(
    system: lti[_Float] | _ToLTIReal,
    U: _ToFloat12D | onp.ToFloat | None,
    T: onp.ToInt1D,
    X0: onp.ToFloat1D | None = None,
    interp: op.CanBool = True,
) -> tuple[_Float1D, _Float1D, _Array012D[_Float]]: ...
@overload
def lsim(
    system: lti | _ToLTIComplex,
    U: _ToFloat12D | onp.ToFloat | None,
    T: onp.ToInt1D,
    X0: onp.ToComplex1D | None = None,
    interp: op.CanBool = True,
) -> tuple[_Float1D, _Inexact1D, _Array012D[_Inexact]]: ...

#
@overload
def dlsim(
    system: StateSpaceDiscrete,
    u: _ToFloat12D | onp.ToFloat | None,
    t: onp.ToFloat1D | None = None,
    x0: onp.ToFloat1D | None = None,
) -> tuple[_Float64_1D, _Float64_2D, _Float64_2D]: ...
@overload
def dlsim(
    system: _ToDLTI, u: _ToFloat12D | onp.ToFloat | None, t: onp.ToFloat1D | None = None, x0: onp.ToFloat1D | None = None
) -> tuple[_Float64_1D, _Float64_1D]: ...

#
@overload
def impulse(
    system: lti[_ZerosT, _PolesT],
    X0: onp.ToComplex1D | None = None,
    T: onp.ToFloat1D | None = None,
    N: onp.ToJustInt | None = None,
) -> tuple[onp.Array1D[_PolesT], onp.Array1D[_ZerosT]]: ...
@overload
def impulse(
    system: _ToLTIReal, X0: onp.ToFloat1D | None = None, T: onp.ToFloat1D | None = None, N: onp.ToJustInt | None = None
) -> tuple[_Float1D, _Float1D]: ...
@overload
def impulse(
    system: _ToLTIComplex, X0: onp.ToComplex1D | None = None, T: onp.ToFloat1D | None = None, N: onp.ToJustInt | None = None
) -> tuple[_Float1D, _Inexact1D]: ...

#
@overload
def dimpulse(
    system: dlti[_ZerosT, _PolesT],
    x0: onp.ToComplex1D | None = None,
    t: onp.ToFloat1D | None = None,
    n: onp.ToJustInt | None = None,
) -> tuple[onp.Array1D[_PolesT], onp.Array1D[_ZerosT]]: ...
@overload
def dimpulse(
    system: _ToDLTI, x0: onp.ToFloat1D | None = None, t: onp.ToFloat1D | None = None, n: onp.ToJustInt | None = None
) -> tuple[_Float1D, _Float1D]: ...

#
@overload
def step(
    system: lti[_ZerosT, _PolesT],
    X0: onp.ToComplex1D | None = None,
    T: onp.ToFloat1D | None = None,
    N: onp.ToJustInt | None = None,
) -> tuple[onp.Array1D[_PolesT], onp.Array1D[_ZerosT]]: ...
@overload
def step(
    system: _ToLTIReal, X0: onp.ToFloat1D | None = None, T: onp.ToFloat1D | None = None, N: onp.ToJustInt | None = None
) -> tuple[_Float1D, _Float1D]: ...
@overload
def step(
    system: _ToLTIComplex, X0: onp.ToComplex1D | None = None, T: onp.ToFloat1D | None = None, N: onp.ToJustInt | None = None
) -> tuple[_Float1D, _Inexact1D]: ...

#
@overload
def dstep(
    system: dlti[_ZerosT, _PolesT],
    x0: onp.ToComplex1D | None = None,
    t: onp.ToFloat1D | None = None,
    n: onp.ToJustInt | None = None,
) -> tuple[onp.Array1D[_PolesT], onp.Array1D[_ZerosT]]: ...
@overload
def dstep(
    system: _ToDLTI, x0: onp.ToFloat1D | None = None, t: onp.ToFloat1D | None = None, n: onp.ToJustInt | None = None
) -> tuple[_Float1D, _Float1D]: ...

#
def bode(
    system: lti | _ToLTIComplex, w: onp.ToFloat1D | None = None, n: onp.ToJustInt = 100
) -> tuple[_Float1D, _Float1D, _Float1D]: ...

#
def dbode(
    system: dlti | _ToDLTI, w: onp.ToFloat1D | None = None, n: onp.ToJustInt = 100
) -> tuple[_Float1D, _Float1D, _Float1D]: ...

#
def freqresp(
    system: lti | _ToLTIComplex, w: onp.ToFloat1D | None = None, n: onp.ToJustInt = 10_000
) -> tuple[_Float1D, _Complex1D]: ...

#
def dfreqresp(
    system: dlti | _ToDLTI, w: onp.ToFloat1D | None = None, n: onp.ToJustInt = 10_000, whole: op.CanBool = False
) -> tuple[_Float1D, _Complex1D]: ...
