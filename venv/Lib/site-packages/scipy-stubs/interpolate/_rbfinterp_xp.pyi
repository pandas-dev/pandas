from _typeshed import Incomplete
from collections.abc import Callable
from types import ModuleType
from typing import Any, Final, Literal

import optype as op

###

type _IntMat = Incomplete

type _KernelFunc = Callable[[Any, ModuleType], Any]
type _KernelName = Literal[
    "linear", "thin_plate_spline", "cubic", "quintic", "multiquadric", "inverse_multiquadric", "inverse_quadratic", "gaussian"
]

###
# scipy-stubs does not yet support array-api, so we assume numpy arrays here.

def _monomial_powers(ndim: int, degree: int, xp: ModuleType) -> _IntMat: ...
def _build_and_solve_system[FloatMatT, FloatVecT](
    y: FloatMatT, d: FloatMatT, smoothing: FloatVecT, kernel: str, epsilon: float, powers: _IntMat, xp: ModuleType
) -> tuple[FloatMatT, FloatVecT, FloatVecT]: ...

#
def linear[T](r: op.CanNeg[T], xp: ModuleType) -> T: ...
def cubic[T](r: op.CanPow2[Literal[3], T], xp: ModuleType) -> T: ...
def quintic[T](r: op.CanPow2[Literal[5], T], xp: ModuleType) -> T: ...
def thin_plate_spline[InexactArrT](r: InexactArrT, xp: ModuleType) -> InexactArrT: ...
def multiquadric[InexactArrT](r: InexactArrT, xp: ModuleType) -> InexactArrT: ...
def inverse_multiquadric[InexactArrT](r: InexactArrT, xp: ModuleType) -> InexactArrT: ...
def inverse_quadratic[InexactArrT](r: InexactArrT, xp: ModuleType) -> InexactArrT: ...
def gaussian[InexactArrT](r: InexactArrT, xp: ModuleType) -> InexactArrT: ...

NAME_TO_FUNC: Final[dict[_KernelName, _KernelFunc]] = ...

def kernel_matrix[FloatMatT, ModuleT: ModuleType, T](
    x: FloatMatT, kernel_func: Callable[[FloatMatT, ModuleT], T], xp: ModuleT
) -> T: ...
def polynomial_matrix[FloatMatT](x: FloatMatT, powers: _IntMat, xp: ModuleType) -> FloatMatT: ...

#
def _build_system[FloatMatT, FloatVecT](
    y: FloatMatT, d: FloatMatT, smoothing: FloatVecT, kernel: _KernelName, epsilon: float, powers: _IntMat, xp: ModuleType
) -> tuple[FloatMatT, FloatMatT, FloatVecT, FloatVecT]: ...
def _build_evaluation_coefficients[FloatMatT, FloatVecT](
    x: FloatMatT,
    y: FloatMatT,
    kernel: _KernelName,
    epsilon: float,
    powers: _IntMat,
    shift: FloatVecT,
    scale: FloatVecT,
    xp: ModuleType,
) -> FloatMatT: ...
def compute_interpolation[FloatMatT, FloatVecT](
    x: FloatMatT,
    y: FloatMatT,
    kernel: _KernelName,
    epsilon: float,
    powers: _IntMat,
    shift: FloatVecT,
    scale: FloatVecT,
    coeffs: FloatMatT,
    xp: ModuleType,
) -> FloatMatT: ...
