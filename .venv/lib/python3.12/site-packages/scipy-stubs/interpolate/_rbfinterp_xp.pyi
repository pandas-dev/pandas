from _typeshed import Incomplete
from collections.abc import Callable
from types import ModuleType
from typing import Any, Final, Literal, TypeAlias, TypeVar

import optype as op

_KernelFunc: TypeAlias = Callable[[Any, ModuleType], Any]

_T = TypeVar("_T")
_ModuleT = TypeVar("_ModuleT", bound=ModuleType)

# TODO: Add upper bounds (a protocol with one or two of the required array-api methods)
# For now, we'll go with a "type-safety by naming convention" approach; type-safety purely on vibes.
_FloatVecT = TypeVar("_FloatVecT")
_FloatMatT = TypeVar("_FloatMatT")
_InexactArrT = TypeVar("_InexactArrT")

_IntMat: TypeAlias = Incomplete

_KernelName: TypeAlias = Literal[
    "linear", "thin_plate_spline", "cubic", "quintic", "multiquadric", "inverse_multiquadric", "inverse_quadratic", "gaussian"
]

###
# scipy-stubs does not yet support array-api, so we assume numpy arrays here.

def _monomial_powers(ndim: int, degree: int, xp: ModuleType) -> _IntMat: ...
def _build_and_solve_system(
    y: _FloatMatT, d: _FloatMatT, smoothing: _FloatVecT, kernel: str, epsilon: float, powers: _IntMat, xp: ModuleType
) -> tuple[_FloatMatT, _FloatVecT, _FloatVecT]: ...

#
def linear(r: op.CanNeg[_T], xp: ModuleType) -> _T: ...
def cubic(r: op.CanPow2[Literal[3], _T], xp: ModuleType) -> _T: ...
def quintic(r: op.CanPow2[Literal[5], _T], xp: ModuleType) -> _T: ...
def thin_plate_spline(r: _InexactArrT, xp: ModuleType) -> _InexactArrT: ...
def multiquadric(r: _InexactArrT, xp: ModuleType) -> _InexactArrT: ...
def inverse_multiquadric(r: _InexactArrT, xp: ModuleType) -> _InexactArrT: ...
def inverse_quadratic(r: _InexactArrT, xp: ModuleType) -> _InexactArrT: ...
def gaussian(r: _InexactArrT, xp: ModuleType) -> _InexactArrT: ...

NAME_TO_FUNC: Final[dict[_KernelName, _KernelFunc]] = ...

def kernel_matrix(x: _FloatMatT, kernel_func: Callable[[_FloatMatT, _ModuleT], _T], xp: _ModuleT) -> _T: ...
def polynomial_matrix(x: _FloatMatT, powers: _IntMat, xp: ModuleType) -> _FloatMatT: ...

#
def _build_system(
    y: _FloatMatT, d: _FloatMatT, smoothing: _FloatVecT, kernel: _KernelName, epsilon: float, powers: _IntMat, xp: ModuleType
) -> tuple[_FloatMatT, _FloatMatT, _FloatVecT, _FloatVecT]: ...
def _build_evaluation_coefficients(
    x: _FloatMatT,
    y: _FloatMatT,
    kernel: _KernelName,
    epsilon: float,
    powers: _IntMat,
    shift: _FloatVecT,
    scale: _FloatVecT,
    xp: ModuleType,
) -> _FloatMatT: ...
def compute_interpolation(
    x: _FloatMatT,
    y: _FloatMatT,
    kernel: _KernelName,
    epsilon: float,
    powers: _IntMat,
    shift: _FloatVecT,
    scale: _FloatVecT,
    coeffs: _FloatMatT,
    xp: ModuleType,
) -> _FloatMatT: ...
