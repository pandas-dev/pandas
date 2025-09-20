from collections.abc import Callable
from typing import Concatenate, Final, Generic, Literal, TypedDict, overload, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp

from scipy._lib._util import _RichResult

__all__ = ["nsum"]

_SCT = TypeVar("_SCT")

@type_check_only
class _Tolerances(TypedDict, total=False):
    rtol: onp.ToFloat
    atol: onp.ToFloat

@type_check_only
class _TanhSinhResult(_RichResult[bool | int | _SCT], Generic[_SCT]):
    success: Final[bool]
    status: Final[Literal[0, -1, -2, -3, -4, 1]]
    integral: _SCT
    error: Final[np.float64]
    maxlevel: Final[int]
    nfev: Final[int]

@type_check_only
class _NSumResult(_RichResult[bool | int | _SCT], Generic[_SCT]):
    success: Final[bool]
    status: Final[Literal[0, -1, -2, -3]]
    sum: _SCT
    error: Final[np.float64]
    maxlevel: Final[int]
    nfev: Final[int]

###

@overload
def tanhsinh(
    f: Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToFloat1D],
    a: onp.ToFloat | onp.ToFloat1D,
    b: onp.ToFloat | onp.ToFloat1D,
    *,
    args: tuple[object, ...] = (),
    log: op.CanBool = False,
    maxlevel: onp.ToInt | None = None,
    minlevel: onp.ToInt | None = 2,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    preserve_shape: bool = False,
    callback: Callable[[_TanhSinhResult[np.float64]], None] | None = None,
) -> _TanhSinhResult[np.float64]: ...
@overload
def tanhsinh(
    f: Callable[Concatenate[onp.Array1D[np.float64 | np.complex128], ...], onp.ToComplex1D],
    a: onp.ToFloat | onp.ToFloat1D,
    b: onp.ToFloat | onp.ToFloat1D,
    *,
    args: tuple[object, ...] = (),
    log: op.CanBool = False,
    maxlevel: onp.ToInt | None = None,
    minlevel: onp.ToInt | None = 2,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    preserve_shape: bool = False,
    callback: Callable[[_TanhSinhResult[np.float64 | np.complex128]], None] | None = None,
) -> _TanhSinhResult[np.float64 | np.complex128]: ...

#
@overload
def nsum(
    f: Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToFloat1D],
    a: onp.ToFloat | onp.ToFloat1D,
    b: onp.ToFloat | onp.ToFloat1D,
    *,
    step: onp.ToFloat | onp.ToFloat1D = 1,
    args: tuple[object, ...] = (),
    log: op.CanBool = False,
    maxterms: onp.ToInt = 0x10_00_00,
    tolerances: _Tolerances | None = None,
) -> _NSumResult[np.float64]: ...
@overload
def nsum(
    f: Callable[Concatenate[onp.Array1D[np.float64 | np.complex128], ...], onp.ToComplex1D],
    a: onp.ToFloat | onp.ToFloat1D,
    b: onp.ToFloat | onp.ToFloat1D,
    *,
    step: onp.ToFloat | onp.ToFloat1D = 1,
    args: tuple[object, ...] = (),
    log: op.CanBool = False,
    maxterms: onp.ToInt = 0x10_00_00,
    tolerances: _Tolerances | None = None,
) -> _NSumResult[np.float64 | np.complex128]: ...
