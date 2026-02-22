from collections.abc import Callable
from typing import Any, Concatenate, Final, Generic, Literal, TypeAlias, TypedDict, overload, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._lib._util import _RichResult

__all__ = ["nsum"]

_ScalarT = TypeVar("_ScalarT", bound=np.generic)
_InT = TypeVar("_InT")
_ResultT = TypeVar("_ResultT")
_ResultT_co = TypeVar("_ResultT_co", covariant=True)
_SuccessT_co = TypeVar("_SuccessT_co", covariant=True)
_StatusT_co = TypeVar("_StatusT_co", covariant=True)
_MaxLevelT_co = TypeVar("_MaxLevelT_co", covariant=True)

_ArgsND: TypeAlias = tuple[onp.ToScalar | onp.ToArrayND, ...]
_Callback: TypeAlias = Callable[[_ResultT], object]  # return value is ignored

_Integrand: TypeAlias = Callable[Concatenate[_InT, ...], onp.ArrayND[_ScalarT]]
_IntegrandReal: TypeAlias = _Integrand[onp.ArrayND[np.float64], npc.floating]
_IntegrandComplex: TypeAlias = _Integrand[onp.ArrayND[np.float64] | onp.ArrayND[np.complex128], npc.complexfloating]

@type_check_only
class _TanhSinhResult(
    _RichResult[_ResultT_co | _SuccessT_co | _StatusT_co | _MaxLevelT_co],
    Generic[_ResultT_co, _SuccessT_co, _StatusT_co, _MaxLevelT_co],
):
    integral: _ResultT_co
    error: _ResultT_co
    success: _SuccessT_co
    status: _StatusT_co
    nfev: _StatusT_co
    maxlevel: _MaxLevelT_co

_TanhSinhResult0: TypeAlias = _TanhSinhResult[_ScalarT, np.bool_, np.int32, np.int64]
_TanhSinhResultN: TypeAlias = _TanhSinhResult[
    onp.ArrayND[_ScalarT],
    onp.ArrayND[np.bool_],
    onp.ArrayND[np.int32],
    onp.ArrayND[np.int64],
]  # fmt: skip
_TanhSinhResultN_: TypeAlias = _TanhSinhResult[
    onp.ArrayND[_ScalarT] | Any,
    onp.ArrayND[np.bool_] | Any,
    onp.ArrayND[np.int32] | Any,
    onp.ArrayND[np.int64] | Any,
]  # fmt: skip

@type_check_only
class _Tolerances(TypedDict, total=False):
    rtol: float
    atol: float

@type_check_only
class _NSumResult0(_RichResult[np.bool_ | np.int32 | np.float64]):
    sum: Final[np.float64]
    error: Final[np.float64]
    success: Final[np.bool_]
    status: Final[np.int32]
    nfev: Final[np.int32]

@type_check_only
class _NSumResultN(_RichResult[onp.ArrayND[np.bool_] | onp.ArrayND[np.int32] | onp.ArrayND[np.float64]]):
    sum: Final[onp.ArrayND[np.float64]]
    error: Final[onp.ArrayND[np.float64]]
    success: Final[onp.ArrayND[np.bool_]]
    status: Final[onp.ArrayND[np.int32]]
    nfev: Final[onp.ArrayND[np.int32]]

###

#
@overload  # real f, scalar a, scalar b, preserve_shape=False
def tanhsinh(
    f: _IntegrandReal,
    a: onp.ToFloat,
    b: onp.ToFloat,
    *,
    args: tuple[onp.ToScalar, ...] = (),
    log: bool = False,
    maxlevel: int | None = None,
    minlevel: int | None = 2,
    atol: float | None = None,
    rtol: float | None = None,
    preserve_shape: Literal[False] = False,
    callback: _Callback[_TanhSinhResult0[np.float64]] | None = None,
) -> _TanhSinhResult0[np.float64]: ...
@overload  # real f, scalar/array a, array b
def tanhsinh(
    f: _IntegrandReal,
    a: onp.ToFloat | onp.ToFloatND,
    b: onp.ToFloatND,
    *,
    args: _ArgsND = (),
    log: bool = False,
    maxlevel: int | None = None,
    minlevel: int | None = 2,
    atol: float | None = None,
    rtol: float | None = None,
    preserve_shape: bool = False,
    callback: _Callback[_TanhSinhResultN[np.float64]] | None = None,
) -> _TanhSinhResultN[np.float64]: ...
@overload  # real f, array a, scalar/array b
def tanhsinh(
    f: _IntegrandReal,
    a: onp.ToFloatND,
    b: onp.ToFloat | onp.ToFloatND,
    *,
    args: _ArgsND = (),
    log: bool = False,
    maxlevel: int | None = None,
    minlevel: int | None = 2,
    atol: float | None = None,
    rtol: float | None = None,
    preserve_shape: bool = False,
    callback: _Callback[_TanhSinhResultN[np.float64]] | None = None,
) -> _TanhSinhResultN[np.float64]: ...
@overload  # real f, fallback
def tanhsinh(
    f: _IntegrandReal,
    a: onp.ToFloat | onp.ToFloatND,
    b: onp.ToFloat | onp.ToFloatND,
    *,
    args: _ArgsND = (),
    log: bool = False,
    maxlevel: int | None = None,
    minlevel: int | None = 2,
    atol: float | None = None,
    rtol: float | None = None,
    preserve_shape: bool = False,
    callback: _Callback[_TanhSinhResultN_[np.float64]] | None = None,
) -> _TanhSinhResultN_[np.float64]: ...
@overload  # complex f, scalar a, scalar b
def tanhsinh(
    f: _IntegrandComplex,
    a: onp.ToFloat,
    b: onp.ToFloat,
    *,
    args: tuple[onp.ToScalar, ...] = (),
    log: bool = False,
    maxlevel: int | None = None,
    minlevel: int | None = 2,
    atol: float | None = None,
    rtol: float | None = None,
    preserve_shape: Literal[False] = False,
    callback: _Callback[_TanhSinhResult0[np.complex128]] | None = None,
) -> _TanhSinhResult0[np.complex128]: ...
@overload  # complex f, scalar/array a, array b
def tanhsinh(
    f: _IntegrandComplex,
    a: onp.ToFloat | onp.ToFloatND,
    b: onp.ToFloatND,
    *,
    args: _ArgsND = (),
    log: bool = False,
    maxlevel: int | None = None,
    minlevel: int | None = 2,
    atol: float | None = None,
    rtol: float | None = None,
    preserve_shape: bool = False,
    callback: _Callback[_TanhSinhResultN[np.complex128]] | None = None,
) -> _TanhSinhResultN[np.complex128]: ...
@overload  # complex f, array a, scalar/array b
def tanhsinh(
    f: _IntegrandComplex,
    a: onp.ToFloatND,
    b: onp.ToFloat | onp.ToFloatND,
    *,
    args: _ArgsND = (),
    log: bool = False,
    maxlevel: int | None = None,
    minlevel: int | None = 2,
    atol: float | None = None,
    rtol: float | None = None,
    preserve_shape: bool = False,
    callback: _Callback[_TanhSinhResultN[np.complex128]] | None = None,
) -> _TanhSinhResultN[np.complex128]: ...
@overload  # complex f, fallback
def tanhsinh(
    f: _IntegrandComplex,
    a: onp.ToFloat | onp.ToFloatND,
    b: onp.ToFloat | onp.ToFloatND,
    *,
    args: _ArgsND = (),
    log: bool = False,
    maxlevel: int | None = None,
    minlevel: int | None = 2,
    atol: float | None = None,
    rtol: float | None = None,
    preserve_shape: bool = False,
    callback: _Callback[_TanhSinhResultN_[np.complex128]] | None = None,
) -> _TanhSinhResultN_[np.complex128]: ...

#
@overload  # scalar a, scalar b, scalar step
def nsum(
    f: _IntegrandReal,
    a: onp.ToFloat,
    b: onp.ToFloat,
    *,
    step: onp.ToFloat = 1,
    args: tuple[onp.ToScalar, ...] = (),
    log: bool = False,
    maxterms: int = 0x10_00_00,
    tolerances: _Tolerances | None = None,
) -> _NSumResult0: ...
@overload  # scalar/array a, array b, scalar/array step
def nsum(
    f: _IntegrandReal,
    a: onp.ToFloat | onp.ToFloatND,
    b: onp.ToFloatND,
    *,
    step: onp.ToFloat | onp.ToFloatND = 1,
    args: _ArgsND = (),
    log: bool = False,
    maxterms: int = 0x10_00_00,
    tolerances: _Tolerances | None = None,
) -> _NSumResultN: ...
@overload  # array a, scalar/array b, array step
def nsum(
    f: _IntegrandReal,
    a: onp.ToFloatND,
    b: onp.ToFloat | onp.ToFloatND,
    *,
    step: onp.ToFloat | onp.ToFloatND = 1,
    args: _ArgsND = (),
    log: bool = False,
    maxterms: int = 0x10_00_00,
    tolerances: _Tolerances | None = None,
) -> _NSumResultN: ...
@overload  # scalar/array a, scalar/array b, array step
def nsum(
    f: _IntegrandReal,
    a: onp.ToFloat | onp.ToFloatND,
    b: onp.ToFloat | onp.ToFloatND,
    *,
    step: onp.ToFloatND,
    args: _ArgsND = (),
    log: bool = False,
    maxterms: int = 0x10_00_00,
    tolerances: _Tolerances | None = None,
) -> _NSumResultN: ...
