from collections.abc import Callable, Mapping
from typing import Any, Concatenate, Final, Generic, TypedDict, final, overload, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._lib._util import _RichResult

_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)

type _Bracket[T] = tuple[T, T]

@type_check_only
class _Tolerances(TypedDict, total=False):
    xatol: onp.ToFloat
    xrtol: onp.ToFloat
    fatol: onp.ToFloat
    frtol: onp.ToFloat

@type_check_only
@final
class _FindResult(_RichResult, Generic[_ShapeT_co]):
    success: onp.ArrayND[np.bool, _ShapeT_co]
    status: onp.ArrayND[np.int32, _ShapeT_co]
    x: onp.ArrayND[np.float64, _ShapeT_co]
    f_x: onp.ArrayND[np.float64, _ShapeT_co]
    nfev: onp.ArrayND[np.int32, _ShapeT_co]
    nit: onp.ArrayND[np.int32, _ShapeT_co]
    bracket: _Bracket[onp.ArrayND[np.float64, _ShapeT_co]]
    f_bracket: _Bracket[onp.ArrayND[np.float64, _ShapeT_co]]
    _order_keys: Final = ["success", "status", "x", "f_x", "nfev", "nit", "bracket", "f_bracket"]

@type_check_only
@final
class _BracketResult(_RichResult, Generic[_ShapeT_co]):
    success: onp.ArrayND[np.bool, _ShapeT_co]
    status: onp.ArrayND[np.int32, _ShapeT_co]
    bracket: _Bracket[onp.ArrayND[np.float64, _ShapeT_co]]
    f_bracket: _Bracket[onp.ArrayND[np.float64, _ShapeT_co]]
    nfev: onp.ArrayND[np.int32, _ShapeT_co]
    nit: onp.ArrayND[np.int32, _ShapeT_co]

###

# TODO(@jorenham): array-api support
@overload
def find_root[ShapeT: tuple[int, ...]](
    f: Callable[[onp.ArrayND[np.float64, ShapeT]], onp.ArrayND[npc.floating]],
    init: _Bracket[onp.ToFloat] | _Bracket[onp.ToFloatND],
    /,
    *,
    args: tuple[()] = (),
    kwargs: None = None,
    tolerances: _Tolerances | None = None,
    maxiter: int | None = None,
    callback: Callable[[_FindResult[ShapeT]], None] | None = None,
) -> _FindResult[ShapeT]: ...
@overload
def find_root[ShapeT: tuple[int, ...]](
    f: Callable[Concatenate[onp.ArrayND[np.float64, ShapeT], ...], onp.ArrayND[npc.floating]],
    init: _Bracket[onp.ToFloat] | _Bracket[onp.ToFloatND],
    /,
    *,
    args: tuple[object, ...],
    kwargs: Mapping[str, object] | None = None,
    tolerances: _Tolerances | None = None,
    maxiter: int | None = None,
    callback: Callable[[_FindResult[ShapeT]], None] | None = None,
) -> _FindResult[ShapeT]: ...
@overload
def find_root[ShapeT: tuple[int, ...]](
    f: Callable[Concatenate[onp.ArrayND[np.float64, ShapeT], ...], onp.ArrayND[npc.floating]],
    init: _Bracket[onp.ToFloat] | _Bracket[onp.ToFloatND],
    /,
    *,
    args: tuple[object, ...] = (),
    kwargs: Mapping[str, object],
    tolerances: _Tolerances | None = None,
    maxiter: int | None = None,
    callback: Callable[[_FindResult[ShapeT]], None] | None = None,
) -> _FindResult[ShapeT]: ...

# TODO(@jorenham): array-api support
@overload
def find_minimum[ShapeT: tuple[int, ...]](
    f: Callable[[onp.ArrayND[np.float64, ShapeT]], onp.ArrayND[npc.floating]],
    init: tuple[onp.ToFloat, onp.ToFloat, onp.ToFloat] | tuple[onp.ToFloatND, onp.ToFloatND, onp.ToFloatND],
    /,
    *,
    args: tuple[()] = (),
    kwargs: None = None,
    tolerances: _Tolerances | None = None,
    maxiter: int = 100,
    callback: Callable[[_FindResult[ShapeT]], None] | None = None,
) -> _FindResult[ShapeT]: ...
@overload
def find_minimum[ShapeT: tuple[int, ...]](
    f: Callable[Concatenate[onp.ArrayND[np.float64, ShapeT], ...], onp.ArrayND[npc.floating]],
    init: tuple[onp.ToFloat, onp.ToFloat, onp.ToFloat] | tuple[onp.ToFloatND, onp.ToFloatND, onp.ToFloatND],
    /,
    *,
    args: tuple[object, ...],
    kwargs: Mapping[str, object] | None = None,
    tolerances: _Tolerances | None = None,
    maxiter: int = 100,
    callback: Callable[[_FindResult[ShapeT]], None] | None = None,
) -> _FindResult[ShapeT]: ...
@overload
def find_minimum[ShapeT: tuple[int, ...]](
    f: Callable[Concatenate[onp.ArrayND[np.float64, ShapeT], ...], onp.ArrayND[npc.floating]],
    init: tuple[onp.ToFloat, onp.ToFloat, onp.ToFloat] | tuple[onp.ToFloatND, onp.ToFloatND, onp.ToFloatND],
    /,
    *,
    args: tuple[object, ...] = (),
    kwargs: Mapping[str, object],
    tolerances: _Tolerances | None = None,
    maxiter: int = 100,
    callback: Callable[[_FindResult[ShapeT]], None] | None = None,
) -> _FindResult[ShapeT]: ...

# TODO(@jorenham): array-api support
@overload
def bracket_root[ShapeT: tuple[int, ...]](
    f: Callable[[onp.ArrayND[np.float64, ShapeT]], onp.ArrayND[npc.floating]],
    xl0: onp.ToFloat | onp.ToFloatND,
    xr0: onp.ToFloat | onp.ToFloatND | None = None,
    *,
    xmin: onp.ToFloat | onp.ToFloatND | None = None,
    xmax: onp.ToFloat | onp.ToFloatND | None = None,
    factor: onp.ToFloat | onp.ToFloatND | None = None,
    args: tuple[()] = (),
    kwargs: None = None,
    maxiter: int = 1_000,
) -> _BracketResult[ShapeT]: ...
@overload
def bracket_root[ShapeT: tuple[int, ...]](
    f: Callable[Concatenate[onp.ArrayND[np.float64, ShapeT], ...], onp.ArrayND[npc.floating]],
    xl0: onp.ToFloat | onp.ToFloatND,
    xr0: onp.ToFloat | onp.ToFloatND | None = None,
    *,
    xmin: onp.ToFloat | onp.ToFloatND | None = None,
    xmax: onp.ToFloat | onp.ToFloatND | None = None,
    factor: onp.ToFloat | onp.ToFloatND | None = None,
    args: tuple[object, ...],
    kwargs: Mapping[str, object] | None = None,
    maxiter: int = 1_000,
) -> _BracketResult[ShapeT]: ...
@overload
def bracket_root[ShapeT: tuple[int, ...]](
    f: Callable[Concatenate[onp.ArrayND[np.float64, ShapeT], ...], onp.ArrayND[npc.floating]],
    xl0: onp.ToFloat | onp.ToFloatND,
    xr0: onp.ToFloat | onp.ToFloatND | None = None,
    *,
    xmin: onp.ToFloat | onp.ToFloatND | None = None,
    xmax: onp.ToFloat | onp.ToFloatND | None = None,
    factor: onp.ToFloat | onp.ToFloatND | None = None,
    args: tuple[object, ...] = (),
    kwargs: Mapping[str, object],
    maxiter: int = 1_000,
) -> _BracketResult[ShapeT]: ...

# TODO(@jorenham): array-api support
@overload
def bracket_minimum[ShapeT: tuple[int, ...]](
    f: Callable[[onp.ArrayND[np.float64, ShapeT]], onp.ArrayND[npc.floating]],
    xm0: onp.ToFloat | onp.ToFloatND,
    *,
    xl0: onp.ToFloat | onp.ToFloatND | None = None,
    xr0: onp.ToFloat | onp.ToFloatND | None = None,
    xmin: onp.ToFloat | onp.ToFloatND | None = None,
    xmax: onp.ToFloat | onp.ToFloatND | None = None,
    factor: onp.ToFloat | onp.ToFloatND | None = None,
    args: tuple[()] = (),
    kwargs: None = None,
    maxiter: int = 1_000,
) -> _BracketResult[ShapeT]: ...
@overload
def bracket_minimum[ShapeT: tuple[int, ...]](
    f: Callable[Concatenate[onp.ArrayND[np.float64, ShapeT], ...], onp.ArrayND[npc.floating]],
    xm0: onp.ToFloat | onp.ToFloatND,
    *,
    xl0: onp.ToFloat | onp.ToFloatND | None = None,
    xr0: onp.ToFloat | onp.ToFloatND | None = None,
    xmin: onp.ToFloat | onp.ToFloatND | None = None,
    xmax: onp.ToFloat | onp.ToFloatND | None = None,
    factor: onp.ToFloat | onp.ToFloatND | None = None,
    args: tuple[object, ...],
    kwargs: Mapping[str, object] | None = None,
    maxiter: int = 1_000,
) -> _BracketResult[ShapeT]: ...
@overload
def bracket_minimum[ShapeT: tuple[int, ...]](
    f: Callable[Concatenate[onp.ArrayND[np.float64, ShapeT], ...], onp.ArrayND[npc.floating]],
    xm0: onp.ToFloat | onp.ToFloatND,
    *,
    xl0: onp.ToFloat | onp.ToFloatND | None = None,
    xr0: onp.ToFloat | onp.ToFloatND | None = None,
    xmin: onp.ToFloat | onp.ToFloatND | None = None,
    xmax: onp.ToFloat | onp.ToFloatND | None = None,
    factor: onp.ToFloat | onp.ToFloatND | None = None,
    args: tuple[object, ...] = (),
    kwargs: Mapping[str, object],
    maxiter: int = 1_000,
) -> _BracketResult[ShapeT]: ...
