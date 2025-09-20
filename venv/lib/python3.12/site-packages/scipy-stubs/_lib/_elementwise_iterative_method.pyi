from collections.abc import Callable, Iterable, Mapping, Sequence
from types import ModuleType
from typing import Any, Concatenate, Final, TypeAlias
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._util import _RichResult

###

_FloatT = TypeVar("_FloatT", bound=npc.floating, default=np.float64)
_ShapeT = TypeVar("_ShapeT", bound=onp.AtLeast1D, default=onp.AtLeast0D[Any])
_FuncRealT = TypeVar("_FuncRealT", bound=Callable[Concatenate[onp.ArrayND[np.float64], ...], object])
_ModuleT = TypeVar("_ModuleT", bound=ModuleType, default=ModuleType)
_WorkT = TypeVar("_WorkT", bound=Mapping[str, Any])
_ResT = TypeVar("_ResT", bound=_RichResult, default=_RichResult)
_ToShapeT = TypeVar("_ToShapeT", bound=op.CanIndex | tuple[op.CanIndex, ...], default=onp.AtLeast0D)

_Ignored: TypeAlias = _ResT

###

_ESIGNERR: Final = -1
_ECONVERR: Final = -2
_EVALUEERR: Final = -3
_ECALLBACK: Final = -4
_EINPUTERR: Final = -5
_ECONVERGED: Final = 0
_EINPROGRESS: Final = 1

# TODO: complex
def _initialize(
    func: _FuncRealT,
    xs: Sequence[onp.ToFloat1D],
    args: tuple[onp.ToFloat1D, ...],
    complex_ok: onp.ToFalse = False,
    preserve_shape: bool | None = None,
    xp: _ModuleT | None = None,
) -> tuple[
    _FuncRealT,  # func
    list[onp.Array1D[_FloatT]],  # xs
    list[onp.Array1D[_FloatT]],  # fs
    list[onp.Array1D[npc.floating]],  # args
    onp.AtLeast1D,  # shape
    _FloatT,  # xfat
    _ModuleT,  # xp
]: ...

# TODO: `_RichResult` subtype
def _loop(
    work: _ResT,
    callback: Callable[[_ResT], _Ignored],
    shape: Sequence[op.CanIndex],
    maxiter: int,
    func: Callable[[onp.Array[_ShapeT, _FloatT]], onp.ToComplexND],
    args: tuple[onp.ArrayND[npc.floating], ...],
    dtype: npc.inexact,
    pre_func_eval: Callable[[_ResT], onp.Array[_ShapeT, _FloatT]],
    post_func_eval: Callable[[onp.Array[_ShapeT, _FloatT], onp.Array[_ShapeT, npc.floating], _ResT], _Ignored],
    check_termination: Callable[[_ResT], onp.Array[_ShapeT, np.bool_]],
    post_termination_check: Callable[[_ResT], _Ignored],
    customize_result: Callable[[_ResT, _ToShapeT], tuple[int, ...]],
    res_work_pairs: Iterable[tuple[str, str]],
    xp: ModuleType,
    preserve_shape: bool | None = False,
) -> _ResT: ...

#
def _check_termination(
    work: _WorkT,
    res: Mapping[str, onp.Array[_ShapeT, _FloatT]],
    res_work_pairs: Iterable[tuple[str, str]],
    active: onp.Array[_ShapeT, npc.integer],
    check_termination: Callable[[_WorkT], onp.Array[_ShapeT, np.bool_]],
    preserve_shape: bool | None,
    xp: ModuleType,
) -> onp.Array1D[np.intp]: ...

#
def _update_active(
    work: Mapping[str, onp.Array[_ShapeT, _FloatT]],
    res: Mapping[str, onp.Array[_ShapeT, _FloatT]],
    res_work_pairs: Iterable[tuple[str, str]],
    active: onp.Array[_ShapeT, npc.integer],
    mask: onp.Array[_ShapeT, np.bool_] | None,
    preserve_shape: bool | None,
    xp: ModuleType,
) -> None: ...

#
def _prepare_result(
    work: Mapping[str, onp.Array[_ShapeT, _FloatT]],
    res: _ResT,
    res_work_pairs: Iterable[tuple[str, str]],
    active: onp.Array[_ShapeT, npc.integer],
    shape: _ToShapeT,
    customize_result: Callable[[_ResT, _ToShapeT], tuple[int, ...]],
    preserve_shape: bool | None,
    xp: ModuleType,
) -> _ResT: ...
