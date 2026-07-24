from _typeshed import Unused
from collections.abc import Callable, Iterable, Mapping, Sequence
from types import ModuleType
from typing import Any, Concatenate, Final, Literal, SupportsIndex

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._util import _RichResult

###

_ESIGNERR: Final = -1
_ECONVERR: Final = -2
_EVALUEERR: Final = -3
_ECALLBACK: Final = -4
_EINPUTERR: Final = -5
_ECONVERGED: Final = 0
_EINPROGRESS: Final = 1

# TODO: complex
def _initialize[FuncT: Callable[Concatenate[onp.ArrayND[np.float64], ...], object], ModuleT: ModuleType, FloatT: npc.floating](
    func: FuncT,
    xs: Sequence[onp.ToFloat1D],
    args: tuple[onp.ToFloat1D, ...],
    kwargs: dict[str, Any] | None = None,
    complex_ok: Literal[False] = False,
    preserve_shape: bool | None = None,
    xp: ModuleT | None = None,
) -> tuple[
    FuncT,  # func
    list[onp.Array1D[FloatT]],  # xs
    list[onp.Array1D[FloatT]],  # fs
    list[onp.Array1D[npc.floating]],  # args
    onp.AtLeast1D,  # shape
    FloatT,  # xfat
    ModuleT,  # xp
]: ...
def _loop[ResT: _RichResult[Any], ShapeT: tuple[int, ...], FloatT: npc.floating](
    work: ResT,
    callback: Callable[[ResT], Unused],
    shape: Sequence[SupportsIndex],
    maxiter: int,
    func: Callable[[onp.Array[ShapeT, FloatT]], onp.ToComplexND],
    args: tuple[onp.ArrayND[npc.floating], ...],
    dtype: npc.inexact,
    pre_func_eval: Callable[[ResT], onp.Array[ShapeT, FloatT]],
    post_func_eval: Callable[[onp.Array[ShapeT, FloatT], onp.Array[ShapeT, npc.floating], ResT], Unused],
    check_termination: Callable[[ResT], onp.Array[ShapeT, np.bool]],
    post_termination_check: Callable[[ResT], Unused],
    customize_result: Callable[[ResT, tuple[int, ...]], tuple[int, ...]],
    res_work_pairs: Iterable[tuple[str, str]],
    xp: ModuleType,
    preserve_shape: bool | None = False,
) -> ResT: ...

#
def _check_termination[WorkT: Mapping[str, Any], ShapeT: tuple[int, ...], FloatT: npc.floating](
    work: WorkT,
    res: Mapping[str, onp.Array[ShapeT, FloatT]],
    res_work_pairs: Iterable[tuple[str, str]],
    active: onp.Array[ShapeT, npc.integer],
    check_termination: Callable[[WorkT], onp.Array[ShapeT, np.bool]],
    preserve_shape: bool | None,
    xp: ModuleType,
) -> onp.Array1D[np.intp]: ...

#
def _update_active[ShapeT: tuple[int, ...], FloatT: npc.floating](
    work: Mapping[str, onp.Array[ShapeT, FloatT]],
    res: Mapping[str, onp.Array[ShapeT, FloatT]],
    res_work_pairs: Iterable[tuple[str, str]],
    active: onp.Array[ShapeT, npc.integer],
    mask: onp.Array[ShapeT, np.bool] | None,
    preserve_shape: bool | None,
    xp: ModuleType,
) -> None: ...

#
def _prepare_result[
    ResT: _RichResult[Any],
    ShapeT: tuple[int, ...],
    ToShapeT: SupportsIndex | tuple[SupportsIndex, ...],
    FloatT: npc.floating,
](
    work: Mapping[str, onp.Array[ShapeT, FloatT]],
    res: ResT,
    res_work_pairs: Iterable[tuple[str, str]],
    active: onp.Array[ShapeT, npc.integer],
    shape: ToShapeT,
    customize_result: Callable[[ResT, ToShapeT], tuple[int, ...]],
    preserve_shape: bool | None,
    xp: ModuleType,
) -> ResT: ...
