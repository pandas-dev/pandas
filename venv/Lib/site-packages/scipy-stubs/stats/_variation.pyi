from collections.abc import Sequence
from typing import Any, Literal, Never, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._typing import NanPolicy

###

# workaround for https://github.com/microsoft/pyright/issues/10232
type _JustAnyShape = tuple[Never, ...]

type _co_integer = npc.integer | np.bool  # noqa: PYI042

###

# NOTE: Pyright reports false positives for some overloads involving gradual shape types on numpy<2.1
# pyright: reportOverlappingOverload=false

# NOTE: Hilariously, mypy also reports false positive `overload-overlap` errors on numpy<2.1,
# but for completely different overloads.
# mypy: disable-error-code=overload-overlap

@overload  # ?d +int: inexact
def variation(
    a: onp.ArrayND[_co_integer, _JustAnyShape],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> onp.ArrayND[np.float64] | Any: ...
@overload  # 1d +float64
def variation(
    a: onp.ToArrayStrict1D[float, _co_integer],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload  # 2d +float64
def variation(
    a: onp.ToArrayStrict2D[float, _co_integer],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> onp.Array1D[np.float64]: ...
@overload  # >1d +float64
def variation(
    a: Sequence[onp.SequenceND[float]],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> onp.ArrayND[np.float64]: ...
@overload  # nd +float64
def variation(
    a: onp.ToArrayND[float, _co_integer],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> onp.ArrayND[np.float64] | Any: ...
@overload  # nd +float64, axis=None
def variation(
    a: onp.ToArrayND[float, _co_integer],
    axis: None,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> np.float64: ...
@overload  # T:nd +float64, keepdims=True
def variation[ShapeT: tuple[int, ...]](
    a: onp.ArrayND[_co_integer, ShapeT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[True],
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # nd +float64, keepdims=True
def variation(
    a: onp.ToArrayND[float, _co_integer],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[True],
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d ~T:floating
def variation[FloatingT: npc.floating](
    a: onp.ArrayND[FloatingT, _JustAnyShape],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> onp.ArrayND[FloatingT] | Any: ...
@overload  # 1d ~T:floating
def variation[FloatingT: npc.floating](
    a: onp.ToArrayStrict1D[FloatingT, FloatingT],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> FloatingT: ...
@overload  # 2d ~T:floating
def variation[FloatingT: npc.floating](
    a: onp.ToArrayStrict2D[FloatingT, FloatingT],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> onp.Array1D[FloatingT]: ...
@overload  # nd ~T:floating
def variation[FloatingT: npc.floating](
    a: onp.ToArrayND[FloatingT, FloatingT],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> onp.ArrayND[FloatingT] | Any: ...
@overload  # nd ~T:floating, axis=None
def variation[FloatingT: npc.floating](
    a: onp.ToArrayND[FloatingT, FloatingT],
    axis: None,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> FloatingT: ...
@overload  # T:nd ~T:floating, keepdims=True
def variation[FloatingT: npc.floating, ShapeT: tuple[int, ...]](
    a: onp.ArrayND[FloatingT, ShapeT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[True],
) -> onp.ArrayND[FloatingT, ShapeT]: ...
@overload  # nd ~T:floating, keepdims=True
def variation[FloatingT: npc.floating](
    a: onp.ToArrayND[FloatingT, FloatingT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[True],
) -> onp.ArrayND[FloatingT]: ...
@overload  # fallback
def variation(
    a: onp.ToFloatND, axis: int | None = 0, nan_policy: NanPolicy = "propagate", ddof: int = 0, *, keepdims: bool = False
) -> onp.ArrayND[np.float64 | Any] | Any: ...
