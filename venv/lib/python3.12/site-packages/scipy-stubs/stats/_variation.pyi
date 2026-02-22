from collections.abc import Sequence
from typing import Any, Literal, Never, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._typing import NanPolicy

_FloatingT = TypeVar("_FloatingT", bound=npc.floating)
_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])

# workaround for https://github.com/microsoft/pyright/issues/10232
_JustAnyShape: TypeAlias = tuple[Never, ...]

_co_integer: TypeAlias = npc.integer | np.bool_  # noqa: PYI042

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
def variation(
    a: onp.ArrayND[_co_integer, _ShapeT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[True],
) -> onp.ArrayND[np.float64, _ShapeT]: ...
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
def variation(
    a: onp.ArrayND[_FloatingT, _JustAnyShape],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> onp.ArrayND[_FloatingT] | Any: ...
@overload  # 1d ~T:floating
def variation(
    a: onp.ToArrayStrict1D[_FloatingT, _FloatingT],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> _FloatingT: ...
@overload  # 2d ~T:floating
def variation(
    a: onp.ToArrayStrict2D[_FloatingT, _FloatingT],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> onp.Array1D[_FloatingT]: ...
@overload  # nd ~T:floating
def variation(
    a: onp.ToArrayND[_FloatingT, _FloatingT],
    axis: int = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> onp.ArrayND[_FloatingT] | Any: ...
@overload  # nd ~T:floating, axis=None
def variation(
    a: onp.ToArrayND[_FloatingT, _FloatingT],
    axis: None,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[False] = False,
) -> _FloatingT: ...
@overload  # T:nd ~T:floating, keepdims=True
def variation(
    a: onp.ArrayND[_FloatingT, _ShapeT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[True],
) -> onp.ArrayND[_FloatingT, _ShapeT]: ...
@overload  # nd ~T:floating, keepdims=True
def variation(
    a: onp.ToArrayND[_FloatingT, _FloatingT],
    axis: int | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: int = 0,
    *,
    keepdims: Literal[True],
) -> onp.ArrayND[_FloatingT]: ...
@overload  # fallback
def variation(
    a: onp.ToFloatND, axis: int | None = 0, nan_policy: NanPolicy = "propagate", ddof: int = 0, *, keepdims: bool = False
) -> onp.ArrayND[np.float64 | Any] | Any: ...
