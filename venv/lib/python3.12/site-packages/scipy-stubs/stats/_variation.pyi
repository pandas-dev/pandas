from collections.abc import Sequence
from typing import overload
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc
from numpy._typing import _ArrayLike, _NestedSequence

from ._typing import NanPolicy

_SCT_fc = TypeVar("_SCT_fc", bound=npc.inexact)

###

# NOTE: This demonstrates the ridiculous complexity that's required to properly annotate this simple function with "array-likes".
# NOTE: Shape-typing hasn't even been included, as that would require even more overloads.

# sequences of `builtins.float`, that implicitly (and inevitably) also cover `builtins.int` and `builtins.bool`
@overload
def variation(
    a: Sequence[float],
    axis: op.CanIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: onp.ToInt = 0,
    *,
    keepdims: onp.ToFalse = False,
) -> np.float64: ...
@overload
def variation(
    a: onp.ToIntND | _NestedSequence[float],
    axis: None,
    nan_policy: NanPolicy = "propagate",
    ddof: onp.ToInt = 0,
    *,
    keepdims: onp.ToFalse = False,
) -> np.float64: ...
@overload
def variation(
    a: onp.ToIntND | _NestedSequence[float],
    axis: op.CanIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: onp.ToInt = 0,
    *,
    keepdims: onp.ToTrue,
) -> onp.ArrayND[np.float64]: ...
@overload
def variation(
    a: onp.ToIntND | _NestedSequence[float],
    axis: op.CanIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: onp.ToInt = 0,
    *,
    keepdims: onp.ToBool = False,
) -> np.float64 | onp.ArrayND[np.float64]: ...

# array-like input with known floating- or complex-floating dtypes
@overload
def variation(
    a: _ArrayLike[_SCT_fc], axis: None, nan_policy: NanPolicy = "propagate", ddof: onp.ToInt = 0, *, keepdims: onp.ToFalse = False
) -> _SCT_fc: ...
@overload
def variation(
    a: _ArrayLike[_SCT_fc],
    axis: op.CanIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: onp.ToInt = 0,
    *,
    keepdims: onp.ToTrue,
) -> onp.ArrayND[_SCT_fc]: ...
@overload
def variation(
    a: _ArrayLike[_SCT_fc],
    axis: op.CanIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: onp.ToInt = 0,
    *,
    keepdims: onp.ToBool = False,
) -> _SCT_fc | onp.ArrayND[_SCT_fc]: ...

# sequences of `builtin.complex`, which behave as if `float <: complex` and therefore "overlaps" with the `builtins.float`
# overloads, hence the `complex128 | float64` returns
@overload
def variation(
    a: Sequence[complex],
    axis: op.CanIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: onp.ToInt = 0,
    *,
    keepdims: onp.ToFalse = False,
) -> np.complex128 | np.float64: ...
@overload
def variation(
    a: _NestedSequence[complex],
    axis: None,
    nan_policy: NanPolicy = "propagate",
    ddof: onp.ToInt = 0,
    *,
    keepdims: onp.ToFalse = False,
) -> np.complex128 | np.float64: ...
@overload
def variation(
    a: _NestedSequence[complex],
    axis: op.CanIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: onp.ToInt = 0,
    *,
    keepdims: onp.ToTrue,
) -> onp.ArrayND[np.complex128 | np.float64]: ...
@overload
def variation(
    a: _NestedSequence[complex],
    axis: op.CanIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: onp.ToInt = 0,
    *,
    keepdims: onp.ToBool = False,
) -> np.complex128 | np.float64 | onp.ArrayND[np.complex128 | np.float64]: ...

# catch-all in case of broad gradual types
@overload
def variation(
    a: onp.ToComplexND,
    axis: op.CanIndex | None = 0,
    nan_policy: NanPolicy = "propagate",
    ddof: onp.ToInt = 0,
    *,
    keepdims: onp.ToBool = False,
) -> npc.inexact | onp.ArrayND[npc.inexact]: ...
