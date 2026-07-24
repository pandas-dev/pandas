from collections.abc import Iterable
from typing import Any, Literal, SupportsIndex, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

###

type _Mode = Literal["nearest", "wrap", "reflect", "grid-mirror", "mirror", "constant", "grid-wrap", "grid-constant"]
type _ModeCode = Literal[0, 1, 2, 3, 4, 5, 6]

###

def _extend_mode_to_code(mode: _Mode, is_filter: bool = False) -> _ModeCode: ...

#
@overload
def _normalize_sequence(input: str, rank: int) -> list[str]: ...
@overload
def _normalize_sequence[T](input: Iterable[T], rank: int) -> list[T]: ...
@overload
def _normalize_sequence[T](input: T, rank: int) -> list[T]: ...

#
@overload
def _get_output(
    output: type[complex | npc.number] | np.dtype[npc.number] | str | onp.ToComplexND | None,
    input: onp.ArrayND[npc.number],
    shape: tuple[int, ...] | None = None,
    complex_output: Literal[False] = False,
) -> onp.ArrayND[npc.number]: ...
@overload
def _get_output(
    output: type[complex | npc.complexfloating] | np.dtype[npc.complexfloating] | str | onp.ToJustComplexND | None,
    input: onp.ArrayND[npc.complexfloating],
    shape: tuple[int, ...] | None = None,
    *,
    complex_output: Literal[True],
) -> onp.ArrayND[npc.complexfloating]: ...
@overload
def _get_output(
    output: type[complex | npc.complexfloating] | np.dtype[npc.complexfloating] | str | onp.ToJustComplexND | None,
    input: onp.ArrayND[npc.complexfloating],
    shape: tuple[int, ...] | None,
    complex_output: Literal[True],
) -> onp.ArrayND[npc.complexfloating]: ...

#
def _check_axes(axes: SupportsIndex | Iterable[SupportsIndex], ndim: int) -> tuple[int, ...]: ...

#
@overload
def _skip_if_dtype(arg: str | np.dtype[Any] | type[np.generic]) -> None: ...
@overload
def _skip_if_dtype[T](arg: T) -> T | None: ...
