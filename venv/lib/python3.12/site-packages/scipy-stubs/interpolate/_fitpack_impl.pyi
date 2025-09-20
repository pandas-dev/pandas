from collections.abc import Sequence
from typing import Literal, LiteralString, TypeAlias, overload

import numpy as np
import optype.numpy as onp

__all__ = ["bisplev", "bisplrep", "insert", "spalde", "splantider", "splder", "splev", "splint", "splprep", "splrep", "sproot"]

_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]
_FloatND: TypeAlias = onp.ArrayND[np.float64]

_Task: TypeAlias = Literal[-1, 0, 1]
_Ext: TypeAlias = Literal[0, 1, 2, 3]

# This coincidentally works for both the univariate and bivariate cases
_ToTCK: TypeAlias = Sequence[onp.ToFloat1D | onp.ToFloat2D | int]

# `(t, c, k)`
_OutTCK1: TypeAlias = tuple[_Float1D, _Float1D, int]
# `[tx, ty, c, kx, ky]`
_OutTCK2: TypeAlias = list[_Float1D | _Float2D | int]
# `([t, c, k], u)`
_OutTCKU1: TypeAlias = tuple[list[_Float1D | list[_Float1D] | int], _Float1D]

###

# NOTE: The docs are incorrect about the return type of `splgrep`
@overload  # full_output: falsy = ...
def splprep(
    x: onp.ToFloat2D,
    w: onp.ToFloat1D | None = None,
    u: onp.ToFloat1D | None = None,
    ub: onp.ToFloat | None = None,
    ue: onp.ToFloat | None = None,
    k: int = 3,
    task: _Task = 0,
    s: onp.ToFloat | None = None,
    t: onp.ToFloat1D | None = None,
    full_output: onp.ToFalse = 0,
    nest: int | None = None,
    per: onp.ToBool = 0,
    quiet: onp.ToBool = 1,
) -> _OutTCKU1: ...
@overload  # full_output: truthy (positional)
def splprep(
    x: onp.ToFloat2D,
    w: onp.ToFloat1D | None,
    u: onp.ToFloat1D | None,
    ub: onp.ToFloat | None,
    ue: onp.ToFloat | None,
    k: int,
    task: _Task,
    s: onp.ToFloat | None,
    t: onp.ToFloat1D | None,
    full_output: onp.ToTrue,
    nest: int | None = None,
    per: onp.ToBool = 0,
    quiet: onp.ToBool = 1,
) -> tuple[_OutTCKU1, float, int, LiteralString]: ...
@overload  # full_output: truthy (keyword)
def splprep(
    x: onp.ToFloat2D,
    w: onp.ToFloat1D | None = None,
    u: onp.ToFloat1D | None = None,
    ub: onp.ToFloat | None = None,
    ue: onp.ToFloat | None = None,
    k: int = 3,
    task: _Task = 0,
    s: onp.ToFloat | None = None,
    t: onp.ToFloat1D | None = None,
    *,
    full_output: onp.ToTrue,
    nest: int | None = None,
    per: onp.ToBool = 0,
    quiet: onp.ToBool = 1,
) -> tuple[_OutTCKU1, float, int, LiteralString]: ...

#
@overload  # full_output: falsy = ...
def splrep(
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    w: onp.ToFloat1D | None = None,
    xb: onp.ToFloat | None = None,
    xe: onp.ToFloat | None = None,
    k: int = 3,
    task: _Task = 0,
    s: onp.ToFloat | None = None,
    t: onp.ToFloat1D | None = None,
    full_output: onp.ToFalse = 0,
    per: onp.ToBool = 0,
    quiet: onp.ToBool = 1,
) -> _OutTCK1: ...
@overload  # full_output: truthy (positional)
def splrep(
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    w: onp.ToFloat1D | None,
    xb: onp.ToFloat | None,
    xe: onp.ToFloat | None,
    k: int,
    task: _Task,
    s: onp.ToFloat | None,
    t: onp.ToFloat1D | None,
    full_output: onp.ToTrue,
    per: onp.ToBool = 0,
    quiet: onp.ToBool = 1,
) -> tuple[_OutTCK1, float, int, LiteralString]: ...
@overload  # full_output: truthy (keyword)
def splrep(
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    w: onp.ToFloat1D | None = None,
    xb: onp.ToFloat | None = None,
    xe: onp.ToFloat | None = None,
    k: int = 3,
    task: _Task = 0,
    s: onp.ToFloat | None = None,
    t: onp.ToFloat1D | None = None,
    *,
    full_output: onp.ToTrue,
    per: onp.ToBool = 0,
    quiet: onp.ToBool = 1,
) -> tuple[_OutTCK1, float, int, LiteralString]: ...

#
def splev(x: onp.ToFloatND, tck: _ToTCK, der: int = 0, ext: _Ext = 0) -> _FloatND: ...

#
@overload  # full_output: falsy
def splint(a: onp.ToFloat, b: onp.ToFloat, tck: _ToTCK, full_output: onp.ToFalse = 0) -> float | list[float]: ...
@overload  # full_output: truthy
def splint(a: onp.ToFloat, b: onp.ToFloat, tck: _ToTCK, full_output: onp.ToTrue) -> tuple[float | list[float], _Float1D]: ...

#
def sproot(tck: _ToTCK, mest: int = 10) -> _Float1D | list[_Float1D]: ...

#
@overload  # x: 1-d
def spalde(x: onp.ToFloatStrict1D, tck: _ToTCK) -> _Float1D: ...
@overload  # x: 2-d
def spalde(x: onp.ToFloatStrict2D, tck: _ToTCK) -> list[_Float1D]: ...
@overload  # x: {1,2}-d
def spalde(x: onp.ToFloat1D | onp.ToFloat2D, tck: _ToTCK) -> _Float1D | list[_Float1D]: ...

#
@overload  # full_output: falsy = ...
def bisplrep(
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    z: onp.ToFloat1D,
    w: onp.ToFloat1D | None = None,
    xb: onp.ToFloat | None = None,
    xe: onp.ToFloat | None = None,
    yb: onp.ToFloat | None = None,
    ye: onp.ToFloat | None = None,
    kx: int = 3,
    ky: int = 3,
    task: _Task = 0,
    s: onp.ToFloat | None = None,
    eps: onp.ToFloat = 1e-16,
    tx: onp.ToFloat1D | None = None,
    ty: onp.ToFloat1D | None = None,
    full_output: onp.ToFalse = 0,
    nxest: onp.ToFloat | None = None,
    nyest: onp.ToFloat | None = None,
    quiet: onp.ToBool = 1,
) -> _OutTCK2: ...
@overload  # full_output: truthy (positional)
def bisplrep(
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    z: onp.ToFloat1D,
    w: onp.ToFloat1D | None,
    xb: onp.ToFloat | None,
    xe: onp.ToFloat | None,
    yb: onp.ToFloat | None,
    ye: onp.ToFloat | None,
    kx: int,
    ky: int,
    task: _Task,
    s: onp.ToFloat | None,
    eps: onp.ToFloat,
    tx: onp.ToFloat1D | None,
    ty: onp.ToFloat1D | None,
    full_output: onp.ToTrue,
    nxest: onp.ToFloat | None = None,
    nyest: onp.ToFloat | None = None,
    quiet: onp.ToBool = 1,
) -> tuple[_OutTCK2, float, int, LiteralString]: ...
@overload  # full_output: truthy (keyword)
def bisplrep(
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    z: onp.ToFloat1D,
    w: onp.ToFloat1D | None = None,
    xb: onp.ToFloat | None = None,
    xe: onp.ToFloat | None = None,
    yb: onp.ToFloat | None = None,
    ye: onp.ToFloat | None = None,
    kx: int = 3,
    ky: int = 3,
    task: _Task = 0,
    s: onp.ToFloat | None = None,
    eps: onp.ToFloat = 1e-16,
    tx: onp.ToFloat1D | None = None,
    ty: onp.ToFloat1D | None = None,
    *,
    full_output: onp.ToTrue,
    nxest: onp.ToFloat | None = None,
    nyest: onp.ToFloat | None = None,
    quiet: onp.ToBool = 1,
) -> tuple[_OutTCK2, float, int, LiteralString]: ...

# requires `len(tck) == 5`
def bisplev(x: onp.ToFloat1D, y: onp.ToFloat1D, tck: _ToTCK, dx: int = 0, dy: int = 0) -> _Float2D: ...
def dblint(xa: onp.ToFloat, xb: onp.ToFloat, ya: onp.ToFloat, yb: onp.ToFloat, tck: _ToTCK) -> float: ...

# requires `len(tck) == 3`
def insert(x: onp.ToFloat, tck: _ToTCK, m: int = 1, per: onp.ToBool = 0) -> _OutTCK1: ...
def splder(tck: _ToTCK, n: int = 1) -> _OutTCK1: ...
def splantider(tck: _ToTCK, n: int = 1) -> _OutTCK1: ...
