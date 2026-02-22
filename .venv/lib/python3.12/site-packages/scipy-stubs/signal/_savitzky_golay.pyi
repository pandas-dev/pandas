from _typeshed import Incomplete
from collections.abc import Sequence
from types import ModuleType
from typing import Any, Literal, TypeAlias, TypeVar, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

###

# we use `tuple[Any, ...]` instead of `tuple[int, ...]` here to work around bugs in both mypy and pyright (sigh...)
_ShapeT = TypeVar("_ShapeT", bound=tuple[Any, ...])

_Mode: TypeAlias = Literal["mirror", "constant", "nearest", "wrap", "interp"]
# savgol_filter casts everything except for `float32` to `float64` (kinda weird, but ok)
_ToFloat64Savgol: TypeAlias = npc.floating64 | npc.floating16 | npc.floating80 | npc.integer | np.bool_

###

@overload
def savgol_coeffs(
    window_length: int,
    polyorder: int,
    deriv: int = 0,
    delta: float = 1.0,
    pos: int | None = None,
    use: Literal["conv", "dot"] = "conv",
    *,
    xp: None = None,
    device: None = None,
) -> onp.Array1D[np.float64]: ...
@overload
def savgol_coeffs(
    window_length: int,
    polyorder: int,
    deriv: int = 0,
    delta: float = 1.0,
    pos: int | None = None,
    use: Literal["conv", "dot"] = "conv",
    *,
    xp: ModuleType,
    device: object | None = None,
) -> Incomplete: ...

# NOTE: On `numpy<2.1` the shape-type of `ndarray` was invariant, and pyright has a bug that incorrectly causes it to report
# the overlapping overloads as incompatible, even though the shape-type defaults to `tuple[Any, ...]` and is therefore compatible
# with all other shape-types. Thus, we're forced to disable pyright's overlapping overload check here.

# pyright: reportOverlappingOverload=false

#
@overload  # f64, shape known
def savgol_filter(
    x: onp.ArrayND[_ToFloat64Savgol, _ShapeT],
    window_length: int,
    polyorder: int,
    deriv: int = 0,
    delta: float = 1.0,
    axis: int = -1,
    mode: _Mode = "interp",
    cval: float = 0.0,
) -> onp.ArrayND[np.float64, _ShapeT]: ...
@overload  # float, 1d
def savgol_filter(
    x: Sequence[float],
    window_length: int,
    polyorder: int,
    deriv: int = 0,
    delta: float = 1.0,
    axis: int = -1,
    mode: _Mode = "interp",
    cval: float = 0.0,
) -> onp.Array1D[np.float64]: ...
@overload  # float, 2d
def savgol_filter(
    x: Sequence[Sequence[float]],
    window_length: int,
    polyorder: int,
    deriv: int = 0,
    delta: float = 1.0,
    axis: int = -1,
    mode: _Mode = "interp",
    cval: float = 0.0,
) -> onp.Array2D[np.float64]: ...
@overload  # float, 3d
def savgol_filter(
    x: Sequence[Sequence[Sequence[float]]],
    window_length: int,
    polyorder: int,
    deriv: int = 0,
    delta: float = 1.0,
    axis: int = -1,
    mode: _Mode = "interp",
    cval: float = 0.0,
) -> onp.Array3D[np.float64]: ...
@overload  # f64, shape unknown
def savgol_filter(
    x: onp.ToArrayND[float, _ToFloat64Savgol],
    window_length: int,
    polyorder: int,
    deriv: int = 0,
    delta: float = 1.0,
    axis: int = -1,
    mode: _Mode = "interp",
    cval: float = 0.0,
) -> onp.ArrayND[np.float64]: ...
@overload  # f32, shape known
def savgol_filter(
    x: onp.ArrayND[npc.floating32, _ShapeT],
    window_length: int,
    polyorder: int,
    deriv: int = 0,
    delta: float = 1.0,
    axis: int = -1,
    mode: _Mode = "interp",
    cval: float = 0.0,
) -> onp.ArrayND[np.float32, _ShapeT]: ...
@overload  # f32, shape unknown
def savgol_filter(
    x: onp.ToJustFloat32_ND,
    window_length: int,
    polyorder: int,
    deriv: int = 0,
    delta: float = 1.0,
    axis: int = -1,
    mode: _Mode = "interp",
    cval: float = 0.0,
) -> onp.ArrayND[np.float32]: ...
@overload  # fallback
def savgol_filter(
    x: onp.ToComplexND,
    window_length: int,
    polyorder: int,
    deriv: int = 0,
    delta: float = 1.0,
    axis: int = -1,
    mode: _Mode = "interp",
    cval: float = 0.0,
) -> onp.ArrayND[np.float64 | np.float32]: ...
