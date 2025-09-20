from collections.abc import Iterable
from typing import Any, Generic, Literal, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._ndbspline import _SolverFunc

__all__ = ["RegularGridInterpolator", "interpn"]

###

_CT_co = TypeVar("_CT_co", bound=np.float64 | np.complex128, default=np.float64, covariant=True)

_MethodReal: TypeAlias = Literal["linear", "nearest", "slinear", "cubic", "quintic"]
_Method: TypeAlias = Literal[_MethodReal, "pchip"]
_ToPoints: TypeAlias = Iterable[onp.ToFloat1D]

###

class RegularGridInterpolator(Generic[_CT_co]):
    grid: tuple[onp.ArrayND[_CT_co], ...]
    values: onp.ArrayND[_CT_co]
    method: _Method
    fill_value: float | None
    bounds_error: bool

    @overload
    def __init__(
        self: RegularGridInterpolator[np.float64],
        /,
        points: _ToPoints,
        values: onp.ToFloatND,
        method: _Method = "linear",
        bounds_error: onp.ToBool = True,
        fill_value: onp.ToFloat | None = ...,  # np.nan
        *,
        solver: _SolverFunc[npc.floating | npc.integer] | None = None,
        solver_args: tuple[object, ...] | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: RegularGridInterpolator[np.complex128],
        /,
        points: _ToPoints,
        values: onp.ToJustComplexND,
        method: _MethodReal = "linear",
        bounds_error: onp.ToBool = True,
        fill_value: onp.ToComplex | None = ...,  # np.nan
        *,
        solver: _SolverFunc[npc.number] | None = None,
        solver_args: tuple[object, ...] | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: RegularGridInterpolator[Any],
        /,
        points: _ToPoints,
        values: onp.ToComplexND,
        method: _MethodReal = "linear",
        bounds_error: onp.ToBool = True,
        fill_value: onp.ToComplex | None = ...,  # np.nan
        *,
        solver: _SolverFunc[npc.number] | None = None,
        solver_args: tuple[object, ...] | None = None,
    ) -> None: ...

    #
    def __call__(
        self, /, xi: onp.ToFloatND, method: _Method | None = None, *, nu: onp.ToJustInt1D | None = None
    ) -> onp.ArrayND[_CT_co]: ...

@overload
def interpn(
    points: _ToPoints,
    values: onp.ToFloatND,
    xi: onp.ToFloatND,
    method: _Method = "linear",
    bounds_error: onp.ToBool = True,
    fill_value: onp.ToFloat = ...,  # np.nan
) -> onp.ArrayND[np.float64]: ...
@overload
def interpn(
    points: _ToPoints,
    values: onp.ToJustComplex1D,
    xi: onp.ToFloatND,
    method: _Method = "linear",
    bounds_error: onp.ToBool = True,
    fill_value: onp.ToComplex = ...,  # np.nan
) -> onp.ArrayND[np.complex128]: ...
@overload
def interpn(
    points: _ToPoints,
    values: onp.ToComplex1D,
    xi: onp.ToFloatND,
    method: _Method = "linear",
    bounds_error: onp.ToBool = True,
    fill_value: onp.ToComplex = ...,  # np.nan
) -> onp.ArrayND[Any]: ...
