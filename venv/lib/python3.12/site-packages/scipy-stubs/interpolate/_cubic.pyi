from typing import Any, Generic, Literal, Never, TypeAlias, overload
from typing_extensions import TypeVar, override

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._interpolate import PPoly

__all__ = ["Akima1DInterpolator", "CubicHermiteSpline", "CubicSpline", "PchipInterpolator", "pchip_interpolate"]

_T = TypeVar("_T")
_CT = TypeVar("_CT", bound=np.float64 | np.complex128)
_CT_co = TypeVar("_CT_co", bound=np.float64 | np.complex128, default=np.float64, covariant=True)
_AxisT = TypeVar("_AxisT", bound=_ToAxis)

_Tuple2: TypeAlias = tuple[_T, _T]
_ToAxis: TypeAlias = int | npc.integer

_Akima1DMethod: TypeAlias = Literal["akima", "makima"]
_Extrapolate: TypeAlias = Literal["periodic"] | bool
_CubicBCName: TypeAlias = Literal["not-a-knot", "clamped", "natural"]
_CubicBCOrder: TypeAlias = Literal[1, 2]
_CubicBCType: TypeAlias = Literal[_CubicBCName, "periodic"] | _Tuple2[_CubicBCName | tuple[_CubicBCOrder, onp.ToComplexND]]

_PreparedInput: TypeAlias = tuple[
    onp.Array1D[np.float64],  # x
    onp.Array1D[np.float64],  # dx
    onp.ArrayND[_CT],  # y
    _AxisT,  # axis
    onp.ArrayND[_CT],  # dydx
]

###

class CubicHermiteSpline(PPoly[_CT_co]):
    @overload
    def __init__(
        self: CubicHermiteSpline[np.float64],
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatND,
        dydx: onp.ToFloatND,
        axis: _ToAxis = 0,
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: CubicHermiteSpline[np.complex128],
        /,
        x: onp.ToFloat1D,
        y: onp.ToJustComplexND,
        dydx: onp.ToComplexND,
        axis: _ToAxis = 0,
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: CubicHermiteSpline[Any],
        /,
        x: onp.ToFloat1D,
        y: onp.ToComplexND,
        dydx: onp.ToComplexND,
        axis: _ToAxis = 0,
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...

class PchipInterpolator(CubicHermiteSpline[np.float64]):
    def __init__(self, /, x: onp.ToFloat1D, y: onp.ToFloatND, axis: _ToAxis = 0, extrapolate: bool | None = None) -> None: ...

class Akima1DInterpolator(CubicHermiteSpline[np.float64]):
    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatND,
        axis: _ToAxis = 0,
        *,
        method: _Akima1DMethod = "akima",
        extrapolate: onp.ToBool | None = None,
    ) -> None: ...

    # the following (class)methods will raise `NotImplementedError` when called
    @override
    def extend(self, /, c: Never, x: Never, right: bool = True) -> Never: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]
    @classmethod
    @override
    def from_spline(cls, tck: Never, extrapolate: None = None) -> Never: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]
    @classmethod
    @override
    def from_bernstein_basis(cls, bp: Never, extrapolate: None = None) -> Never: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]

class CubicSpline(CubicHermiteSpline[_CT_co], Generic[_CT_co]):
    @overload
    def __init__(
        self: CubicSpline[np.float64],
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatND,
        axis: _ToAxis = 0,
        bc_type: _CubicBCType = "not-a-knot",
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: CubicSpline[np.complex128],
        /,
        x: onp.ToFloat1D,
        y: onp.ToJustComplexND,
        axis: _ToAxis = 0,
        bc_type: _CubicBCType = "not-a-knot",
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self: CubicSpline[Any],
        /,
        x: onp.ToFloat1D,
        y: onp.ToComplexND,
        axis: _ToAxis = 0,
        bc_type: _CubicBCType = "not-a-knot",
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...

@overload
def pchip_interpolate(
    xi: onp.ToFloat1D, yi: onp.ToFloat1D, x: onp.ToFloat, der: onp.ToInt = 0, axis: _ToAxis = 0
) -> np.float64: ...
@overload
def pchip_interpolate(
    xi: onp.ToFloat1D, yi: onp.ToFloat1D, x: onp.ToFloat1D, der: onp.ToInt | onp.ToInt1D = 0, axis: _ToAxis = 0
) -> onp.ArrayND[np.float64]: ...

# undocumented
@overload
def prepare_input(
    x: onp.ToFloat1D, y: onp.ToFloatND, axis: _AxisT, dydx: onp.ToFloatND | None = None
) -> _PreparedInput[np.float64, _AxisT]: ...
@overload
def prepare_input(
    x: onp.ToFloat1D, y: onp.ToJustComplexND, axis: _AxisT, dydx: onp.ToComplexND | None = None
) -> _PreparedInput[np.complex128, _AxisT]: ...
@overload
def prepare_input(
    x: onp.ToFloat1D, y: onp.ToComplexND, axis: _AxisT, dydx: onp.ToComplexND | None = None
) -> _PreparedInput[Any, _AxisT]: ...
