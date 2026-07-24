from _typeshed import Incomplete
from types import ModuleType
from typing import Any, ClassVar, Generic, Literal, Never, overload, override
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._interpolate import PPoly

__all__ = ["Akima1DInterpolator", "CubicHermiteSpline", "CubicSpline", "PchipInterpolator", "pchip_interpolate"]

###

type _Tuple2[T] = tuple[T, T]
type _ToAxis = int | npc.integer

type _Akima1DMethod = Literal["akima", "makima"]
type _Extrapolate = Literal["periodic"] | bool
type _CubicBCName = Literal["not-a-knot", "clamped", "natural"]
type _CubicBCOrder = Literal[1, 2]
type _CubicBCType = Literal[_CubicBCName, "periodic"] | _Tuple2[_CubicBCName | tuple[_CubicBCOrder, onp.ToComplexND]]

type _PreparedInput[CT: np.float64 | np.complex128, AxisT: _ToAxis] = tuple[
    onp.Array1D[np.float64],  # x
    onp.Array1D[np.float64],  # dx
    onp.ArrayND[CT],  # y
    AxisT,  # axis
    onp.ArrayND[CT],  # dydx
]

# workaround for https://github.com/microsoft/pyright/issues/10232
type _JustAnyShape = tuple[Never, Never, Never]

_CT_co = TypeVar("_CT_co", bound=np.float64 | np.complex128, default=np.float64, covariant=True)
_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)

###

class CubicHermiteSpline(PPoly[_CT_co, _ShapeT_co], Generic[_CT_co, _ShapeT_co]):
    @overload  # ?d real
    def __init__(
        self: CubicHermiteSpline[np.float64],
        /,
        x: onp.ToFloat1D,
        y: onp.ArrayND[npc.floating | npc.integer | np.bool, _JustAnyShape],
        dydx: onp.ToFloatND,
        axis: _ToAxis = 0,
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # ?d complex
    def __init__(
        self: CubicHermiteSpline[np.complex128],
        /,
        x: onp.ToFloat1D,
        y: onp.ArrayND[npc.complexfloating, _JustAnyShape],
        dydx: onp.ToComplexND,
        axis: _ToAxis = 0,
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # 0d real
    def __init__(
        self: CubicHermiteSpline[np.float64, tuple[()]],
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatStrict1D,
        dydx: onp.ToFloat1D,
        axis: _ToAxis = 0,
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # 0d complex
    def __init__(
        self: CubicHermiteSpline[np.complex128, tuple[()]],
        /,
        x: onp.ToFloat1D,
        y: onp.ToJustComplexStrict1D,
        dydx: onp.ToComplex1D,
        axis: _ToAxis = 0,
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # 1d real
    def __init__(
        self: CubicHermiteSpline[np.float64, tuple[int]],
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatStrict2D,
        dydx: onp.ToFloat2D,
        axis: _ToAxis = 0,
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # 1d complex
    def __init__(
        self: CubicHermiteSpline[np.complex128, tuple[int]],
        /,
        x: onp.ToFloat1D,
        y: onp.ToJustComplexStrict2D,
        dydx: onp.ToComplex2D,
        axis: _ToAxis = 0,
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # Nd real
    def __init__(
        self: CubicHermiteSpline[np.float64],
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatND,
        dydx: onp.ToFloatND,
        axis: _ToAxis = 0,
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # Nd complex
    def __init__(
        self: CubicHermiteSpline[np.complex128],
        /,
        x: onp.ToFloat1D,
        y: onp.ToJustComplexND,
        dydx: onp.ToComplexND,
        axis: _ToAxis = 0,
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # fallback
    def __init__(
        self: CubicHermiteSpline[Any],
        /,
        x: onp.ToFloat1D,
        y: onp.ToComplexND,
        dydx: onp.ToComplexND,
        axis: _ToAxis = 0,
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...

class PchipInterpolator(CubicHermiteSpline[np.float64, _ShapeT_co], Generic[_ShapeT_co]):
    # pyrefly: ignore [bad-override]
    __class_getitem__: ClassVar[None] = None  # type:ignore[assignment]  # pyright:ignore[reportIncompatibleMethodOverride]

    @overload  # ?d
    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        y: onp.ArrayND[npc.floating | npc.integer | np.bool, _JustAnyShape],
        axis: _ToAxis = 0,
        extrapolate: bool | None = None,
    ) -> None: ...
    @overload  # 0d
    def __init__(
        self: PchipInterpolator[tuple[()]],
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatStrict1D,
        axis: _ToAxis = 0,
        extrapolate: bool | None = None,
    ) -> None: ...
    @overload  # 1d
    def __init__(
        self: PchipInterpolator[tuple[int]],
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatStrict2D,
        axis: _ToAxis = 0,
        extrapolate: bool | None = None,
    ) -> None: ...
    @overload  # Nd
    def __init__(
        self: PchipInterpolator[tuple[Any, ...]],
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatND,
        axis: _ToAxis = 0,
        extrapolate: bool | None = None,
    ) -> None: ...

class Akima1DInterpolator(CubicHermiteSpline[np.float64, _ShapeT_co], Generic[_ShapeT_co]):
    # pyrefly: ignore [bad-override]
    __class_getitem__: ClassVar[None] = None  # type:ignore[assignment]  # pyright:ignore[reportIncompatibleMethodOverride]

    @overload  # ?d
    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        y: onp.ArrayND[npc.floating | npc.integer | np.bool, _JustAnyShape],
        axis: _ToAxis = 0,
        *,
        method: _Akima1DMethod = "akima",
        extrapolate: bool | None = None,
    ) -> None: ...
    @overload  # 0d
    def __init__(
        self: Akima1DInterpolator[tuple[()]],
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatStrict1D,
        axis: _ToAxis = 0,
        *,
        method: _Akima1DMethod = "akima",
        extrapolate: bool | None = None,
    ) -> None: ...
    @overload  # 1d
    def __init__(
        self: Akima1DInterpolator[tuple[int]],
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatStrict2D,
        axis: _ToAxis = 0,
        *,
        method: _Akima1DMethod = "akima",
        extrapolate: bool | None = None,
    ) -> None: ...
    @overload  # Nd
    def __init__(
        self: Akima1DInterpolator[tuple[Any, ...]],
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatND,
        axis: _ToAxis = 0,
        *,
        method: _Akima1DMethod = "akima",
        extrapolate: bool | None = None,
    ) -> None: ...

    # the following (class)methods will raise `NotImplementedError` when called
    @override
    # pyrefly: ignore [bad-override]
    def extend(self, /, c: Never, x: Never, right: bool = True) -> Never: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]  # ty: ignore[invalid-method-override]
    @classmethod
    @override
    # pyrefly: ignore [bad-override]
    def from_spline(cls, tck: Never, extrapolate: None = None) -> Never: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]  # ty: ignore[invalid-method-override]
    @classmethod
    @override
    # pyrefly: ignore [bad-override]
    def from_bernstein_basis(cls, bp: Never, extrapolate: None = None) -> Never: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]  # ty: ignore[invalid-method-override]

class CubicSpline(CubicHermiteSpline[_CT_co, _ShapeT_co], Generic[_CT_co, _ShapeT_co]):
    @overload  # ?d real
    def __init__(
        self: CubicSpline[np.float64],
        /,
        x: onp.ToFloat1D,
        y: onp.ArrayND[npc.floating | npc.integer | np.bool, _JustAnyShape],
        axis: _ToAxis = 0,
        bc_type: _CubicBCType = "not-a-knot",
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # ?d complex
    def __init__(
        self: CubicSpline[np.complex128],
        /,
        x: onp.ToFloat1D,
        y: onp.ArrayND[npc.complexfloating, _JustAnyShape],
        axis: _ToAxis = 0,
        bc_type: _CubicBCType = "not-a-knot",
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # 0d real
    def __init__(
        self: CubicSpline[np.float64, tuple[()]],
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatStrict1D,
        axis: _ToAxis = 0,
        bc_type: _CubicBCType = "not-a-knot",
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # 0d complex
    def __init__(
        self: CubicSpline[np.complex128, tuple[()]],
        /,
        x: onp.ToFloat1D,
        y: onp.ToJustComplexStrict1D,
        axis: _ToAxis = 0,
        bc_type: _CubicBCType = "not-a-knot",
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # 1d real
    def __init__(
        self: CubicSpline[np.float64, tuple[int]],
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatStrict2D,
        axis: _ToAxis = 0,
        bc_type: _CubicBCType = "not-a-knot",
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # 1d complex
    def __init__(
        self: CubicSpline[np.complex128, tuple[int]],
        /,
        x: onp.ToFloat1D,
        y: onp.ToJustComplexStrict2D,
        axis: _ToAxis = 0,
        bc_type: _CubicBCType = "not-a-knot",
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # Nd real
    def __init__(
        self: CubicSpline[np.float64],
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloatND,
        axis: _ToAxis = 0,
        bc_type: _CubicBCType = "not-a-knot",
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # Nd complex
    def __init__(
        self: CubicSpline[np.complex128],
        /,
        x: onp.ToFloat1D,
        y: onp.ToJustComplexND,
        axis: _ToAxis = 0,
        bc_type: _CubicBCType = "not-a-knot",
        extrapolate: _Extrapolate | None = None,
    ) -> None: ...
    @overload  # fallback
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
def prepare_input[AxisT: _ToAxis](
    x: onp.ToFloat1D, y: onp.ToFloatND, axis: AxisT, dydx: onp.ToFloatND | None = None, xp: None = None
) -> _PreparedInput[np.float64, AxisT]: ...
@overload
def prepare_input[AxisT: _ToAxis](
    x: onp.ToFloat1D, y: onp.ToJustComplexND, axis: AxisT, dydx: onp.ToComplexND | None = None, xp: None = None
) -> _PreparedInput[np.complex128, AxisT]: ...
@overload
def prepare_input[AxisT: _ToAxis](
    x: onp.ToFloat1D, y: onp.ToComplexND, axis: AxisT, dydx: onp.ToComplexND | None = None, xp: None = None
) -> _PreparedInput[Any, AxisT]: ...
@overload
def prepare_input[AxisT: _ToAxis](
    x: onp.ToFloat1D, y: onp.ToComplexND, axis: AxisT, dydx: onp.ToComplexND | None = None, *, xp: ModuleType
) -> tuple[Incomplete, Incomplete, Incomplete, AxisT, Incomplete]: ...
