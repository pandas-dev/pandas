from collections.abc import Sequence
from typing import Literal, Never, TypeAlias
from typing_extensions import override

import numpy as np
import optype.numpy as onp

__all__ = [
    "BivariateSpline",
    "InterpolatedUnivariateSpline",
    "LSQBivariateSpline",
    "LSQSphereBivariateSpline",
    "LSQUnivariateSpline",
    "RectBivariateSpline",
    "RectSphereBivariateSpline",
    "SmoothBivariateSpline",
    "SmoothSphereBivariateSpline",
    "UnivariateSpline",
]

_Float1D: TypeAlias = onp.Array1D[np.float64]
_FloatND: TypeAlias = onp.Array2D[np.float64]

_Degree: TypeAlias = Literal[1, 2, 3, 4, 5]

_ExtInt: TypeAlias = Literal[0, 1, 2, 3]
_ExtStr: TypeAlias = Literal["extrapolate", "zeroes", "raise", "const"]
_Ext: TypeAlias = _ExtInt | _ExtStr

_BBox: TypeAlias = onp.Array[tuple[Literal[2]], np.float64]
_ToBBox: TypeAlias = Sequence[onp.ToFloat | None]

###

class UnivariateSpline:
    @staticmethod
    def validate_input(
        x: onp.ToFloat1D,
        y: onp.ToFloat1D,
        w: onp.ToFloat1D,
        bbox: _ToBBox,
        k: _Degree,
        s: onp.ToFloat | None,
        ext: _Ext,
        check_finite: onp.ToBool,
    ) -> tuple[_Float1D, _Float1D, _Float1D, _BBox, _ExtInt]: ...

    # at runtime the `__init__` might change the `__class__` attribute...
    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloat1D,
        w: onp.ToFloat1D | None = None,
        bbox: _ToBBox = [None, None],  # size 2
        k: _Degree = 3,
        s: onp.ToFloat | None = None,
        ext: _Ext = 0,
        check_finite: onp.ToBool = False,
    ) -> None: ...
    def __call__(self, /, x: onp.ToFloat1D, nu: int = 0, ext: _Ext | None = None) -> _Float1D: ...

    #
    def get_knots(self, /) -> _Float1D: ...
    def get_coeffs(self, /) -> _Float1D: ...
    def get_residual(self, /) -> _Float1D: ...
    def set_smoothing_factor(self, /, s: onp.ToFloat) -> None: ...

    #
    def roots(self, /) -> _Float1D: ...  # requires `self.k == 3`
    def derivatives(self, /, x: onp.ToFloat) -> _Float1D: ...
    def derivative(self, /, n: int = 1) -> UnivariateSpline: ...
    def antiderivative(self, /, n: int = 1) -> UnivariateSpline: ...
    def integral(self, /, a: onp.ToFloat, b: onp.ToFloat) -> float: ...

class InterpolatedUnivariateSpline(UnivariateSpline):
    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloat1D,
        w: onp.ToFloat1D | None = None,
        bbox: _ToBBox = [None, None],  # size 2
        k: _Degree = 3,
        ext: _Ext = 0,
        check_finite: onp.ToBool = False,
    ) -> None: ...

class LSQUnivariateSpline(UnivariateSpline):
    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloat1D,
        t: onp.ToFloat1D,
        w: onp.ToFloat1D | None = None,
        bbox: _ToBBox = [None, None],  # size 2
        k: _Degree = 3,
        ext: _Ext = 0,
        check_finite: onp.ToBool = False,
    ) -> None: ...

class _BivariateSplineBase:  # undocumented
    def __call__(self, /, x: onp.ToFloatND, y: onp.ToFloatND, dx: int = 0, dy: int = 0, grid: onp.ToBool = True) -> _Float1D: ...
    def get_residual(self, /) -> _Float1D: ...
    def get_knots(self, /) -> tuple[_Float1D, _Float1D]: ...
    def get_coeffs(self, /) -> _Float1D: ...
    def partial_derivative(self, /, dx: int, dy: int) -> _DerivedBivariateSpline: ...

class BivariateSpline(_BivariateSplineBase):
    def ev(self, /, xi: onp.ToFloatND, yi: onp.ToFloatND, dx: int = 0, dy: int = 0) -> _FloatND: ...
    def integral(self, /, xa: onp.ToFloat, xb: onp.ToFloat, ya: onp.ToFloat, yb: onp.ToFloat) -> float: ...

class _DerivedBivariateSpline(_BivariateSplineBase):  # undocumented
    @property
    def fp(self, /) -> Never: ...

class SmoothBivariateSpline(BivariateSpline):
    fp: float
    tck: tuple[_Float1D, _Float1D, int]
    degrees: tuple[int, int]

    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloat1D,
        z: onp.ToFloat1D,
        w: onp.ToFloat1D | None = None,
        bbox: _ToBBox = [None, None, None, None],
        kx: int = 3,
        ky: int = 3,
        s: onp.ToFloat | None = None,
        eps: onp.ToFloat = 1e-16,
    ) -> None: ...

class LSQBivariateSpline(BivariateSpline):
    fp: float
    tck: tuple[_Float1D, _Float1D, int]
    degrees: tuple[int, int]

    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloat1D,
        z: onp.ToFloat1D,
        tx: onp.ToFloat1D,
        ty: onp.ToFloat1D,
        w: onp.ToFloat1D | None = None,
        bbox: _ToBBox = [None, None, None, None],
        kx: int = 3,
        ky: int = 3,
        eps: onp.ToFloat | None = None,
    ) -> None: ...

class RectBivariateSpline(BivariateSpline):
    fp: float
    tck: tuple[_Float1D, _Float1D, int]
    degrees: tuple[int, int]

    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloat1D,
        z: onp.ToFloat2D,
        bbox: _ToBBox = [None, None, None, None],
        kx: int = 3,
        ky: int = 3,
        s: onp.ToFloat = 0,
        maxit: int = 20,
    ) -> None: ...

class SphereBivariateSpline(_BivariateSplineBase):
    @override
    def __call__(  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]
        self, /, theta: onp.ToFloat1D, phi: onp.ToFloat1D, dtheta: int = 0, dphi: int = 0, grid: onp.ToBool = True
    ) -> _FloatND: ...
    def ev(self, /, theta: onp.ToFloatND, phi: onp.ToFloatND, dtheta: int = 0, dphi: int = 0) -> _FloatND: ...

class SmoothSphereBivariateSpline(SphereBivariateSpline):
    fp: float
    tck: tuple[_Float1D, _Float1D, int]
    degrees: tuple[int, int]

    def __init__(
        self,
        /,
        theta: onp.ToFloat1D,
        phi: onp.ToFloat1D,
        r: onp.ToFloat1D,
        w: onp.ToFloat1D | None = None,
        s: onp.ToFloat = 0.0,
        eps: onp.ToFloat = 1e-16,
    ) -> None: ...

class LSQSphereBivariateSpline(SphereBivariateSpline):
    fp: float
    tck: tuple[_Float1D, _Float1D, int]
    degrees: tuple[int, int]

    def __init__(
        self,
        /,
        theta: onp.ToFloat1D,
        phi: onp.ToFloat1D,
        r: onp.ToFloat1D,
        tt: onp.ToFloat1D,
        tp: onp.ToFloat1D,
        w: onp.ToFloat1D | None = None,
        eps: onp.ToFloat = 1e-16,
    ) -> None: ...

class RectSphereBivariateSpline(SphereBivariateSpline):
    fp: float
    tck: tuple[_Float1D, _Float1D, int]
    degrees: tuple[int, int]
    v0: np.float64

    def __init__(
        self,
        /,
        u: onp.ToFloat1D,
        v: onp.ToFloat1D,
        r: onp.ToFloat2D,
        s: onp.ToFloat = 0.0,
        pole_continuity: onp.ToBool | tuple[onp.ToBool, onp.ToBool] = False,
        pole_values: onp.ToFloat | tuple[onp.ToFloat, onp.ToFloat] | None = None,
        pole_exact: onp.ToBool | tuple[onp.ToBool, onp.ToBool] = False,
        pole_flat: onp.ToBool | tuple[onp.ToBool, onp.ToBool] = False,
    ) -> None: ...
