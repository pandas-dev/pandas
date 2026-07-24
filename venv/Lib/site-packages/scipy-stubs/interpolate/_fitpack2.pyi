from collections.abc import Sequence
from typing import Final, Literal, Never, overload, override

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

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

type _Float1D = onp.Array1D[np.float64]
type _FloatND = onp.ArrayND[np.float64]

type _Degree = Literal[1, 2, 3, 4, 5]

type _ExtInt = Literal[0, 1, 2, 3]
type _ExtStr = Literal["extrapolate", "zeros", "raise", "const"]
type _Ext = _ExtInt | _ExtStr

type _ToBBox = Sequence[onp.ToFloat | None]

###

# pyright reports a false positive `reportOverlappingOverload` error on numpy<2.1
# pyright: reportOverlappingOverload=false

class UnivariateSpline:
    ext: int

    # at runtime the `__init__` might change the `__class__` attribute...
    def __init__(
        self,
        /,
        x: onp.ToFloat1D,
        y: onp.ToFloat1D,
        w: onp.ToFloat1D | None = None,
        bbox: _ToBBox = [None, None],  # size 2
        k: _Degree = 3,
        s: float | None = None,
        ext: _Ext = 0,
        check_finite: bool = False,
    ) -> None: ...

    #
    @overload  # 1d
    def __call__(self, /, x: onp.ToFloatStrict1D, nu: int = 0, ext: _Ext | None = None) -> onp.Array1D[np.float64]: ...
    @overload  # 2d
    def __call__(self, /, x: onp.ToFloatStrict2D, nu: int = 0, ext: _Ext | None = None) -> onp.Array2D[np.float64]: ...
    @overload  # 3d
    def __call__(self, /, x: onp.ToFloatStrict3D, nu: int = 0, ext: _Ext | None = None) -> onp.Array3D[np.float64]: ...
    @overload  # known shape
    def __call__[ShapeT: tuple[int, ...]](
        self, /, x: onp.ArrayND[npc.floating | npc.integer, ShapeT], nu: int = 0, ext: _Ext | None = None
    ) -> onp.ArrayND[np.float64, ShapeT]: ...
    @overload  # ?d
    def __call__(self, /, x: onp.ToFloatND, nu: int = 0, ext: _Ext | None = None) -> onp.ArrayND[np.float64]: ...

    #
    def get_knots(self, /) -> _Float1D: ...
    def get_coeffs(self, /) -> _Float1D: ...
    def get_residual(self, /) -> float: ...
    def set_smoothing_factor(self, /, s: float) -> None: ...

    #
    def roots(self, /) -> _Float1D: ...  # requires `self.k == 3`
    def derivatives(self, /, x: onp.ToFloat) -> _Float1D: ...
    def derivative(self, /, n: int = 1) -> UnivariateSpline: ...
    def antiderivative(self, /, n: int = 1) -> UnivariateSpline: ...
    def integral(self, /, a: onp.ToFloat, b: onp.ToFloat) -> float: ...

    #
    @staticmethod
    def validate_input(
        x: onp.ToFloat1D,
        y: onp.ToFloat1D,
        w: onp.ToFloat1D,
        bbox: _ToBBox,
        k: _Degree,
        s: float | None,
        ext: _Ext,
        check_finite: bool,
    ) -> tuple[_Float1D, _Float1D, _Float1D, _Float1D, _ExtInt]: ...

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
        check_finite: bool = False,
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
        check_finite: bool = False,
    ) -> None: ...

class _BivariateSplineBase:  # undocumented
    @overload  # grid=True  (default)
    def __call__(
        self, /, x: onp.ToFloatND, y: onp.ToFloatND, dx: int = 0, dy: int = 0, grid: Literal[True] = True
    ) -> onp.Array2D[np.float64]: ...
    @overload  # grid=False
    def __call__(
        self, /, x: onp.ToFloatND, y: onp.ToFloatND, dx: int = 0, dy: int = 0, *, grid: Literal[False]
    ) -> onp.Array1D[np.float64]: ...

    #
    def get_residual(self, /) -> float: ...
    def get_knots(self, /) -> tuple[_Float1D, _Float1D]: ...
    def get_coeffs(self, /) -> _Float1D: ...
    def partial_derivative(self, /, dx: int, dy: int) -> _DerivedBivariateSpline: ...

class BivariateSpline(_BivariateSplineBase):
    def ev(self, /, xi: onp.ToFloatND, yi: onp.ToFloatND, dx: int = 0, dy: int = 0) -> _FloatND: ...
    def integral(self, /, xa: onp.ToFloat, xb: onp.ToFloat, ya: onp.ToFloat, yb: onp.ToFloat) -> float: ...

class _DerivedBivariateSpline(_BivariateSplineBase):  # undocumented
    @property
    def fp(self, /) -> Never: ...
    @override
    def get_residual(self, /) -> Never: ...  # raises AttributeError

class SmoothBivariateSpline(BivariateSpline):
    fp: Final[float]
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
        s: float | None = None,
        eps: float = 1e-16,
    ) -> None: ...

class LSQBivariateSpline(BivariateSpline):
    fp: Final[float]
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
        eps: float | None = None,
    ) -> None: ...

class RectBivariateSpline(BivariateSpline):
    fp: Final[float]
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
        s: float = 0,
        maxit: int = 20,
    ) -> None: ...

class SphereBivariateSpline(_BivariateSplineBase):
    @override
    def __call__(  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]  # ty: ignore[invalid-method-override]
        self, /, theta: onp.ToFloat1D, phi: onp.ToFloat1D, dtheta: int = 0, dphi: int = 0, grid: bool = True
    ) -> _FloatND: ...
    def ev(self, /, theta: onp.ToFloatND, phi: onp.ToFloatND, dtheta: int = 0, dphi: int = 0) -> _FloatND: ...

class SmoothSphereBivariateSpline(SphereBivariateSpline):
    fp: Final[float]
    tck: tuple[_Float1D, _Float1D, int]
    degrees: tuple[int, int]

    def __init__(
        self,
        /,
        theta: onp.ToFloat1D,
        phi: onp.ToFloat1D,
        r: onp.ToFloat1D,
        w: onp.ToFloat1D | None = None,
        s: float = 0.0,
        eps: float = 1e-16,
    ) -> None: ...

class LSQSphereBivariateSpline(SphereBivariateSpline):
    fp: Final[float]
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
        eps: float = 1e-16,
    ) -> None: ...

class RectSphereBivariateSpline(SphereBivariateSpline):
    fp: Final[float]
    tck: tuple[_Float1D, _Float1D, int]
    degrees: tuple[int, int]
    v0: np.float64

    def __init__(
        self,
        /,
        u: onp.ToFloat1D,
        v: onp.ToFloat1D,
        r: onp.ToFloat2D,
        s: float = 0.0,
        pole_continuity: bool | tuple[bool, bool] = False,
        pole_values: onp.ToFloat | tuple[onp.ToFloat, onp.ToFloat] | None = None,
        pole_exact: bool | tuple[bool, bool] = False,
        pole_flat: bool | tuple[bool, bool] = False,
    ) -> None: ...
