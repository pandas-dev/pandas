from typing import Any, Generic, Literal, Self, TypeAlias, TypeVar, overload

import numpy as np
import optype as op
import optype.numpy as onp

from scipy.interpolate import CubicSpline
from scipy.sparse import csr_array

__all__ = ["BSpline", "make_interp_spline", "make_lsq_spline", "make_smoothing_spline"]

###

_CT_co = TypeVar("_CT_co", bound=np.float64 | np.complex128, default=np.float64, covariant=True)

_Extrapolate: TypeAlias = Literal["periodic"] | bool
_BCType: TypeAlias = Literal["not-a-knot", "natural", "clamped", "periodic"]
_ToBCType: TypeAlias = tuple[onp.ToFloat, onp.ToFloat] | _BCType
_LSQMethod: TypeAlias = Literal["qr", "norm-eq"]

###

class BSpline(Generic[_CT_co]):
    t: onp.Array1D[np.float64]
    c: onp.Array[onp.AtLeast1D, _CT_co]
    k: int
    axis: int
    extrapolate: _Extrapolate

    @property
    def tck(self, /) -> tuple[onp.Array1D[np.float64], onp.Array[onp.AtLeast1D, _CT_co], int]: ...

    #
    @overload
    def __init__(
        self: BSpline[np.float64],
        /,
        t: onp.ToFloat1D,
        c: onp.ToFloatND,
        k: op.CanIndex,
        extrapolate: _Extrapolate = True,
        axis: op.CanIndex = 0,
    ) -> None: ...
    @overload
    def __init__(
        self: BSpline[np.complex128],
        /,
        t: onp.ToFloat1D,
        c: onp.ToJustComplexND,
        k: op.CanIndex,
        extrapolate: _Extrapolate = True,
        axis: op.CanIndex = 0,
    ) -> None: ...
    @overload
    def __init__(
        self: BSpline[Any],
        /,
        t: onp.ToFloat1D,
        c: onp.ToComplex1D,
        k: op.CanIndex,
        extrapolate: _Extrapolate = True,
        axis: op.CanIndex = 0,
    ) -> None: ...

    # NOTE: Complex `x` will unsafely be cast to `float64`, even if the coefficients are complex
    def __call__(
        self, /, x: onp.ToComplex | onp.ToComplexND, nu: int = 0, extrapolate: _Extrapolate | None = None
    ) -> onp.ArrayND[_CT_co]: ...

    #
    def derivative(self, /, nu: int = 1) -> Self: ...
    def antiderivative(self, /, nu: int = 1) -> Self: ...

    #
    def insert_knot(self, /, x: onp.ToFloat, m: op.CanIndex = 1) -> Self: ...

    # NOTE: `integrate` will raise a (cryptic) `ValueError` for complex coefficients
    def integrate(
        self: BSpline[np.float64], /, a: onp.ToFloat, b: onp.ToFloat, extrapolate: _Extrapolate | None = None
    ) -> onp.ArrayND[np.float64]: ...

    #
    @overload
    @classmethod
    def basis_element(cls, t: onp.ToFloatND, extrapolate: _Extrapolate = True) -> BSpline[np.float64]: ...
    @overload
    @classmethod
    def basis_element(cls, t: onp.ToJustComplexND, extrapolate: _Extrapolate = True) -> BSpline[np.complex128]: ...
    @overload
    @classmethod
    def basis_element(cls, t: onp.ToComplexND, extrapolate: _Extrapolate = True) -> BSpline[Any]: ...

    #
    @classmethod
    def from_power_basis(cls, pp: CubicSpline[_CT_co], bc_type: _BCType = "not-a-knot") -> Self: ...

    #
    @classmethod
    def construct_fast(
        cls, t: onp.ArrayND[np.float64], c: onp.ArrayND[_CT_co], k: int, extrapolate: _Extrapolate = True, axis: int = 0
    ) -> Self: ...

    #
    @classmethod
    def design_matrix(
        cls, x: onp.ToFloat1D, t: onp.ToFloat1D, k: op.CanIndex, extrapolate: _Extrapolate = False
    ) -> csr_array[np.float64, tuple[int, int]]: ...

#
@overload
def make_interp_spline(
    x: onp.ToFloat1D,
    y: onp.ToFloatND,
    k: op.CanIndex = 3,
    t: onp.ToFloat1D | None = None,
    bc_type: _ToBCType | None = None,
    axis: op.CanIndex = 0,
    check_finite: onp.ToBool = True,
) -> BSpline[np.float64]: ...
@overload
def make_interp_spline(
    x: onp.ToFloat1D,
    y: onp.ToJustComplexND,
    k: op.CanIndex = 3,
    t: onp.ToFloat1D | None = None,
    bc_type: _ToBCType | None = None,
    axis: op.CanIndex = 0,
    check_finite: onp.ToBool = True,
) -> BSpline[np.complex128]: ...
@overload
def make_interp_spline(
    x: onp.ToFloat1D,
    y: onp.ToComplexND,
    k: op.CanIndex = 3,
    t: onp.ToFloat1D | None = None,
    bc_type: _ToBCType | None = None,
    axis: op.CanIndex = 0,
    check_finite: onp.ToBool = True,
) -> BSpline[Any]: ...

#
@overload
def make_lsq_spline(
    x: onp.ToFloat1D,
    y: onp.ToFloatND,
    t: onp.ToFloat1D,
    k: op.CanIndex = 3,
    w: onp.ToFloat1D | None = None,
    axis: op.CanIndex = 0,
    check_finite: onp.ToBool = True,
    *,
    method: _LSQMethod = "qr",
) -> BSpline[np.float64]: ...
@overload
def make_lsq_spline(
    x: onp.ToFloat1D,
    y: onp.ToJustComplexND,
    t: onp.ToFloat1D,
    k: op.CanIndex = 3,
    w: onp.ToFloat1D | None = None,
    axis: op.CanIndex = 0,
    check_finite: onp.ToBool = True,
    *,
    method: _LSQMethod = "qr",
) -> BSpline[np.complex128]: ...
@overload
def make_lsq_spline(
    x: onp.ToFloat1D,
    y: onp.ToComplexND,
    t: onp.ToFloat1D,
    k: op.CanIndex = 3,
    w: onp.ToFloat1D | None = None,
    axis: op.CanIndex = 0,
    check_finite: onp.ToBool = True,
    *,
    method: _LSQMethod = "qr",
) -> BSpline[Any]: ...

# NOTE: will unsafely cast complex input to float64
def make_smoothing_spline(
    x: onp.ToComplex1D,
    y: onp.ToComplexND,
    w: onp.ToFloat1D | None = None,
    lam: onp.ToFloat | None = None,
    *,
    axis: op.CanIndex = 0,
) -> BSpline[np.float64]: ...

#
def fpcheck(x: onp.ToFloat1D, t: onp.ToFloat1D, k: onp.ToJustInt) -> None: ...  # undocumented
