from _typeshed import Incomplete
from collections.abc import Callable, Sequence
from typing import Any, ClassVar, Final, Generic, Protocol, Self, SupportsIndex, TypeAlias, final, overload, type_check_only
from typing_extensions import TypeIs, TypeVar, override

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = [
    "BarycentricInterpolator",
    "KroghInterpolator",
    "approximate_taylor_polynomial",
    "barycentric_interpolate",
    "krogh_interpolate",
]

_XT = TypeVar("_XT", bound=npc.integer | npc.floating)
_XT_co = TypeVar("_XT_co", bound=npc.integer | npc.floating, default=Any, covariant=True)
_YT = TypeVar("_YT", bound=np.float64 | np.complex128)
_YT_co = TypeVar("_YT_co", bound=np.float64 | np.complex128, default=np.float64 | np.complex128, covariant=True)

_MultiIndex: TypeAlias = SupportsIndex | tuple[SupportsIndex, ...]

@type_check_only
class _HasShape0(Protocol):
    @property
    def shape(self) -> tuple[()]: ...

###

class _Interpolator1D(Generic[_YT_co]):  # undocumented
    __slots__: ClassVar[tuple[str, ...]] = "_y_axis", "_y_extra_shape", "dtype"

    _y_axis: int | None
    _y_extra_shape: tuple[int, ...] | None
    dtype: type[_YT_co] | None

    @overload  # no yi (default); unknown dtype
    def __init__(self, /, xi: onp.ToFloatND | None = None, yi: None = None, axis: int | None = None) -> None: ...
    @overload  # floating yi (positional)
    def __init__(
        self: _Interpolator1D[np.float64], /, xi: onp.ToFloatND | None, yi: onp.ToFloatND, axis: int | None = None
    ) -> None: ...
    @overload  # floating yi (keyword)
    def __init__(
        self: _Interpolator1D[np.float64], /, xi: onp.ToFloatND | None = None, *, yi: onp.ToFloatND, axis: int | None = None
    ) -> None: ...
    @overload  # complex yi (positional)
    def __init__(
        self: _Interpolator1D[np.complex128],
        /,
        xi: onp.ToFloatND | None,
        yi: onp.ToJustComplexND | None = None,
        axis: int | None = None,
    ) -> None: ...
    @overload  # complex yi (keyword)
    def __init__(
        self: _Interpolator1D[np.complex128],
        /,
        xi: onp.ToFloatND | None = None,
        *,
        yi: onp.ToJustComplexND,
        axis: int | None = None,
    ) -> None: ...

    #
    def __call__(self, /, x: onp.ToFloat | onp.ToFloatND) -> onp.ArrayND[_YT_co]: ...

    #
    def _evaluate(self, /, x: onp.ToFloat1D) -> onp.ArrayND[_YT_co]: ...  # undocumented, not implemented
    @final
    def _prepare_x(
        self, /, x: onp.ToFloat | onp.ToFloatND
    ) -> tuple[onp.Array1D[Incomplete], tuple[int, ...]]: ...  # undocumented
    @final
    def _finish_y(self, /, y: onp.Array2D[_YT_co], x_shape: tuple[int, ...]) -> onp.ArrayND[_YT_co]: ...  # undocumented

    #
    @overload
    def _reshape_yi(
        self: _Interpolator1D[np.float64], /, yi: onp.ToFloatND, check: bool = False
    ) -> onp.Array2D[_YT_co]: ...  # undocumented
    @overload
    def _reshape_yi(
        self: _Interpolator1D[np.complex128], /, yi: onp.ToJustComplexND, check: bool = False
    ) -> onp.Array2D[_YT_co]: ...  # undocumented

    #
    @overload
    def _set_yi(
        self: _Interpolator1D[np.float64], /, yi: onp.ToFloatND, xi: onp.ToFloatND | None = None, axis: int | None = None
    ) -> None: ...  # undocumented
    @overload
    def _set_yi(
        self: _Interpolator1D[np.complex128], /, yi: onp.ToJustComplexND, xi: onp.ToFloatND | None = None, axis: int | None = None
    ) -> None: ...  # undocumented

    #
    @final
    def _set_dtype(
        self: _Interpolator1D[_YT], /, dtype: np.dtype[_YT] | type[_YT], union: bool = False
    ) -> None: ...  # undocumented

class _Interpolator1DWithDerivatives(_Interpolator1D[_YT_co], Generic[_YT_co]):  # undocumented
    def derivatives(self, /, x: onp.ToFloat | onp.ToFloatND, der: _MultiIndex | None = None) -> onp.ArrayND[_YT_co]: ...
    def derivative(self, /, x: onp.ToFloat | onp.ToFloatND, der: SupportsIndex = 1) -> onp.ArrayND[_YT_co]: ...
    def _evaluate_derivatives(
        self, /, x: onp.ToFloat1D, der: int | None = None
    ) -> onp.ArrayND[_YT_co]: ...  # undocumented, not implemented

# NOTE: `KroghInterpolator` is not generic at runtime (`scipy<1.17`):
# https://github.com/scipy/scipy-stubs/issues/653
class KroghInterpolator(_Interpolator1DWithDerivatives[_YT_co], Generic[_YT_co, _XT_co]):
    xi: onp.Array1D[_XT_co]
    yi: onp.Array2D[_YT_co]
    c: onp.Array2D[_YT_co]  # undocumented
    n: Final[int]  # undocumented
    r: Final[int]  # undocumented

    @overload  # xi: T, yi: f64
    def __init__(
        self: KroghInterpolator[np.float64, _XT],
        /,
        xi: onp.CanArray[Any, np.dtype[_XT]] | Sequence[_XT],
        yi: onp.ToFloatND,
        axis: int = 0,
    ) -> None: ...
    @overload  # xi: T, yi: c128
    def __init__(
        self: KroghInterpolator[np.complex128, _XT],
        /,
        xi: onp.CanArray[Any, np.dtype[_XT]] | Sequence[_XT],
        yi: onp.ToJustComplexND,
        axis: int = 0,
    ) -> None: ...
    @overload  # xi: i64, yi: f64
    def __init__(
        self: KroghInterpolator[np.float64, np.float64], /, xi: onp.ToJustInt64_1D, yi: onp.ToFloatND, axis: int = 0
    ) -> None: ...
    @overload  # xi: i64, yi: c128
    def __init__(
        self: KroghInterpolator[np.complex128, np.float64], /, xi: onp.ToJustInt64_1D, yi: onp.ToJustComplexND, axis: int = 0
    ) -> None: ...
    @overload  # xi: f64, yi: f64
    def __init__(
        self: KroghInterpolator[np.float64, np.float64], /, xi: onp.ToJustFloat64_1D, yi: onp.ToFloatND, axis: int = 0
    ) -> None: ...
    @overload  # xi: f64, yi: c128
    def __init__(
        self: KroghInterpolator[np.complex128, np.float64], /, xi: onp.ToJustFloat64_1D, yi: onp.ToJustComplexND, axis: int = 0
    ) -> None: ...

    #
    @override
    def _evaluate(self, /, x: onp.ToFloat1D) -> onp.Array2D[_YT_co]: ...  # undocumented
    @override
    def _evaluate_derivatives(self, /, x: onp.ToFloat1D, der: int | None = None) -> onp.Array3D[_YT_co]: ...  # undocumented

# NOTE: `BarycentricInterpolator` is not generic at runtime (`scipy<1.17`):
# https://github.com/scipy/scipy-stubs/issues/653
class BarycentricInterpolator(_Interpolator1DWithDerivatives[_YT_co], Generic[_YT_co]):
    xi: onp.Array1D[np.float64]
    yi: onp.Array2D[_YT_co] | None
    wi: onp.Array1D[npc.floating] | None

    n: Final[int]  # undocumented
    r: Final[int]  # undocumented
    _diff_baryint: Self | None  # undocumented

    @overload  # yi: f64
    def __init__(
        self: BarycentricInterpolator[np.float64],
        /,
        xi: onp.ToFloat1D,
        yi: onp.ToFloatND | None = None,
        axis: int = 0,
        *,
        wi: onp.ArrayND[npc.floating] | None = None,
        rng: onp.random.ToRNG | None = None,
        random_state: onp.random.ToRNG | None = None,  # legacy
    ) -> None: ...
    @overload  # yi: c128
    def __init__(
        self: BarycentricInterpolator[np.complex128],
        /,
        xi: onp.ToFloat1D,
        yi: onp.ToJustComplexND,
        axis: int = 0,
        *,
        wi: onp.ArrayND[npc.floating] | None = None,
        rng: onp.random.ToRNG | None = None,
        random_state: onp.random.ToRNG | None = None,  # legacy
    ) -> None: ...

    #
    @overload
    def set_yi(self: BarycentricInterpolator[np.float64], /, yi: onp.ToFloatND, axis: int | None = None) -> None: ...
    @overload
    def set_yi(self: BarycentricInterpolator[np.complex128], /, yi: onp.ToJustComplexND, axis: int | None = None) -> None: ...

    #
    @overload
    def add_xi(self: BarycentricInterpolator[np.float64], /, xi: onp.ToFloat1D, yi: onp.ToFloat2D | None = None) -> None: ...
    @overload
    def add_xi(self: BarycentricInterpolator[np.complex128], /, xi: onp.ToFloat1D, yi: onp.ToComplex2D | None = None) -> None: ...

    #
    @override
    def _evaluate_derivatives(
        self, /, x: onp.ToFloat1D, der: int | None = None, all_lower: bool = True
    ) -> onp.Array2D[_YT_co] | onp.Array3D[_YT_co]: ...  # undocumented

#
def _isscalar(x: object) -> TypeIs[np.generic | complex | str | bytes | _HasShape0]: ...  # undocumented

#
@overload
def krogh_interpolate(
    xi: onp.ToFloat1D, yi: onp.ToFloatND, x: onp.ToFloat | onp.ToFloat1D, der: int | list[int] | None = 0, axis: int = 0
) -> onp.ArrayND[np.float64]: ...
@overload
def krogh_interpolate(
    xi: onp.ToFloat1D, yi: onp.ToJustComplexND, x: onp.ToFloat | onp.ToFloat1D, der: int | list[int] | None = 0, axis: int = 0
) -> onp.ArrayND[np.complex128]: ...
@overload
def krogh_interpolate(
    xi: onp.ToFloat1D, yi: onp.ToComplexND, x: onp.ToFloat | onp.ToFloat1D, der: int | list[int] | None = 0, axis: int = 0
) -> onp.ArrayND[np.float64 | np.complex128]: ...

#
def approximate_taylor_polynomial(
    f: Callable[[onp.Array1D[np.float64]], onp.ToComplexND] | np.ufunc,
    x: onp.ToFloat,
    degree: int,
    scale: onp.ToFloat,
    order: int | None = None,
) -> np.poly1d: ...

#
@overload  # 0d f64
def barycentric_interpolate(
    xi: onp.ToFloat1D,
    yi: onp.ToFloatND,
    x: onp.ToFloat,
    axis: int = 0,
    *,
    der: int | list[int] | None = 0,
    rng: onp.random.ToRNG | None = None,
) -> np.float64: ...
@overload  # 1d f64
def barycentric_interpolate(
    xi: onp.ToFloat1D,
    yi: onp.ToFloatND,
    x: onp.ToFloatStrict1D,
    axis: int = 0,
    *,
    der: int | list[int] | None = 0,
    rng: onp.random.ToRNG | None = None,
) -> onp.Array1D[np.float64]: ...
@overload  # nd f64
def barycentric_interpolate(
    xi: onp.ToFloat1D,
    yi: onp.ToFloatND,
    x: onp.ToFloatND,
    axis: int = 0,
    *,
    der: int | list[int] | None = 0,
    rng: onp.random.ToRNG | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # 0d c128
def barycentric_interpolate(
    xi: onp.ToFloat1D,
    yi: onp.ToJustComplexND,
    x: onp.ToFloat,
    axis: int = 0,
    *,
    der: int | list[int] | None = 0,
    rng: onp.random.ToRNG | None = None,
) -> np.complex128: ...
@overload  # 1d c128
def barycentric_interpolate(
    xi: onp.ToFloat1D,
    yi: onp.ToJustComplexND,
    x: onp.ToFloatStrict1D,
    axis: int = 0,
    *,
    der: int | list[int] | None = 0,
    rng: onp.random.ToRNG | None = None,
) -> onp.Array1D[np.complex128]: ...
@overload  # nd c128
def barycentric_interpolate(
    xi: onp.ToFloat1D,
    yi: onp.ToJustComplexND,
    x: onp.ToFloatND,
    axis: int = 0,
    *,
    der: int | list[int] | None = 0,
    rng: onp.random.ToRNG | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload  # 0d f64 or c128
def barycentric_interpolate(
    xi: onp.ToFloat1D,
    yi: onp.ToComplexND,
    x: onp.ToFloat,
    axis: int = 0,
    *,
    der: int | list[int] | None = 0,
    rng: onp.random.ToRNG | None = None,
) -> np.float64 | np.complex128: ...
@overload  # 1d f64 or c128
def barycentric_interpolate(
    xi: onp.ToFloat1D,
    yi: onp.ToComplexND,
    x: onp.ToFloatStrict1D,
    axis: int = 0,
    *,
    der: int | list[int] | None = 0,
    rng: onp.random.ToRNG | None = None,
) -> onp.Array1D[np.float64 | np.complex128]: ...
@overload  # nd f64 or c128
def barycentric_interpolate(
    xi: onp.ToFloat1D,
    yi: onp.ToComplexND,
    x: onp.ToFloatND,
    axis: int = 0,
    *,
    der: int | list[int] | None = 0,
    rng: onp.random.ToRNG | None = None,
) -> onp.ArrayND[np.float64 | np.complex128]: ...
