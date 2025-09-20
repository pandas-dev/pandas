# scipy/interpolate/interpnd.pyx

from typing import Any, Generic, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp

from scipy.spatial._qhull import Delaunay, DelaunayInfo_t

###

_CT_co = TypeVar("_CT_co", bound=np.float64 | np.complex128, default=np.float64, covariant=True)

###

class GradientEstimationWarning(Warning): ...

class NDInterpolatorBase(Generic[_CT_co]):
    points: onp.Array2D[np.float64]
    values: onp.ArrayND[_CT_co] | None
    is_complex: bool
    scale: onp.Array1D[np.float64] | None
    offset: onp.Array1D[np.float64]  # only if rescale=True

    @overload
    def __init__(
        self: NDInterpolatorBase[np.float64],
        /,
        points: onp.ToFloat2D | Delaunay,
        values: onp.ToFloatND,
        fill_value: onp.ToFloat = ...,  # np.nan
        ndim: int | None = None,
        rescale: bool = False,
        need_contiguous: bool = True,
        need_values: bool = True,
    ) -> None: ...
    @overload
    def __init__(
        self: NDInterpolatorBase[np.complex128],
        /,
        points: onp.ToFloat2D | Delaunay,
        values: onp.ToJustComplexND | None,
        fill_value: onp.ToComplex = ...,  # np.nan
        ndim: int | None = None,
        rescale: bool = False,
        need_contiguous: bool = True,
        need_values: bool = True,
    ) -> None: ...
    @overload
    def __init__(
        self: NDInterpolatorBase[Any],
        /,
        points: onp.ToFloat2D | Delaunay,
        values: onp.ToComplexND | None,
        fill_value: onp.ToComplex = ...,  # np.nan
        ndim: int | None = None,
        rescale: bool = False,
        need_contiguous: bool = True,
        need_values: bool = True,
    ) -> None: ...

    #
    def __call__(self, /, *args: onp.ToFloatND) -> onp.ArrayND[_CT_co]: ...

class LinearNDInterpolator(NDInterpolatorBase[_CT_co], Generic[_CT_co]):
    @overload
    def __init__(
        self: LinearNDInterpolator[np.float64],
        /,
        points: onp.ToFloat2D | Delaunay,
        values: onp.ToFloatND,
        fill_value: onp.ToFloat = ...,  # np.nan
        rescale: bool = False,
    ) -> None: ...
    @overload
    def __init__(
        self: LinearNDInterpolator[np.complex128],
        /,
        points: onp.ToFloat2D | Delaunay,
        values: onp.ToJustComplexND,
        fill_value: onp.ToComplex = ...,  # np.nan
        rescale: bool = False,
    ) -> None: ...
    @overload
    def __init__(
        self: LinearNDInterpolator[Any],
        /,
        points: onp.ToFloat2D | Delaunay,
        values: onp.ToComplexND,
        fill_value: onp.ToComplex = ...,  # np.nan
        rescale: bool = False,
    ) -> None: ...

class CloughTocher2DInterpolator(NDInterpolatorBase[_CT_co], Generic[_CT_co]):
    @overload
    def __init__(
        self: CloughTocher2DInterpolator[np.float64],
        /,
        points: onp.ToFloat2D | Delaunay,
        values: onp.ToFloatND,
        fill_value: onp.ToFloat = ...,  # np.nan
        tol: onp.ToFloat = 1e-06,
        maxiter: int = 400,
        rescale: bool = False,
    ) -> None: ...
    @overload
    def __init__(
        self: CloughTocher2DInterpolator[np.complex128],
        /,
        points: onp.ToFloat2D | Delaunay,
        values: onp.ToJustComplexND,
        fill_value: onp.ToComplex = ...,  # np.nan
        tol: onp.ToFloat = 1e-06,
        maxiter: int = 400,
        rescale: bool = False,
    ) -> None: ...
    @overload
    def __init__(
        self: CloughTocher2DInterpolator[Any],
        /,
        points: onp.ToFloat2D | Delaunay,
        values: onp.ToComplexND,
        fill_value: onp.ToComplex = ...,  # np.nan
        tol: onp.ToFloat = 1e-06,
        maxiter: int = 400,
        rescale: bool = False,
    ) -> None: ...

@overload
def estimate_gradients_2d_global(
    tri: DelaunayInfo_t, y: onp.ToFloat1D | onp.ToFloat2D, maxiter: onp.ToJustInt = 400, tol: float = 1e-6
) -> onp.Array3D[np.float64]: ...
@overload
def estimate_gradients_2d_global(
    tri: DelaunayInfo_t, y: onp.ToJustComplex1D | onp.ToJustComplex2D, maxiter: onp.ToJustInt = 400, tol: float = 1e-6
) -> onp.Array3D[np.complex128]: ...
@overload
def estimate_gradients_2d_global(
    tri: DelaunayInfo_t, y: onp.ToComplex1D | onp.ToComplex2D, maxiter: onp.ToJustInt = 400, tol: float = 1e-6
) -> onp.Array3D[np.float64] | onp.Array3D[np.complex128]: ...
