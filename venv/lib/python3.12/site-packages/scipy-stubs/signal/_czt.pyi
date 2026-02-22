from typing import Final, TypeAlias, overload

import numpy as np
import optype as op
import optype.numpy as onp

__all__ = ["CZT", "ZoomFFT", "czt", "czt_points", "zoom_fft"]

_Complex: TypeAlias = np.complex128 | np.clongdouble

###

class CZT:
    w: Final[complex | np.complex128]
    a: Final[complex | np.complex128]
    m: Final[int]
    n: Final[int]

    def __init__(
        self, /, n: int, m: int | None = None, w: complex | np.complex128 | None = None, a: complex | np.complex128 = 1 + 0j
    ) -> None: ...
    @overload
    def __call__(self, /, x: onp.ToComplexStrict1D, *, axis: op.CanIndex = -1) -> onp.Array1D[_Complex]: ...
    @overload
    def __call__(self, /, x: onp.ToComplexStrict2D, *, axis: op.CanIndex = -1) -> onp.Array2D[_Complex]: ...
    @overload
    def __call__(self, /, x: onp.ToComplexStrict3D, *, axis: op.CanIndex = -1) -> onp.Array3D[_Complex]: ...
    @overload
    def __call__(self, /, x: onp.ToComplexND, *, axis: op.CanIndex = -1) -> onp.ArrayND[_Complex]: ...
    def points(self, /) -> onp.Array1D[np.complex128]: ...

class ZoomFFT(CZT):
    f1: onp.ToFloat
    f2: onp.ToFloat
    fs: float | np.float64

    def __init__(
        self,
        /,
        n: int,
        fn: float | np.float64 | onp.ToFloat1D,
        m: int | None = None,
        *,
        fs: float | np.float64 = 2,
        endpoint: onp.ToBool = False,
    ) -> None: ...

#
def _validate_sizes(n: int, m: int | None) -> int: ...

#
def czt_points(
    m: int, w: complex | np.complex128 | None = None, a: complex | np.complex128 = 1 + 0j
) -> onp.Array1D[np.complex128]: ...

#
@overload
def czt(
    x: onp.ToComplexStrict1D,
    m: int | None = None,
    w: complex | np.complex128 | None = None,
    a: complex | np.complex128 = 1 + 0j,
    *,
    axis: op.CanIndex = -1,
) -> onp.Array1D[_Complex]: ...
@overload
def czt(
    x: onp.ToComplexStrict2D,
    m: int | None = None,
    w: complex | np.complex128 | None = None,
    a: complex | np.complex128 = 1 + 0j,
    *,
    axis: op.CanIndex = -1,
) -> onp.Array2D[_Complex]: ...
@overload
def czt(
    x: onp.ToComplexStrict3D,
    m: int | None = None,
    w: complex | np.complex128 | None = None,
    a: complex | np.complex128 = 1 + 0j,
    *,
    axis: op.CanIndex = -1,
) -> onp.Array3D[_Complex]: ...
@overload
def czt(
    x: onp.ToComplexND,
    m: int | None = None,
    w: complex | np.complex128 | None = None,
    a: complex | np.complex128 = 1 + 0j,
    *,
    axis: op.CanIndex = -1,
) -> onp.ArrayND[_Complex]: ...

#
@overload
def zoom_fft(
    x: onp.ToComplexStrict1D,
    fn: float | np.float64 | onp.ToFloat1D,
    m: int | None = None,
    *,
    fs: float | np.float64 = 2,
    endpoint: onp.ToBool = False,
    axis: op.CanIndex = -1,
) -> onp.Array1D[_Complex]: ...
@overload
def zoom_fft(
    x: onp.ToComplexStrict2D,
    fn: float | np.float64 | onp.ToFloat1D,
    m: int | None = None,
    *,
    fs: float | np.float64 = 2,
    endpoint: onp.ToBool = False,
    axis: op.CanIndex = -1,
) -> onp.Array2D[_Complex]: ...
@overload
def zoom_fft(
    x: onp.ToComplexStrict3D,
    fn: float | np.float64 | onp.ToFloat1D,
    m: int | None = None,
    *,
    fs: float | np.float64 = 2,
    endpoint: onp.ToBool = False,
    axis: op.CanIndex = -1,
) -> onp.Array3D[_Complex]: ...
@overload
def zoom_fft(
    x: onp.ToComplexND,
    fn: float | np.float64 | onp.ToFloat1D,
    m: int | None = None,
    *,
    fs: float | np.float64 = 2,
    endpoint: onp.ToBool = False,
    axis: op.CanIndex = -1,
) -> onp.ArrayND[_Complex]: ...
