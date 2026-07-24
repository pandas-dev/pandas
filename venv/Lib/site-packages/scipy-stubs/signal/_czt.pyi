from typing import Final, Never, SupportsIndex, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["CZT", "ZoomFFT", "czt", "czt_points", "zoom_fft"]

###

type _Complex = np.complex128 | np.clongdouble

# workaround for non-overload-spec-compliant type-checkers
type _JustAnyShape = tuple[Never, Never, Never, Never]
type _ToComplexStrictND = onp.ArrayND[npc.number | np.bool, _JustAnyShape]

###

class CZT:
    w: Final[complex]
    a: Final[complex]
    m: Final[int]
    n: Final[int]

    def __init__(self, /, n: int, m: int | None = None, w: complex | None = None, a: complex = 1 + 0j) -> None: ...

    #
    @overload
    def __call__(self, /, x: _ToComplexStrictND, *, axis: SupportsIndex = -1) -> onp.ArrayND[_Complex]: ...
    @overload
    def __call__(self, /, x: onp.ToComplexStrict1D, *, axis: SupportsIndex = -1) -> onp.Array1D[_Complex]: ...
    @overload
    def __call__(self, /, x: onp.ToComplexStrict2D, *, axis: SupportsIndex = -1) -> onp.Array2D[_Complex]: ...
    @overload
    def __call__(self, /, x: onp.ToComplexStrict3D, *, axis: SupportsIndex = -1) -> onp.Array3D[_Complex]: ...
    @overload
    def __call__(self, /, x: onp.ToComplexND, *, axis: SupportsIndex = -1) -> onp.ArrayND[_Complex]: ...

    #
    def points(self, /) -> onp.Array1D[np.complex128]: ...

class ZoomFFT(CZT):
    f1: onp.ToFloat
    f2: onp.ToFloat
    fs: float

    def __init__(
        self, /, n: int, fn: float | onp.ToFloat1D, m: int | None = None, *, fs: float = 2, endpoint: bool = False
    ) -> None: ...

#
def _validate_sizes(n: int, m: int | None) -> int: ...

#
def czt_points(m: int, w: complex | None = None, a: complex = 1 + 0j) -> onp.Array1D[np.complex128]: ...

#
@overload
def czt(
    x: _ToComplexStrictND, m: int | None = None, w: complex | None = None, a: complex = 1 + 0j, *, axis: SupportsIndex = -1
) -> onp.ArrayND[_Complex]: ...
@overload
def czt(
    x: onp.ToComplexStrict1D, m: int | None = None, w: complex | None = None, a: complex = 1 + 0j, *, axis: SupportsIndex = -1
) -> onp.Array1D[_Complex]: ...
@overload
def czt(
    x: onp.ToComplexStrict2D, m: int | None = None, w: complex | None = None, a: complex = 1 + 0j, *, axis: SupportsIndex = -1
) -> onp.Array2D[_Complex]: ...
@overload
def czt(
    x: onp.ToComplexStrict3D, m: int | None = None, w: complex | None = None, a: complex = 1 + 0j, *, axis: SupportsIndex = -1
) -> onp.Array3D[_Complex]: ...
@overload
def czt(
    x: onp.ToComplexND, m: int | None = None, w: complex | None = None, a: complex = 1 + 0j, *, axis: SupportsIndex = -1
) -> onp.ArrayND[_Complex]: ...

#
@overload
def zoom_fft(
    x: _ToComplexStrictND,
    fn: float | onp.ToFloat1D,
    m: int | None = None,
    *,
    fs: float = 2,
    endpoint: bool = False,
    axis: SupportsIndex = -1,
) -> onp.ArrayND[_Complex]: ...
@overload
def zoom_fft(
    x: onp.ToComplexStrict1D,
    fn: float | onp.ToFloat1D,
    m: int | None = None,
    *,
    fs: float = 2,
    endpoint: bool = False,
    axis: SupportsIndex = -1,
) -> onp.Array1D[_Complex]: ...
@overload
def zoom_fft(
    x: onp.ToComplexStrict2D,
    fn: float | onp.ToFloat1D,
    m: int | None = None,
    *,
    fs: float = 2,
    endpoint: bool = False,
    axis: SupportsIndex = -1,
) -> onp.Array2D[_Complex]: ...
@overload
def zoom_fft(
    x: onp.ToComplexStrict3D,
    fn: float | onp.ToFloat1D,
    m: int | None = None,
    *,
    fs: float = 2,
    endpoint: bool = False,
    axis: SupportsIndex = -1,
) -> onp.Array3D[_Complex]: ...
@overload
def zoom_fft(
    x: onp.ToComplexND,
    fn: float | onp.ToFloat1D,
    m: int | None = None,
    *,
    fs: float = 2,
    endpoint: bool = False,
    axis: SupportsIndex = -1,
) -> onp.ArrayND[_Complex]: ...
