# NOTE: Due to the many false positive overlap mypy errors in this file, we disable the error code, and instead rely on pyright to
# catch incompatible overlapping overloads.

# mypy: disable-error-code=overload-overlap

from collections.abc import Callable, Sequence
from typing import Any, Literal as L, TypeAlias, TypeVar, TypedDict, overload, type_check_only

import numpy as np
import numpy_typing_compat as nptc
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._ltisys import dlti
from .windows._windows import _ToWindow
from scipy._typing import AnyShape

__all__ = [
    "choose_conv_method",
    "convolve",
    "convolve2d",
    "correlate",
    "correlate2d",
    "correlation_lags",
    "decimate",
    "deconvolve",
    "detrend",
    "envelope",
    "fftconvolve",
    "filtfilt",
    "hilbert",
    "hilbert2",
    "invres",
    "invresz",
    "lfilter",
    "lfilter_zi",
    "lfiltic",
    "medfilt",
    "medfilt2d",
    "oaconvolve",
    "order_filter",
    "resample",
    "resample_poly",
    "residue",
    "residuez",
    "sosfilt",
    "sosfilt_zi",
    "sosfiltfilt",
    "unique_roots",
    "vectorstrength",
    "wiener",
]

###

_T = TypeVar("_T")
_InexactT = TypeVar("_InexactT", bound=npc.inexact)
_NumericT = TypeVar("_NumericT", bound=npc.number | np.bool_)
_InexactT2 = TypeVar("_InexactT2", bound=np.float32 | np.float64 | np.complex64 | np.complex128)
_InexactT3 = TypeVar("_InexactT3", bound=np.float32 | np.float64 | npc.floating80 | npc.complexfloating)
_CoFloat64T = TypeVar("_CoFloat64T", bound=np.float64 | np.float32 | npc.integer)
_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])

_AnyShapeT = TypeVar("_AnyShapeT", tuple[int], tuple[int, int], tuple[int, int, int], tuple[Any, ...])
_AnyInexact64T = TypeVar("_AnyInexact64T", np.float64, np.complex128)

_Tuple2: TypeAlias = tuple[_T, _T]

_ConvMethod: TypeAlias = L["direct", "fft"]
_ToConvMethod: TypeAlias = L["auto", _ConvMethod]
_BoundaryConditions: TypeAlias = L["fill", "wrap", "symm"]
_ResidueType: TypeAlias = L["avg", "mean", "min", "minimum", "max", "maximum"]
_Domain: TypeAlias = L["time", "freq"]
_TrendType: TypeAlias = L["linear", "constant"]
_PadType: TypeAlias = L["constant", "line", "mean", "median", "maximum", "minimum", "symmetric", "reflect", "edge", "wrap"]
_FiltFiltPadType: TypeAlias = L["odd", "even", "constant"] | None
_FiltFiltMethod: TypeAlias = L["pad", "gust"]
_ResidualKind: TypeAlias = L["lowpass", "all"]
_FilterType: TypeAlias = L["iir", "fir"] | dlti

_F16_64: TypeAlias = np.float64 | np.float32 | np.float16
_C64_128: TypeAlias = np.complex128 | np.complex64

_ToResampleWindow: TypeAlias = Callable[[onp.Array1D[_InexactT]], onp.ToFloat1D] | onp.ToFloat1D | _ToWindow

# workaround for a strange bug in pyright's overlapping overload detection with `numpy<2.1`
_WorkaroundForPyright: TypeAlias = tuple[int] | tuple[Any, ...]

@type_check_only
class _ConvMeasureDict(TypedDict):
    direct: float
    fft: float

###

# undocumented
@overload
def _valfrommode(mode: L["valid"]) -> L[0]: ...
@overload
def _valfrommode(mode: L["same"]) -> L[1]: ...
@overload
def _valfrommode(mode: L["full"]) -> L[2]: ...
@overload
def _valfrommode(mode: str) -> int: ...

# undocumented
@overload
def _bvalfromboundary(boundary: L["fill", "pad"]) -> L[0]: ...
@overload
def _bvalfromboundary(boundary: L["symm", "symmetric"]) -> L[4]: ...
@overload
def _bvalfromboundary(boundary: L["wrap", "circular"]) -> L[8]: ...
@overload
def _bvalfromboundary(boundary: L["reflect"]) -> L[16]: ...
@overload
def _bvalfromboundary(boundary: str) -> int: ...

#
@overload
def choose_conv_method(
    in1: onp.ToIntND, in2: onp.ToIntND, mode: onp.ConvolveMode = "full", measure: onp.ToFalse = False
) -> L["direct"]: ...
@overload
def choose_conv_method(
    in1: onp.ToComplexND, in2: onp.ToComplexND, mode: onp.ConvolveMode = "full", measure: onp.ToFalse = False
) -> _ConvMethod: ...
@overload
def choose_conv_method(
    in1: onp.ToComplexND, in2: onp.ToComplexND, mode: onp.ConvolveMode, measure: onp.ToTrue
) -> tuple[_ConvMethod, _ConvMeasureDict]: ...
@overload
def choose_conv_method(
    in1: onp.ToComplexND, in2: onp.ToComplexND, mode: onp.ConvolveMode = "full", *, measure: onp.ToTrue
) -> tuple[_ConvMethod, _ConvMeasureDict]: ...

# NOTE: keep in sync with `correlate`
@overload  # ~bool, ~bool
def convolve(
    in1: onp.ToJustBoolND, in2: onp.ToJustBoolND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[np.bool_]: ...
@overload  # generic
def convolve(
    in1: onp.CanArray[_AnyShapeT, np.dtype[_NumericT]],
    in2: onp.CanArray[_AnyShapeT, np.dtype[_NumericT]],
    mode: onp.ConvolveMode = "full",
    method: _ToConvMethod = "auto",
) -> onp.ArrayND[_NumericT, _AnyShapeT]: ...
@overload  # ~int64, +int64
def convolve(
    in1: onp.ToJustInt64_ND, in2: onp.ToIntND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[np.int64]: ...
@overload  # +int64, ~int64
def convolve(
    in1: onp.ToIntND, in2: onp.ToJustInt64_ND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[np.int64]: ...
@overload  # ~float64, +float64
def convolve(
    in1: onp.ToJustFloat64_ND, in2: onp.ToFloat64_ND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[np.float64]: ...
@overload  # +float64, ~float64
def convolve(
    in1: onp.ToFloat64_ND, in2: onp.ToJustFloat64_ND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[np.float64]: ...
@overload  # ~complex128, +complex128
def convolve(
    in1: onp.ToJustComplex128_ND, in2: onp.ToComplex128_ND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[np.complex128]: ...
@overload  # +complex128, ~complex128
def convolve(
    in1: onp.ToComplex128_ND, in2: onp.ToJustComplex128_ND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[np.complex128]: ...
@overload  # fallback
def convolve(
    in1: onp.ToComplexND, in2: onp.ToComplexND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[Any]: ...

# NOTE: keep in sync with `convolve`
@overload  # ~bool, ~bool
def correlate(
    in1: onp.ToJustBoolND, in2: onp.ToJustBoolND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[np.bool_]: ...
@overload  # generic
def correlate(
    in1: onp.CanArray[_AnyShapeT, np.dtype[_NumericT]],
    in2: onp.CanArray[_AnyShapeT, np.dtype[_NumericT]],
    mode: onp.ConvolveMode = "full",
    method: _ToConvMethod = "auto",
) -> onp.ArrayND[_NumericT, _AnyShapeT]: ...
@overload  # ~int64, +int64
def correlate(
    in1: onp.ToJustInt64_ND, in2: onp.ToIntND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[np.int64]: ...
@overload  # +int64, ~int64
def correlate(
    in1: onp.ToIntND, in2: onp.ToJustInt64_ND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[np.int64]: ...
@overload  # ~float64, +float64
def correlate(
    in1: onp.ToJustFloat64_ND, in2: onp.ToFloat64_ND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[np.float64]: ...
@overload  # +float64, ~float64
def correlate(
    in1: onp.ToFloat64_ND, in2: onp.ToJustFloat64_ND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[np.float64]: ...
@overload  # ~complex128, +complex128
def correlate(
    in1: onp.ToJustComplex128_ND, in2: onp.ToComplex128_ND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[np.complex128]: ...
@overload  # +complex128, ~complex128
def correlate(
    in1: onp.ToComplex128_ND, in2: onp.ToJustComplex128_ND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[np.complex128]: ...
@overload  # fallback
def correlate(
    in1: onp.ToComplexND, in2: onp.ToComplexND, mode: onp.ConvolveMode = "full", method: _ToConvMethod = "auto"
) -> onp.ArrayND[Any]: ...

# NOTE: keep in sync with `correlate2d`
@overload  # generic
def convolve2d(
    in1: onp.CanArrayND[_NumericT],
    in2: onp.CanArrayND[_NumericT],
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToComplex = 0,
) -> onp.Array2D[_NumericT]: ...
@overload  # ~int64, +int64
def convolve2d(
    in1: onp.ToJustInt64_ND,
    in2: onp.ToIntND,
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToFloat = 0,
) -> onp.Array2D[np.int64]: ...
@overload  # +int64, ~int64
def convolve2d(
    in1: onp.ToIntND,
    in2: onp.ToJustInt64_ND,
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToFloat = 0,
) -> onp.Array2D[np.int64]: ...
@overload  # ~float64, +float64
def convolve2d(
    in1: onp.ToJustFloat64_ND,
    in2: onp.ToFloat64_ND,
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToFloat = 0,
) -> onp.Array2D[np.float64]: ...
@overload  # +float64, ~float64
def convolve2d(
    in1: onp.ToFloat64_ND,
    in2: onp.ToJustFloat64_ND,
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToFloat = 0,
) -> onp.Array2D[np.float64]: ...
@overload  # ~complex128, +complex128
def convolve2d(
    in1: onp.ToJustComplex128_ND,
    in2: onp.ToComplex128_ND,
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToComplex = 0,
) -> onp.Array2D[np.complex128]: ...
@overload  # +complex128, ~complex128
def convolve2d(
    in1: onp.ToComplex128_ND,
    in2: onp.ToJustComplex128_ND,
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToComplex = 0,
) -> onp.Array2D[np.complex128]: ...
@overload  # fallback
def convolve2d(
    in1: onp.ToComplex2D,
    in2: onp.ToComplex2D,
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToComplex = 0,
) -> onp.Array2D[Any]: ...

# NOTE: keep in sync with `convolve2d`
@overload  # generic
def correlate2d(
    in1: onp.CanArrayND[_NumericT],
    in2: onp.CanArrayND[_NumericT],
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToComplex = 0,
) -> onp.Array2D[_NumericT]: ...
@overload  # ~int64, +int64
def correlate2d(
    in1: onp.ToJustInt64_ND,
    in2: onp.ToIntND,
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToFloat = 0,
) -> onp.Array2D[np.int64]: ...
@overload  # +int64, ~int64
def correlate2d(
    in1: onp.ToIntND,
    in2: onp.ToJustInt64_ND,
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToFloat = 0,
) -> onp.Array2D[np.int64]: ...
@overload  # ~float64, +float64
def correlate2d(
    in1: onp.ToJustFloat64_ND,
    in2: onp.ToFloat64_ND,
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToFloat = 0,
) -> onp.Array2D[np.float64]: ...
@overload  # +float64, ~float64
def correlate2d(
    in1: onp.ToFloat64_ND,
    in2: onp.ToJustFloat64_ND,
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToFloat = 0,
) -> onp.Array2D[np.float64]: ...
@overload  # ~complex128, +complex128
def correlate2d(
    in1: onp.ToJustComplex128_ND,
    in2: onp.ToComplex128_ND,
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToComplex = 0,
) -> onp.Array2D[np.complex128]: ...
@overload  # +complex128, ~complex128
def correlate2d(
    in1: onp.ToComplex128_ND,
    in2: onp.ToJustComplex128_ND,
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToComplex = 0,
) -> onp.Array2D[np.complex128]: ...
@overload  # fallback
def correlate2d(
    in1: onp.ToComplex2D,
    in2: onp.ToComplex2D,
    mode: onp.ConvolveMode = "full",
    boundary: _BoundaryConditions = "fill",
    fillvalue: onp.ToComplex = 0,
) -> onp.Array2D[Any]: ...

# NOTE: keep in sync with `oaconvolve`
@overload  # float32 | float16, float32 | float16, generic shape
def fftconvolve(
    in1: onp.ArrayND[np.float16 | np.float32, _AnyShapeT],
    in2: onp.ArrayND[np.float16 | np.float32, _AnyShapeT],
    mode: onp.ConvolveMode = "full",
    axes: None = None,
) -> onp.ArrayND[np.float32, _AnyShapeT]: ...
@overload  # generic dtype, generic, shape
def fftconvolve(
    in1: onp.ArrayND[_InexactT3, _AnyShapeT],
    in2: onp.ArrayND[_InexactT3, _AnyShapeT],
    mode: onp.ConvolveMode = "full",
    axes: None = None,
) -> onp.ArrayND[_InexactT3, _AnyShapeT]: ...
@overload  # ~float64, +float64
def fftconvolve(
    in1: onp.ToJustFloat64_ND, in2: onp.ToFloat64_ND, mode: onp.ConvolveMode = "full", axes: AnyShape | None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # ~float64, +float64
def fftconvolve(
    in1: onp.ToFloat64_ND, in2: onp.ToJustFloat64_ND, mode: onp.ConvolveMode = "full", axes: AnyShape | None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # ~complex128, +complex128
def fftconvolve(
    in1: onp.ToJustComplex128_ND, in2: onp.ToComplex128_ND, mode: onp.ConvolveMode = "full", axes: AnyShape | None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # ~complex128, +complex128
def fftconvolve(
    in1: onp.ToComplex128_ND, in2: onp.ToJustComplex128_ND, mode: onp.ConvolveMode = "full", axes: AnyShape | None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # fallback
def fftconvolve(
    in1: onp.ToComplexND, in2: onp.ToComplexND, mode: onp.ConvolveMode = "full", axes: AnyShape | None = None
) -> onp.ArrayND[Any, _WorkaroundForPyright]: ...

# NOTE: keep in sync with `fftconvolve`
@overload  # float32 | float16, float32 | float16, generic shape
def oaconvolve(
    in1: onp.ArrayND[np.float16 | np.float32, _AnyShapeT],
    in2: onp.ArrayND[np.float16 | np.float32, _AnyShapeT],
    mode: onp.ConvolveMode = "full",
    axes: None = None,
) -> onp.ArrayND[np.float32, _AnyShapeT]: ...
@overload  # generic dtype, generic, shape
def oaconvolve(
    in1: onp.ArrayND[_InexactT3, _AnyShapeT],
    in2: onp.ArrayND[_InexactT3, _AnyShapeT],
    mode: onp.ConvolveMode = "full",
    axes: None = None,
) -> onp.ArrayND[_InexactT3, _AnyShapeT]: ...
@overload  # ~float64, +float64
def oaconvolve(
    in1: onp.ToJustFloat64_ND, in2: onp.ToFloat64_ND, mode: onp.ConvolveMode = "full", axes: AnyShape | None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # ~float64, +float64
def oaconvolve(
    in1: onp.ToFloat64_ND, in2: onp.ToJustFloat64_ND, mode: onp.ConvolveMode = "full", axes: AnyShape | None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # ~complex128, +complex128
def oaconvolve(
    in1: onp.ToJustComplex128_ND, in2: onp.ToComplex128_ND, mode: onp.ConvolveMode = "full", axes: AnyShape | None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # ~complex128, +complex128
def oaconvolve(
    in1: onp.ToComplex128_ND, in2: onp.ToJustComplex128_ND, mode: onp.ConvolveMode = "full", axes: AnyShape | None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # fallback
def oaconvolve(
    in1: onp.ToComplexND, in2: onp.ToComplexND, mode: onp.ConvolveMode = "full", axes: AnyShape | None = None
) -> onp.ArrayND[Any, _WorkaroundForPyright]: ...

#
@overload  # +float64, +float64
def deconvolve(signal: onp.ToFloat64_1D, divisor: onp.ToFloat64_1D) -> _Tuple2[onp.Array1D[np.float64]]: ...
@overload  # ~longdouble, +floating
def deconvolve(signal: onp.ToJustLongDouble1D, divisor: onp.ToFloat1D) -> _Tuple2[onp.Array1D[npc.floating80]]: ...
@overload  # +floating, ~longdouble
def deconvolve(signal: onp.ToFloat1D, divisor: onp.ToJustLongDouble1D) -> _Tuple2[onp.Array1D[npc.floating80]]: ...
@overload  # ~complex128 | ~complex64, +complex128
def deconvolve(
    signal: onp.ToArray1D[op.JustComplex, _C64_128], divisor: onp.ToComplex128_1D
) -> _Tuple2[onp.Array1D[np.complex128]]: ...
@overload  # +complex128, ~complex128 | ~complex64
def deconvolve(
    signal: onp.ToComplex128_1D, divisor: onp.ToArray1D[op.JustComplex, _C64_128]
) -> _Tuple2[onp.Array1D[np.complex128]]: ...
@overload  # ~clongdouble, +complexfloating
def deconvolve(signal: onp.ToJustCLongDouble1D, divisor: onp.ToComplex1D) -> _Tuple2[onp.Array1D[np.clongdouble]]: ...
@overload  # +complexfloating, ~clongdouble
def deconvolve(signal: onp.ToComplex1D, divisor: onp.ToJustCLongDouble1D) -> _Tuple2[onp.Array1D[np.clongdouble]]: ...
@overload  # fallback
def deconvolve(signal: onp.ToComplex1D, divisor: onp.ToComplex1D) -> _Tuple2[onp.Array1D[Any]]: ...

#
def correlation_lags(in1_len: int, in2_len: int, mode: onp.ConvolveMode = "full") -> onp.Array1D[np.int_]: ...

###

#
@overload  # +float64, ~float64
def lfilter_zi(b: onp.ToJustFloat64_1D, a: onp.ToFloat64_1D) -> onp.Array1D[np.float64]: ...
@overload  # ~float64, +float64
def lfilter_zi(b: onp.ToFloat64_1D, a: onp.ToJustFloat64_1D) -> onp.Array1D[np.float64]: ...
@overload  # ~float32, +float32
def lfilter_zi(b: onp.ToJustFloat32_1D, a: onp.ToFloat32_1D) -> onp.Array1D[np.float32]: ...
@overload  # +float32, ~float32
def lfilter_zi(b: onp.ToFloat32_1D, a: onp.ToJustFloat32_1D) -> onp.Array1D[np.float32]: ...
@overload  # +complex128, ~complex128
def lfilter_zi(b: onp.ToJustComplex128_1D, a: onp.ToComplex128_1D) -> onp.Array1D[np.complex128]: ...
@overload  # ~complex128, +complex128
def lfilter_zi(b: onp.ToComplex128_1D, a: onp.ToJustComplex128_1D) -> onp.Array1D[np.complex128]: ...
@overload  # +complex64, ~complex64
def lfilter_zi(b: onp.ToJustComplex64_1D, a: onp.ToComplex64_1D) -> onp.Array1D[np.complex64]: ...
@overload  # ~complex64, +complex64
def lfilter_zi(b: onp.ToComplex64_1D, a: onp.ToJustComplex64_1D) -> onp.Array1D[np.complex64]: ...
@overload  # fallback
def lfilter_zi(b: onp.ToComplex1D, a: onp.ToComplex1D) -> onp.Array1D[Any]: ...

#
@overload  # ~float64 | integer, +float64, +float64, +float64 | None
def lfiltic(
    b: onp.ToArray1D[float, np.float64 | npc.integer], a: onp.ToFloat64_1D, y: onp.ToFloat64_1D, x: onp.ToFloat64_1D | None = None
) -> onp.Array1D[np.float64]: ...
@overload  # +float64, ~float64 | integer, +float64, +float64 | None
def lfiltic(
    b: onp.ToFloat64_1D, a: onp.ToArray1D[float, np.float64 | npc.integer], y: onp.ToFloat64_1D, x: onp.ToFloat64_1D | None = None
) -> onp.Array1D[np.float64]: ...
@overload  # +float64, +float64, ~float64 | integer, +float64 | None
def lfiltic(
    b: onp.ToFloat64_1D, a: onp.ToFloat64_1D, y: onp.ToArray1D[float, np.float64 | npc.integer], x: onp.ToFloat64_1D | None = None
) -> onp.Array1D[np.float64]: ...
@overload  # ~float32, +float32, +float32, +float32 | None
def lfiltic(
    b: onp.ToJustFloat32_1D, a: onp.ToFloat32_1D, y: onp.ToFloat32_1D, x: onp.ToFloat32_1D | None = None
) -> onp.Array1D[np.float32]: ...
@overload  # +float32, ~float32, +float32, +float32 | None
def lfiltic(
    b: onp.ToFloat32_1D, a: onp.ToJustFloat32_1D, y: onp.ToFloat32_1D, x: onp.ToFloat32_1D | None = None
) -> onp.Array1D[np.float32]: ...
@overload  # +float32, +float32, ~float32, +float32 | None
def lfiltic(
    b: onp.ToFloat32_1D, a: onp.ToFloat32_1D, y: onp.ToJustFloat32_1D, x: onp.ToFloat32_1D | None = None
) -> onp.Array1D[np.float32]: ...
@overload  # ~complex128, +complex128, +complex128, +complex128 | None
def lfiltic(
    b: onp.ToJustComplex128_1D, a: onp.ToComplex128_1D, y: onp.ToComplex128_1D, x: onp.ToComplex128_1D | None = None
) -> onp.Array1D[np.complex128]: ...
@overload  # +complex128, ~complex128, +complex128, +complex128 | None
def lfiltic(
    b: onp.ToComplex128_1D, a: onp.ToJustComplex128_1D, y: onp.ToComplex128_1D, x: onp.ToComplex128_1D | None = None
) -> onp.Array1D[np.complex128]: ...
@overload  # +complex128, +complex128, ~complex128, +complex128 | None
def lfiltic(
    b: onp.ToComplex128_1D, a: onp.ToComplex128_1D, y: onp.ToJustComplex128_1D, x: onp.ToComplex128_1D | None = None
) -> onp.Array1D[np.complex128]: ...
@overload  # ~complex64, +complex64, +complex64, +complex64 | None
def lfiltic(
    b: onp.ToJustComplex64_1D, a: onp.ToComplex64_1D, y: onp.ToComplex64_1D, x: onp.ToComplex64_1D | None = None
) -> onp.Array1D[np.complex64]: ...
@overload  # +complex64, ~complex64, +complex64, +complex64 | None
def lfiltic(
    b: onp.ToComplex64_1D, a: onp.ToJustComplex64_1D, y: onp.ToComplex64_1D, x: onp.ToComplex64_1D | None = None
) -> onp.Array1D[np.complex64]: ...
@overload  # +complex64, +complex64, ~complex64, +complex64 | None
def lfiltic(
    b: onp.ToComplex64_1D, a: onp.ToComplex64_1D, y: onp.ToJustComplex64_1D, x: onp.ToComplex64_1D | None = None
) -> onp.Array1D[np.complex64]: ...
@overload  # fallback
def lfiltic(b: onp.ToComplex1D, a: onp.ToComplex1D, y: onp.ToComplex1D, x: onp.ToComplex1D | None = None) -> onp.Array1D[Any]: ...

#
@overload  # ~float64 | integer, +float64, +float64, zi: None (default)
def lfilter(
    b: onp.ToArray1D[float, np.float64 | npc.integer], a: onp.ToFloat64_1D, x: onp.ToFloat64_ND, axis: int = -1, zi: None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # ~float64 | integer, +float64, +float64, *, zi: +float64
def lfilter(
    b: onp.ToArray1D[float, np.float64 | npc.integer],
    a: onp.ToFloat64_1D,
    x: onp.ToFloat64_ND,
    axis: int = -1,
    *,
    zi: onp.ToFloat64_ND,
) -> _Tuple2[onp.ArrayND[np.float64]]: ...
@overload  # +float64, ~float64 | integer, +float64, zi: None (default)
def lfilter(
    b: onp.ToFloat64_1D, a: onp.ToArray1D[float, np.float64 | npc.integer], x: onp.ToFloat64_ND, axis: int = -1, zi: None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # +float64, ~float64 | integer, +float64, zi: *, zi: +float64
def lfilter(
    b: onp.ToFloat64_1D,
    a: onp.ToArray1D[float, np.float64 | npc.integer],
    x: onp.ToFloat64_ND,
    axis: int = -1,
    *,
    zi: onp.ToFloat64_ND,
) -> _Tuple2[onp.ArrayND[np.float64]]: ...
@overload  # +float64, +float64, ~float64 | integer, zi: None (default)
def lfilter(
    b: onp.ToFloat64_1D, a: onp.ToFloat64_1D, x: onp.ToArrayND[float, np.float64 | npc.integer], axis: int = -1, zi: None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # +float64, +float64, ~float64 | integer, *, zi: +float64
def lfilter(
    b: onp.ToFloat64_1D,
    a: onp.ToFloat64_1D,
    x: onp.ToArrayND[float, np.float64 | npc.integer],
    axis: int = -1,
    *,
    zi: onp.ToFloat64_ND,
) -> _Tuple2[onp.ArrayND[np.float64]]: ...
@overload  # ~float32, +float32, +float32, zi: None (default)
def lfilter(
    b: onp.ToJustFloat32_1D, a: onp.ToFloat32_1D, x: onp.ToFloat32_ND, axis: int = -1, zi: None = None
) -> onp.ArrayND[np.float32]: ...
@overload  # ~float32, +float32, +float32, *, zi: +float32
def lfilter(
    b: onp.ToJustFloat32_1D, a: onp.ToFloat32_1D, x: onp.ToFloat32_ND, axis: int = -1, *, zi: onp.ToFloat32_ND
) -> _Tuple2[onp.ArrayND[np.float32]]: ...
@overload  # +float32, ~float32, +float32, zi: None (default)
def lfilter(
    b: onp.ToFloat32_1D, a: onp.ToJustFloat32_1D, x: onp.ToFloat32_ND, axis: int = -1, zi: None = None
) -> onp.ArrayND[np.float32]: ...
@overload  # +float32, ~float32, +float32, *, zi: +float32
def lfilter(
    b: onp.ToFloat32_1D, a: onp.ToJustFloat32_1D, x: onp.ToFloat32_ND, axis: int = -1, *, zi: onp.ToFloat32_ND
) -> _Tuple2[onp.ArrayND[np.float32]]: ...
@overload  # +float32, +float32, ~float32, zi: None (default)
def lfilter(
    b: onp.ToFloat32_1D, a: onp.ToFloat32_1D, x: onp.ToJustFloat32_ND, axis: int = -1, zi: None = None
) -> onp.ArrayND[np.float32]: ...
@overload  # +float32, +float32, ~float32, *, zi: +float32
def lfilter(
    b: onp.ToFloat32_1D, a: onp.ToFloat32_1D, x: onp.ToJustFloat32_ND, axis: int = -1, *, zi: onp.ToFloat32_ND
) -> _Tuple2[onp.ArrayND[np.float32]]: ...
@overload  # ~complex128, +complex128, +complex128, zi: None (default)
def lfilter(
    b: onp.ToJustComplex128_1D, a: onp.ToComplex128_1D, x: onp.ToComplex128_ND, axis: int = -1, zi: None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # ~complex128, +complex128, +complex128, *, zi: +complex128
def lfilter(
    b: onp.ToJustComplex128_1D, a: onp.ToComplex128_1D, x: onp.ToComplex128_ND, axis: int = -1, *, zi: onp.ToComplex128_ND
) -> _Tuple2[onp.ArrayND[np.complex128]]: ...
@overload  # +complex128, ~complex128, +complex128, zi: None (default)
def lfilter(
    b: onp.ToComplex128_1D, a: onp.ToJustComplex128_1D, x: onp.ToComplex128_ND, axis: int = -1, zi: None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # +complex128, ~complex128, +complex128, *, zi: +complex128
def lfilter(
    b: onp.ToComplex128_1D, a: onp.ToJustComplex128_1D, x: onp.ToComplex128_ND, axis: int = -1, *, zi: onp.ToComplex128_ND
) -> _Tuple2[onp.ArrayND[np.complex128]]: ...
@overload  # +complex128, +complex128, ~complex128, zi: None (default)
def lfilter(
    b: onp.ToComplex128_1D, a: onp.ToComplex128_1D, x: onp.ToJustComplex128_ND, axis: int = -1, zi: None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # +complex128, +complex128, ~complex128, *, zi: +complex128
def lfilter(
    b: onp.ToComplex128_1D, a: onp.ToComplex128_1D, x: onp.ToJustComplex128_ND, axis: int = -1, *, zi: onp.ToComplex128_ND
) -> _Tuple2[onp.ArrayND[np.complex128]]: ...
@overload  # ~complex64, +complex64, +complex64, zi: None (default)
def lfilter(
    b: onp.ToJustComplex64_1D, a: onp.ToComplex64_1D, x: onp.ToComplex64_ND, axis: int = -1, zi: None = None
) -> onp.ArrayND[np.complex64]: ...
@overload  # ~complex64, +complex64, +complex64, *, zi: +complex64
def lfilter(
    b: onp.ToJustComplex64_1D, a: onp.ToComplex64_1D, x: onp.ToComplex64_ND, axis: int = -1, *, zi: onp.ToComplex64_ND
) -> _Tuple2[onp.ArrayND[np.complex64]]: ...
@overload  # +complex64, ~complex64, +complex64, zi: None (default)
def lfilter(
    b: onp.ToComplex64_1D, a: onp.ToJustComplex64_1D, x: onp.ToComplex64_ND, axis: int = -1, zi: None = None
) -> onp.ArrayND[np.complex64]: ...
@overload  # +complex64, ~complex64, +complex64, *, zi: +complex64
def lfilter(
    b: onp.ToComplex64_1D, a: onp.ToJustComplex64_1D, x: onp.ToComplex64_ND, axis: int = -1, *, zi: onp.ToComplex64_ND
) -> _Tuple2[onp.ArrayND[np.complex64]]: ...
@overload  # +complex64, +complex64, ~complex64, zi: None (default)
def lfilter(
    b: onp.ToComplex64_1D, a: onp.ToComplex64_1D, x: onp.ToJustComplex64_ND, axis: int = -1, zi: None = None
) -> onp.ArrayND[np.complex64]: ...
@overload  # +complex64, +complex64, ~complex64, *, zi: +complex64
def lfilter(
    b: onp.ToComplex64_1D, a: onp.ToComplex64_1D, x: onp.ToJustComplex64_ND, axis: int = -1, *, zi: onp.ToComplex64_ND
) -> _Tuple2[onp.ArrayND[np.complex64]]: ...
@overload  # fallback, zi: None (default)
def lfilter(b: onp.ToComplex1D, a: onp.ToComplex1D, x: onp.ToComplexND, axis: int = -1, zi: None = None) -> onp.ArrayND[Any]: ...
@overload  # fallback, *, zi: +complex
def lfilter(
    b: onp.ToComplex1D, a: onp.ToComplex1D, x: onp.ToComplexND, axis: int = -1, *, zi: onp.ToComplexND
) -> _Tuple2[onp.ArrayND[Any]]: ...

#
@overload  # ~float64 | integer, +float64, +float64
def filtfilt(
    b: onp.ToArray1D[float, np.float64 | npc.integer],
    a: onp.ToFloat64_1D,
    x: onp.ToFloat64_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
    method: _FiltFiltMethod = "pad",
    irlen: int | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # +float64, ~float64 | integer, +float64
def filtfilt(
    b: onp.ToFloat64_1D,
    a: onp.ToArray1D[float, np.float64 | npc.integer],
    x: onp.ToFloat64_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
    method: _FiltFiltMethod = "pad",
    irlen: int | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # +float64, +float64, ~float64 | integer
def filtfilt(
    b: onp.ToFloat64_1D,
    a: onp.ToFloat64_1D,
    x: onp.ToArrayND[float, np.float64 | npc.integer],
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
    method: _FiltFiltMethod = "pad",
    irlen: int | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # ~float32, +float32, +float32
def filtfilt(
    b: onp.ToJustFloat32_1D,
    a: onp.ToFloat32_1D,
    x: onp.ToFloat32_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
    method: _FiltFiltMethod = "pad",
    irlen: int | None = None,
) -> onp.ArrayND[np.float32]: ...
@overload  # +float32, ~float32, +float32
def filtfilt(
    b: onp.ToFloat32_1D,
    a: onp.ToJustFloat32_1D,
    x: onp.ToFloat32_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
    method: _FiltFiltMethod = "pad",
    irlen: int | None = None,
) -> onp.ArrayND[np.float32]: ...
@overload  # +float32, +float32, ~float32
def filtfilt(
    b: onp.ToFloat32_1D,
    a: onp.ToFloat32_1D,
    x: onp.ToJustFloat32_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
    method: _FiltFiltMethod = "pad",
    irlen: int | None = None,
) -> onp.ArrayND[np.float32]: ...
@overload  # ~complex128, +complex128, +complex128
def filtfilt(
    b: onp.ToJustComplex128_1D,
    a: onp.ToComplex128_1D,
    x: onp.ToComplex128_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
    method: _FiltFiltMethod = "pad",
    irlen: int | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload  # +complex128, ~complex128, +complex128
def filtfilt(
    b: onp.ToComplex128_1D,
    a: onp.ToJustComplex128_1D,
    x: onp.ToComplex128_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
    method: _FiltFiltMethod = "pad",
    irlen: int | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload  # +complex128, +complex128, ~complex128
def filtfilt(
    b: onp.ToComplex128_1D,
    a: onp.ToComplex128_1D,
    x: onp.ToJustComplex128_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
    method: _FiltFiltMethod = "pad",
    irlen: int | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload  # ~complex64, +complex64, +complex64
def filtfilt(
    b: onp.ToJustComplex64_1D,
    a: onp.ToComplex64_1D,
    x: onp.ToComplex64_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
    method: _FiltFiltMethod = "pad",
    irlen: int | None = None,
) -> onp.ArrayND[np.complex64]: ...
@overload  # +complex64, ~complex64, +complex64
def filtfilt(
    b: onp.ToComplex64_1D,
    a: onp.ToJustComplex64_1D,
    x: onp.ToComplex64_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
    method: _FiltFiltMethod = "pad",
    irlen: int | None = None,
) -> onp.ArrayND[np.complex64]: ...
@overload  # +complex64, +complex64, ~complex64
def filtfilt(
    b: onp.ToComplex64_1D,
    a: onp.ToComplex64_1D,
    x: onp.ToJustComplex64_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
    method: _FiltFiltMethod = "pad",
    irlen: int | None = None,
) -> onp.ArrayND[np.complex64]: ...
@overload  # fallback
def filtfilt(
    b: onp.ToComplex1D,
    a: onp.ToComplex1D,
    x: onp.ToComplexND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
    method: _FiltFiltMethod = "pad",
    irlen: int | None = None,
) -> onp.ArrayND[Any]: ...

#
@overload  # T -> T
def sosfilt_zi(sos: onp.CanArrayND[_InexactT]) -> onp.Array2D[_InexactT]: ...
@overload  # +float -> float64
def sosfilt_zi(sos: onp.ToArray2D[float, np.float64 | npc.integer]) -> onp.Array2D[np.float64]: ...
@overload  # ~complex -> complex128
def sosfilt_zi(sos: onp.ToJustComplex128_2D) -> onp.Array2D[np.complex128]: ...
@overload  # fallback
def sosfilt_zi(sos: onp.ToComplex2D) -> onp.Array2D[Any]: ...

#
@overload  # ~float64 | integer, +float64, zi: None (default)
def sosfilt(
    sos: onp.ToArray2D[float, np.float64 | npc.integer], x: onp.ToFloat64_ND, axis: int = -1, zi: None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # ~float64 | integer, +float64, *, zi: +float64
def sosfilt(
    sos: onp.ToArray2D[float, np.float64 | npc.integer], x: onp.ToFloat64_ND, *, axis: int = -1, zi: onp.ToFloat64_ND
) -> _Tuple2[onp.ArrayND[np.float64]]: ...
@overload  # +float64, ~float64 | integer, zi: None (default)
def sosfilt(
    sos: onp.ToFloat64_2D, x: onp.ToArrayND[float, np.float64 | npc.integer], axis: int = -1, zi: None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # +float64, ~float64 | integer, *, zi: +float64
def sosfilt(
    sos: onp.ToFloat64_2D, x: onp.ToArrayND[float, np.float64 | npc.integer], *, axis: int = -1, zi: onp.ToFloat64_ND
) -> _Tuple2[onp.ArrayND[np.float64]]: ...
@overload  # ~float32, +float32, zi: None (default)
def sosfilt(sos: onp.ToJustFloat32_2D, x: onp.ToFloat32_ND, axis: int = -1, zi: None = None) -> onp.ArrayND[np.float32]: ...
@overload  # ~float32, +float32, *, zi: +float32
def sosfilt(
    sos: onp.ToJustFloat32_2D, x: onp.ToFloat32_ND, *, axis: int = -1, zi: onp.ToFloat32_ND
) -> _Tuple2[onp.ArrayND[np.float32]]: ...
@overload  # +float32, ~float32, zi: None (default)
def sosfilt(sos: onp.ToFloat32_2D, x: onp.ToJustFloat32_ND, axis: int = -1, zi: None = None) -> onp.ArrayND[np.float32]: ...
@overload  # +float32, ~float32, *, zi: +float32
def sosfilt(
    sos: onp.ToFloat32_2D, x: onp.ToJustFloat32_ND, *, axis: int = -1, zi: onp.ToFloat32_ND
) -> _Tuple2[onp.ArrayND[np.float32]]: ...
@overload  # ~complex128, +complex128, zi: None (default)
def sosfilt(
    sos: onp.ToJustComplex128_2D, x: onp.ToComplex128_ND, axis: int = -1, zi: None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # ~complex128, +complex128, *, zi: +complex128
def sosfilt(
    sos: onp.ToJustComplex128_2D, x: onp.ToComplex128_ND, *, axis: int = -1, zi: onp.ToComplex128_ND
) -> _Tuple2[onp.ArrayND[np.complex128]]: ...
@overload  # +complex128, ~complex128, zi: None (default)
def sosfilt(
    sos: onp.ToComplex128_2D, x: onp.ToJustComplex128_ND, axis: int = -1, zi: None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # +complex128, ~complex128, *, zi: +complex128
def sosfilt(
    sos: onp.ToComplex128_2D, x: onp.ToJustComplex128_ND, *, axis: int = -1, zi: onp.ToComplex128_ND
) -> _Tuple2[onp.ArrayND[np.complex128]]: ...
@overload  # ~complex64, +complex64, zi: None (default)
def sosfilt(sos: onp.ToJustComplex64_2D, x: onp.ToComplex64_ND, axis: int = -1, zi: None = None) -> onp.ArrayND[np.complex64]: ...
@overload  # ~complex64, +complex64, *, zi: +complex64
def sosfilt(
    sos: onp.ToJustComplex64_2D, x: onp.ToComplex64_ND, *, axis: int = -1, zi: onp.ToComplex64_ND
) -> _Tuple2[onp.ArrayND[np.complex64]]: ...
@overload  # +complex64, ~complex64, zi: None (default)
def sosfilt(sos: onp.ToComplex64_2D, x: onp.ToJustComplex64_ND, axis: int = -1, zi: None = None) -> onp.ArrayND[np.complex64]: ...
@overload  # +complex64, ~complex64, *, zi: +complex64
def sosfilt(
    sos: onp.ToComplex64_2D, x: onp.ToJustComplex64_ND, *, axis: int = -1, zi: onp.ToComplex64_ND
) -> _Tuple2[onp.ArrayND[np.complex64]]: ...

#
@overload  # ~float64 | integer, +float64
def sosfiltfilt(
    sos: onp.ToArray2D[float, np.float64 | npc.integer],
    x: onp.ToFloat64_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # +float64, ~float64 | integer
def sosfiltfilt(
    sos: onp.ToFloat64_2D,
    x: onp.ToArrayND[float, np.float64 | npc.integer],
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # ~float32, +float32
def sosfiltfilt(
    sos: onp.ToJustFloat32_2D, x: onp.ToFloat32_ND, axis: int = -1, padtype: _FiltFiltPadType = "odd", padlen: int | None = None
) -> onp.ArrayND[np.float32]: ...
@overload  # +float32, ~float32
def sosfiltfilt(
    sos: onp.ToFloat32_2D, x: onp.ToJustFloat32_ND, axis: int = -1, padtype: _FiltFiltPadType = "odd", padlen: int | None = None
) -> onp.ArrayND[np.float32]: ...
@overload  # ~complex128, +complex128
def sosfiltfilt(
    sos: onp.ToJustComplex128_2D,
    x: onp.ToComplex128_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload  # +complex128, ~complex128
def sosfiltfilt(
    sos: onp.ToComplex128_2D,
    x: onp.ToJustComplex128_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload  # ~complex64, +complex64
def sosfiltfilt(
    sos: onp.ToJustComplex64_2D,
    x: onp.ToComplex64_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
) -> onp.ArrayND[np.complex64]: ...
@overload  # +complex64, ~complex64
def sosfiltfilt(
    sos: onp.ToComplex64_2D,
    x: onp.ToJustComplex64_ND,
    axis: int = -1,
    padtype: _FiltFiltPadType = "odd",
    padlen: int | None = None,
) -> onp.ArrayND[np.complex64]: ...
@overload  # fallback
def sosfiltfilt(
    sos: onp.ToComplex2D, x: onp.ToComplexND, axis: int = -1, padtype: _FiltFiltPadType = "odd", padlen: int | None = None
) -> onp.ArrayND[Any]: ...

#
@overload
def order_filter(a: onp.ToJustInt64_ND, domain: int | onp.ToIntND, rank: int) -> onp.ArrayND[np.int_]: ...
@overload
def order_filter(a: onp.ToJustFloat64_ND, domain: int | onp.ToIntND, rank: int) -> onp.ArrayND[np.float64]: ...
@overload
def order_filter(
    a: onp.CanArray[_ShapeT, np.dtype[_CoFloat64T]], domain: int | onp.ToIntND, rank: int
) -> onp.ArrayND[_CoFloat64T, _ShapeT]: ...
@overload
def order_filter(a: onp.ToFloatND, domain: int | onp.ToIntND, rank: int) -> onp.ArrayND[Any]: ...

#
@overload
def medfilt(volume: onp.ToJustInt64_ND, kernel_size: int | onp.ToInt1D | None = None) -> onp.ArrayND[np.int_]: ...
@overload
def medfilt(volume: onp.ToJustFloat64_ND, kernel_size: int | onp.ToInt1D | None = None) -> onp.ArrayND[np.float64]: ...
@overload
def medfilt(
    volume: onp.CanArray[_ShapeT, np.dtype[_CoFloat64T]], kernel_size: int | onp.ToInt1D | None = None
) -> onp.ArrayND[_CoFloat64T, _ShapeT]: ...
@overload
def medfilt(volume: onp.ToFloatND, kernel_size: int | onp.ToInt1D | None = None) -> onp.ArrayND[Any]: ...

#
@overload
def medfilt2d(input: onp.ToJustInt64_2D, kernel_size: int | onp.ToInt1D = 3) -> onp.Array2D[np.int_]: ...
@overload
def medfilt2d(input: onp.ToJustFloat64_2D, kernel_size: int | onp.ToInt1D = 3) -> onp.Array2D[np.float64]: ...
@overload
def medfilt2d(input: onp.CanArrayND[_CoFloat64T], kernel_size: int | onp.ToInt1D = 3) -> onp.Array2D[_CoFloat64T]: ...
@overload
def medfilt2d(input: onp.ToFloat2D, kernel_size: int | onp.ToInt1D = 3) -> onp.Array2D[Any]: ...

#
@overload  # +float64
def wiener(
    im: onp.ToFloat64_ND, mysize: int | onp.ToInt1D | None = None, noise: float | None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # ~longdouble
def wiener(
    im: onp.ToJustLongDoubleND, mysize: int | onp.ToInt1D | None = None, noise: float | None = None
) -> onp.ArrayND[npc.floating80]: ...
@overload  # ~complex128 | ~complex64
def wiener(
    im: onp.ToArrayND[op.JustComplex, np.complex128 | np.complex64],
    mysize: int | onp.ToInt1D | None = None,
    noise: float | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload  # ~clongdouble
def wiener(
    im: onp.ToJustCLongDoubleND, mysize: int | onp.ToInt1D | None = None, noise: float | None = None
) -> onp.ArrayND[np.clongdouble]: ...
@overload  # fallback
def wiener(im: onp.ToComplexND, mysize: int | onp.ToInt1D | None = None, noise: float | None = None) -> onp.ArrayND[Any]: ...

#
@overload  # float64 | integer, known shape
def hilbert(
    x: onp.CanArray[_ShapeT, np.dtype[np.float64 | npc.integer]], N: int | None = None, axis: int = -1
) -> onp.ArrayND[np.complex128, _ShapeT]: ...
@overload  # float32 | float16, known shape
def hilbert(
    x: onp.CanArray[_ShapeT, np.dtype[np.float32 | np.float16]], N: int | None = None, axis: int = -1
) -> onp.ArrayND[np.complex64, _ShapeT]: ...
@overload  # longdouble, known shape
def hilbert(
    x: onp.CanArray[_ShapeT, np.dtype[npc.floating80]], N: int | None = None, axis: int = -1
) -> onp.ArrayND[np.clongdouble, _ShapeT]: ...
@overload  # float64 | integer, unknown shape
def hilbert(
    x: onp.ToArrayND[float, np.float64 | npc.integer], N: int | None = None, axis: int = -1
) -> onp.ArrayND[np.complex128]: ...
@overload  # fallback
def hilbert(x: onp.ToFloatND, N: int | None = None, axis: int = -1) -> onp.ArrayND[npc.complexfloating]: ...

#
@overload  # float64 | integer
def hilbert2(
    x: onp.ToArray2D[float, np.float64 | npc.integer], N: int | tuple[int, int] | None = None
) -> onp.Array2D[np.complex128]: ...
@overload  # float32 | float16
def hilbert2(
    x: onp.ToJustFloat32_2D | onp.ToJustFloat16_2D, N: int | tuple[int, int] | None = None
) -> onp.Array2D[np.complex64]: ...
@overload  # longdouble
def hilbert2(x: onp.ToJustLongDouble2D, N: int | tuple[int, int] | None = None) -> onp.Array2D[np.clongdouble]: ...
@overload  # fallback
def hilbert2(x: onp.ToFloat2D, N: int | tuple[int, int] | None = None) -> onp.Array2D[npc.complexfloating]: ...

###

# TODO(jorenham): improve
@overload
def detrend(
    data: onp.ToFloatND, axis: int = -1, type: _TrendType = "linear", bp: int | onp.ToJustIntND = 0, overwrite_data: bool = False
) -> onp.ArrayND[_F16_64]: ...
@overload
def detrend(
    data: onp.ToComplexND,
    axis: int = -1,
    type: _TrendType = "linear",
    bp: int | onp.ToJustIntND = 0,
    overwrite_data: bool = False,
) -> onp.ArrayND[_C64_128 | _F16_64]: ...

# TODO(jorenham): improve
@overload
def unique_roots(
    p: onp.ToFloat1D, tol: float = 0.001, rtype: _ResidueType = "min"
) -> tuple[onp.Array1D[npc.floating], onp.Array1D[np.int_]]: ...
@overload
def unique_roots(
    p: onp.ToComplex1D, tol: float = 0.001, rtype: _ResidueType = "min"
) -> tuple[onp.Array1D[npc.inexact], onp.Array1D[np.int_]]: ...

#
def residue(
    b: onp.ToComplex1D, a: onp.ToComplex1D, tol: float = 0.001, rtype: _ResidueType = "avg"
) -> tuple[onp.Array1D[np.complex128], onp.Array1D[np.complex128], onp.Array1D[np.float64]]: ...
def residuez(
    b: onp.ToComplex1D, a: onp.ToComplex1D, tol: float = 0.001, rtype: _ResidueType = "avg"
) -> tuple[onp.Array1D[np.complex128], onp.Array1D[np.complex128], onp.Array1D[np.float64]]: ...

#
def invres(
    r: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToFloat1D, tol: float = 0.001, rtype: _ResidueType = "avg"
) -> tuple[onp.Array1D[np.complex128], onp.Array1D[np.complex128]]: ...
def invresz(
    r: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToFloat1D, tol: float = 0.001, rtype: _ResidueType = "avg"
) -> tuple[onp.Array1D[np.complex128], onp.Array1D[np.complex128]]: ...

# NOTE: We use `_AnyInexact64T` as "free" type parameter, which behaves exactly as
# the (hypothetical) `AnyOf[np.float64, np.complex128]` gradual type.
@overload  # known dtype, known shape, t=None (default)
def resample(
    x: nptc.CanArray[_ShapeT, np.dtype[_InexactT3]],
    num: int,
    t: None = None,
    axis: int = 0,
    window: _ToResampleWindow[_AnyInexact64T] | None = None,
    domain: _Domain = "time",
) -> onp.ArrayND[_InexactT3, _ShapeT]: ...
@overload  # known dtype, known shape, t=<given>
def resample(
    x: nptc.CanArray[_ShapeT, np.dtype[_InexactT3]],
    num: int,
    t: onp.ToFloat1D,
    axis: int = 0,
    window: _ToResampleWindow[_AnyInexact64T] | None = None,
    domain: _Domain = "time",
) -> tuple[onp.ArrayND[_InexactT3, _ShapeT], onp.Array1D[np.float64]]: ...
@overload  # +integer, known shape, t=None (default)
def resample(
    x: nptc.CanArray[_ShapeT, np.dtype[npc.integer | np.bool_]],
    num: int,
    t: None = None,
    axis: int = 0,
    window: _ToResampleWindow[np.float64] | None = None,
    domain: _Domain = "time",
) -> onp.ArrayND[np.float64, _ShapeT]: ...
@overload  # +integer, known shape, t=<given>
def resample(
    x: nptc.CanArray[_ShapeT, np.dtype[npc.integer | np.bool_]],
    num: int,
    t: onp.ToFloat1D,
    axis: int = 0,
    window: _ToResampleWindow[np.float64] | None = None,
    domain: _Domain = "time",
) -> tuple[onp.ArrayND[np.float64, _ShapeT], onp.Array1D[np.float64]]: ...
@overload  # ~float16, known shape, t=None (default)
def resample(
    x: nptc.CanArray[_ShapeT, np.dtype[np.float16]],
    num: int,
    t: None = None,
    axis: int = 0,
    window: _ToResampleWindow[np.float64] | None = None,
    domain: _Domain = "time",
) -> onp.ArrayND[np.float32, _ShapeT]: ...
@overload  # ~float16, known shape, t=<given>
def resample(
    x: nptc.CanArray[_ShapeT, np.dtype[np.float16]],
    num: int,
    t: onp.ToFloat1D,
    axis: int = 0,
    window: _ToResampleWindow[np.float64] | None = None,
    domain: _Domain = "time",
) -> tuple[onp.ArrayND[np.float32, _ShapeT], onp.Array1D[np.float64]]: ...
@overload  # +float, unknown shape, t=None (default)
def resample(
    x: onp.SequenceND[float],
    num: int,
    t: None = None,
    axis: int = 0,
    window: _ToResampleWindow[np.float64] | None = None,
    domain: _Domain = "time",
) -> onp.ArrayND[np.float64]: ...
@overload  # +float, unknown shape, t=<given>
def resample(
    x: onp.SequenceND[float],
    num: int,
    t: onp.ToFloat1D,
    axis: int = 0,
    window: _ToResampleWindow[np.float64] | None = None,
    domain: _Domain = "time",
) -> tuple[onp.ArrayND[np.float64], onp.Array1D[np.float64]]: ...
@overload  # ~complex, unknown shape, t=None (default)
def resample(
    x: onp.SequenceND[op.JustComplex | np.complex128],
    num: int,
    t: None = None,
    axis: int = 0,
    window: _ToResampleWindow[np.complex128] | None = None,
    domain: _Domain = "time",
) -> onp.ArrayND[np.complex128]: ...
@overload  # ~complex, unknown shape, t=<given>
def resample(
    x: onp.SequenceND[op.JustComplex | np.complex128],
    num: int,
    t: onp.ToFloat1D,
    axis: int = 0,
    window: _ToResampleWindow[np.complex128] | None = None,
    domain: _Domain = "time",
) -> tuple[onp.ArrayND[np.complex128], onp.Array1D[np.float64]]: ...
@overload  # unknown dtype, unknown shape, t=None (default)
def resample(
    x: onp.ToComplexND,
    num: int,
    t: None = None,
    axis: int = 0,
    window: _ToResampleWindow[_AnyInexact64T] | None = None,
    domain: _Domain = "time",
) -> onp.ArrayND[Any, _WorkaroundForPyright]: ...
@overload  # unknown dtype, unknown shape, t=<given>
def resample(
    x: onp.ToComplexND,
    num: int,
    t: onp.ToFloat1D,
    axis: int = 0,
    window: _ToResampleWindow[_AnyInexact64T] | None = None,
    domain: _Domain = "time",
) -> tuple[onp.ArrayND[Any, _WorkaroundForPyright], onp.Array1D[np.float64]]: ...

# NOTE: This does not support the (useless) `up == down` case, which at runtime can return ANY dtype.
@overload  # known dtype, known shape
def resample_poly(
    x: nptc.CanArray[_ShapeT, np.dtype[_InexactT2]],
    up: int,
    down: int,
    axis: int = 0,
    window: _ToWindow = ("kaiser", 5.0),
    padtype: _PadType = "constant",
    cval: float | None = None,
) -> onp.ArrayND[_InexactT2, _ShapeT]: ...
@overload  # +integer, known shape
def resample_poly(
    x: nptc.CanArray[_ShapeT, np.dtype[npc.integer | np.bool_]],
    up: int,
    down: int,
    axis: int = 0,
    window: _ToWindow = ("kaiser", 5.0),
    padtype: _PadType = "constant",
    cval: float | None = None,
) -> onp.ArrayND[np.float64, _ShapeT]: ...
@overload  # ~float16, known shape
def resample_poly(
    x: nptc.CanArray[_ShapeT, np.dtype[np.float16]],
    up: int,
    down: int,
    axis: int = 0,
    window: _ToWindow = ("kaiser", 5.0),
    padtype: _PadType = "constant",
    cval: float | None = None,
) -> onp.ArrayND[np.float32, _ShapeT]: ...
@overload  # +float, 1d
def resample_poly(
    x: Sequence[float],
    up: int,
    down: int,
    axis: int = 0,
    window: _ToWindow = ("kaiser", 5.0),
    padtype: _PadType = "constant",
    cval: float | None = None,
) -> onp.Array1D[np.float64]: ...
@overload  # +float, unknown shape
def resample_poly(
    x: onp.SequenceND[float],
    up: int,
    down: int,
    axis: int = 0,
    window: _ToWindow = ("kaiser", 5.0),
    padtype: _PadType = "constant",
    cval: float | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # ~complex, 1d
def resample_poly(
    x: Sequence[op.JustComplex | np.complex128],
    up: int,
    down: int,
    axis: int = 0,
    window: _ToWindow = ("kaiser", 5.0),
    padtype: _PadType = "constant",
    cval: float | None = None,
) -> onp.Array1D[np.complex128]: ...
@overload  # ~complex, unknown shape
def resample_poly(
    x: onp.SequenceND[op.JustComplex | np.complex128],
    up: int,
    down: int,
    axis: int = 0,
    window: _ToWindow = ("kaiser", 5.0),
    padtype: _PadType = "constant",
    cval: float | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload  # unknown dtype, unknown shape
def resample_poly(
    x: onp.ToComplex128_ND,
    up: int,
    down: int,
    axis: int = 0,
    window: _ToWindow = ("kaiser", 5.0),
    padtype: _PadType = "constant",
    cval: float | None = None,
) -> onp.ArrayND[Any, _WorkaroundForPyright]: ...

# TODO(jorenham): improve
@overload
def decimate(
    x: onp.ToFloatND, q: int, n: int | None = None, ftype: _FilterType = "iir", axis: int = -1, zero_phase: bool = True
) -> onp.ArrayND[npc.floating]: ...
@overload
def decimate(
    x: onp.ToComplexND, q: int, n: int | None = None, ftype: _FilterType = "iir", axis: int = -1, zero_phase: bool = True
) -> onp.ArrayND[npc.inexact]: ...

# TODO(jorenham): improve
@overload
def vectorstrength(events: onp.ToFloat1D, period: onp.ToFloat) -> _Tuple2[npc.floating]: ...
@overload
def vectorstrength(events: onp.ToFloat1D, period: onp.ToFloat1D) -> _Tuple2[onp.Array1D[npc.floating]]: ...
@overload
def vectorstrength(events: onp.ToComplex1D, period: onp.ToComplex) -> _Tuple2[npc.inexact]: ...
@overload
def vectorstrength(events: onp.ToComplex1D, period: onp.ToComplex1D) -> _Tuple2[onp.Array1D[npc.inexact]]: ...

#
@overload
def envelope(
    z: onp.Array1D[np.float16],
    bp_in: tuple[int | None, int | None] = (1, None),
    *,
    n_out: int | None = None,
    squared: bool = False,
    residual: _ResidualKind | None = "lowpass",
    axis: int = -1,
) -> onp.Array2D[np.float32]: ...
@overload
def envelope(
    z: onp.Array2D[np.float16],
    bp_in: tuple[int | None, int | None] = (1, None),
    *,
    n_out: int | None = None,
    squared: bool = False,
    residual: _ResidualKind | None = "lowpass",
    axis: int = -1,
) -> onp.Array3D[np.float32]: ...
@overload
def envelope(
    z: onp.ArrayND[np.float16],
    bp_in: tuple[int | None, int | None] = (1, None),
    *,
    n_out: int | None = None,
    squared: bool = False,
    residual: _ResidualKind | None = "lowpass",
    axis: int = -1,
) -> onp.ArrayND[np.float32]: ...
@overload
def envelope(
    z: onp.Array1D[_InexactT3],
    bp_in: tuple[int | None, int | None] = (1, None),
    *,
    n_out: int | None = None,
    squared: bool = False,
    residual: _ResidualKind | None = "lowpass",
    axis: int = -1,
) -> onp.Array2D[_InexactT3]: ...
@overload
def envelope(
    z: onp.Array2D[_InexactT3],
    bp_in: tuple[int | None, int | None] = (1, None),
    *,
    n_out: int | None = None,
    squared: bool = False,
    residual: _ResidualKind | None = "lowpass",
    axis: int = -1,
) -> onp.Array3D[_InexactT3]: ...
@overload
def envelope(
    z: onp.ArrayND[_InexactT3],
    bp_in: tuple[int | None, int | None] = (1, None),
    *,
    n_out: int | None = None,
    squared: bool = False,
    residual: _ResidualKind | None = "lowpass",
    axis: int = -1,
) -> onp.ArrayND[_InexactT3]: ...
