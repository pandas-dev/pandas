from collections.abc import Callable, Sequence
from typing import Any, Concatenate, Literal, TypeAlias, TypedDict, overload, type_check_only
from typing_extensions import TypeVar, Unpack

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._lib._ccallback import LowLevelCallable

__all__ = [
    "convolve",
    "convolve1d",
    "correlate",
    "correlate1d",
    "gaussian_filter",
    "gaussian_filter1d",
    "gaussian_gradient_magnitude",
    "gaussian_laplace",
    "generic_filter",
    "generic_filter1d",
    "generic_gradient_magnitude",
    "generic_laplace",
    "laplace",
    "maximum_filter",
    "maximum_filter1d",
    "median_filter",
    "minimum_filter",
    "minimum_filter1d",
    "percentile_filter",
    "prewitt",
    "rank_filter",
    "sobel",
    "uniform_filter",
    "uniform_filter1d",
    "vectorized_filter",
]

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])
_DTypeT = TypeVar("_DTypeT", bound=np.dtype[np.bool_ | npc.number])
_ScalarT = TypeVar("_ScalarT", bound=np.bool_ | npc.number, default=Any)

_Ignored: TypeAlias = object

_Mode: TypeAlias = Literal["reflect", "constant", "nearest", "mirror", "wrap", "grid-constant", "grid-mirror", "grid-wrap"]
_Modes: TypeAlias = _Mode | Sequence[_Mode]
_Ints: TypeAlias = int | Sequence[int]
_AnyOutput: TypeAlias = onp.ArrayND[npc.number | np.bool_] | onp.AnyDType

_FilterFunc1D: TypeAlias = Callable[Concatenate[onp.Array1D[np.float64], onp.Array1D[np.float64], ...], _Ignored]
_FilterFuncND: TypeAlias = Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToComplex | onp.ToComplexND]
_Derivative: TypeAlias = Callable[
    # (input, axis, output, mode, cval, *extra_arguments, **extra_keywords)
    Concatenate[np.ndarray[Any, _DTypeT], int, onp.Array[Any, _ScalarT] | np.dtype[_ScalarT], _Mode, onp.ToComplex, ...],
    onp.ArrayND[Any],
]

@type_check_only
class _GaussianKwargs(TypedDict, total=False):
    truncate: float
    radius: _Ints

###

@overload
def vectorized_filter(
    input: onp.ToComplexND,
    function: Callable[..., _ScalarT],
    *,
    size: op.CanIndex | tuple[op.CanIndex, ...] | None = None,
    footprint: onp.Array | None = None,
    output: None = None,
    mode: Literal["reflect", "constant", "nearest", "mirror", "wrap"] = "reflect",
    cval: onp.ToFloat | None = None,
    origin: onp.ToInt | onp.ToInt1D | None = None,
    axes: tuple[op.CanIndex, ...] | None = None,
    batch_memory: int = 1_073_741_824,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def vectorized_filter(
    input: onp.ToComplexND,
    function: Callable[..., onp.ToComplex],
    *,
    size: op.CanIndex | tuple[op.CanIndex, ...] | None = None,
    footprint: onp.Array | None = None,
    output: onp.ArrayND[_ScalarT] | np.dtype[_ScalarT] | type[_ScalarT] | None = None,
    mode: Literal["reflect", "constant", "nearest", "mirror", "wrap"] = "reflect",
    cval: onp.ToFloat | None = None,
    origin: onp.ToInt | onp.ToInt1D | None = None,
    axes: tuple[op.CanIndex, ...] | None = None,
    batch_memory: int = 1_073_741_824,
) -> onp.ArrayND[_ScalarT]: ...

# keep roughly in sync with sobel
@overload
def correlate1d(
    input: onp.CanArray[_ShapeT, _DTypeT],
    weights: onp.ToFloat1D,
    axis: int = -1,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def correlate1d(
    input: onp.ToArrayND[op.JustInt, np.intp],
    weights: onp.ToFloat1D,
    axis: int = -1,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.intp]: ...
@overload
def correlate1d(
    input: onp.ToJustFloat64_ND,
    weights: onp.ToFloat1D,
    axis: int = -1,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.float64]: ...
@overload
def correlate1d(
    input: onp.ToJustComplex128_ND,
    weights: onp.ToFloat1D,
    axis: int = -1,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.complex128]: ...
@overload
def correlate1d(
    input: onp.ToComplexND,
    weights: onp.ToFloat1D,
    axis: int = -1,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def correlate1d(
    input: onp.ToComplexND,
    weights: onp.ToFloat1D,
    axis: int = -1,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> onp.ArrayND[Any]: ...

# keep in sync with correlate1d
@overload
def correlate(
    input: onp.CanArray[_ShapeT, _DTypeT],
    weights: onp.ToFloatND,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
    *,
    axes: _ShapeT | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def correlate(
    input: onp.ToArrayND[op.JustInt, np.intp],
    weights: onp.ToFloatND,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: int = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def correlate(
    input: onp.ToJustFloat64_ND,
    weights: onp.ToFloatND,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: int = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def correlate(
    input: onp.ToJustComplex128_ND,
    weights: onp.ToFloatND,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    origin: int = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def correlate(
    input: onp.ToComplexND,
    weights: onp.ToFloatND,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def correlate(
    input: onp.ToComplexND,
    weights: onp.ToFloatND,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[Any]: ...

# keep in sync with correlate1d
@overload
def convolve1d(
    input: onp.CanArray[_ShapeT, _DTypeT],
    weights: onp.ToFloat1D,
    axis: int = -1,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def convolve1d(
    input: onp.ToArrayND[op.JustInt, np.intp],
    weights: onp.ToFloat1D,
    axis: int = -1,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.intp]: ...
@overload
def convolve1d(
    input: onp.ToJustFloat64_ND,
    weights: onp.ToFloat1D,
    axis: int = -1,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.float64]: ...
@overload
def convolve1d(
    input: onp.ToJustComplex128_ND,
    weights: onp.ToFloat1D,
    axis: int = -1,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.complex128]: ...
@overload
def convolve1d(
    input: onp.ToComplexND,
    weights: onp.ToFloat1D,
    axis: int = -1,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def convolve1d(
    input: onp.ToComplexND,
    weights: onp.ToFloat1D,
    axis: int = -1,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> onp.ArrayND[Any]: ...

# keep in sync with correlate
@overload
def convolve(
    input: onp.CanArray[_ShapeT, _DTypeT],
    weights: onp.ToFloatND,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
    *,
    axes: _ShapeT | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def convolve(
    input: onp.ToArrayND[op.JustInt, np.intp],
    weights: onp.ToFloatND,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: int = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def convolve(
    input: onp.ToJustFloat64_ND,
    weights: onp.ToFloatND,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: int = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def convolve(
    input: onp.ToJustComplex128_ND,
    weights: onp.ToFloatND,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    origin: int = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def convolve(
    input: onp.ToComplexND,
    weights: onp.ToFloatND,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def convolve(
    input: onp.ToComplexND,
    weights: onp.ToFloatND,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[Any]: ...

# keep in sync with sobel
@overload
def prewitt(
    input: onp.CanArray[_ShapeT, _DTypeT],
    axis: int = -1,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def prewitt(
    input: onp.ToArrayND[op.JustInt, np.intp],
    axis: int = -1,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
) -> onp.ArrayND[np.intp]: ...
@overload
def prewitt(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    axis: int = -1,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
) -> onp.ArrayND[np.float64]: ...
@overload
def prewitt(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    axis: int = -1,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
) -> onp.ArrayND[np.complex128]: ...
@overload
def prewitt(
    input: onp.ToComplex | onp.ToComplexND,
    axis: int = -1,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def prewitt(
    input: onp.ToComplex | onp.ToComplexND,
    axis: int = -1,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
) -> onp.ArrayND[Any]: ...

# keep in sync with prewitt
@overload
def sobel(
    input: onp.CanArray[_ShapeT, _DTypeT],
    axis: int = -1,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def sobel(
    input: onp.ToArrayND[op.JustInt, np.intp],
    axis: int = -1,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
) -> onp.ArrayND[np.intp]: ...
@overload
def sobel(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    axis: int = -1,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
) -> onp.ArrayND[np.float64]: ...
@overload
def sobel(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    axis: int = -1,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
) -> onp.ArrayND[np.complex128]: ...
@overload
def sobel(
    input: onp.ToComplex | onp.ToComplexND,
    axis: int = -1,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def sobel(
    input: onp.ToComplex | onp.ToComplexND,
    axis: int = -1,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
) -> onp.ArrayND[Any]: ...

#
@overload
def laplace(
    input: onp.CanArray[_ShapeT, _DTypeT],
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    *,
    axes: _ShapeT | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def laplace(
    input: onp.ToArrayND[op.JustInt, np.intp],
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def laplace(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def laplace(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def laplace(
    input: onp.ToComplex | onp.ToComplexND,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def laplace(
    input: onp.ToComplex | onp.ToComplexND,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[Any]: ...

# keep in sync with laplace
@overload
def gaussian_laplace(
    input: onp.CanArray[_ShapeT, _DTypeT],
    sigma: onp.ToFloat | onp.ToFloatND,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    *,
    axes: _ShapeT | None = None,
    **kwargs: Unpack[_GaussianKwargs],
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def gaussian_laplace(
    input: onp.ToArrayND[op.JustInt, np.intp],
    sigma: onp.ToFloat | onp.ToFloatND,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
    **kwargs: Unpack[_GaussianKwargs],
) -> onp.ArrayND[np.intp]: ...
@overload
def gaussian_laplace(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    sigma: onp.ToFloat | onp.ToFloatND,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
    **kwargs: Unpack[_GaussianKwargs],
) -> onp.ArrayND[np.float64]: ...
@overload
def gaussian_laplace(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    sigma: onp.ToFloat | onp.ToFloatND,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
    **kwargs: Unpack[_GaussianKwargs],
) -> onp.ArrayND[np.complex128]: ...
@overload
def gaussian_laplace(
    input: onp.ToComplex | onp.ToComplexND,
    sigma: onp.ToFloat | onp.ToFloatND,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
    **kwargs: Unpack[_GaussianKwargs],
) -> onp.ArrayND[_ScalarT]: ...
@overload
def gaussian_laplace(
    input: onp.ToComplex | onp.ToComplexND,
    sigma: onp.ToFloat | onp.ToFloatND,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
    **kwargs: Unpack[_GaussianKwargs],
) -> onp.ArrayND[Any]: ...

# keep in sync with laplace
@overload
def generic_laplace(
    input: onp.CanArray[_ShapeT, _DTypeT],
    derivative2: _Derivative[_DTypeT, Any],
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: _ShapeT | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def generic_laplace(
    input: onp.ToArrayND[op.JustInt, np.intp],
    derivative2: _Derivative[np.dtype[np.intp], np.intp],
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def generic_laplace(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    derivative2: _Derivative[np.dtype[np.float64], np.float64],
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def generic_laplace(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    derivative2: _Derivative[np.dtype[np.complex128], np.complex128],
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def generic_laplace(
    input: onp.ToComplex | onp.ToComplexND,
    derivative2: _Derivative[Any, _ScalarT],
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def generic_laplace(
    input: onp.ToComplex | onp.ToComplexND,
    derivative2: _Derivative[Any, Any],
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[Any]: ...

# keep in sync with gaussian_laplace
@overload
def gaussian_gradient_magnitude(
    input: onp.CanArray[_ShapeT, _DTypeT],
    sigma: onp.ToFloat | onp.ToFloatND,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    *,
    axes: _ShapeT | None = None,
    **kwargs: Unpack[_GaussianKwargs],
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def gaussian_gradient_magnitude(
    input: onp.ToArrayND[op.JustInt, np.intp],
    sigma: onp.ToFloat | onp.ToFloatND,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
    **kwargs: Unpack[_GaussianKwargs],
) -> onp.ArrayND[np.intp]: ...
@overload
def gaussian_gradient_magnitude(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    sigma: onp.ToFloat | onp.ToFloatND,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
    **kwargs: Unpack[_GaussianKwargs],
) -> onp.ArrayND[np.float64]: ...
@overload
def gaussian_gradient_magnitude(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    sigma: onp.ToFloat | onp.ToFloatND,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
    **kwargs: Unpack[_GaussianKwargs],
) -> onp.ArrayND[np.complex128]: ...
@overload
def gaussian_gradient_magnitude(
    input: onp.ToComplex | onp.ToComplexND,
    sigma: onp.ToFloat | onp.ToFloatND,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
    **kwargs: Unpack[_GaussianKwargs],
) -> onp.ArrayND[_ScalarT]: ...
@overload
def gaussian_gradient_magnitude(
    input: onp.ToComplex | onp.ToComplexND,
    sigma: onp.ToFloat | onp.ToFloatND,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    *,
    axes: tuple[int, ...] | None = None,
    **kwargs: Unpack[_GaussianKwargs],
) -> onp.ArrayND[Any]: ...

#
@overload
def generic_gradient_magnitude(
    input: onp.CanArray[_ShapeT, _DTypeT],
    derivative: _Derivative[_DTypeT, Any],
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: _ShapeT | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def generic_gradient_magnitude(
    input: onp.ToArrayND[op.JustInt, np.intp],
    derivative: _Derivative[np.dtype[np.intp], np.intp],
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def generic_gradient_magnitude(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    derivative: _Derivative[np.dtype[np.float64], np.float64],
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def generic_gradient_magnitude(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    derivative: _Derivative[np.dtype[np.complex128], np.complex128],
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def generic_gradient_magnitude(
    input: onp.ToComplex | onp.ToComplexND,
    derivative: _Derivative[Any, _ScalarT],
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def generic_gradient_magnitude(
    input: onp.ToComplex | onp.ToComplexND,
    derivative: _Derivative[Any, Any],
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[Any]: ...

#
@overload
def gaussian_filter1d(
    input: onp.CanArray[_ShapeT, _DTypeT],
    sigma: onp.ToFloat,
    axis: int = -1,
    order: int = 0,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToComplex = 0.0,
    truncate: float = 4.0,
    *,
    radius: int | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def gaussian_filter1d(
    input: onp.ToArrayND[op.JustInt, np.intp],
    sigma: onp.ToFloat,
    axis: int = -1,
    order: int = 0,
    output: onp.AnyIntPDType | None = None,
    mode: _Mode = "reflect",
    cval: float | onp.ToInt = 0.0,
    truncate: float = 4.0,
    *,
    radius: int | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def gaussian_filter1d(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    sigma: onp.ToFloat,
    axis: int = -1,
    order: int = 0,
    output: onp.AnyFloat64DType | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToFloat64 = 0.0,
    truncate: float = 4.0,
    *,
    radius: int | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def gaussian_filter1d(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    sigma: onp.ToFloat,
    axis: int = -1,
    order: int = 0,
    output: onp.AnyComplex128DType | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToComplex128 = 0.0,
    truncate: float = 4.0,
    *,
    radius: int | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def gaussian_filter1d(
    input: onp.ToComplex | onp.ToComplexND,
    sigma: onp.ToFloat,
    axis: int = -1,
    order: int = 0,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToComplex = 0.0,
    truncate: float = 4.0,
    *,
    radius: int | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def gaussian_filter1d(
    input: onp.ToComplex | onp.ToComplexND,
    sigma: onp.ToFloat,
    axis: int = -1,
    order: int = 0,
    output: _AnyOutput | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToComplex = 0.0,
    truncate: float = 4.0,
    *,
    radius: int | None = None,
) -> onp.ArrayND[Any]: ...

# keep in sync with generic_laplace
@overload
def gaussian_filter(
    input: onp.CanArray[_ShapeT, _DTypeT],
    sigma: onp.ToFloat | onp.ToFloatND,
    order: _Ints = 0,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    truncate: float = 4.0,
    *,
    radius: _Ints | None = None,
    axes: tuple[int, ...] | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def gaussian_filter(
    input: onp.ToArrayND[op.JustInt, np.intp],
    sigma: onp.ToFloat | onp.ToFloatND,
    order: _Ints = 0,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    truncate: float = 4.0,
    *,
    radius: _Ints | None = None,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def gaussian_filter(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    sigma: onp.ToFloat | onp.ToFloatND,
    order: _Ints = 0,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    truncate: float = 4.0,
    *,
    radius: _Ints | None = None,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def gaussian_filter(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    sigma: onp.ToFloat | onp.ToFloatND,
    order: _Ints = 0,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    truncate: float = 4.0,
    *,
    radius: _Ints | None = None,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def gaussian_filter(
    input: onp.ToComplex | onp.ToComplexND,
    sigma: onp.ToFloat | onp.ToFloatND,
    order: _Ints = 0,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    truncate: float = 4.0,
    *,
    radius: _Ints | None = None,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def gaussian_filter(
    input: onp.ToComplex | onp.ToComplexND,
    sigma: onp.ToFloat | onp.ToFloatND,
    order: _Ints = 0,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    truncate: float = 4.0,
    *,
    radius: _Ints | None = None,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[Any]: ...

#
@overload
def generic_filter1d(
    input: onp.CanArray[_ShapeT, _DTypeT],
    function: _FilterFunc1D | LowLevelCallable,
    filter_size: float,
    axis: int = -1,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToFloat = 0.0,
    origin: int = 0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def generic_filter1d(
    input: onp.ToArrayND[op.JustInt, np.intp],
    function: _FilterFunc1D | LowLevelCallable,
    filter_size: float,
    axis: int = -1,
    output: onp.AnyIntPDType | None = None,
    mode: _Mode = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: int = 0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def generic_filter1d(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    function: _FilterFunc1D | LowLevelCallable,
    filter_size: float,
    axis: int = -1,
    output: onp.AnyFloat64DType | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: int = 0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def generic_filter1d(
    input: onp.ToFloat | onp.ToFloatND,
    function: _FilterFunc1D | LowLevelCallable,
    filter_size: float,
    axis: int = -1,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToFloat = 0.0,
    origin: int = 0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def generic_filter1d(
    input: onp.ToFloat | onp.ToFloatND,
    function: _FilterFunc1D | LowLevelCallable,
    filter_size: float,
    axis: int = -1,
    output: _AnyOutput | None = None,
    mode: _Mode = "reflect",
    cval: onp.ToFloat = 0.0,
    origin: int = 0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
) -> onp.ArrayND[Any]: ...

# keep in sync with generic_filter1d
@overload
def generic_filter(
    input: onp.CanArray[_ShapeT, _DTypeT],
    function: _FilterFuncND | LowLevelCallable,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat = 0.0,
    origin: _Ints = 0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def generic_filter(
    input: onp.ToArrayND[op.JustInt, np.intp],
    function: _FilterFuncND | LowLevelCallable,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: _Ints = 0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def generic_filter(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    function: _FilterFuncND | LowLevelCallable,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: _Ints = 0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def generic_filter(
    input: onp.ToFloat | onp.ToFloatND,
    function: _FilterFuncND | LowLevelCallable,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat = 0.0,
    origin: _Ints = 0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def generic_filter(
    input: onp.ToFloat | onp.ToFloatND,
    function: _FilterFuncND | LowLevelCallable,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat = 0.0,
    origin: _Ints = 0,
    extra_arguments: tuple[object, ...] = (),
    extra_keywords: dict[str, Any] | None = None,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[Any]: ...

#
@overload
def uniform_filter1d(
    input: onp.CanArray[_ShapeT, _DTypeT],
    size: int,
    axis: int = -1,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def uniform_filter1d(
    input: onp.ToArrayND[op.JustInt, np.intp],
    size: int,
    axis: int = -1,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.intp]: ...
@overload
def uniform_filter1d(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    size: int,
    axis: int = -1,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.float64]: ...
@overload
def uniform_filter1d(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    size: int,
    axis: int = -1,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.complex128]: ...
@overload
def uniform_filter1d(
    input: onp.ToComplex | onp.ToComplexND,
    size: int,
    axis: int = -1,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def uniform_filter1d(
    input: onp.ToComplex | onp.ToComplexND,
    size: int,
    axis: int = -1,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> onp.ArrayND[Any]: ...

# keep roughly in sync with uniform_filter1d
@overload
def uniform_filter(
    input: onp.CanArray[_ShapeT, _DTypeT],
    size: _Ints = 3,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def uniform_filter(
    input: onp.ToArrayND[op.JustInt, np.intp],
    size: _Ints = 3,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def uniform_filter(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    size: _Ints = 3,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def uniform_filter(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    size: _Ints = 3,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def uniform_filter(
    input: onp.ToComplex | onp.ToComplexND,
    size: _Ints = 3,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def uniform_filter(
    input: onp.ToComplex | onp.ToComplexND,
    size: _Ints = 3,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[Any]: ...

# keep in sync with uniform_filter1d
@overload
def minimum_filter1d(
    input: onp.CanArray[_ShapeT, _DTypeT],
    size: int,
    axis: int = -1,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def minimum_filter1d(
    input: onp.ToArrayND[op.JustInt, np.intp],
    size: int,
    axis: int = -1,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.intp]: ...
@overload
def minimum_filter1d(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    size: int,
    axis: int = -1,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.float64]: ...
@overload
def minimum_filter1d(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    size: int,
    axis: int = -1,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.complex128]: ...
@overload
def minimum_filter1d(
    input: onp.ToComplex | onp.ToComplexND,
    size: int,
    axis: int = -1,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def minimum_filter1d(
    input: onp.ToComplex | onp.ToComplexND,
    size: int,
    axis: int = -1,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> onp.ArrayND[Any]: ...

# keep in sync with uniform_filter (plus `footprint`)
@overload
def minimum_filter(
    input: onp.CanArray[_ShapeT, _DTypeT],
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def minimum_filter(
    input: onp.ToArrayND[op.JustInt, np.intp],
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def minimum_filter(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def minimum_filter(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def minimum_filter(
    input: onp.ToComplex | onp.ToComplexND,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def minimum_filter(
    input: onp.ToComplex | onp.ToComplexND,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[Any]: ...

# keep in sync with uniform_filter1d
@overload
def maximum_filter1d(
    input: onp.CanArray[_ShapeT, _DTypeT],
    size: int,
    axis: int = -1,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def maximum_filter1d(
    input: onp.ToArrayND[op.JustInt, np.intp],
    size: int,
    axis: int = -1,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.intp]: ...
@overload
def maximum_filter1d(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    size: int,
    axis: int = -1,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.float64]: ...
@overload
def maximum_filter1d(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    size: int,
    axis: int = -1,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    origin: int = 0,
) -> onp.ArrayND[np.complex128]: ...
@overload
def maximum_filter1d(
    input: onp.ToComplex | onp.ToComplexND,
    size: int,
    axis: int = -1,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def maximum_filter1d(
    input: onp.ToComplex | onp.ToComplexND,
    size: int,
    axis: int = -1,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: int = 0,
) -> onp.ArrayND[Any]: ...

# keep in sync with minimum_filter
@overload
def maximum_filter(
    input: onp.CanArray[_ShapeT, _DTypeT],
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def maximum_filter(
    input: onp.ToArrayND[op.JustInt, np.intp],
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def maximum_filter(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def maximum_filter(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def maximum_filter(
    input: onp.ToComplex | onp.ToComplexND,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def maximum_filter(
    input: onp.ToComplex | onp.ToComplexND,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[Any]: ...

# keep in sync with maximum_filter
@overload
def median_filter(
    input: onp.CanArray[_ShapeT, _DTypeT],
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def median_filter(
    input: onp.ToArrayND[op.JustInt, np.intp],
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def median_filter(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def median_filter(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def median_filter(
    input: onp.ToComplex | onp.ToComplexND,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def median_filter(
    input: onp.ToComplex | onp.ToComplexND,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[Any]: ...

# keep in sync with median_filter (plus `rank`)
@overload
def rank_filter(
    input: onp.CanArray[_ShapeT, _DTypeT],
    rank: int,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def rank_filter(
    input: onp.ToArrayND[op.JustInt, np.intp],
    rank: int,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def rank_filter(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    rank: int,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def rank_filter(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    rank: int,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def rank_filter(
    input: onp.ToComplex | onp.ToComplexND,
    rank: int,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def rank_filter(
    input: onp.ToComplex | onp.ToComplexND,
    rank: int,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[Any]: ...

# keep in sync with median_filter (plus `percentile` and minus `rank`)
@overload
def percentile_filter(
    input: onp.CanArray[_ShapeT, _DTypeT],
    percentile: onp.ToFloat,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: np.ndarray[_ShapeT, _DTypeT] | _DTypeT | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def percentile_filter(
    input: onp.ToArrayND[op.JustInt, np.intp],
    percentile: onp.ToFloat,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyIntPDType | None = None,
    mode: _Modes = "reflect",
    cval: float | onp.ToInt = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.intp]: ...
@overload
def percentile_filter(
    input: onp.ToJustFloat64 | onp.ToJustFloat64_ND,
    percentile: onp.ToFloat,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyFloat64DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToFloat64 = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def percentile_filter(
    input: onp.ToJustComplex128 | onp.ToJustComplex128_ND,
    percentile: onp.ToFloat,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.AnyComplex128DType | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex128 = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def percentile_filter(
    input: onp.ToComplex | onp.ToComplexND,
    percentile: onp.ToFloat,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: onp.ArrayND[_ScalarT] | onp.ToDType[_ScalarT] | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[_ScalarT]: ...
@overload
def percentile_filter(
    input: onp.ToComplex | onp.ToComplexND,
    percentile: onp.ToFloat,
    size: _Ints | None = None,
    footprint: onp.ToInt | onp.ToIntND | None = None,
    output: _AnyOutput | None = None,
    mode: _Modes = "reflect",
    cval: onp.ToComplex = 0.0,
    origin: _Ints = 0,
    *,
    axes: tuple[int, ...] | None = None,
) -> onp.ArrayND[Any]: ...
