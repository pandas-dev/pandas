from typing import Any, SupportsFloat as CanFloat, SupportsIndex as CanIndex, TypeAlias, TypeVar, overload

import numpy as np
import optype.numpy as onp

__all__ = ["fourier_ellipsoid", "fourier_gaussian", "fourier_shift", "fourier_uniform"]

###

_ScalarComplex: TypeAlias = np.complex64 | np.complex128
_OutputScalarComplexT = TypeVar("_OutputScalarComplexT", bound=_ScalarComplex)
_OutputArrayComplexT = TypeVar("_OutputArrayComplexT", bound=onp.ArrayND[_ScalarComplex])

_Scalar: TypeAlias = np.float32 | np.float64 | _ScalarComplex
_OutputScalarT = TypeVar("_OutputScalarT", bound=_Scalar)
_OutputArrayT = TypeVar("_OutputArrayT", bound=onp.ArrayND[_Scalar])

_Sigma: TypeAlias = CanFloat | onp.ToFloat1D
_AsF32: TypeAlias = np.float16 | np.float32
_InputF64: TypeAlias = onp.ToJustFloat64_ND | onp.ToIntND
_InputC128: TypeAlias = onp.ToJustComplex128_ND  # pro forma
# these *should* be equivalent to `onp.ArrayND[Never, {}]`, but both mypy and pyright currently have bugs with `Never` unions
_InputF32: TypeAlias = onp.CanArrayND[_AsF32] | onp.SequenceND[onp.CanArray[Any, np.dtype[_AsF32]]]
_InputC64: TypeAlias = onp.CanArrayND[np.complex64] | onp.SequenceND[onp.CanArray[Any, np.dtype[np.complex64]]]

###
# NOTE: The gaussian, uniform, and ellipsoid function signatures are equivalent (except for the 2nd *name*): Keep them in sync!
# NOTE: The [overload-overlap] mypy errors are false positives (probably a union/join thing).

# undocumented
@overload  # output: <T: ndarray>
def _get_output_fourier(output: _OutputArrayT, input: onp.ToComplex128_ND) -> _OutputArrayT: ...
@overload  # output: <T: scalar type>
def _get_output_fourier(output: type[_OutputScalarT], input: onp.ToComplex128_ND) -> onp.ArrayND[_OutputScalarT]: ...
@overload  # +float32
def _get_output_fourier(output: None, input: _InputF32) -> onp.ArrayND[np.float32]: ...
@overload  # +float64
def _get_output_fourier(output: None, input: _InputF64) -> onp.ArrayND[np.float64]: ...
@overload  # ~complex64
def _get_output_fourier(output: None, input: _InputC64) -> onp.ArrayND[np.complex64]: ...
@overload  # ~complex128
def _get_output_fourier(output: None, input: _InputC128) -> onp.ArrayND[np.complex128]: ...
@overload  # fallback
def _get_output_fourier(output: None, input: onp.ToComplex128_ND) -> onp.ArrayND[_Scalar]: ...

# undocumented
@overload  # output: complex64 | complex128 array or scalar-type
def _get_output_fourier_complex(
    output: onp.ArrayND[_OutputScalarComplexT] | type[_OutputScalarComplexT], input: onp.ToComplex128_ND
) -> onp.ArrayND[_OutputScalarComplexT]: ...
@overload  # ~complex64
def _get_output_fourier_complex(output: None, input: _InputC64) -> onp.ArrayND[np.complex64]: ...
@overload  # fallback
def _get_output_fourier_complex(output: None, input: onp.ToComplex128_ND) -> onp.ArrayND[_ScalarComplex]: ...

# NOTE: Keep in sync with `fourier_uniform` and `fourier_ellipsoid` (but note the different 2nd parameter names)
@overload  # output: <T: ndarray> (positional)
def fourier_gaussian(
    input: onp.ToComplex128_ND, sigma: _Sigma, n: CanIndex, axis: int, output: _OutputArrayT
) -> _OutputArrayT: ...
@overload  # output: <T: ndarray> (keyword)
def fourier_gaussian(
    input: onp.ToComplex128_ND, sigma: _Sigma, n: CanIndex = -1, axis: int = -1, *, output: _OutputArrayT
) -> _OutputArrayT: ...
@overload  # output: <T: scalar type> (positional)
def fourier_gaussian(
    input: onp.ToComplex128_ND, sigma: _Sigma, n: CanIndex, axis: int, output: type[_OutputScalarT]
) -> onp.ArrayND[_OutputScalarT]: ...
@overload  # output: <T: scalar type> (keyword)
def fourier_gaussian(
    input: onp.ToComplex128_ND, sigma: _Sigma, n: CanIndex = -1, axis: int = -1, *, output: type[_OutputScalarT]
) -> onp.ArrayND[_OutputScalarT]: ...
@overload  # +float32
def fourier_gaussian(
    input: _InputF32, sigma: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.float32]: ...
@overload  # +float64
def fourier_gaussian(
    input: _InputF64, sigma: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # ~complex64
def fourier_gaussian(
    input: _InputC64, sigma: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex64]: ...
@overload  # ~complex128
def fourier_gaussian(
    input: _InputC128, sigma: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # fallback
def fourier_gaussian(
    input: onp.ToComplex128_ND, sigma: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[_Scalar]: ...

# NOTE: Keep in sync with `fourier_ellipsoid` and `fourier_gaussian` (but note the different 2nd parameter name)
@overload  # output: <T: ndarray> (positional)
def fourier_uniform(input: onp.ToComplex128_ND, size: _Sigma, n: CanIndex, axis: int, output: _OutputArrayT) -> _OutputArrayT: ...
@overload  # output: <T: ndarray> (keyword)
def fourier_uniform(
    input: onp.ToComplex128_ND, size: _Sigma, n: CanIndex = -1, axis: int = -1, *, output: _OutputArrayT
) -> _OutputArrayT: ...
@overload  # output: <T: scalar type> (positional)
def fourier_uniform(
    input: onp.ToComplex128_ND, size: _Sigma, n: CanIndex, axis: int, output: type[_OutputScalarT]
) -> onp.ArrayND[_OutputScalarT]: ...
@overload  # output: <T: scalar type> (keyword)
def fourier_uniform(
    input: onp.ToComplex128_ND, size: _Sigma, n: CanIndex = -1, axis: int = -1, *, output: type[_OutputScalarT]
) -> onp.ArrayND[_OutputScalarT]: ...
@overload  # +float32
def fourier_uniform(
    input: _InputF32, size: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.float32]: ...
@overload  # +float64
def fourier_uniform(
    input: _InputF64, size: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # ~complex64
def fourier_uniform(
    input: _InputC64, size: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex64]: ...
@overload  # ~complex128
def fourier_uniform(
    input: _InputC128, size: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # fallback
def fourier_uniform(
    input: onp.ToComplex128_ND, size: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[_Scalar]: ...

# NOTE: Keep in sync with `fourier_uniform` and `fourier_gaussian` (but note the different 2nd parameter name)
@overload  # output: <T: ndarray> (positional)
def fourier_ellipsoid(
    input: onp.ToComplex128_ND, size: _Sigma, n: CanIndex, axis: int, output: _OutputArrayT
) -> _OutputArrayT: ...
@overload  # output: <T: ndarray> (keyword)
def fourier_ellipsoid(
    input: onp.ToComplex128_ND, size: _Sigma, n: CanIndex = -1, axis: int = -1, *, output: _OutputArrayT
) -> _OutputArrayT: ...
@overload  # output: <T: scalar type> (positional)
def fourier_ellipsoid(
    input: onp.ToComplex128_ND, size: _Sigma, n: CanIndex, axis: int, output: type[_OutputScalarT]
) -> onp.ArrayND[_OutputScalarT]: ...
@overload  # output: <T: scalar type> (keyword)
def fourier_ellipsoid(
    input: onp.ToComplex128_ND, size: _Sigma, n: CanIndex = -1, axis: int = -1, *, output: type[_OutputScalarT]
) -> onp.ArrayND[_OutputScalarT]: ...
@overload  # +float32
def fourier_ellipsoid(
    input: _InputF32, size: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.float32]: ...
@overload  # +float64
def fourier_ellipsoid(
    input: _InputF64, size: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # ~complex64
def fourier_ellipsoid(
    input: _InputC64, size: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex64]: ...
@overload  # ~complex128
def fourier_ellipsoid(
    input: _InputC128, size: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # fallback
def fourier_ellipsoid(
    input: onp.ToComplex128_ND, size: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[_Scalar]: ...

# NOTE: Unlike the other three functions, this always returns complex output
@overload  # output: <T: ndarray> (positional)
def fourier_shift(
    input: onp.ToComplex128_ND, shift: _Sigma, n: CanIndex, axis: int, output: _OutputArrayComplexT
) -> _OutputArrayComplexT: ...
@overload  # output: <T: ndarray> (keyword)
def fourier_shift(
    input: onp.ToComplex128_ND, shift: _Sigma, n: CanIndex = -1, axis: int = -1, *, output: _OutputArrayComplexT
) -> _OutputArrayComplexT: ...
@overload  # output: <T: scalar type> (positional)
def fourier_shift(
    input: onp.ToComplex128_ND, shift: _Sigma, n: CanIndex, axis: int, output: type[_OutputScalarComplexT]
) -> onp.ArrayND[_OutputScalarComplexT]: ...
@overload  # output: <T: scalar type> (keyword)
def fourier_shift(
    input: onp.ToComplex128_ND, shift: _Sigma, n: CanIndex = -1, axis: int = -1, *, output: type[_OutputScalarComplexT]
) -> onp.ArrayND[_OutputScalarComplexT]: ...
@overload  # ~complex64
def fourier_shift(
    input: _InputC64, shift: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex64]: ...
@overload  # ~complex128 | +floating
def fourier_shift(
    input: _InputC128 | onp.ToFloat64_ND, shift: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # fallback
def fourier_shift(
    input: onp.ToComplex128_ND, shift: _Sigma, n: CanIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[_ScalarComplex]: ...
