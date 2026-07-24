from typing import SupportsFloat, SupportsIndex, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["fourier_ellipsoid", "fourier_gaussian", "fourier_shift", "fourier_uniform"]

###

type _ScalarComplex = np.complex64 | np.complex128
type _Scalar = np.float32 | np.float64 | _ScalarComplex

type _Sigma = SupportsFloat | onp.ToFloat1D
type _InputF64 = onp.ToJustFloat64_ND | onp.ToIntND
type _InputF32 = onp.ToJustFloat32_ND
type _Input = onp.ToArrayND[complex, npc.number64 | npc.number32 | npc.integer | np.bool]

###

# NOTE: The gaussian, uniform, and ellipsoid function signatures are equivalent (except for the 2nd *name*): Keep them in sync!
# NOTE: The [overload-overlap] mypy errors are false positives (probably a union/join thing).

# undocumented
@overload  # output: <T: ndarray>
def _get_output_fourier[OutputArrayT: onp.ArrayND[_Scalar]](output: OutputArrayT, input: onp.ToComplex128_ND) -> OutputArrayT: ...
@overload  # output: <T: scalar type>
def _get_output_fourier[OutputScalarT: _Scalar](
    output: type[OutputScalarT], input: onp.ToComplex128_ND
) -> onp.ArrayND[OutputScalarT]: ...
@overload  # +float64
def _get_output_fourier(output: None, input: _InputF64) -> onp.ArrayND[np.float64]: ...
@overload  # +float32
def _get_output_fourier(output: None, input: _InputF32) -> onp.ArrayND[np.float32]: ...
@overload  # ~complex128
def _get_output_fourier(output: None, input: onp.ToJustComplex128_ND) -> onp.ArrayND[np.complex128]: ...
@overload  # ~complex64
def _get_output_fourier(output: None, input: onp.ToJustComplex64_ND) -> onp.ArrayND[np.complex64]: ...
@overload  # fallback
def _get_output_fourier(output: None, input: onp.ToComplex128_ND) -> onp.ArrayND[_Scalar]: ...

# undocumented
@overload  # output: complex64 | complex128 array or scalar-type
def _get_output_fourier_complex[OutputScalarComplexT: _ScalarComplex](
    output: onp.ArrayND[OutputScalarComplexT] | type[OutputScalarComplexT], input: _Input
) -> onp.ArrayND[OutputScalarComplexT]: ...
@overload  # ~complex64
def _get_output_fourier_complex(output: None, input: onp.ToJustComplex64_ND) -> onp.ArrayND[np.complex64]: ...
@overload  # fallback
def _get_output_fourier_complex(output: None, input: _Input) -> onp.ArrayND[_ScalarComplex]: ...

# NOTE: Keep in sync with `fourier_uniform` and `fourier_ellipsoid` (but note the different 2nd parameter names)
@overload  # output: <T: ndarray> (positional)
def fourier_gaussian[OutputArrayT: onp.ArrayND[_Scalar]](
    input: _Input, sigma: _Sigma, n: SupportsIndex, axis: int, output: OutputArrayT
) -> OutputArrayT: ...
@overload  # output: <T: ndarray> (keyword)
def fourier_gaussian[OutputArrayT: onp.ArrayND[_Scalar]](
    input: _Input, sigma: _Sigma, n: SupportsIndex = -1, axis: int = -1, *, output: OutputArrayT
) -> OutputArrayT: ...
@overload  # output: <T: scalar type> (positional)
def fourier_gaussian[OutputScalarT: _Scalar](
    input: _Input, sigma: _Sigma, n: SupportsIndex, axis: int, output: type[OutputScalarT]
) -> onp.ArrayND[OutputScalarT]: ...
@overload  # output: <T: scalar type> (keyword)
def fourier_gaussian[OutputScalarT: _Scalar](
    input: _Input, sigma: _Sigma, n: SupportsIndex = -1, axis: int = -1, *, output: type[OutputScalarT]
) -> onp.ArrayND[OutputScalarT]: ...
@overload  # +float64
def fourier_gaussian(
    input: _InputF64, sigma: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # +float32
def fourier_gaussian(
    input: _InputF32, sigma: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.float32]: ...
@overload  # ~complex128
def fourier_gaussian(
    input: onp.ToJustComplex128_ND, sigma: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # ~complex64
def fourier_gaussian(
    input: onp.ToJustComplex64_ND, sigma: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex64]: ...
@overload  # fallback
def fourier_gaussian(
    input: _Input, sigma: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[_Scalar]: ...

# NOTE: Keep in sync with `fourier_ellipsoid` and `fourier_gaussian` (but note the different 2nd parameter name)
@overload  # output: <T: ndarray> (positional)
def fourier_uniform[OutputArrayT: onp.ArrayND[_Scalar]](
    input: _Input, size: _Sigma, n: SupportsIndex, axis: int, output: OutputArrayT
) -> OutputArrayT: ...
@overload  # output: <T: ndarray> (keyword)
def fourier_uniform[OutputArrayT: onp.ArrayND[_Scalar]](
    input: _Input, size: _Sigma, n: SupportsIndex = -1, axis: int = -1, *, output: OutputArrayT
) -> OutputArrayT: ...
@overload  # output: <T: scalar type> (positional)
def fourier_uniform[OutputScalarT: _Scalar](
    input: _Input, size: _Sigma, n: SupportsIndex, axis: int, output: type[OutputScalarT]
) -> onp.ArrayND[OutputScalarT]: ...
@overload  # output: <T: scalar type> (keyword)
def fourier_uniform[OutputScalarT: _Scalar](
    input: _Input, size: _Sigma, n: SupportsIndex = -1, axis: int = -1, *, output: type[OutputScalarT]
) -> onp.ArrayND[OutputScalarT]: ...
@overload  # +float64
def fourier_uniform(
    input: _InputF64, size: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # +float32
def fourier_uniform(
    input: _InputF32, size: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.float32]: ...
@overload  # ~complex128
def fourier_uniform(
    input: onp.ToJustComplex128_ND, size: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # ~complex64
def fourier_uniform(
    input: onp.ToJustComplex64_ND, size: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex64]: ...
@overload  # fallback
def fourier_uniform(
    input: _Input, size: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[_Scalar]: ...

# NOTE: Keep in sync with `fourier_uniform` and `fourier_gaussian` (but note the different 2nd parameter name)
@overload  # output: <T: ndarray> (positional)
def fourier_ellipsoid[OutputArrayT: onp.ArrayND[_Scalar]](
    input: _Input, size: _Sigma, n: SupportsIndex, axis: int, output: OutputArrayT
) -> OutputArrayT: ...
@overload  # output: <T: ndarray> (keyword)
def fourier_ellipsoid[OutputArrayT: onp.ArrayND[_Scalar]](
    input: _Input, size: _Sigma, n: SupportsIndex = -1, axis: int = -1, *, output: OutputArrayT
) -> OutputArrayT: ...
@overload  # output: <T: scalar type> (positional)
def fourier_ellipsoid[OutputScalarT: _Scalar](
    input: _Input, size: _Sigma, n: SupportsIndex, axis: int, output: type[OutputScalarT]
) -> onp.ArrayND[OutputScalarT]: ...
@overload  # output: <T: scalar type> (keyword)
def fourier_ellipsoid[OutputScalarT: _Scalar](
    input: _Input, size: _Sigma, n: SupportsIndex = -1, axis: int = -1, *, output: type[OutputScalarT]
) -> onp.ArrayND[OutputScalarT]: ...
@overload  # +float64
def fourier_ellipsoid(
    input: _InputF64, size: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # +float32
def fourier_ellipsoid(
    input: _InputF32, size: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.float32]: ...
@overload  # ~complex128
def fourier_ellipsoid(
    input: onp.ToJustComplex128_ND, size: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # ~complex64
def fourier_ellipsoid(
    input: onp.ToJustComplex64_ND, size: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex64]: ...
@overload  # fallback
def fourier_ellipsoid(
    input: _Input, size: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[_Scalar]: ...

# NOTE: Unlike the other three functions, this always returns complex output
@overload  # output: <T: ndarray> (positional)
def fourier_shift[OutputArrayComplexT: onp.ArrayND[_ScalarComplex]](
    input: _Input, shift: _Sigma, n: SupportsIndex, axis: int, output: OutputArrayComplexT
) -> OutputArrayComplexT: ...
@overload  # output: <T: ndarray> (keyword)
def fourier_shift[OutputArrayComplexT: onp.ArrayND[_ScalarComplex]](
    input: _Input, shift: _Sigma, n: SupportsIndex = -1, axis: int = -1, *, output: OutputArrayComplexT
) -> OutputArrayComplexT: ...
@overload  # output: <T: scalar type> (positional)
def fourier_shift[OutputScalarComplexT: _ScalarComplex](
    input: _Input, shift: _Sigma, n: SupportsIndex, axis: int, output: type[OutputScalarComplexT]
) -> onp.ArrayND[OutputScalarComplexT]: ...
@overload  # output: <T: scalar type> (keyword)
def fourier_shift[OutputScalarComplexT: _ScalarComplex](
    input: _Input, shift: _Sigma, n: SupportsIndex = -1, axis: int = -1, *, output: type[OutputScalarComplexT]
) -> onp.ArrayND[OutputScalarComplexT]: ...
@overload  # ~complex128 | +floating
def fourier_shift(
    input: onp.ToJustComplex128_ND | _InputF64 | _InputF32,
    shift: _Sigma,
    n: SupportsIndex = -1,
    axis: int = -1,
    output: None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload  # ~complex64
def fourier_shift(
    input: onp.ToJustComplex64_ND, shift: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[np.complex64]: ...
@overload  # fallback
def fourier_shift(
    input: _Input, shift: _Sigma, n: SupportsIndex = -1, axis: int = -1, output: None = None
) -> onp.ArrayND[_ScalarComplex]: ...
