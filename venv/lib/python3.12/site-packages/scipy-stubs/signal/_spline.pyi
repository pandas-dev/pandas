from typing import TypeAlias, overload

import numpy as np
import optype.numpy as onp

_Float: TypeAlias = np.float32 | np.float64
_Complex: TypeAlias = np.complex64 | np.complex128

_ToFloat: TypeAlias = float | _Float

_FloatND: TypeAlias = onp.ArrayND[_Float]
_ComplexND: TypeAlias = onp.ArrayND[_Complex]
_InexactND: TypeAlias = onp.ArrayND[_Float | _Complex]

###

@overload
def sepfir2d(input: _FloatND, hrow: _FloatND, hcol: _FloatND) -> _FloatND: ...
@overload
def sepfir2d(input: _ComplexND, hrow: _ComplexND, hcol: _ComplexND) -> _ComplexND: ...
@overload
def sepfir2d(input: _InexactND, hrow: _InexactND, hcol: _InexactND) -> _InexactND: ...

#
@overload
def symiirorder1_ic(signal: _FloatND, c0: _ToFloat, z1: _ToFloat, precision: _ToFloat) -> _FloatND: ...
@overload
def symiirorder1_ic(signal: _ComplexND, c0: _ToFloat, z1: _ToFloat, precision: _ToFloat) -> _ComplexND: ...
@overload
def symiirorder1_ic(signal: _InexactND, c0: _ToFloat, z1: _ToFloat, precision: _ToFloat) -> _InexactND: ...

#
def symiirorder2_ic_fwd(signal: _FloatND, r: _ToFloat, omega: _ToFloat, precision: _ToFloat) -> _FloatND: ...
def symiirorder2_ic_bwd(signal: _FloatND, r: _ToFloat, omega: _ToFloat, precision: _ToFloat) -> _FloatND: ...
