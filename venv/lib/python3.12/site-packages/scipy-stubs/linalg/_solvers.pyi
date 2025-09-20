from _typeshed import Incomplete
from typing import Final, Literal, TypeAlias, TypeVar, overload

import numpy as np
import optype as op
import optype.numpy as onp

__all__ = [
    "solve_continuous_are",
    "solve_continuous_lyapunov",
    "solve_discrete_are",
    "solve_discrete_lyapunov",
    "solve_lyapunov",
    "solve_sylvester",
]

_InexactT = TypeVar("_InexactT", bound=np.float32 | np.float64 | np.complex64 | np.complex128)

_FloatND: TypeAlias = onp.ArrayND[np.float32 | np.float64]
_ComplexND: TypeAlias = onp.ArrayND[np.complex64 | np.complex128]

_DiscreteMethod: TypeAlias = Literal["direct", "bilinear"]

###

@overload  # real
def solve_sylvester(a: onp.ToFloatND, b: onp.ToFloatND, q: onp.ToFloatND) -> _FloatND: ...
@overload  # ~complex, +complex, +complex
def solve_sylvester(a: onp.ToJustComplexND, b: onp.ToComplexND, q: onp.ToComplexND) -> _ComplexND: ...
@overload  # +complex, ~complex, +complex
def solve_sylvester(a: onp.ToComplexND, b: onp.ToJustComplexND, q: onp.ToComplexND) -> _ComplexND: ...
@overload  # +complex, +complex, ~complex
def solve_sylvester(a: onp.ToComplexND, b: onp.ToComplexND, q: onp.ToJustComplexND) -> _ComplexND: ...

#
@overload  # real
def solve_continuous_lyapunov(a: onp.ToFloatND, q: onp.ToFloatND) -> _FloatND: ...
@overload  # ~complex, +complex
def solve_continuous_lyapunov(a: onp.ToJustComplexND, q: onp.ToComplexND) -> _ComplexND: ...
@overload  # +complex, ~complex
def solve_continuous_lyapunov(a: onp.ToComplexND, q: onp.ToJustComplexND) -> _ComplexND: ...

#
solve_lyapunov: Final = solve_continuous_lyapunov

#
def _solve_discrete_lyapunov_direct(a: onp.Array2D[_InexactT], q: onp.Array2D[_InexactT]) -> onp.Array2D[_InexactT]: ...
def _solve_discrete_lyapunov_bilinear(a: onp.Array2D[_InexactT], q: onp.Array2D[_InexactT]) -> onp.Array2D[_InexactT]: ...

#
@overload  # real
def solve_discrete_lyapunov(a: onp.ToFloatND, q: onp.ToFloatND, method: _DiscreteMethod | None = None) -> _FloatND: ...
@overload  # ~complex, +complex
def solve_discrete_lyapunov(a: onp.ToJustComplexND, q: onp.ToComplexND, method: _DiscreteMethod | None = None) -> _ComplexND: ...
@overload  # +complex, ~complex
def solve_discrete_lyapunov(a: onp.ToComplexND, q: onp.ToJustComplexND, method: _DiscreteMethod | None = None) -> _ComplexND: ...

#
@overload  # real
def solve_continuous_are(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    q: onp.ToFloatND,
    r: onp.ToFloatND,
    e: onp.ToFloatND | None = None,
    s: onp.ToFloatND | None = None,
    balanced: op.CanBool = True,
) -> _FloatND: ...
@overload  # ~complex, +complex, +complex, +complex, +complex?, +complex?
def solve_continuous_are(
    a: onp.ToJustComplexND,
    b: onp.ToComplexND,
    q: onp.ToComplexND,
    r: onp.ToComplexND,
    e: onp.ToComplexND | None = None,
    s: onp.ToComplexND | None = None,
    balanced: op.CanBool = True,
) -> _ComplexND: ...
@overload  # +complex, ~complex, +complex, +complex, +complex?, +complex?
def solve_continuous_are(
    a: onp.ToComplexND,
    b: onp.ToJustComplexND,
    q: onp.ToComplexND,
    r: onp.ToComplexND,
    e: onp.ToComplexND | None = None,
    s: onp.ToComplexND | None = None,
    balanced: op.CanBool = True,
) -> _ComplexND: ...
@overload  # +complex, +complex, ~complex, +complex, +complex?, +complex?
def solve_continuous_are(
    a: onp.ToComplexND,
    b: onp.ToComplexND,
    q: onp.ToJustComplexND,
    r: onp.ToComplexND,
    e: onp.ToComplexND | None = None,
    s: onp.ToComplexND | None = None,
    balanced: op.CanBool = True,
) -> _ComplexND: ...
@overload  # +complex, +complex, +complex, ~complex, +complex?, +complex?
def solve_continuous_are(
    a: onp.ToComplexND,
    b: onp.ToComplexND,
    q: onp.ToComplexND,
    r: onp.ToJustComplexND,
    e: onp.ToComplexND | None = None,
    s: onp.ToComplexND | None = None,
    balanced: op.CanBool = True,
) -> _ComplexND: ...
@overload  # +complex, +complex, +complex, +complex, ~complex, +complex?
def solve_continuous_are(
    a: onp.ToComplexND,
    b: onp.ToComplexND,
    q: onp.ToComplexND,
    r: onp.ToComplexND,
    e: onp.ToJustComplexND,
    s: onp.ToComplexND | None = None,
    balanced: op.CanBool = True,
) -> _ComplexND: ...
@overload  # +complex, +complex, +complex, +complex, +complex, ~complex
def solve_continuous_are(
    a: onp.ToComplexND,
    b: onp.ToComplexND,
    q: onp.ToComplexND,
    r: onp.ToComplexND,
    e: onp.ToComplexND | None,
    s: onp.ToJustComplexND,
    balanced: op.CanBool = True,
) -> _ComplexND: ...
@overload  # +complex, +complex, +complex, +complex, +complex?, *, ~complex
def solve_continuous_are(
    a: onp.ToComplexND,
    b: onp.ToComplexND,
    q: onp.ToComplexND,
    r: onp.ToComplexND,
    e: onp.ToComplexND | None = None,
    *,
    s: onp.ToJustComplexND,
    balanced: op.CanBool = True,
) -> _ComplexND: ...

#
@overload  # real
def solve_discrete_are(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    q: onp.ToFloatND,
    r: onp.ToFloatND,
    e: onp.ToFloatND | None = None,
    s: onp.ToFloatND | None = None,
    balanced: op.CanBool = True,
) -> _FloatND: ...
@overload  # ~complex, +complex, +complex, +complex, +complex?, +complex?
def solve_discrete_are(
    a: onp.ToJustComplexND,
    b: onp.ToComplexND,
    q: onp.ToComplexND,
    r: onp.ToComplexND,
    e: onp.ToComplexND | None = None,
    s: onp.ToComplexND | None = None,
    balanced: op.CanBool = True,
) -> _ComplexND: ...
@overload  # +complex, ~complex, +complex, +complex, +complex?, +complex?
def solve_discrete_are(
    a: onp.ToComplexND,
    b: onp.ToJustComplexND,
    q: onp.ToComplexND,
    r: onp.ToComplexND,
    e: onp.ToComplexND | None = None,
    s: onp.ToComplexND | None = None,
    balanced: op.CanBool = True,
) -> _ComplexND: ...
@overload  # +complex, +complex, ~complex, +complex, +complex?, +complex?
def solve_discrete_are(
    a: onp.ToComplexND,
    b: onp.ToComplexND,
    q: onp.ToJustComplexND,
    r: onp.ToComplexND,
    e: onp.ToComplexND | None = None,
    s: onp.ToComplexND | None = None,
    balanced: op.CanBool = True,
) -> _ComplexND: ...
@overload  # +complex, +complex, +complex, ~complex, +complex?, +complex?
def solve_discrete_are(
    a: onp.ToComplexND,
    b: onp.ToComplexND,
    q: onp.ToComplexND,
    r: onp.ToJustComplexND,
    e: onp.ToComplexND | None = None,
    s: onp.ToComplexND | None = None,
    balanced: op.CanBool = True,
) -> _ComplexND: ...
@overload  # +complex, +complex, +complex, +complex, ~complex, +complex?
def solve_discrete_are(
    a: onp.ToComplexND,
    b: onp.ToComplexND,
    q: onp.ToComplexND,
    r: onp.ToComplexND,
    e: onp.ToJustComplexND,
    s: onp.ToComplexND | None = None,
    balanced: op.CanBool = True,
) -> _ComplexND: ...
@overload  # +complex, +complex, +complex, +complex, +complex, ~complex
def solve_discrete_are(
    a: onp.ToComplexND,
    b: onp.ToComplexND,
    q: onp.ToComplexND,
    r: onp.ToComplexND,
    e: onp.ToComplexND | None,
    s: onp.ToJustComplexND,
    balanced: op.CanBool = True,
) -> _ComplexND: ...
@overload  # +complex, +complex, +complex, +complex, +complex?, *, ~complex
def solve_discrete_are(
    a: onp.ToComplexND,
    b: onp.ToComplexND,
    q: onp.ToComplexND,
    r: onp.ToComplexND,
    e: onp.ToComplexND | None = None,
    *,
    s: onp.ToJustComplexND,
    balanced: op.CanBool = True,
) -> _ComplexND: ...

#
def _are_validate_args(
    a: onp.ToComplexND,
    b: onp.ToComplexND,
    q: onp.ToComplexND,
    r: onp.ToComplexND,
    e: onp.ToComplexND | None,
    s: onp.ToComplexND | None,
    eq_type: Literal["care", "dare"] = "care",
) -> tuple[
    onp.ArrayND[Incomplete],  # a
    onp.ArrayND[Incomplete],  # b
    onp.ArrayND[Incomplete],  # q
    onp.ArrayND[Incomplete],  # r
    onp.ArrayND[Incomplete],  # e
    onp.ArrayND[Incomplete],  # s
    int,  # m
    int,  # n
    type[float | complex],  # r_or_c
    bool,  # gen_or_not
]: ...
