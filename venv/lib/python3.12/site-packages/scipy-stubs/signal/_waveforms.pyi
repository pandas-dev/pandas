from collections.abc import Iterable
from typing import Literal, TypeAlias, TypeVar, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc
from numpy._typing import _DTypeLike

from scipy._typing import AnyShape

__all__ = ["chirp", "gausspulse", "sawtooth", "square", "sweep_poly", "unit_impulse"]

_SCT = TypeVar("_SCT", bound=np.generic)

_ToFloat0ND: TypeAlias = onp.ToFloat | onp.ToFloatND
_FloatingND: TypeAlias = onp.ArrayND[npc.floating]
_FloatND: TypeAlias = onp.ArrayND[np.float64]

_ChirpMethod: TypeAlias = Literal["linear", "quadratic", "logarithmic", "hyperbolic"]

###

def sawtooth(t: _ToFloat0ND, width: _ToFloat0ND = 1) -> _FloatND: ...
def square(t: _ToFloat0ND, duty: _ToFloat0ND = 0.5) -> _FloatND: ...

#
@overload  # real scalars
def chirp(
    t: onp.ToFloat,
    f0: onp.ToFloat,
    t1: onp.ToFloat,
    f1: onp.ToFloat,
    method: _ChirpMethod = "linear",
    phi: onp.ToFloat = 0,
    vertex_zero: op.CanBool = True,
    *,
    complex: onp.ToFalse = False,
) -> npc.floating: ...
@overload  # real arrays
def chirp(
    t: onp.ToFloatND,
    f0: onp.ToFloat,
    t1: onp.ToFloat,
    f1: onp.ToFloat,
    method: _ChirpMethod = "linear",
    phi: onp.ToFloat = 0,
    vertex_zero: op.CanBool = True,
    *,
    complex: onp.ToFalse = False,
) -> _FloatingND: ...
@overload  # complex scalars
def chirp(
    t: onp.ToFloat,
    f0: onp.ToFloat,
    t1: onp.ToFloat,
    f1: onp.ToFloat,
    method: _ChirpMethod = "linear",
    phi: onp.ToFloat = 0,
    vertex_zero: op.CanBool = True,
    *,
    complex: onp.ToTrue,
) -> npc.complexfloating: ...
@overload  # complex arrays
def chirp(
    t: onp.ToFloatND,
    f0: onp.ToFloat,
    t1: onp.ToFloat,
    f1: onp.ToFloat,
    method: _ChirpMethod = "linear",
    phi: onp.ToFloat = 0,
    vertex_zero: op.CanBool = True,
    *,
    complex: onp.ToTrue,
) -> onp.ArrayND[npc.complexfloating]: ...

#
def sweep_poly(t: _ToFloat0ND, poly: onp.ToFloatND | np.poly1d, phi: onp.ToFloat = 0) -> _FloatND: ...

#
@overload  # dtype is not given
def unit_impulse(
    shape: AnyShape, idx: op.CanIndex | Iterable[op.CanIndex] | Literal["mid"] | None = None, dtype: type[float] = ...
) -> _FloatND: ...
@overload  # dtype is given
def unit_impulse(
    shape: AnyShape, idx: op.CanIndex | Iterable[op.CanIndex] | Literal["mid"] | None, dtype: _DTypeLike[_SCT]
) -> onp.ArrayND[_SCT]: ...

# Overloads for gausspulse when `t` is `"cutoff"`
@overload  # retquad: False = ..., retenv: False = ...
def gausspulse(
    t: Literal["cutoff"],
    fc: onp.ToFloat = 1000,
    bw: onp.ToFloat = 0.5,
    bwr: onp.ToFloat = -6,
    tpr: onp.ToFloat = -60,
    retquad: op.CanBool = False,
    retenv: op.CanBool = False,
) -> np.float64: ...

# Overloads for gausspulse when `t` is scalar
@overload  # retquad: False = ..., retenv: False = ...
def gausspulse(
    t: onp.ToFloat,
    fc: onp.ToFloat = 1000,
    bw: onp.ToFloat = 0.5,
    bwr: onp.ToFloat = -6,
    tpr: onp.ToFloat = -60,
    retquad: onp.ToFalse = False,
    retenv: onp.ToFalse = False,
) -> np.float64: ...
@overload  # retquad: False = ..., retenv: True (keyword)
def gausspulse(
    t: onp.ToFloat,
    fc: onp.ToFloat = 1000,
    bw: onp.ToFloat = 0.5,
    bwr: onp.ToFloat = -6,
    tpr: onp.ToFloat = -60,
    retquad: onp.ToFalse = False,
    *,
    retenv: onp.ToTrue,
) -> tuple[np.float64, np.float64]: ...
@overload  # retquad: False (positional), retenv: False (positional)
def gausspulse(
    t: onp.ToFloat, fc: onp.ToFloat, bw: onp.ToFloat, bwr: onp.ToFloat, tpr: onp.ToFloat, retquad: onp.ToFalse, retenv: onp.ToTrue
) -> tuple[np.float64, np.float64]: ...
@overload  # retquad: True (positional), retenv: False = ...
def gausspulse(
    t: onp.ToFloat,
    fc: onp.ToFloat,
    bw: onp.ToFloat,
    bwr: onp.ToFloat,
    tpr: onp.ToFloat,
    retquad: onp.ToTrue,
    retenv: onp.ToFalse = False,
) -> tuple[np.float64, np.float64]: ...
@overload  # retquad: True (keyword), retenv: False = ...
def gausspulse(
    t: onp.ToFloat,
    fc: onp.ToFloat = 1000,
    bw: onp.ToFloat = 0.5,
    bwr: onp.ToFloat = -6,
    tpr: onp.ToFloat = -60,
    *,
    retquad: onp.ToTrue,
    retenv: onp.ToFalse = False,
) -> tuple[np.float64, np.float64]: ...
@overload  # retquad: True (positional), retenv: True (positional/keyword)
def gausspulse(
    t: onp.ToFloat, fc: onp.ToFloat, bw: onp.ToFloat, bwr: onp.ToFloat, tpr: onp.ToFloat, retquad: onp.ToTrue, retenv: onp.ToTrue
) -> tuple[np.float64, np.float64, np.float64]: ...
@overload  # retquad: True (keyword), retenv: True
def gausspulse(
    t: onp.ToFloat,
    fc: onp.ToFloat = 1000,
    bw: onp.ToFloat = 0.5,
    bwr: onp.ToFloat = -6,
    tpr: onp.ToFloat = -60,
    *,
    retquad: onp.ToTrue,
    retenv: onp.ToTrue,
) -> tuple[np.float64, np.float64, np.float64]: ...

# Overloads for `gausspulse` when `t` is a non-scalar array like
@overload  # retquad: False = ..., retenv: False = ...
def gausspulse(
    t: onp.ToFloatND,
    fc: onp.ToFloat = 1000,
    bw: onp.ToFloat = 0.5,
    bwr: onp.ToFloat = -6,
    tpr: onp.ToFloat = -60,
    retquad: onp.ToFalse = False,
    retenv: onp.ToFalse = False,
) -> _FloatND: ...
@overload  # retquad: False = ..., retenv: True (keyword)
def gausspulse(
    t: onp.ToFloatND,
    fc: onp.ToFloat = 1000,
    bw: onp.ToFloat = 0.5,
    bwr: onp.ToFloat = -6,
    tpr: onp.ToFloat = -60,
    retquad: onp.ToFalse = False,
    *,
    retenv: onp.ToTrue,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # retquad: False (positional), retenv: False (positional)
def gausspulse(
    t: onp.ToFloatND,
    fc: onp.ToFloat,
    bw: onp.ToFloat,
    bwr: onp.ToFloat,
    tpr: onp.ToFloat,
    retquad: onp.ToFalse,
    retenv: onp.ToTrue,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # retquad: True (positional), retenv: False = ...
def gausspulse(
    t: onp.ToFloatND,
    fc: onp.ToFloat,
    bw: onp.ToFloat,
    bwr: onp.ToFloat,
    tpr: onp.ToFloat,
    retquad: onp.ToTrue,
    retenv: onp.ToFalse = False,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # retquad: True (keyword), retenv: False = ...
def gausspulse(
    t: onp.ToFloatND,
    fc: onp.ToFloat = 1000,
    bw: onp.ToFloat = 0.5,
    bwr: onp.ToFloat = -6,
    tpr: onp.ToFloat = -60,
    *,
    retquad: onp.ToTrue,
    retenv: onp.ToFalse = False,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # retquad: True (positional), retenv: True (positional/keyword)
def gausspulse(
    t: onp.ToFloatND,
    fc: onp.ToFloat,
    bw: onp.ToFloat,
    bwr: onp.ToFloat,
    tpr: onp.ToFloat,
    retquad: onp.ToTrue,
    retenv: onp.ToTrue,
) -> tuple[_FloatND, _FloatND, _FloatND]: ...
@overload  # retquad: True (keyword), retenv: True
def gausspulse(
    t: onp.ToFloatND,
    fc: onp.ToFloat = 1000,
    bw: onp.ToFloat = 0.5,
    bwr: onp.ToFloat = -6,
    tpr: onp.ToFloat = -60,
    *,
    retquad: onp.ToTrue,
    retenv: onp.ToTrue,
) -> tuple[_FloatND, _FloatND, _FloatND]: ...
