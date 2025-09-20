from _typeshed import Incomplete
from types import ModuleType
from typing import Literal, TypeVar, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._typing import AnyShape

_InexactT = TypeVar("_InexactT", bound=npc.inexact)

def next_fast_len(target: op.CanIndex, real: op.CanBool = False) -> int: ...
def prev_fast_len(target: op.CanIndex, real: op.CanBool = False) -> int: ...

#
@overload  # xp: None -> np.fft.fftfreq
def fftfreq(
    n: int | npc.integer, d: onp.ToFloat = 1.0, *, xp: None = None, device: Literal["cpu"] | None = None
) -> onp.Array1D[np.float64]: ...
@overload  # xp: ModuleType -> xp.fft.fftfreq
def fftfreq(n: int, d: float = 1.0, *, xp: ModuleType, device: object | None = None) -> Incomplete: ...

#
@overload  # np.fft.rfftfreq
def rfftfreq(
    n: int | npc.integer, d: onp.ToFloat = 1.0, *, xp: None = None, device: Literal["cpu"] | None = None
) -> onp.Array1D[np.float64]: ...
@overload
def rfftfreq(n: int, d: float = 1.0, *, xp: ModuleType, device: object | None = None) -> Incomplete: ...

#
@overload
def fftshift(x: onp.ToIntND | onp.ToJustFloat64_ND, axes: AnyShape | None = None) -> onp.ArrayND[np.float64]: ...
@overload
def fftshift(x: onp.ToJustComplex128_ND, axes: AnyShape | None = None) -> onp.ArrayND[np.complex128]: ...
@overload
def fftshift(x: onp.ToArrayND[_InexactT, _InexactT], axes: AnyShape | None = None) -> onp.ArrayND[_InexactT]: ...
@overload
def fftshift(x: onp.ToComplexND, axes: AnyShape | None = None) -> onp.ArrayND[npc.inexact]: ...

#
@overload
def ifftshift(x: onp.ToIntND | onp.ToJustFloat64_ND, axes: AnyShape | None = None) -> onp.ArrayND[np.float64]: ...
@overload
def ifftshift(x: onp.ToJustComplex128_ND, axes: AnyShape | None = None) -> onp.ArrayND[np.complex128]: ...
@overload
def ifftshift(x: onp.ToArrayND[_InexactT, _InexactT], axes: AnyShape | None = None) -> onp.ArrayND[_InexactT]: ...
@overload
def ifftshift(x: onp.ToComplexND, axes: AnyShape | None = None) -> onp.ArrayND[npc.inexact]: ...
