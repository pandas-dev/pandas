from _typeshed import Incomplete
from types import ModuleType
from typing import Any, Literal, SupportsIndex, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._typing import AnyShape

def next_fast_len(target: SupportsIndex, real: bool = False) -> int: ...
def prev_fast_len(target: SupportsIndex, real: bool = False) -> int: ...

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
def fftshift[ScalarT: np.generic, ShapeT: tuple[int, ...]](
    x: onp.ArrayND[ScalarT, ShapeT], axes: AnyShape | None = None
) -> onp.ArrayND[ScalarT, ShapeT]: ...
@overload
def fftshift(x: onp.ToJustBoolND, axes: AnyShape | None = None) -> onp.ArrayND[np.bool]: ...
@overload
def fftshift(x: onp.ToJustInt64_ND, axes: AnyShape | None = None) -> onp.ArrayND[np.int64]: ...
@overload
def fftshift(x: onp.ToJustFloat64_ND, axes: AnyShape | None = None) -> onp.ArrayND[np.float64]: ...
@overload
def fftshift(x: onp.ToJustComplex128_ND, axes: AnyShape | None = None) -> onp.ArrayND[np.complex128]: ...
@overload
def fftshift[InexactT: npc.inexact](
    x: onp.ToArrayND[InexactT, InexactT], axes: AnyShape | None = None
) -> onp.ArrayND[InexactT]: ...
@overload
def fftshift(x: onp.ToComplexND, axes: AnyShape | None = None) -> onp.ArrayND[Any]: ...

#
@overload
def ifftshift[ScalarT: np.generic, ShapeT: tuple[int, ...]](
    x: onp.ArrayND[ScalarT, ShapeT], axes: AnyShape | None = None
) -> onp.ArrayND[ScalarT, ShapeT]: ...
@overload
def ifftshift(x: onp.ToJustBoolND, axes: AnyShape | None = None) -> onp.ArrayND[np.bool]: ...
@overload
def ifftshift(x: onp.ToJustInt64_ND, axes: AnyShape | None = None) -> onp.ArrayND[np.int64]: ...
@overload
def ifftshift(x: onp.ToJustFloat64_ND, axes: AnyShape | None = None) -> onp.ArrayND[np.float64]: ...
@overload
def ifftshift(x: onp.ToJustComplex128_ND, axes: AnyShape | None = None) -> onp.ArrayND[np.complex128]: ...
@overload
def ifftshift[InexactT: npc.inexact](
    x: onp.ToArrayND[InexactT, InexactT], axes: AnyShape | None = None
) -> onp.ArrayND[InexactT]: ...
@overload
def ifftshift(x: onp.ToComplexND, axes: AnyShape | None = None) -> onp.ArrayND[Any]: ...
