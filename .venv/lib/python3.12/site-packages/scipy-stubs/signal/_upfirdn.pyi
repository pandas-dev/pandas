from _typeshed import Incomplete
from typing import Literal, TypeAlias, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["_output_len", "upfirdn"]

###

_FIRMode: TypeAlias = Literal["constant", "symmetric", "reflect", "wrap"]
_int64_t: TypeAlias = int | np.int64  # noqa: PYI042

###

class _UpFIRDn:  # undocumented
    def __init__(self, /, h: onp.ArrayND[npc.floating], x_dtype: np.dtype[npc.floating], up: int, down: int) -> None: ...
    def apply_filter(
        self, /, x: onp.ArrayND[npc.number], axis: int = -1, mode: _FIRMode = "constant", cval: int = 0
    ) -> onp.ArrayND[npc.floating]: ...

# The mypy `overload-overlap` errors are false positives
@overload  # ~f64, ~f64
def upfirdn(  # type: ignore[overload-overlap]
    h: onp.ToJustFloat64_1D | onp.ToInt1D,
    x: onp.ToJustFloat64_ND | onp.ToIntND,
    up: int = 1,
    down: int = 1,
    axis: int = -1,
    mode: _FIRMode = "constant",
    cval: float = 0,
) -> onp.ArrayND[np.float64]: ...
@overload  # +f32, ~f32
def upfirdn(
    h: onp.ToFloat32_1D,
    x: onp.ToJustFloat32_1D | onp.ToJustFloat16_1D,
    up: int = 1,
    down: int = 1,
    axis: int = -1,
    mode: _FIRMode = "constant",
    cval: float = 0,
) -> onp.ArrayND[np.float32]: ...
@overload  # ~f32, +f32
def upfirdn(
    h: onp.ToJustFloat32_1D | onp.ToJustFloat16_1D,
    x: onp.ToFloat32_ND,
    up: int = 1,
    down: int = 1,
    axis: int = -1,
    mode: _FIRMode = "constant",
    cval: float = 0,
) -> onp.ArrayND[np.float32]: ...
@overload  # +c128, ~c128
def upfirdn(  # type: ignore[overload-overlap]
    h: onp.ToComplex128_1D,
    x: onp.ToJustComplex128_ND,
    up: int = 1,
    down: int = 1,
    axis: int = -1,
    mode: _FIRMode = "constant",
    cval: float = 0,
) -> onp.ArrayND[np.complex128]: ...
@overload  # ~c128, +c128
def upfirdn(  # type: ignore[overload-overlap]
    h: onp.ToJustComplex128_1D,
    x: onp.ToComplex128_ND,
    up: int = 1,
    down: int = 1,
    axis: int = -1,
    mode: _FIRMode = "constant",
    cval: float = 0,
) -> onp.ArrayND[np.complex128]: ...
@overload  # +c64, ~c64
def upfirdn(
    h: onp.ToComplex64_1D,
    x: onp.ToJustComplex64_ND,
    up: int = 1,
    down: int = 1,
    axis: int = -1,
    mode: _FIRMode = "constant",
    cval: float = 0,
) -> onp.ArrayND[np.complex64]: ...
@overload  # ~c64, +c64
def upfirdn(
    h: onp.ToJustComplex64_1D,
    x: onp.ToComplex64_ND,
    up: int = 1,
    down: int = 1,
    axis: int = -1,
    mode: _FIRMode = "constant",
    cval: float = 0,
) -> onp.ArrayND[np.complex64]: ...
@overload  # fallback
def upfirdn(
    h: onp.ToComplex1D,
    x: onp.ToComplexND,
    up: int = 1,
    down: int = 1,
    axis: int = -1,
    mode: _FIRMode = "constant",
    cval: float = 0,
) -> onp.ArrayND[Incomplete]: ...

# originally defined in `scipy/signal/_upfirdn_apply.pyx` (as `(((in_len - 1) * up + len_h) - 1) // down + 1`)
def _output_len(len_h: _int64_t, in_len: _int64_t, up: _int64_t, down: _int64_t) -> int: ...  # undocumented
