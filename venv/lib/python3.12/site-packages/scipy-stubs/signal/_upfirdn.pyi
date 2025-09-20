from typing import Literal, TypeAlias

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["_output_len", "upfirdn"]

_FIRMode: TypeAlias = Literal["constant", "symmetric", "reflect", "wrap"]
_int64_t: TypeAlias = int | np.int64  # noqa: PYI042

class _UpFIRDn:
    def __init__(self, /, h: onp.ArrayND[npc.floating], x_dtype: np.dtype[npc.floating], up: int, down: int) -> None: ...
    def apply_filter(
        self, /, x: onp.ArrayND[npc.number], axis: int = -1, mode: _FIRMode = "constant", cval: int = 0
    ) -> onp.ArrayND[npc.floating]: ...

def upfirdn(
    h: onp.AnyFloatingArray,
    x: onp.AnyIntegerArray | onp.AnyFloatingArray,
    up: int = 1,
    down: int = 1,
    axis: int = -1,
    mode: _FIRMode = "constant",
    cval: float = 0,
) -> onp.ArrayND[npc.floating]: ...

# originally defined in `scipy/signal/_upfirdn_apply.pyx` (as `(((in_len - 1) * up + len_h) - 1) // down + 1`)
def _output_len(len_h: _int64_t, in_len: _int64_t, up: _int64_t, down: _int64_t) -> _int64_t: ...
