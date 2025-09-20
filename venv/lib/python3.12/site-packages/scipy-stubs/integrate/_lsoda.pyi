from collections.abc import Callable
from typing import Concatenate, Final, Protocol, TypeAlias, final, type_check_only

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._dop import types as types  # same signature, so no need to redefine

_VecI32: TypeAlias = onp.Array1D[np.int32]
_VecF64: TypeAlias = onp.Array1D[np.float64]

# (t, y) -> ydot
_FnCallback: TypeAlias = Callable[Concatenate[float, _VecF64, ...], onp.ArrayND[npc.floating]]

@type_check_only
@final
class _fortran_lsoda(Protocol):
    def __call__(
        self,
        /,
        f: _FnCallback,
        y: _VecF64,
        t: float,
        tout: float,
        rtol: _VecF64,
        atol: _VecF64,
        itask: int,
        istate: int,
        rwork: _VecF64,
        iwork: _VecI32,
        jac: _FnCallback,
        jt: int,
        f_extra_args: tuple[object, ...] = (),
        overwrite_y: onp.ToBool = 0,
        jac_extra_args: tuple[object, ...] = (),
    ) -> tuple[_VecF64, float, int]: ...  # (y, t, istate)

###

__f2py_numpy_version__: Final[str] = ...
__version__: Final[str] = ...

lsoda: _fortran_lsoda = ...
