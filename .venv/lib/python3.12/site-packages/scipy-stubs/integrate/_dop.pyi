from collections.abc import Callable
from typing import Concatenate, Final, Protocol, TypeAlias, final, type_check_only

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

_VecI32: TypeAlias = onp.Array1D[np.int32]
_VecF64: TypeAlias = onp.Array1D[np.float64]

# (x, y) -> f
_FnCallback: TypeAlias = Callable[Concatenate[float, _VecF64, ...], onp.ArrayND[npc.floating]]
# (nr, xold, x, y, con, icomp[, nd]) -> irtn
_FnSolOut: TypeAlias = Callable[Concatenate[int, float, float, _VecF64, _VecF64, _VecI32, ...], int]

@type_check_only
@final
class _fortran_dop(Protocol):
    def __call__(
        self,
        /,
        fcn: _FnCallback,
        x: float,
        y: _VecF64,
        xend: float,
        rtol: _VecF64,
        atol: _VecF64,
        solout: _FnSolOut,
        iout: int,
        work: _VecF64,
        iwork: _VecI32,
        fcn_extra_args: tuple[object, ...] = (),
        overwrite_y: onp.ToBool = 0,
        solout_extra_args: tuple[object, ...] = (),
    ) -> tuple[float, _VecF64, _VecI32, int]: ...  # (x, y, iwork, idid)

###

__version__: Final[str] = ...

dopri853: _fortran_dop = ...
dopri5: _fortran_dop = ...
