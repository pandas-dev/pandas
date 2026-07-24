from collections.abc import Callable
from typing import Concatenate, Final, Literal, Protocol, final, type_check_only

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

# (t, y) -> ydot or jac
type _FnCallbackD = Callable[Concatenate[float, onp.Array1D[np.float64], ...], onp.ArrayND[npc.floating]]
type _FnCallbackZ = Callable[Concatenate[float, onp.Array1D[np.complex128], ...], onp.ArrayND[npc.complexfloating]]

@type_check_only
@final
class _fortran_dvode(Protocol):
    def __call__(
        self,
        /,
        f: _FnCallbackD,
        jac: _FnCallbackD,
        y: onp.Array1D[np.float64],
        t: float,
        tout: float,
        rtol: onp.Array1D[np.float64],
        atol: onp.Array1D[np.float64],
        itask: int,
        istate: int,
        rwork: onp.Array1D[np.float64],
        iwork: onp.Array1D[np.int32],
        mf: int,
        f_extra_args: tuple[object, ...] = (),
        jac_extra_args: tuple[object, ...] = (),
        overwrite_y: bool | Literal[0, 1] = 0,
    ) -> tuple[onp.Array1D[np.float64], float, int]: ...  # (y, t, istate)

@type_check_only
@final
class _fortran_zvode(Protocol):
    def __call__(
        self,
        /,
        f: _FnCallbackZ,
        jac: _FnCallbackZ,
        y: onp.Array1D[np.complex128],
        t: float,
        tout: float,
        rtol: onp.Array1D[np.float64],
        atol: onp.Array1D[np.float64],
        itask: int,
        istate: int,
        zwork: onp.Array1D[np.complex128],
        rwork: onp.Array1D[np.float64],
        iwork: onp.Array1D[np.int32],
        mf: int,
        f_extra_args: tuple[object, ...] = (),
        jac_extra_args: tuple[object, ...] = (),
        overwrite_y: bool | Literal[0, 1] = 0,
    ) -> tuple[onp.Array1D[np.float64], float, int]: ...  # (y, t, istate)

###

class error(Exception): ...

dvode: Final[_fortran_dvode] = ...
zvode: Final[_fortran_zvode] = ...
