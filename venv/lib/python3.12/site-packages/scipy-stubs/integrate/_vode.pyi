from collections.abc import Callable
from typing import Concatenate, Final, Protocol, TypeAlias, final, type_check_only

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._dop import types as types  # same signature, so no need to redefine

_VecI32: TypeAlias = onp.Array1D[np.int32]
_VecF64: TypeAlias = onp.Array1D[np.float64]
_VecC128: TypeAlias = onp.Array1D[np.complex128]

# (t, y) -> ydot or jac
_FnCallbackD: TypeAlias = Callable[Concatenate[float, _VecF64, ...], onp.ArrayND[npc.floating]]
_FnCallbackZ: TypeAlias = Callable[Concatenate[float, _VecC128, ...], onp.ArrayND[npc.complexfloating]]

@type_check_only
@final
class _fortran_dvode(Protocol):
    def __call__(
        self,
        /,
        f: _FnCallbackD,
        jac: _FnCallbackD,
        y: _VecF64,
        t: float,
        tout: float,
        rtol: _VecF64,
        atol: _VecF64,
        itask: int,
        istate: int,
        rwork: _VecF64,
        iwork: _VecI32,
        mf: int,
        f_extra_args: tuple[object, ...] = (),
        jac_extra_args: tuple[object, ...] = (),
        overwrite_y: onp.ToBool = 0,
    ) -> tuple[_VecF64, float, int]: ...  # (y, t, istate)

@type_check_only
@final
class _fortran_zvode(Protocol):
    def __call__(
        self,
        /,
        f: _FnCallbackZ,
        jac: _FnCallbackZ,
        y: _VecC128,
        t: float,
        tout: float,
        rtol: _VecF64,
        atol: _VecF64,
        itask: int,
        istate: int,
        zwork: _VecC128,
        rwork: _VecF64,
        iwork: _VecI32,
        mf: int,
        f_extra_args: tuple[object, ...] = (),
        jac_extra_args: tuple[object, ...] = (),
        overwrite_y: onp.ToBool = 0,
    ) -> tuple[_VecF64, float, int]: ...  # (y, t, istate)

###

__f2py_numpy_version__: Final[str] = ...
__version__: Final[str] = ...

dvode: _fortran_dvode = ...
zvode: _fortran_zvode = ...
