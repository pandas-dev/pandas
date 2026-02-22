import enum
from typing import TypeAlias

import numpy as np
import optype.numpy as onp

###

_Int1D: TypeAlias = onp.Array1D[np.int32]
_Double1D: TypeAlias = onp.Array1D[np.float64]
_CDouble1D: TypeAlias = onp.Array1D[np.complex128]

###

class Status(enum.IntEnum):  # undocumented
    OK = 0
    NoMemory = 1
    Other = 2

def airyzo(nt: int, kf: int) -> tuple[_Double1D, _Double1D, _Double1D, _Double1D]: ...  # undocumented
def bernob(n: int) -> _Double1D: ...  # undocumented
def cerzo(nt: int) -> _CDouble1D: ...  # undocumented
def cpbdn(n: int, z: complex) -> tuple[_CDouble1D, _CDouble1D]: ...  # undocumented
def cyzo(nt: int, kf: int, kc: int) -> tuple[_CDouble1D, _CDouble1D]: ...  # undocumented
def eulerb(n: int) -> _Double1D: ...  # undocumented
def fcoef(kd: int, m: int, q: float, a: float) -> _Double1D: ...  # undocumented
def fcszo(kf: int, nt: int) -> _CDouble1D: ...  # undocumented
def jdzo(nt: int) -> tuple[_Int1D, _Int1D, _Int1D, _Double1D]: ...  # undocumented
def jyzo(n: int, nt: int) -> tuple[_Double1D, _Double1D, _Double1D, _Double1D]: ...  # undocumented
def klvnzo(nt: int, kd: int) -> _Double1D: ...  # undocumented
def lamn(n: int, x: float) -> tuple[int, _Double1D, _Double1D]: ...  # undocumented
def lamv(v: float, x: float) -> tuple[float, _Double1D, _Double1D]: ...  # undocumented
def pbdv(v: float, x: float) -> tuple[_Double1D, _Double1D, float, float]: ...  # undocumented
def pbvv(v: float, x: float) -> tuple[_Double1D, _Double1D, float, float]: ...  # undocumented
def sdmn(m: int, n: int, c: float, cv: float, kd: int) -> _Double1D: ...  # undocumented
def segv(m: int, n: int, c: float, kd: int) -> tuple[float, _Double1D]: ...  # undocumented
