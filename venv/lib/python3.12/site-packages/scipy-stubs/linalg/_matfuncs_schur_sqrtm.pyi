from typing import Literal, TypeVar

import numpy as np
import optype.numpy as onp

_InexactT = TypeVar("_InexactT", bound=np.float32 | np.float64 | np.complex64 | np.complex128)
_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])

###

class error(Exception): ...

def recursive_schur_sqrtm(
    A: onp.ArrayND[_InexactT, _ShapeT],
) -> tuple[onp.ArrayND[_InexactT, _ShapeT], Literal[0, 1], Literal[0, 1], int]: ...
