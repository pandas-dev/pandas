# scipy/linalg/_matfuncsmodule.c

from typing import Literal

import numpy as np
import optype.numpy as onp

###

type _Inexact = np.float32 | np.float64 | np.complex64 | np.complex128

###

class error(Exception): ...

#
def recursive_schur_sqrtm[InexactT: _Inexact, ShapeT: onp.AtLeast2D](
    A: onp.ArrayND[InexactT, ShapeT], /
) -> tuple[onp.ArrayND[InexactT, ShapeT], Literal[0, 1], Literal[0, 1], int]: ...

#
def matrix_exponential[InexactT: _Inexact, ShapeT: onp.AtLeast2D](
    A: onp.ArrayND[InexactT, ShapeT], /
) -> tuple[onp.ArrayND[InexactT, ShapeT], int]: ...
