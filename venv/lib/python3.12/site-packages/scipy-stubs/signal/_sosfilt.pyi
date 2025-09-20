# defined in scipy/signal/_sosfilt.pyx

from typing import TypeVar

import numpy as np
import optype.numpy as onp

# emulate `ctypedef fused DTYPE_t`
_DTypeT = TypeVar("_DTypeT", np.float32, np.float64, np.longdouble, np.complex64, np.complex128, np.clongdouble, np.object_)

###

def _sosfilt_object(sos: onp.Array2D[np.object_], x: onp.Array2D[np.object_], zi: onp.Array3D[np.object_]) -> None: ...
def _sosfilt(sos: onp.Array2D[_DTypeT], x: onp.Array2D[_DTypeT], zi: onp.Array3D[_DTypeT]) -> None: ...
