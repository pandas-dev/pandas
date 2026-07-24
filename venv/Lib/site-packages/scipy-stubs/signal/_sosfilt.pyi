# defined in scipy/signal/_sosfilt.pyx

import numpy as np
import optype.numpy as onp

def _sosfilt_object(sos: onp.Array2D[np.object_], x: onp.Array2D[np.object_], zi: onp.Array3D[np.object_]) -> None: ...
def _sosfilt[
    # emulate `ctypedef fused DTYPE_t`
    DTypeT: (np.float32, np.float64, np.longdouble, np.complex64, np.complex128, np.clongdouble, np.object_)
](sos: onp.Array2D[DTypeT], x: onp.Array2D[DTypeT], zi: onp.Array3D[DTypeT]) -> None: ...
