from typing import overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

# undocumented
@overload
def tukeylambda_variance(lam: onp.ToFloat) -> onp.Array0D[np.float64]: ...
@overload
def tukeylambda_variance[ShapeT: tuple[int, ...]](
    lam: onp.CanArray[ShapeT, np.dtype[npc.floating | npc.integer | np.bool]],
) -> onp.Array[ShapeT, np.float64]: ...
@overload
def tukeylambda_variance(lam: onp.ToFloatND) -> onp.ArrayND[np.float64]: ...

# undocumented
@overload
def tukeylambda_kurtosis(lam: onp.ToFloat) -> onp.Array0D[np.float64]: ...
@overload
def tukeylambda_kurtosis[ShapeT: tuple[int, ...]](
    lam: onp.CanArray[ShapeT, np.dtype[npc.floating | npc.integer | np.bool]],
) -> onp.Array[ShapeT, np.float64]: ...
@overload
def tukeylambda_kurtosis(lam: onp.ToFloatND) -> onp.ArrayND[np.float64]: ...
