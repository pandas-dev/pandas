from typing import TypeVar, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])

@overload
def tukeylambda_variance(lam: onp.ToFloat) -> onp.Array0D[np.float64]: ...
@overload
def tukeylambda_variance(
    lam: onp.CanArray[_ShapeT, np.dtype[npc.floating | npc.integer | np.bool_]],
) -> onp.Array[_ShapeT, np.float64]: ...
@overload
def tukeylambda_variance(lam: onp.ToFloatND) -> onp.ArrayND[np.float64]: ...

#
@overload
def tukeylambda_kurtosis(lam: onp.ToFloat) -> onp.Array0D[np.float64]: ...
@overload
def tukeylambda_kurtosis(
    lam: onp.CanArray[_ShapeT, np.dtype[npc.floating | npc.integer | np.bool_]],
) -> onp.Array[_ShapeT, np.float64]: ...
@overload
def tukeylambda_kurtosis(lam: onp.ToFloatND) -> onp.ArrayND[np.float64]: ...
