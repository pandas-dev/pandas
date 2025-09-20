from typing import TypeVar

import numpy as np
import optype.numpy as onp

_FloatT = TypeVar("_FloatT", np.float32, np.float64)
_ShapeT = TypeVar("_ShapeT", bound=tuple[int] | tuple[int, int])

def update_cluster_means(
    obs: onp.ArrayND[_FloatT, _ShapeT], labels: onp.Array1D[np.int32], nc: int
) -> tuple[onp.ArrayND[_FloatT, _ShapeT], onp.Array1D[np.int32]]: ...
def vq(
    obs: onp.ArrayND[_FloatT, tuple[int] | tuple[int, int]], codes: onp.Array2D[_FloatT]
) -> tuple[onp.Array1D[np.int32], onp.Array1D[_FloatT]]: ...
