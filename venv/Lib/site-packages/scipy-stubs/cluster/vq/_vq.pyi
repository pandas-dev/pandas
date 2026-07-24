import numpy as np
import optype.numpy as onp

#
def update_cluster_means[FloatT: np.float32 | np.float64, ShapeT: tuple[int] | tuple[int, int]](
    obs: onp.ArrayND[FloatT, ShapeT], labels: onp.Array1D[np.int32], nc: int
) -> tuple[onp.ArrayND[FloatT, ShapeT], onp.Array1D[np.int32]]: ...

#
def vq[FloatT: (np.float32, np.float64)](
    obs: onp.ArrayND[FloatT, tuple[int] | tuple[int, int]], codes: onp.Array2D[FloatT]
) -> tuple[onp.Array1D[np.int32], onp.Array1D[FloatT]]: ...
