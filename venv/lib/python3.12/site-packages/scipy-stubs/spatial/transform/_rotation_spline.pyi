from typing import ClassVar, Final, Literal, overload

import numpy as np
import optype.numpy as onp

from ._rotation import Rotation
from scipy.interpolate import PPoly

class RotationSpline:
    MAX_ITER: ClassVar[int] = 10
    TOL: ClassVar[float] = 1e-9

    times: Final[onp.Array1D[np.int32 | np.int64 | np.float32 | np.float64]]
    rotations: Final[Rotation]
    interpolator: Final[PPoly]

    def __init__(self, /, times: onp.ToFloat1D, rotations: Rotation) -> None: ...
    #
    @overload
    def __call__(self, /, times: onp.ToFloat1D, order: Literal[0] = 0) -> Rotation | onp.ArrayND[np.float64]: ...
    @overload
    def __call__(self, /, times: onp.ToFloat1D, order: Literal[1, 2]) -> onp.ArrayND[np.float64]: ...
