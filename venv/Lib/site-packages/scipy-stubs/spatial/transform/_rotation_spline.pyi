from typing import ClassVar, Final, Literal, overload

import numpy as np
import optype.numpy as onp

from ._rotation import Rotation
from scipy.interpolate import PPoly

def _create_skew_matrix(x: onp.Array2D[np.float64]) -> onp.Array3D[np.float64]: ...
def _matrix_vector_product_of_stacks(A: onp.Array3D[np.float64], b: onp.Array2D[np.float64]) -> onp.Array2D[np.float64]: ...
def _angular_rate_to_rotvec_dot_matrix(rotvecs: onp.Array2D[np.float64]) -> onp.Array3D[np.float64]: ...
def _rotvec_dot_to_angular_rate_matrix(rotvecs: onp.Array2D[np.float64]) -> onp.Array3D[np.float64]: ...
def _angular_acceleration_nonlinear_term(
    rotvecs: onp.Array2D[np.float64], rotvecs_dot: onp.Array2D[np.float64]
) -> onp.Array2D[np.float64]: ...
def _compute_angular_rate(rotvecs: onp.Array2D[np.float64], rotvecs_dot: onp.Array2D[np.float64]) -> onp.Array2D[np.float64]: ...
def _compute_angular_acceleration(
    rotvecs: onp.Array2D[np.float64], rotvecs_dot: onp.Array2D[np.float64], rotvecs_dot_dot: onp.Array2D[np.float64]
) -> onp.Array2D[np.float64]: ...
def _create_block_3_diagonal_matrix(
    A: onp.Array3D[np.float64], B: onp.Array3D[np.float64], d: onp.Array1D[np.float64]
) -> onp.Array2D[np.float64]: ...

class RotationSpline:
    MAX_ITER: ClassVar[Literal[10]] = 10
    TOL: ClassVar[float] = 1e-9

    times: Final[onp.Array1D[np.float64]]
    rotations: Final[Rotation[tuple[int]]]
    interpolator: Final[PPoly[np.float64]]

    #
    def _solve_for_angular_rates(
        self, /, dt: onp.Array1D[np.float64], angular_rates: onp.Array2D[np.float64], rotvecs: onp.Array2D[np.float64]
    ) -> tuple[onp.Array2D[np.float64], onp.Array2D[np.float64]]: ...

    #
    def __init__(self, /, times: onp.ToFloat1D, rotations: Rotation[tuple[int]]) -> None: ...

    #
    @overload
    def __call__(self, /, times: onp.ToFloat, order: Literal[0] = 0) -> Rotation[tuple[()]]: ...
    @overload
    def __call__(self, /, times: onp.ToFloat, order: Literal[1, 2]) -> onp.Array1D[np.float64]: ...
    @overload
    def __call__(self, /, times: onp.ToFloat1D, order: Literal[0] = 0) -> Rotation[tuple[int]]: ...
    @overload
    def __call__(self, /, times: onp.ToFloat1D, order: Literal[1, 2]) -> onp.Array2D[np.float64]: ...
