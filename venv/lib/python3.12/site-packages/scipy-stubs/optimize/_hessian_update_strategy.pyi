from typing import Final, Literal, Protocol, TypeAlias, overload

import numpy as np
import optype.numpy as onp

__all__ = ["BFGS", "SR1", "HessianUpdateStrategy"]

_HessApprox: TypeAlias = Literal["hess", "inv_hess"]
_ExceptionStrategy: TypeAlias = Literal["skip_update", "damp_update"]
_ToInitScale: TypeAlias = onp.ToFloat | onp.ToFloat2D | Literal["auto"]

###

# NOTE: this isn't a Protocol in practice, but it avoids having to deal with`abc`
class HessianUpdateStrategy(Protocol):
    #
    @overload
    def __matmul__(self, p: onp.ToFloat1D, /) -> onp.Array1D[np.float64]: ...
    @overload
    def __matmul__(self, p: onp.ToComplex1D, /) -> onp.Array1D[np.float64 | np.complex128]: ...
    #
    @overload
    def dot(self, /, p: onp.ToFloat1D) -> onp.Array1D[np.float64]: ...
    @overload
    def dot(self, /, p: onp.ToComplex1D) -> onp.Array1D[np.float64 | np.complex128]: ...
    #
    def get_matrix(self, /) -> onp.Array2D[np.float64]: ...
    def initialize(self, /, n: int, approx_type: _HessApprox) -> None: ...
    def update(self, /, delta_x: onp.Array1D[np.float64], delta_grad: onp.Array1D[np.float64]) -> None: ...

class FullHessianUpdateStrategy(HessianUpdateStrategy):
    init_scale: Final[_ToInitScale]

    B: onp.Array2D[np.float64] | None
    H: onp.Array2D[np.float64] | None
    first_iteration: bool | None
    approx_type: _HessApprox | None
    n: int  # set in `initialize`

    def __init__(self, /, init_scale: _ToInitScale = "auto") -> None: ...

class BFGS(FullHessianUpdateStrategy):
    exception_strategy: Final[_ExceptionStrategy]
    min_curvature: Final[float]

    def __init__(
        self,
        /,
        exception_strategy: _ExceptionStrategy = "skip_update",
        min_curvature: float | None = None,
        init_scale: _ToInitScale = "auto",
    ) -> None: ...

class SR1(FullHessianUpdateStrategy):
    min_denominator: Final[float]

    def __init__(self, /, min_denominator: float = 1e-08, init_scale: _ToInitScale = "auto") -> None: ...
