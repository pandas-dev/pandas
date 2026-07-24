from collections.abc import Callable
from typing import Concatenate, overload

import numpy as np
import optype.numpy as onp

__all__ = ["convolve", "convolve_z", "destroy_convolve_cache", "init_convolution_kernel"]

# NOTE: this doesn't do anything; nothing is cached
# undocumented
def destroy_convolve_cache() -> None: ...

#
def convolve(
    inout: onp.ToFloat1D, omega: onp.ToFloat1D, swap_real_imag: bool = False, overwrite_x: bool = False
) -> onp.Array1D[np.float64]: ...

# undocumented
def convolve_z(
    inout: onp.ToFloat1D, omega_real: onp.ToFloat1D, omega_imag: onp.ToFloat1D, overwrite_x: bool = False
) -> onp.Array1D[np.float64]: ...

# undocumented
@overload
def init_convolution_kernel(
    n: int,
    kernel_func: Callable[[int], float],
    d: int = 0,
    zero_nyquist: int | None = None,
    kernel_func_extra_args: tuple[()] = (),
) -> onp.Array1D[np.float64]: ...
@overload
def init_convolution_kernel(
    n: int,
    kernel_func: Callable[Concatenate[int, ...], float],
    d: int,
    zero_nyquist: int | None,
    kernel_func_extra_args: tuple[object, ...],
) -> onp.Array1D[np.float64]: ...
@overload
def init_convolution_kernel(
    n: int,
    kernel_func: Callable[Concatenate[int, ...], float],
    d: int = 0,
    zero_nyquist: int | None = None,
    *,
    kernel_func_extra_args: tuple[object, ...],
) -> onp.Array1D[np.float64]: ...
