from collections.abc import Callable
from typing import Final, Literal

import numpy as np
import optype.numpy as onp

from scipy.stats.distributions import rv_continuous

__all__ = ["levy_stable", "levy_stable_gen", "pdf_from_cf_with_fft"]

Cotes: Final[onp.Array2D[np.float64]] = ...
Cotes_table: Final[onp.Array1D[np.object_]] = ...
levy_stable: Final[levy_stable_gen] = ...

class levy_stable_gen(rv_continuous):
    parameterization: Literal["S0", "S1"]
    pdf_default_method: Literal["piecewise", "best", "zolotarev", "dni", "quadrature", "fft-simpson"]
    cdf_default_method: Literal["piecewise", "fft-simpson"]
    quad_eps: float
    piecewise_x_tol_near_zeta: float
    piecewise_alpha_tol_near_one: float
    pdf_fft_min_points_threshold: float | None
    pdf_fft_grid_spacing: float
    pdf_fft_n_points_two_power: float | None
    pdf_fft_interpolation_level: int
    pdf_fft_interpolation_degree: int

def pdf_from_cf_with_fft(
    cf: Callable[[float], complex], h: float = 0.01, q: int = 9, level: int = 3
) -> tuple[onp.Array1D[np.float64], onp.Array1D[np.float64]]: ...
