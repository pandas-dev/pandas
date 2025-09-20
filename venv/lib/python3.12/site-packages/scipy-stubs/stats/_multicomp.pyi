from dataclasses import dataclass

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._common import ConfidenceInterval
from ._typing import Alternative

__all__ = ["dunnett"]

@dataclass
class DunnettResult:
    statistic: onp.Array1D[np.float64]
    pvalue: onp.Array1D[np.float64]

    _alternative: Alternative
    _rho: onp.Array2D[np.float64]
    _df: int
    _std: np.float64
    _mean_samples: onp.Array1D[np.float64]
    _mean_control: np.float64  # incorrectly annotated as `ndarray` at runtime
    _n_samples: onp.Array1D[np.int_]
    _n_control: int
    _rng: np.random.Generator | np.random.RandomState

    _ci: ConfidenceInterval | None = None
    _ci_cl: float | npc.floating | None = None

    def confidence_interval(self, /, confidence_level: float | npc.floating = 0.95) -> ConfidenceInterval: ...

def dunnett(
    *samples: onp.ToFloat1D,
    control: onp.ToFloat1D,
    alternative: Alternative = "two-sided",
    rng: onp.random.ToRNG | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> DunnettResult: ...
