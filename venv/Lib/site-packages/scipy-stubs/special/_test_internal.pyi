from typing import Literal

import numpy as np
import optype.numpy as onp

def have_fenv() -> bool: ...  # undocumented
def random_double(size: int, rng: np.random.RandomState) -> onp.Array1D[np.float64]: ...  # undocumented
def test_add_round(size: int, mode: Literal["up", "down"], rng: np.random.RandomState) -> None: ...  # undocumented
def _dd_exp(xhi: float, xlo: float) -> tuple[float, float]: ...  # undocumented
def _dd_log(xhi: float, xlo: float) -> tuple[float, float]: ...  # undocumented
