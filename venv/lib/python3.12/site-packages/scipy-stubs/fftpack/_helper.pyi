# NOTE: this ignore is required for `numpy==1.23.5` compat
# pyright: reportUnknownVariableType=false

import numpy as np
import optype.numpy as onp
from numpy.fft import fftfreq, fftshift, ifftshift

__all__ = ["fftfreq", "fftshift", "ifftshift", "next_fast_len", "rfftfreq"]

def rfftfreq(n: onp.ToInt, d: onp.ToFloat = 1.0) -> onp.Array1D[np.float64]: ...
def next_fast_len(target: onp.ToInt) -> int: ...
