from typing_extensions import deprecated

import numpy as np
import optype.numpy as onp

__all__ = ["pade"]

@deprecated("This function is deprecated and will be removed in SciPy 1.20.0. Use `mpmath.pade` instead.")
def pade(an: onp.ToComplex1D, m: onp.ToJustInt, n: onp.ToJustInt | None = None) -> tuple[np.poly1d, np.poly1d]: ...
