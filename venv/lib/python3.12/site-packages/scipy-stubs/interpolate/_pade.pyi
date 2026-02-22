import numpy as np
import optype.numpy as onp

__all__ = ["pade"]

def pade(an: onp.ToComplex1D, m: onp.ToJustInt, n: onp.ToJustInt | None = None) -> tuple[np.poly1d, np.poly1d]: ...
