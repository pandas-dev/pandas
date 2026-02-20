import numpy as np

__all__ = ["replace"]


def replace(a, old, new):
    "Slow replace (inplace) used for unaccelerated dtypes."
    if type(a) is not np.ndarray:
        raise TypeError("`a` must be a numpy array.")
    if not issubclass(a.dtype.type, np.inexact):
        if old != old:
            # int arrays do not contain NaN
            return
        if int(old) != old:
            raise ValueError("Cannot safely cast `old` to int.")
        if int(new) != new:
            raise ValueError("Cannot safely cast `new` to int.")
    if old != old:
        mask = np.isnan(a)
    else:
        mask = a == old
    np.putmask(a, mask, new)
