import numpy as np

REVERSED_NAMES: dict[str, str]
UFUNC_ALIASES: dict[str, str]
DISPATCHED_UFUNCS: set[str]

def maybe_dispatch_ufunc_to_dunder_op(
    self, ufunc: np.ufunc, method: str, *inputs, **kwargs
): ...
