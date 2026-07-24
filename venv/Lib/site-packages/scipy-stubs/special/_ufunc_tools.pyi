from collections.abc import Callable, Iterable

import numpy as np
import optype.numpy as onp

def _parse_core_ndims(signature: str) -> tuple[int, ...]: ...
def _with_cache_optimization(
    *, name: str, arg_names: Iterable[str], docstring: str, ufunc: np.ufunc, cache_arg_indices: Iterable[int]
) -> Callable[..., onp.Array]: ...
