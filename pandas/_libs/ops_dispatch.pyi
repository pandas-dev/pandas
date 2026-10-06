from typing import Any

import numpy as np

from pandas._typing import AnyArrayLike

def maybe_dispatch_ufunc_to_dunder_op(
    self: AnyArrayLike,
    ufunc: np.ufunc,
    method: str,
    *inputs: AnyArrayLike,
    **kwargs: Any,
) -> Any: ...
