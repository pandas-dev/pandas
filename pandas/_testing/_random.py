from __future__ import annotations

import string
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from pandas._typing import NpDtype
RANDS_CHARS = np.array(list(string.ascii_letters + string.digits), dtype=(np.str_, 1))


def rands_array(
    nchars, size: int, dtype: NpDtype = "O", replace: bool = True
) -> np.ndarray:
    """
    Generate an array of byte strings.
    """
    retval = (
        np.random.default_rng(2)
        .choice(RANDS_CHARS, size=nchars * np.prod(size), replace=replace)
        .view((np.str_, nchars))
        .reshape(size)
    )
    return retval.astype(dtype)
