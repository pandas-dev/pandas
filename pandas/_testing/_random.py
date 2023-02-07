import string

import numpy as np

from pandas._typing import NpDtype

RANDS_CHARS = np.array(list(string.ascii_letters + string.digits), dtype=(np.str_, 1))


def rands_array(nchars, size, dtype: NpDtype = "O", replace: bool = True) -> np.ndarray:
    """
    Generate an array of byte strings.
    """
    retval = (
        np.random.choice(RANDS_CHARS, size=nchars * np.prod(size), replace=replace)
        .view((np.str_, nchars))
        .reshape(size)
    )
    return retval.astype(dtype)


def rands(nchars) -> str:
    """
    Generate one random byte string.

    See `rands_array` if you want to create an array of random strings.

    """
    return "".join(np.random.choice(RANDS_CHARS, nchars))
