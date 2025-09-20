# NOTE: keep in sync with scipy/stubs/linalg/_decomp_lu_cython.pyi

from typing import Any, TypeVar

import numpy as np
import numpy.typing as npt

# this mimicks the `ctypedef fused lapack_t`
_LapackT = TypeVar("_LapackT", np.float32, np.float64, np.complex64, np.complex128)

def lu_dispatcher(
    a: npt.NDArray[_LapackT], u: npt.NDArray[_LapackT], piv: npt.NDArray[np.integer[Any]], permute_l: bool
) -> None: ...
