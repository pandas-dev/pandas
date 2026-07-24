import numpy as np
import optype.numpy as onp

__pythran__: tuple[str, str] = ...

def _max_len_seq_inner(
    taps: onp.Array1D[np.intp],
    state: onp.Array1D[np.int8],
    nbits: int | np.intp,
    length: int | np.intp,
    out: onp.Array1D[np.int8],
) -> onp.Array1D[np.int8]: ...
