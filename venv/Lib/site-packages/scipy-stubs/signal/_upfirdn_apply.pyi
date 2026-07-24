# defined in scipy/signal/_upfirdn_apply.pyx

from typing import Literal

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

###

type _Mode = Literal["constant", "symmetric", "edge", "smooth", "wrap", "reflect", "antisymmetric", "antireflect", "line"]
type _ModeCode = Literal[0, 1, 2, 3, 4, 5, 6, 7, 8]

###

def _output_len(len_h: np.int64 | int, in_len: np.int64 | int, up: np.int64 | int, down: np.int64 | int) -> int: ...
def mode_enum(mode: _Mode) -> _ModeCode: ...
def _pad_test[DTypeT: (np.float32, np.float64, np.complex64, np.complex128)](
    data: onp.ArrayND[DTypeT], npre: np.intp | int = 0, npost: np.intp | int = 0, mode: _ModeCode = 0
) -> onp.Array1D[DTypeT]: ...
def _apply[DTypeT: (np.float32, np.float64, np.complex64, np.complex128)](
    data: onp.ArrayND[npc.number],
    h_trans_flip: onp.Array1D[DTypeT],
    out: onp.ArrayND[npc.number],
    up: np.intp | int,
    down: np.intp | int,
    axis: np.intp | int,
    mode: np.intp | int,
    cval: DTypeT | complex,
) -> onp.ArrayND[DTypeT]: ...
