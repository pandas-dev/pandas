# defined in scipy/ndimage/src/_rank_filter_1d.cpp

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

def rank_filter(
    in_arr_obj: onp.ArrayND[np.float32 | np.float64 | np.int64],
    rank: int,
    win_len: int,
    out_arr_obj: onp.ArrayND[np.float32 | np.float64 | np.int64],
    mode: int,
    cval_obj: onp.ArrayND[npc.floating | npc.integer],
    origin: int,
    /,
) -> None: ...  # undocumented
