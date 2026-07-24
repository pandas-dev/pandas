from typing import Final

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import csr_array
from scipy.sparse._base import _spbase

###

DTYPE: Final[type[np.float64]] = ...
ITYPE: Final[type[np.int32]] = ...

def minimum_spanning_tree(
    csgraph: onp.ToFloat2D | _spbase[npc.integer | npc.floating, tuple[int, int]], overwrite: bool = False
) -> csr_array[np.float64, tuple[int, int]]: ...
