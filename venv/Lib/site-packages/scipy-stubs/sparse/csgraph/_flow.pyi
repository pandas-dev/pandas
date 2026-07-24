from typing import Final, Literal

import numpy as np
import optype.numpy.compat as npc

from scipy.sparse import csr_array

###

DTYPE: Final[type[np.float64]] = ...
ITYPE: Final[type[np.int32]] = ...

class MaximumFlowResult:
    flow_value: Final[np.int_]
    flow: csr_array[np.int32, tuple[int, int]]

    def __init__(self, /, flow_value: np.int_, flow: csr_array[np.int32, tuple[int, int]]) -> None: ...

def maximum_flow(
    csgraph: csr_array[npc.integer, tuple[int, int]],
    source: int,
    sink: int,
    *,
    method: Literal["edmonds_karp", "dinic"] = "dinic",
) -> MaximumFlowResult: ...
