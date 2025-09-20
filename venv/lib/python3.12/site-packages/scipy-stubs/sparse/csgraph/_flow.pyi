from typing import Final, Literal

import numpy as np

from scipy.sparse import csr_array

###

DTYPE: Final[type[np.float64]] = ...
ITYPE: Final[type[np.int32]] = ...

class MaximumFlowResult:
    flow_value: Final[int | np.int32 | np.int64]
    flow: csr_array[np.float64, tuple[int, int]]

    def __init__(self, /, flow_value: int | np.int32 | np.int64, flow: csr_array[np.float64, tuple[int, int]]) -> None: ...

def maximum_flow(
    csgraph: csr_array, source: int, sink: int, *, method: Literal["edmonds_karp", "dinic"] = "dinic"
) -> MaximumFlowResult: ...
