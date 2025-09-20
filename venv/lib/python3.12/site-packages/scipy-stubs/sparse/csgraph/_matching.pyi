from typing import Final, Literal, TypeAlias

import numpy as np
import optype.numpy as onp

from scipy.sparse import csr_array, csr_matrix
from scipy.sparse._coo import coo_array, coo_matrix
from scipy.sparse._csc import csc_array, csc_matrix

_CXXArray: TypeAlias = csr_array | csr_matrix | csc_array | csc_matrix | coo_array | coo_matrix
_IntVector: TypeAlias = onp.Array1D[np.int32 | np.intp]

DTYPE: Final[type[np.float64]] = ...
ITYPE: Final[type[np.int32]] = ...
BTYPE: Final[type[np.uint8]] = ...

def maximum_bipartite_matching(graph: _CXXArray, perm_type: Literal["row", "column"] = "row") -> _IntVector: ...
def min_weight_full_bipartite_matching(biadjacency: _CXXArray, maximize: bool = False) -> tuple[_IntVector, _IntVector]: ...
