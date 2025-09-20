from typing import Final, TypeAlias

import numpy as np
import optype.numpy as onp

from scipy.sparse import coo_array, coo_matrix, csc_array, csc_matrix, csr_array, csr_matrix

_CSXArray: TypeAlias = csr_array | csr_matrix | csc_array | csc_matrix
_CXXArray: TypeAlias = _CSXArray | coo_array | coo_matrix

DTYPE: Final[type[np.float64]] = ...
ITYPE: Final[type[np.int32]] = ...

def reverse_cuthill_mckee(graph: _CSXArray, symmetric_mode: bool = False) -> onp.Array1D[np.int32]: ...
def structural_rank(graph: _CXXArray) -> np.intp: ...
