from typing import Any, Final, Generic, Literal, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import numpy.typing as npt
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import bsr_array, coo_array, csc_array, csr_array, dia_array, dok_array, lil_array
from scipy.sparse.linalg import LinearOperator

__all__ = ["LaplacianNd"]

_SCT = TypeVar("_SCT", bound=npc.integer | npc.floating, default=Any)

# because `scipy.sparse.sparray` does not implement anything :(
_SpArray: TypeAlias = (
    bsr_array[_SCT] | coo_array[_SCT] | csc_array[_SCT] | csr_array[_SCT] | dia_array[_SCT] | dok_array[_SCT] | lil_array[_SCT]
)

_BCs: TypeAlias = Literal["dirichlet", "neumann", "periodic"]

###

class LaplacianNd(LinearOperator[_SCT], Generic[_SCT]):
    grid_shape: Final[onp.AtLeast1D]
    boundary_conditions: Final[_BCs]

    @overload  # default dtype (int8)
    def __init__(self: LaplacianNd[np.int8], /, grid_shape: onp.AtLeast1D, *, boundary_conditions: _BCs = "neumann") -> None: ...
    @overload  # know dtype
    def __init__(
        self, /, grid_shape: onp.AtLeast1D, *, boundary_conditions: _BCs = "neumann", dtype: onp.ToDType[_SCT]
    ) -> None: ...
    @overload  # unknow dtype
    def __init__(self, /, grid_shape: onp.AtLeast1D, *, boundary_conditions: _BCs = "neumann", dtype: npt.DTypeLike) -> None: ...

    #
    def eigenvalues(self, /, m: onp.ToJustInt | None = None) -> onp.Array1D[np.float64]: ...
    def eigenvectors(self, /, m: onp.ToJustInt | None = None) -> onp.Array2D[np.float64]: ...
    def toarray(self, /) -> onp.Array2D[_SCT]: ...
    def tosparse(self, /) -> _SpArray[_SCT]: ...
