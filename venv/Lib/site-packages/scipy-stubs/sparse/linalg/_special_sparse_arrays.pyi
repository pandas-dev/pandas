from typing import Any, Final, Generic, Literal, overload, override
from typing_extensions import TypeVar

import numpy as np
import numpy.typing as npt
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import bsr_array, coo_array, csc_array, csr_array, dia_array, dok_array, lil_array
from scipy.sparse.linalg import LinearOperator

__all__ = ["LaplacianNd"]

###

# because `scipy.sparse.sparray` does not implement anything :(
type _SpArray[ScalarT: npc.integer | npc.floating] = (
    bsr_array[ScalarT]
    | coo_array[ScalarT]
    | csc_array[ScalarT]
    | csr_array[ScalarT]
    | dia_array[ScalarT]
    | dok_array[ScalarT]
    | lil_array[ScalarT]
)

type _BCs = Literal["dirichlet", "neumann", "periodic"]

_SCT = TypeVar("_SCT", bound=npc.integer | npc.floating, default=Any)

###

class LaplacianNd(LinearOperator[_SCT], Generic[_SCT]):
    grid_shape: Final[onp.AtLeast1D[Any]]
    boundary_conditions: Final[_BCs]

    @override
    @overload  # default dtype (int8)
    def __new__(cls, /, grid_shape: onp.AtLeast1D, *, boundary_conditions: _BCs = "neumann") -> LaplacianNd[np.int8]: ...  # pyrefly:ignore[bad-override]
    @overload  # know dtype
    def __new__(
        cls, /, grid_shape: onp.AtLeast1D, *, boundary_conditions: _BCs = "neumann", dtype: onp.ToDType[_SCT]
    ) -> LaplacianNd[_SCT]: ...
    @overload  # unknow dtype
    def __new__(
        cls, /, grid_shape: onp.AtLeast1D, *, boundary_conditions: _BCs = "neumann", dtype: npt.DTypeLike
    ) -> LaplacianNd[Any]: ...

    #
    @override
    @overload  # default dtype (int8)
    def __init__(self: LaplacianNd[np.int8], /, grid_shape: onp.AtLeast1D, *, boundary_conditions: _BCs = "neumann") -> None: ...  # pyrefly:ignore[bad-override]
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
