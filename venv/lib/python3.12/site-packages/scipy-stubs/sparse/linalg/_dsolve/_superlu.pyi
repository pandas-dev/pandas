from collections.abc import Callable, Mapping
from typing import Any, Final, Generic, Literal, TypeAlias, final, overload
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import csc_array, csc_matrix, csr_matrix

_InexactT_co = TypeVar("_InexactT_co", bound=np.float32 | np.float64 | np.complex64 | np.complex128, default=Any, covariant=True)

_Int1D: TypeAlias = onp.Array1D[np.int32]
_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]
_Complex1D: TypeAlias = onp.Array1D[np.complex128]
_Complex2D: TypeAlias = onp.Array2D[np.complex128]
_Inexact2D: TypeAlias = onp.Array2D[np.float32 | np.float64 | np.complex64 | np.complex128]

_Real: TypeAlias = npc.integer | npc.floating

###

@final
class SuperLU(Generic[_InexactT_co]):
    shape: Final[tuple[int, int]]
    nnz: Final[int]
    perm_r: Final[onp.Array1D[np.intp]]
    perm_c: Final[onp.Array1D[np.intp]]
    L: csc_array[_InexactT_co]  # readonly
    U: csc_array[_InexactT_co]  # readonly

    @overload
    def solve(self, /, rhs: onp.Array1D[_Real]) -> _Float1D: ...
    @overload
    def solve(self, /, rhs: onp.Array1D[npc.complexfloating]) -> _Complex1D: ...
    @overload
    def solve(self, /, rhs: onp.Array2D[_Real]) -> _Float2D: ...
    @overload
    def solve(self, /, rhs: onp.Array2D[npc.complexfloating]) -> _Complex2D: ...
    @overload
    def solve(self, /, rhs: onp.ArrayND[_Real]) -> onp.ArrayND[np.float64]: ...
    @overload
    def solve(self, /, rhs: onp.ArrayND[npc.complexfloating]) -> onp.ArrayND[np.complex128]: ...
    @overload
    def solve(self, /, rhs: onp.ArrayND[npc.number]) -> onp.ArrayND[np.float64 | np.complex128]: ...

def gssv(
    N: op.CanIndex,
    nnz: op.CanIndex,
    nzvals: _Inexact2D,
    colind: _Int1D,
    rowptr: _Int1D,
    B: _Inexact2D,
    csc: onp.ToBool = 0,
    options: Mapping[str, object] = ...,
) -> tuple[csc_matrix | csr_matrix, int]: ...

#
def gstrf(
    N: op.CanIndex,
    nnz: op.CanIndex,
    nzvals: _Inexact2D,
    colind: _Int1D,
    rowptr: _Int1D,
    csc_construct_func: type[csc_array] | Callable[..., csc_array],
    ilu: onp.ToBool = 0,
    options: Mapping[str, object] = ...,
) -> SuperLU: ...

#
def gstrs(
    trans: Literal["N", "T"],
    L_n: op.CanIndex,
    L_nnz: op.CanIndex,
    L_nzvals: _Inexact2D,
    L_rowind: _Int1D,
    L_colptr: _Int1D,
    U_n: op.CanIndex,
    U_nnz: op.CanIndex,
    U_nzvals: _Inexact2D,
    U_rowind: _Int1D,
    U_colptr: _Int1D,
    B: _Inexact2D,
) -> tuple[onp.ArrayND[np.float64 | np.complex128], int]: ...
