from typing import Literal, TypeAlias, overload

import numpy.typing as npt
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

_LaplacianMatrix: TypeAlias = onp.Array2D[npc.number] | _spbase | LinearOperator
_LaplacianDiag: TypeAlias = onp.Array1D[npc.number]
_ToCSGraph: TypeAlias = onp.ToComplex2D | _spbase
_Form: TypeAlias = Literal["array", "function", "lo"]

###

@overload
def laplacian(
    csgraph: _ToCSGraph,
    normed: bool = False,
    return_diag: onp.ToFalse = False,
    use_out_degree: bool = False,
    *,
    copy: bool = True,
    form: _Form = "array",
    dtype: npt.DTypeLike | None = None,
    symmetrized: bool = False,
) -> _LaplacianMatrix: ...
@overload
def laplacian(
    csgraph: _ToCSGraph,
    normed: bool,
    return_diag: onp.ToTrue,
    use_out_degree: bool = False,
    *,
    copy: bool = True,
    form: _Form = "array",
    dtype: npt.DTypeLike | None = None,
    symmetrized: bool = False,
) -> tuple[_LaplacianMatrix, _LaplacianDiag]: ...
@overload
def laplacian(
    csgraph: _ToCSGraph,
    normed: bool = False,
    *,
    return_diag: onp.ToTrue,
    use_out_degree: bool = False,
    copy: bool = True,
    form: _Form = "array",
    dtype: npt.DTypeLike | None = None,
    symmetrized: bool = False,
) -> tuple[_LaplacianMatrix, _LaplacianDiag]: ...
