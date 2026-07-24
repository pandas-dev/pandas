from collections.abc import Callable
from typing import Literal, overload

import numpy.typing as npt
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

###

type _LaplacianFunction = Callable[[onp.ToComplex2D], onp.Array2D[npc.number]]
type _LaplacianMatrix = onp.Array2D[npc.number] | _spbase | LinearOperator
type _LaplacianDiag = onp.Array1D[npc.number]
type _ToCSGraph = onp.ToComplex2D | _spbase
type _Form = Literal["array", "lo"]

###

@overload
def laplacian(
    csgraph: _ToCSGraph,
    normed: bool = False,
    return_diag: onp.ToFalse = False,
    use_out_degree: bool = False,
    *,
    copy: bool = True,
    form: Literal["function"],
    dtype: npt.DTypeLike | None = None,
    symmetrized: bool = False,
) -> _LaplacianFunction: ...
@overload
def laplacian(
    csgraph: _ToCSGraph,
    normed: bool = False,
    *,
    return_diag: onp.ToTrue,
    use_out_degree: bool = False,
    copy: bool = True,
    form: Literal["function"],
    dtype: npt.DTypeLike | None = None,
    symmetrized: bool = False,
) -> tuple[_LaplacianFunction, _LaplacianDiag]: ...
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
