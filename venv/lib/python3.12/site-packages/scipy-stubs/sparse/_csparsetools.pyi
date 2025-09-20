from collections.abc import Iterable
from typing import TypeAlias, TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

_Index: TypeAlias = int | np.intp
_Scalar: TypeAlias = npc.number | np.bool_

_FusedScalarT = TypeVar(
    "_FusedScalarT",
    np.bool_,
    np.uint8,
    np.int8,
    np.uint16,
    np.int16,
    np.uint32,
    np.int32,
    np.uint64,
    np.int64,
    np.float32,
    np.float64,
    npc.floating80,  # longdouble
    np.complex64,
    np.complex128,
    npc.complexfloating160,  # clongdouble
)

###

#
def lil_get1(
    M: _Index, N: _Index, rows: list[list[_FusedScalarT]], datas: list[list[_FusedScalarT]], i: int, j: int
) -> _FusedScalarT: ...
def lil_insert(
    M: _Index, N: _Index, rows: list[list[_FusedScalarT]], datas: list[list[_FusedScalarT]], i: int, j: int, x: _FusedScalarT
) -> None: ...
def lil_get_lengths(input: list[list[_FusedScalarT]], output: onp.Array1D[npc.integer]) -> None: ...
def lil_flatten_to_array(input: onp.Array1D[np.object_ | np.float64], output: onp.Array1D[_Scalar]) -> None: ...
def lil_fancy_get(
    M: _Index,
    N: _Index,
    rows: list[list[_FusedScalarT]],
    datas: list[list[_FusedScalarT]],
    new_rows: list[list[_FusedScalarT]],
    new_datas: list[list[_FusedScalarT]],
    i_idx: onp.Array1D[np.intp],
    j_idx: onp.Array1D[np.intp],
) -> None: ...
def lil_fancy_set(
    M: _Index,
    N: _Index,
    rows: list[list[_FusedScalarT]],
    data: list[list[_FusedScalarT]],
    i_idx: onp.Array1D[np.intp],
    j_idx: onp.Array1D[np.intp],
    values: onp.Array2D[_Scalar],
) -> None: ...
def lil_get_row_ranges(
    M: _Index,
    N: _Index,
    rows: list[list[_FusedScalarT]],
    datas: list[list[_FusedScalarT]],
    new_rows: list[list[_FusedScalarT]],
    new_datas: list[list[_FusedScalarT]],
    irows: Iterable[_Index],
    j_start: _Index,
    j_stop: _Index,
    j_stride: _Index,
    nj: _Index,
) -> None: ...
