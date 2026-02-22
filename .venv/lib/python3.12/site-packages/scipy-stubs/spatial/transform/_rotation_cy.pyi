from types import EllipsisType
from typing import Literal, Never, TypeAlias, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

_bint: TypeAlias = bool | Literal[0, 1]  # noqa: PYI042

_bool_1d: TypeAlias = onp.Array1D[np.bool_]  # noqa: PYI042
_f64_1d: TypeAlias = onp.Array1D[np.float64]  # noqa: PYI042
_f64_2d: TypeAlias = onp.Array2D[np.float64]  # noqa: PYI042
_f64_3d: TypeAlias = onp.Array3D[np.float64]  # noqa: PYI042
_f64_nd: TypeAlias = onp.ArrayND[np.float64]  # noqa: PYI042
_f64_xd: TypeAlias = onp.Array[tuple[Never, ...], np.float64]  # noqa: PYI042
_isize_2d: TypeAlias = onp.Array2D[np.intp]  # noqa: PYI042

_Indexer: TypeAlias = int | slice | EllipsisType | None
_Order: TypeAlias = Literal["e", "extrinsic", "i", "intrinsic"]

###

def compose_quat(p: _f64_2d, q: _f64_2d) -> _f64_2d: ...
def from_quat(quat: _f64_2d, normalize: _bint = True, copy: _bint = True, scalar_first: _bint = False) -> _f64_2d: ...

#
@overload  # workaround for https://github.com/microsoft/pyright/issues/10232
def from_euler(seq: str, angles: onp.Array[tuple[Never, ...], npc.floating], degrees: _bint = False) -> _f64_nd: ...
@overload
def from_euler(seq: str, angles: onp.ToFloat | onp.ToFloatStrict1D, degrees: _bint = False) -> _f64_1d: ...
@overload
def from_euler(seq: str, angles: onp.ToFloatStrict2D, degrees: _bint = False) -> _f64_2d: ...

#
@overload  # workaround for https://github.com/microsoft/pyright/issues/10232
def from_matrix(matrix: onp.Array[tuple[Never, ...], npc.floating], assume_valid: _bint = False) -> _f64_nd: ...
@overload
def from_matrix(matrix: onp.ToFloatStrict2D, assume_valid: _bint = False) -> _f64_1d: ...
@overload
def from_matrix(matrix: onp.ToFloatStrict3D, assume_valid: _bint = False) -> _f64_2d: ...

#
@overload  # workaround for https://github.com/microsoft/pyright/issues/10232
def _from_matrix_orthogonal(matrix: _f64_xd) -> _f64_nd: ...
@overload
def _from_matrix_orthogonal(matrix: _f64_2d) -> _f64_1d: ...
@overload
def _from_matrix_orthogonal(matrix: _f64_3d) -> _f64_2d: ...

#
def from_rotvec(rotvec: onp.ToFloatStrict1D | onp.ToFloatStrict2D, degrees: _bint = False) -> _f64_2d: ...
def from_mrp(mrp: onp.ToFloatStrict1D | onp.ToFloatStrict2D) -> _f64_2d: ...

#
def as_quat(quat: _f64_2d, canonical: _bint = False, *, scalar_first: _bint = False) -> _f64_2d: ...
def as_matrix(quat: _f64_2d) -> _f64_3d: ...
def as_rotvec(quat: _f64_2d, degrees: _bint = False) -> _f64_2d: ...
def as_mrp(quat: _f64_2d) -> _f64_2d: ...
def as_euler(quat: _f64_2d, seq: str, degrees: _bint = False, suppress_warnings: _bint = False) -> _f64_2d: ...
def as_davenport(
    quat: _f64_2d, axes: _f64_2d, order: _Order, degrees: _bint = False, suppress_warnings: _bint = False
) -> _f64_2d: ...

#
def inv(quat: _f64_2d) -> _f64_2d: ...

#
@overload
def random(
    num: int, rng: onp.random.ToRNG | None = None, shape: None = None, *, random_state: onp.random.ToRNG | None = None
) -> _f64_2d: ...
@overload
def random(
    num: None = None,
    rng: onp.random.ToRNG | None = None,
    shape: tuple[()] | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> _f64_1d: ...
@overload
def random(
    num: None = None,
    rng: onp.random.ToRNG | None = None,
    *,
    shape: int | tuple[int],
    random_state: onp.random.ToRNG | None = None,
) -> _f64_2d: ...
@overload
def random(
    num: None = None, rng: onp.random.ToRNG | None = None, *, shape: tuple[int, int], random_state: onp.random.ToRNG | None = None
) -> _f64_3d: ...
@overload
def random(
    num: None = None, rng: onp.random.ToRNG | None = None, *, shape: onp.AtLeast3D, random_state: onp.random.ToRNG | None = None
) -> _f64_nd: ...

#
@overload
def identity(num: int, shape: None = None) -> _f64_2d: ...
@overload
def identity(num: None, shape: tuple[()] | None = None) -> _f64_1d: ...
@overload
def identity(num: None, shape: int | tuple[int]) -> _f64_2d: ...
@overload
def identity(num: None, shape: tuple[int, int]) -> _f64_3d: ...
@overload
def identity(num: None, shape: onp.AtLeast3D) -> _f64_nd: ...

#
def magnitude(quat: _f64_2d) -> _f64_2d: ...
def approx_equal(quat: _f64_2d, other: _f64_2d, atol: float | None = None, degrees: _bint = False) -> _bool_1d: ...

#
@overload  # actually returns `_cyutility._memoryviewslice` if `axis=()`, but we cannot express that here
def mean(quat: _f64_2d, weights: _f64_1d | None = None, *, axis: tuple[()]) -> _f64_2d: ...
@overload
def mean(quat: _f64_2d, weights: _f64_1d | None = None, axis: Literal[0] | tuple[Literal[0]] | None = None) -> _f64_1d: ...

#
@overload  # left: None, right: None
def reduce(quat: _f64_2d, left: None = None, right: None = None) -> tuple[_f64_2d, None, None]: ...
@overload  # left: f64[:, :]
def reduce(quat: _f64_2d, left: _f64_2d, right: _f64_2d | None = None) -> tuple[_f64_2d, _isize_2d, _isize_2d]: ...
@overload  # right: f64[:, :] (positional)
def reduce(quat: _f64_2d, left: _f64_2d | None, right: _f64_2d) -> tuple[_f64_2d, _isize_2d, _isize_2d]: ...
@overload  # right: f64[:, :] (keyword)
def reduce(quat: _f64_2d, left: _f64_2d | None = None, *, right: _f64_2d) -> tuple[_f64_2d, _isize_2d, _isize_2d]: ...

#
def apply(quat: _f64_2d, vectors: _f64_2d, inverse: _bint = False) -> _f64_2d: ...
def setitem(quat: _f64_2d, value: _f64_2d, indexer: _Indexer) -> _f64_2d: ...

#
@overload
def align_vectors(
    a: onp.ToFloat1D | onp.ToFloat2D,
    b: onp.ToFloat1D | onp.ToFloat2D,
    weights: onp.ToFloatStrict1D | None = None,
    return_sensitivity: Literal[False, 0] = False,
) -> tuple[_f64_2d, float, None]: ...
@overload
def align_vectors(
    a: onp.ToFloat1D | onp.ToFloat2D,
    b: onp.ToFloat1D | onp.ToFloat2D,
    weights: onp.ToFloatStrict1D | None = None,
    *,
    return_sensitivity: Literal[True, 1],
) -> tuple[_f64_2d, float, _f64_2d]: ...

#
def pow(quat: _f64_2d, n: int) -> _f64_2d: ...

#
@overload  # workaround for https://github.com/microsoft/pyright/issues/10232
def from_davenport(
    axes: onp.ToIntStrict1D | onp.ToIntStrict2D,
    order: _Order,
    angles: onp.Array[tuple[Never, Never, Never], npc.floating],
    degrees: _bint = False,
) -> _f64_1d: ...
@overload
def from_davenport(
    axes: onp.ToIntStrict1D | onp.ToIntStrict2D, order: _Order, angles: onp.ToFloatStrict1D, degrees: _bint = False
) -> _f64_1d: ...
@overload
def from_davenport(
    axes: onp.ToIntStrict1D | onp.ToIntStrict2D, order: _Order, angles: onp.ToFloatStrict2D, degrees: _bint = False
) -> _f64_2d: ...
